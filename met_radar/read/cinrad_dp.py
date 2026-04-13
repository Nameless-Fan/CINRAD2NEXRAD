# -*- coding: utf-8 -*-
# @Author: Hao Huang
# @Date:   2016-04-19 12:12:57
# @Last Modified by:   Hao Huang
# @Last Modified time: 2021-7-5
'''
pywxrad.io.cinrad_dp
===========================

Utilities for reading cinrad dp BaseData files.

.. autosummary::
    :toctree: generated/

    read_cinrad_dp
'''

import bz2
import gzip
import warnings
from datetime import datetime
import numpy as np
from pyart.config import FileMetadata, get_fillvalue, get_field_name
from pyart.core.radar import Radar
from pyart.io.common import make_time_unit_str
from pyart.io.nexrad_archive import _interpolate_scan
from scipy import constants

PT = {0: 'horizontal', 1: 'horizontal', 2: 'vertical', 3: 'hv_sim', 4: 'hv_alt'}
SCAN_MODE_NAMES = ['azimuth_surveillance', 'ppi', 'rhi',
                   'sector', 'sector', 'elevation_surveillance', 'manual_ppi']
SCAN_TPYE_NAMES = ['ppi', 'ppi', 'rhi', 'sector', 'sector', 'rhi', 'ppi']

MOMENT_DATA_DICT = {
    1: get_field_name('total_power'),
    2: get_field_name('reflectivity'),
    3: get_field_name('velocity'),
    4: get_field_name('spectrum_width'),
    5: get_field_name('normalized_coherent_power'),
    7: get_field_name('differential_reflectivity'),
    8: get_field_name('linear_depolarization_ratio'),
    9: get_field_name('cross_correlation_ratio'),
    10: get_field_name('differential_phase'),
    11: get_field_name('specific_differential_phase'),
    16: get_field_name('signal_to_noise_ratio'),
}


def read_cinrad_dp(filename):
    '''Read a cinrad dp basedata file.

    Parameters
    ==========

    filename : str
        Name of cinrad dualpol Basedata file.

    Returns
    =======

    radar : Radar
        Radar object containing data from cinrad dualpol_level2 file.

    '''
    # create metadata retrival object
    filemetadata = FileMetadata('cfradial')

    # read data
    cdp = CinradDP(filename)
    cdp.parse()

    # time
    time = filemetadata('time')
    tdata = np.array(cdp.time)
    min_time = np.floor(tdata.min())
    time['data'] = (tdata - min_time).astype(np.float64)
    # Return a time unit string from a datetime object
    time['units'] = make_time_unit_str(datetime.utcfromtimestamp(min_time))

    # range
    ngates = int(np.round((cdp.max_last_gate - cdp.min_first_gate) / cdp.min_gate_spacing)) + 1
    if not np.isclose(cdp.min_gate_spacing, (cdp.max_last_gate - cdp.min_first_gate) / (ngates - 1)):
        raise RuntimeError('cdp.max_last_gate - cdp.min_first_gate is not int times of cdp.min_first_gate')

    _range = filemetadata('range')
    _range['data'] = np.linspace(cdp.min_first_gate, cdp.max_last_gate, ngates).astype(np.float32)
    _range['meters_between_gates'] = cdp.min_gate_spacing
    _range['meters_to_center_of_first_gate'] = 0.0

    sitename = cdp.site_name
    if isinstance(sitename, bytes):
        sitename = sitename.decode('utf-8', 'backslashreplace')
    for i in range(len(sitename)):
        if not sitename[i].isalnum():
            sitename = sitename[:i]
            break

    # metadata
    metadata = filemetadata('metadata')
    metadata['instrument_name'] = 'RADAR'
    metadata['instrument_type'] = 'radar'
    metadata['site_name'] = sitename
    metadata['original_container'] = 'CINRAD DP'
    # Where the original data were produced
    metadata['institution'] = 'Nanjing University'
    metadata['platform_is_mobile'] = 'false'
    metadata['n_gates_vary'] = 'false'
    metadata['ray_times_increase'] = 'true'

    # longitude, latitude, altitude
    latitude, longitude, altitude = [filemetadata(_) for _ in ['latitude', 'longitude', 'altitude']]

    latitude['data'] = np.array([cdp.radar_info['latitude']], np.float64)
    longitude['data'] = np.array([cdp.radar_info['longitude']], np.float64)
    altitude['data'] = np.array([cdp.radar_info['altitude']], np.float64)

    # sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
    # sweep_end_ray_index
    sweep_number = filemetadata('sweep_number')
    sweep_mode = filemetadata('sweep_mode')
    fixed_angle = filemetadata('fixed_angle')
    sweep_start_ray_index = filemetadata('sweep_start_ray_index')
    sweep_end_ray_index = filemetadata('sweep_end_ray_index')
    sweep_number['data'] = np.arange(cdp.num_sweeps, dtype=np.int32)
    sweep_mode['data'] = np.array([cdp.scan_mode] * cdp.num_sweeps)
    fixed_angle['data'] = np.array(cdp.fixed_angle, dtype=np.float32)
    sweep_start_ray_index['data'] = np.array(cdp.sweep_start_ray_index, dtype=np.int32)
    sweep_end_ray_index['data'] = np.array(cdp.sweep_end_ray_index, dtype=np.int32)

    # azimuth, elevation
    azimuth, elevation = [filemetadata(_) for _ in ['azimuth', 'elevation']]
    azimuth['data'] = np.array(cdp.azimuth, dtype=np.float32)
    elevation['data'] = np.array(cdp.elevation, dtype=np.float32)

    # instrument parameters
    instrument_parameters = {}
    for key in cdp.instrument_parameters.keys():
        instrument_parameters[key] = filemetadata(key)
        instrument_parameters[key]['data'] = cdp.instrument_parameters[key]

    # fields
    fields = {}
    for k in cdp.datatype_set:
        field_name = MOMENT_DATA_DICT[k]
        field_dic = filemetadata(field_name)
        field_data = np.ma.masked_all((cdp.nrays, ngates), dtype=np.float32)
        for iray in range(cdp.nrays):
            if k not in cdp.moment_radials[iray]:
                continue
            tmp = cdp.moment_radials[iray][k]
            field_data[iray, :tmp.size] = tmp
        field_data = np.ma.masked_invalid(field_data)
        for iswp in range(cdp.num_sweeps):
            if k in [3, 4]:
                if np.isclose(cdp.doppler_gate_spacing[iswp], 4 * cdp.min_gate_spacing):
                    warnings.warn("Gate spacing is not constant, interpolating data (4 times) in " +
                                  "scans %s for moment %s." % (iswp, field_name), UserWarning)
                    _interpolate_scan(field_data, cdp.sweep_start_ray_index[iswp], cdp.sweep_end_ray_index[iswp],
                                      cdp.doppler_gates[k][iswp], '4', linear_interp=False)
                elif np.isclose(cdp.doppler_gate_spacing[iswp], 2 * cdp.min_gate_spacing):
                    warnings.warn("Gate spacing is not constant, interpolating data (2 times) in " +
                                  "scans %s for moment %s." % (iswp, field_name), UserWarning)
                    _interpolate_scan(field_data, cdp.sweep_start_ray_index[iswp], cdp.sweep_end_ray_index[iswp],
                                      cdp.doppler_gates[k][iswp], '2', linear_interp=False)
                elif np.isclose(cdp.doppler_gate_spacing[iswp], cdp.min_gate_spacing):
                    pass
                else:
                    raise IOError('doppler_gate_spacing at sweep{} should be 1, 2, or 4 times dr '.format(iswp))
            else:
                if np.isclose(cdp.log_gate_spacing[iswp], 4 * cdp.min_gate_spacing):
                    warnings.warn("Gate spacing is not constant, interpolating data (4 times) in " +
                                  "scans %s for moment %s." % (iswp, field_name), UserWarning)
                    _interpolate_scan(field_data, cdp.sweep_start_ray_index[iswp], cdp.sweep_end_ray_index[iswp],
                                      cdp.log_gates[k][iswp], '4', linear_interp=True)
                elif np.isclose(cdp.log_gate_spacing[iswp], 2 * cdp.min_gate_spacing):
                    warnings.warn("Gate spacing is not constant, interpolating data (2 times) in " +
                                  "scans %s for moment %s." % (iswp, field_name), UserWarning)
                    _interpolate_scan(field_data, cdp.sweep_start_ray_index[iswp], cdp.sweep_end_ray_index[iswp],
                                      cdp.log_gates[k][iswp], '2', linear_interp=True)
                elif np.isclose(cdp.log_gate_spacing[iswp], cdp.min_gate_spacing):
                    pass
                else:
                    raise IOError('log_gate_spacing at sweep{} should be 1, 2, or 4 times dr '.format(iswp))

        np.ma.set_fill_value(field_data, get_fillvalue())
        field_dic['data'] = field_data
        field_dic['_FillValue'] = get_fillvalue()
        fields[field_name] = field_dic

    if get_field_name('differential_phase') in fields.keys():
        # convert PHIDP to be within valid_min ~ valid_max
        phidp = fields[get_field_name('differential_phase')]
        valid_max = 180.
        valid_min = -180.
        phi_range = valid_max - valid_min

        while np.ma.any(phidp['data'] < valid_min):
            phidp['data'][np.ma.filled(phidp['data'] < valid_min, False)] += phi_range
        while np.ma.any(phidp['data'] > valid_max):
            phidp['data'][np.ma.filled(phidp['data'] > valid_max, False)] -= phi_range
        fields[get_field_name('differential_phase')] = phidp

    if get_field_name('cross_correlation_ratio') in fields.keys():
        rhohv = fields[get_field_name('cross_correlation_ratio')]
        rhohv['data'][rhohv['data'] > 1] = 1.
        fields[get_field_name('cross_correlation_ratio')] = rhohv

    return Radar(
        time, _range, fields, metadata, cdp.scan_types,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index, azimuth, elevation,
        instrument_parameters=instrument_parameters)


class CinradDP(object):
    '''A file object for cinrad dp level2 data.

    Parameters
    ----------
    filename : str or file-like.
        Name of cinrad dp file to read or a file-like object pointing
        to the beginning of such a file.

    Attributes
    ----------
    num_sweeps : int
        Number of sweeps in the volume.
    time : list of ints
        Time in seconds in epoch for each ray in the volume.
    azimuth : list of floats
        Azimuth angle for each ray in the volume in degrees.
    elevation : list of floats
        Elevation angle for each ray in the volume in degrees.
    fixed_angle : list of floats
        Fixed angles for each sweep.
    scan_types : list of ints
        scan type for each sweep.
    radar_info : dict
        Radar information recorded in the file.

    '''

    def __init__(self, filename):
        self.filename = filename

        self.site_name = ''
        self.nrays = None
        self.num_sweeps = None

        self.min_gate_spacing = None
        self.min_first_gate = None
        self.max_last_gate = None
        self.log_gate_spacing = None
        self.doppler_gate_spacing = None
        self.log_gates = None
        self.doppler_gates = None

        self.time = []
        self.azimuth = []
        self.elevation = []

        self.fixed_angle = None
        self.scan_types = None
        self.scan_mode = None
        self.sweep_start_ray_index = []
        self.sweep_end_ray_index = []

        self.radar_info = {}
        self.instrument_parameters = {}

        self.moment_radials = []
        self.datatype_set = None

        # private attributes

        self._fp = 0
        self._buf_size = None
        self._buf = None

    def _read_data(self, file):
        if hasattr(file, 'read'):
            buf = file.read()
            file.close()
            if buf[:4] != b'RSTM':
                raise IOError('input file is not a cinrad DP BaseData')
            else:
                return buf

        with open(file, 'rb') as f:
            begin = f.read(12)
        if begin[:4] == b'RSTM':
            with open(file, 'rb') as f:
                buf = f.read()
        elif begin[:3] == b'BZh':
            f = bz2.BZ2File(file, 'rb')
            buf = self._read_data(f)
        elif begin[:2] == b'\x1f\x8b':
            f = gzip.open(file, 'rb')
            buf = self._read_data(f)
        return buf

    def _init_file_and_func(self):
        self._buf = self._read_data(self.filename)

        self._buf_size = len(self._buf)
        self._fp = 0

    def _unpack_dtype_s(self, data):
        return {key: data[key] for key in data.dtype.fields.keys()}

    def _np_unpack_buf(self, datatype, count=1):
        temp = np.frombuffer(self._buf[self._fp:self._fp + datatype.itemsize * count], dtype=datatype, count=count)
        self._fp += datatype.itemsize * count
        return [self._unpack_dtype_s(temp[i]) for i in range(count)]

    def _skip_buf(self, datatype, count=1):
        self._fp += datatype.itemsize * count

    def _np_unpack_data_buf(self, datatype, count=1):
        temp = np.frombuffer(self._buf[self._fp:self._fp + datatype.itemsize * count], dtype=datatype, count=count)
        self._fp += datatype.itemsize * count
        return temp

    def parse(self):
        self._init_file_and_func()

        unpack_dtype = self._np_unpack_buf
        skip_dtype = self._skip_buf
        unpack_data = self._np_unpack_data_buf
        dtypes = define_cinrad_dp_dtype()
        valid_moment_types = list(MOMENT_DATA_DICT.keys())

        # 1.1.  Read 32 Bytes Generic Header Block
        # common_header = unpack_dtype(dtypes['header'])[0]
        skip_dtype(dtypes['header'])

        # 1.2.Read 128 Bytes SITE CONFIG
        site_info = unpack_dtype(dtypes['site'])[0]
        # 1.3. Read 256 Bytes TASK CONFIG
        task_info = unpack_dtype(dtypes['task'])[0]
        # 1.4. Read 256*N Bytes CUT CONFIG
        cut_config = unpack_dtype(dtypes['cut'], count=task_info['CutNumber'])

        iray, iswp = 0, 0
        datatype_set = set([])
        self.doppler_gates = {_: np.zeros((task_info['CutNumber'],), dtype=int) for _ in valid_moment_types if
                              _ in [3, 4]}
        self.log_gates = {_: np.zeros((task_info['CutNumber'],), dtype=int) for _ in valid_moment_types if
                          _ not in [3, 4]}
        while True:
            radial = unpack_dtype(dtypes['radial_header'])[0]
            if radial['MomentNumber'] == 0:
                continue

            self.time.append(radial['Seconds'] + radial['microseconds'] / 1e6)
            self.azimuth.append(radial['Azimuth'])
            self.elevation.append(radial['Elevation'])

            moments = {}
            for jj in range(radial['MomentNumber']):
                # read moments
                moment_header = unpack_dtype(dtypes['radial_moment_header'])[0]
                tempdata = unpack_data(np.dtype(np.uint8), count=moment_header['Length'])

                # avoid further process if not in our dict
                if moment_header['DataType'] not in valid_moment_types:
                    continue
                datatype_set.add(moment_header['DataType'])

                # resort
                if moment_header['BinLength'] == 2:
                    tempdata2 = np.float32(tempdata[1::2] * 256. + tempdata[0::2])
                else:
                    tempdata2 = np.float32(tempdata)
                # tempdata = np.ma.masked_less(tempdata, 5.) # cost much more time than np.nan
                tempdata2[tempdata2 < 5] = np.nan
                # moment_data = np.float32(
                #     (tempdata2 - np.float64(moment_header['Offset'])) / np.float64(moment_header['Scale']))
                moment_data = (tempdata2 - moment_header['Offset']) / np.float32(moment_header['Scale'])
                moments.update({moment_header['DataType']: moment_data})

            self.moment_radials.append(moments)

            if radial['RadialState'] in [0, 3]:
                self.sweep_start_ray_index.append(iray)
                for k in moments.keys():
                    if k in [3, 4]:  # velocity
                        self.doppler_gates[k][iswp] = moments[k].size
                    else:
                        self.log_gates[k][iswp] = moments[k].size
                iswp += 1
            # Cut END or #Volume END
            if radial['RadialState'] in [2, 4]:
                self.sweep_end_ray_index.append(iray)
            iray += 1

            if task_info['ScanType'] == 2 and radial['RadialState'] in [2, 4]:
                break
            if task_info['ScanType'] == 3 and radial['RadialState'] in [2, 4]:
                break
            if self._fp >= self._buf_size:
                break

        # check sweep_index
        if self.sweep_start_ray_index[0] != 0 or self.sweep_end_ray_index[-1] != iray - 1:
            raise RuntimeError('cut start/end index error')
        self.nrays = iray
        self.num_sweeps = iswp
        self.datatype_set = datatype_set

        self.scan_types = SCAN_TPYE_NAMES[task_info['ScanType']]
        self.scan_mode = SCAN_MODE_NAMES[task_info['ScanType']]

        # resolutions are consistent in each cut, the gate numbers would also be consistent
        self.log_gate_spacing = [cut_config[_]['LogResolution'] for _ in range(self.num_sweeps)]
        self.doppler_gate_spacing = [cut_config[_]['DopplerResolution'] for _ in range(self.num_sweeps)]
        # to keep consistency, we do not use this
        # start_range = [cut_config[_]['StartRange'] for _ in range(self.num_sweeps)]

        self.min_gate_spacing = np.min(self.log_gate_spacing + self.doppler_gate_spacing)
        self.min_first_gate = self.min_gate_spacing
        max_last_gate = 0
        for k in valid_moment_types:
            if k in [3, 4]:
                max_last_gate = max(max_last_gate, np.max(
                    [self.doppler_gates[k][i] * self.doppler_gate_spacing[i] for i in range(self.num_sweeps)]))
            else:
                max_last_gate = max(max_last_gate, np.max(
                    [self.log_gates[k][i] * self.log_gate_spacing[i] for i in range(self.num_sweeps)]))
        self.max_last_gate = max_last_gate

        if self.scan_types in ["ppi", "sector"]:
            self.fixed_angle = np.array([cut_config[_]['Elevation'] for _ in range(self.num_sweeps)])
        elif self.scan_types == "rhi":
            self.fixed_angle = np.array([cut_config[_]['Azimuth'] for _ in range(self.num_sweeps)])

        # location attributions
        self.radar_info.update({'latitude': site_info['Latitude'],
                                'longitude': site_info['Longitude'],
                                'altitude': site_info['Height']})

        # instrument_parameters
        inst_param = {}
        # Pulsing mode Options are: “fixed”, “staggered”, “dual” Assumed “fixed” if missing.
        prt_mode = []
        for _ in range(self.num_sweeps):
            if cut_config[_]['WaveForm'] == 5:
                prt_mode.append('dual')
            else:
                prt_mode.append('fixed')
        inst_param['prt_mode'] = np.array(prt_mode)  # generate string array

        # Polarization Type: Type  1 - Horizontal 2 - Vertical 3 - Simultaneously 4 - Alternation
        inst_param['polarization_mode'] = np.array([PT[task_info['PolarizationType']]] * self.num_sweeps)
        inst_param['pulse_width'] = np.ones((self.nrays,), dtype=np.float32) * task_info['PulseWidth'] * 1e-9
        inst_param['prt'] = np.empty((self.nrays,), dtype=np.float32)
        inst_param['prt_ratio'] = np.empty((self.nrays,), dtype=np.float32)
        inst_param['nyquist_velocity'] = np.empty((self.nrays,), dtype=np.float32)
        inst_param['n_samples'] = np.empty(len(self.time), dtype=np.int32)

        for kk in range(self.num_sweeps):
            s = slice(self.sweep_start_ray_index[kk], self.sweep_end_ray_index[kk] + 1)
            inst_param['prt_ratio'][s] = np.float64(cut_config[kk]['PRF2']) / np.float64(cut_config[kk]['PRF1'])
            inst_param['nyquist_velocity'][s] = cut_config[kk]['NyquistSpeed']
            inst_param['n_samples'][s] = cut_config[kk]['Sample1']
        inst_param['unambiguous_range'] = constants.c / 2. * inst_param['prt']
        inst_param['frequency'] = np.array([site_info['Frequency'] * 1e6], np.float32)  # frequency MHZ->HZ
        inst_param['radar_beam_width_h'] = np.array([site_info['BeamWidthHori']], np.float32)
        inst_param['radar_beam_width_v'] = inst_param['radar_beam_width_h']  # something wrong with beam width

        self.instrument_parameters = inst_param


def define_cinrad_dp_dtype():
    # 1.1.  Read 32 Bytes Generic Header Block
    header = np.dtype([
        ('MagicWord', np.int32),  # 0x4D545352 , Magic word for product
        ('MajorVersion', np.int16),  # Major Version
        ('MinorVersion', np.int16),  # Minor Version
        ('GenericType', np.int32),  # Type of data, see Table 2-3
        ('ProductType', np.int32),  # Type  of Product, not used
        ('Reserved', 'S16'),  # Reserved
    ])

    # 1.2.Read 128 Bytes SITE CONFIG
    site = np.dtype([
        ('Code', 'S8'),  # Site Code  in characters
        ('Name', 'S32'),  # Site Name or  description  in characters
        ('Latitude', np.float32),  # Latitude of Radar Site
        ('Longitude', np.float32),  # Longitude of Radar Site
        ('Height', np.int32),  # Height of  antenna in meters
        ('Ground', np.int32),  # Height  of  ground  in meters
        ('Frequency', np.float32),  # Radar operation frequency in MHz
        ('BeamWidthHori', np.float32),  # Antenna  Beam  Width Hori
        ('BeamWidthVert', np.float32),  # Antenna  Beam  Width Vert
        ('Reserved', 'S60'),  # Reserved
    ])

    # 1.3. Read 256 Bytes TASK CONFIG
    task = np.dtype([
        ('Name', 'S32'),  # Name of the Task   Configuration
        ('Description', 'S128'),  # Description of Task
        # Polarization Type: Type  1 - Horizontal 2 - Vertical
        # 3 - Simultaneously 4 - Alternation
        ('PolarizationType', np.int32),
        # Volume Scan Type 0 - PPI Volume Scan 1- Single PPI 2 - Single RHI
        # 3 -Single Sector 4-SectorVolume Scan 5-RHI Volume Scan 6-Manual scan
        ('ScanType', np.int32),
        ('PulseWidth', np.int32),  # Pulse Width ,Nanoseconds
        # Start time  of volume scan , UTC , seconds
        ('VolumeStartTime', np.int32),
        # Number of Elevation or Azimuth cuts in the task
        ('CutNumber', np.int32),
        # Noise level of horizontal channel, dBm
        ('HorizontalNoise', np.float32),
        ('VerticalNoise', np.float32),  # Noise  level  of  vertical channe
        # System Reflectivity Calibration Const for horizontal channel.
        ('HorizontalCalibration', np.float32),
        # System Reflectivity Calibration Const for vertical channel.
        ('VerticalCalibration', np.float32),
        # System Noise Temperature Const for horizontal channel.
        ('HorizontalNoiseTemperature', np.float32),
        # System Noise Temperature Const for vertical channel.
        ('VerticalNoiseTemperature', np.float32),
        # Reflectivity calibration difference of horizontal and vertical channel
        ('ZdrCalibration', np.float32),
        # Phase calibration  difference of horizontal and vertical channel
        ('PhaseCalibration', np.float32),
        # LDR calibration difference of horizontal and vertical channel
        ('LDRCalibration', np.float32),
        ('Reserved', 'S40'),  # Reserved
    ])

    # 1.4. Read 256*N Bytes CUT CONFIG
    cut = np.dtype([
        # Main  processing  mode  of signal processing algorithm. 1 - PPP 2 - FFT
        ('ProcessMode', np.int32),
        # WSR-88D  defined  wave form  0:CS  1:CD  2:CDX  3:Rx Test 4:BATCH
        # 5:Dual PRF 6:Random Phase 7:SZ
        ('WaveForm', np.int32),
        # Pulse  Repetition  Frequency 1. For wave form Batch  and Dual PRF
        # mode, it is the high PRF, for other modes it is the only PRF.
        ('PRF1', np.float32),
        # Pulse  Repetition  Frequency 2. For wave form Batch  and Dual PRF
        # mode, it is the low PRF, for other modes it is not used
        ('PRF2', np.float32),
        # Dual PRF mode 1:Single PRF  2: 3:2 mode   3: 4:3 mode    4: 5:4 mode
        ('UnfoldMode', np.int32),
        ('Azimuth', np.float32),  # Azimuth degree for RHI scan mode,
        ('Elevation', np.float32),  # Elevation degree for PPI scan mode
        # Start azimuth angle for  PPI  Sector mode. Start (High) Elevation
        # for RHI mode.
        ('StartAngle', np.float32),
        # Stop azimuth angle for PPI Sector mode. Stop (Low) Elevation for RHI mode.
        ('EndAngle', np.float32),
        # Azimuth  resolution  for  PPI  scan, Elevation resolution for RHI mode.
        ('AngleResolution', np.float32),
        # Azimuth scan speed for PPI scan, Elevation scan speed for RHI mode
        ('ScanSpeed', np.float32),
        # Range  bin  resolution  for surveillance data, reflectivity and ZDR,etc
        ('LogResolution', np.int32),
        # Range bin resolution for  Doppler data, velocity  and spectrum, etc
        ('DopplerResolution', np.int32),
        ('MaximumRange', np.int32),  # Maximum range of scan
        ('MaximumRange2', np.int32),  # Maximum range of scan
        ('StartRange', np.int32),  # Start range of scan
        # Pulse sampling number #1. For wave form Batch  and Dual PRF mode,
        # it is for high PRF, for other modes it!s for only PRF.
        ('Sample1', np.int32),
        # Pulse sampling number #2. For wave form Batch  and Dual PRF mode,
        # it is for low PRF, for other modes it!s not used
        ('Sample2', np.int32),
        # Phase modulation mode. 1:Fixed Phase  2:Random Phase  3:SZ Phase
        ('PhaseMode', np.int32),
        # two-way  atmospheric  attenuation factor   dB/km
        ('AtmosphericLoss', np.float32),
        ('NyquistSpeed', np.float32),  # m/s
        # uint64  Bit  mask  indicates  which moments are involved in the
        # scan. See Table 2-7
        ('MomentsMask', np.uint64),
        # uint64 Bit  mask indicates range  length for  moment data in
        # Table 2-7. 0 for 1 byte, 1 for 2 bytes
        ('MomentsSizeMask', np.uint64),
        # uint32 注意：这里有问题，感觉要多读4个字节才行，而且MomentsMask
        # 和MomentsSizeMask感觉都不对
        ('FilterMask', np.uint32),
        ('SQIThreshold', np.float32),  # SQI Threshold for the scan
        ('SIGThreshold', np.float32),  # SIG Threshold for the scan
        ('CSRThreshold', np.float32),  # CSR Threshold for the scan
        ('LOGThreshold', np.float32),  # LOG Threshold for the scan
        ('CPAThreshold', np.float32),  # CPA Threshold for the scan
        ('PMIThreshold', np.float32),  # PMI Threshold for the scan
        ('ThresholdsReserved', 'S8'),  # Reserved
        # Thresholds used for total reflectivity data. Bits  mask start
        # from "SQI Threshold" take is as LSB.
        ('dBTMask', np.int32),
        # Thresholds used for  reflectivity data. Bits  mask start from
        # "SQI Threshold" take is as LSB.
        ('dBZMask', np.int32),
        # Thresholds used for  velocity data. Bits  mask start from "SQI
        # Threshold" take is as LSB.
        ('VelocityMask', np.int32),
        # Thresholds used for  Spectrum  Width data. Bits  mask start from
        # "SQI Threshold" take is as LSB.
        ('SpectrumWidthMask', np.int32),
        # Thresholds used for  ZDR data. Bits  mask start from "SQI
        # Threshold" take is as LSB.
        ('ZDRMask', np.int32),
        ('MaskReserved', 'S12'),  # Reserved for mask
        ('ScanSync', np.int32),  # Reserved
        # Antenna rotate direction, 1= clockwise, 2=counter clockwise
        ('Direction', np.int32),
        # 1 - All data is passed   2 - No data is passed 3 Use Real Time GC
        # Classifier  4 use bypassmap
        ('GroundClutterClassifierType', np.int16),
        # 0- none   1 -Adaptive FFT  4 - IIR
        ('GroundClutterFilterType', np.int16),
        ('GroundClutterFilterNotchWidth', np.int16),  # Scaled by 10
        # -1-none 0 - rect 1-  Hamming 2-  Blackman 3-  Adaptive
        ('GroundClutterFilterWindow', np.int16),
        ('Spare', 'S72'),  # Reserved
    ])
    # 2.1 Basedata Radial Header Block 64bytes
    radial_header = np.dtype([
        # 0= Cut Start 1=Intermediate Data 2=Cut End 3=Volume Start 4=Volume End
        ('RadialState', np.int32),
        ('SpotBlank', np.int32),  # 0=Normal 1=Spot Blank
        ('SequenceNumber', np.int32),  # Sequence Number
        ('RadialNumber', np.int32),  # Radial Number for each cut
        ('ElevationNumber', np.int32),  # Elevation Number
        ('Azimuth', np.float32),  # Azimuth Angle
        ('Elevation', np.float32),  # Elevation Angle
        ('Seconds', np.int32),  # Radial data time in second
        # Radial data time in   microsecond (expect   seconds)
        ('microseconds', np.int32),
        # Length of data  in  this radial, this header is excluded
        ('Lengthofdata', np.int32),
        ('MomentNumber', np.int32),  # Moments available in this radial
        ('Reserved', 'S20'),  # Reserved
    ])

    # 2.2 Base Data Moment Header Block
    radial_moment_header = np.dtype([
        ('DataType', np.int32),  # Moment data type, See Table 2-7
        ('Scale', np.int32),  # Data coding scale Code = value*scale+offset
        ('Offset', np.int32),  # Data coding offset Code=value*scale+offset
        ('BinLength', np.int16),  # Bytes to save each bin of data
        ('Flags', np.int16),  # Bit Mask of flags for data. Reserved now
        # Length of data  of  current moment, this header is excluded.
        ('Length', np.int32),
        ('Reserved', 'S12'),  # Reserved
    ])
    dtypes = {'header': header, 'site': site, 'task': task, 'cut': cut,
              'radial_header': radial_header, 'radial_moment_header': radial_moment_header}

    return dtypes


