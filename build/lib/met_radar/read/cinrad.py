# Module: Read CINRAD radar binary file(.bin/.bin.bz2)
# Author: bo_fan@qq.com(Thanks for any bug report)  
# Version: 20230910 alpha, D:\Seafile\Share\notebooks\20230601ReadFile\CINRAD_test_read.ipynb
#          20230911 beta, D:\Seafile\Share\notebooks\20230601ReadFile\CINRAD_read.py
#          20260406 debug 
import numpy as np
import struct
import warnings

class cinrad(object):   
    """
    Class for accessing all info and data in a CINRAD (WSR-98D) binary file.

    CINRAD files [1]_, also know as WSR-98D weather radar base data. 
    This class supports encoding the whole format of a CINRAD file.
    Files with uncompressed messages and compressed messages are supported. 

    Parameters
    ----------
    inpath : str
        Path of CINRAD file to read.

    Attributes
    ----------
    info: list
        Comman blocks(generic, site, task and cuts) in the file.
    rec : list
        Radial blocks, including header and data(moment_header, moment_data), in the file.

    References
    ----------
    .. [1] https://www.cma.gov.cn/zfxxgk/gknr/flfgbz/bz/202307/t20230712_5642881.html
    .. [2] pyart.io.nexrad_level2
    .. [3] cinrad.io.StandardData

    """
    def __init__(self,inpath):
        with _prepare_file(inpath) as file:
            buf=file.read()
        pos=0
        header=_unpack_from_buf(buf,pos,GENERIC_HEADER)
        pos=pos+_structure_size(GENERIC_HEADER)
        site=_unpack_from_buf(buf,pos,SITE_CONFIG)
        pos=pos+_structure_size(SITE_CONFIG)
        task=_unpack_from_buf(buf,pos,TASK_CONFIG)
        pos=pos+_structure_size(TASK_CONFIG)
        cut=[]
        for i in range(task['cut_number']):
            cut.append(_unpack_from_buf(buf,pos,CUT_CONFIG))
            pos=pos+_structure_size(CUT_CONFIG)

        info={'header':header,'site':site,'task':task,'cut':cut}
        rec={'header':[],'data':[]}
        i=0

        while pos<len(buf):
            rec['header'].append(_unpack_from_buf(buf,pos,RADIAL_HEADER))
            pos=pos+_structure_size(RADIAL_HEADER)
            rec['data'].append({'moment_header':[],'moment_data':[]})

            for j in range(rec['header'][i]['moment_number']):
                moment_header=_unpack_from_buf(buf,pos,MOMENT_HEADER)
                pos=pos+_structure_size(MOMENT_HEADER)

                ngates=int(moment_header['block_length']/moment_header['bin_length'])
                if moment_header['bin_length'] == 2:
                    moment_data=np.frombuffer(buf[pos: pos + ngates * 2], 'u2')
                    pos=pos+ ngates * 2
                elif moment_header['bin_length'] == 1:
                    moment_data=np.frombuffer(buf[pos: pos + ngates], 'u1')
                    pos=pos+ ngates
                else:
                    warnings.warn(
                        'Unsupported bit size: %s. Returning array dtype "B"' % moment_header['bin_length'])
                rec['data'][i]['moment_header'].append(moment_header)
                rec['data'][i]['moment_data'].append(moment_data)
            i=i+1
        print('Reading: %s'%inpath)
        print('First radial num: %d \nLast radial num: %d'%(0,i-1))
        print('All info & data in dict: "info" and "rec"\n')
        
        # raw 
        self.info=info
        self.rec=rec

        # spherical
        self.data,self.gate,self.azimuth,self.elevation=self.get_data()

        # cartesian
        self.spherical2cartesian=spherical2cartesian

    def get_data(self):
        '''
        input: cinrad object
        output: data, gate, azimuth, elevation
        '''
        import numpy as np
        nrays=len(self.rec['data'])
        nswps=self.info['task']['cut_number']
        self.nswps=nswps
        # data={'elevation':[],'azimuth':[]}
        rays_each_swp=[[] for _ in range(nswps)]
        for nray in range(nrays):
            swp=self.rec['header'][nray]['elevation_number']
            rays_each_swp[swp-1].append(nray)
            

        data=[[] for _ in range(nswps)]
        gate=[[] for _ in range(nswps)]
        elevation=[[] for _ in range(nswps)]
        azimuth=[[] for _ in range(nswps)]

        for nswp in range(nswps):
            field_list=[]
            for nray in rays_each_swp[nswp]:
                ele=self.rec['header'][nray]['elevation']
                azi=self.rec['header'][nray]['azimuth']
                # print('nswp=%d nray=%d'%(nswp,nray),'ele=%.1f azi=%.1f'%(ele,azi))
                elevation[nswp].append(ele)
                azimuth[nswp].append(azi)
                if nray==rays_each_swp[nswp][0]:
                    data[nswp]={}
                    gate[nswp]={}
                for nvar in range(len(self.rec['data'][nray]['moment_header'])):
                    field=data_type[self.rec['data'][nray]['moment_header'][nvar]['data_type']]
                    tmp_data=self.rec['data'][nray]['moment_data'][nvar].copy()
                    scale=self.rec['data'][nray]['moment_header'][nvar]['scale']
                    offset=self.rec['data'][nray]['moment_header'][nvar]['offset']
                    mask=tmp_data< 5 # int 0-4 with special meanings
                    tmp_data = np.ma.array((tmp_data - offset) / scale, mask=mask)    
                    if self.info['cut'][nswp]['log_reso']!=self.info['cut'][nswp]['log_reso']:
                        print('Error in range reso!')
                    ngates=len(self.rec['data'][nray]['moment_data'][nvar])
                    tmp_gate=np.array(range(ngates))*self.info['cut'][nswp]['log_reso']+\
                            self.info['cut'][nswp]['start_range']
                    if field not in field_list:
                        gate[nswp][field]=[]
                        data[nswp][field]=[]
                        field_list.append(field)
                    data[nswp][field].append(tmp_data)
                    if len(tmp_gate)> len(gate[nswp][field]):
                        gate[nswp][field]=tmp_gate
                # print(field_list)
            for field in field_list:
                data[nswp][field]=np.ma.vstack(data[nswp][field])
        return data,gate,azimuth,elevation
    
    def get_cut(self,field='REF',nswp=0): # sp for class cinrad
        '''
        input: field name, cut number
        output: range, azimuth, elevation, field data
        '''
        ran=self.gate[nswp][field]
        azi=self.azimuth[nswp]
        ele=self.elevation[nswp]
        field_data=self.data[nswp][field]
        return ran,azi,ele,field_data

def get_cut(data,gate,elevation,azimuth,field='REF',nswp=0):
    '''
    input: data, gate, elevation, azimuth, field name, cut number
    output: range, azimuth, elevation, field data
    '''
    ran=gate[nswp][field]
    azi=azimuth[nswp]
    ele=elevation[nswp]
    field_data=data[nswp][field]
    return ran,azi,ele,field_data

def spherical2cartesian(ran,azi,ele,field_data=None):
    '''
    output: range, azimuth, elevation, field data
    output: xx(km), yy(km), cc
    '''
    from pyart.core.transforms import antenna_vectors_to_cartesian
    x,y,z=antenna_vectors_to_cartesian(ran,azi,ele)
    xx,yy,cc=x/1000,y/1000,field_data
    return xx,yy,cc

 
def _prepare_file(inpath):
    # function from cinrad.io.base

    import bz2
    import gzip
    if hasattr(inpath, "read"):
        return inpath
    f = open(inpath, "rb")
    magic = f.read(3)
    f.close()
    if magic.startswith(b"\x1f\x8b"):
        return gzip.GzipFile(inpath, "rb")
    if magic.startswith(b"BZh"):
        return bz2.BZ2File(inpath, "rb")
    return open(inpath, "rb") 

def _structure_size(structure):
    """ Find the size of a structure in bytes. """
    fmt = '=' + ''.join([i[1] for i in structure])   # CINRAD is small-endian # debug in 2026-Apr-6
    return struct.calcsize(fmt)

def _unpack_structure(string, structure):
    """ Unpack a structure from a string. """
    fmt = '=' + ''.join([i[1] for i in structure])  # CINRAD is small-endian # debug in 2026-Apr-6
    lst = struct.unpack(fmt, string)
    return dict(zip([i[0] for i in structure], lst))

def _unpack_from_buf(buf, pos, structure):
    """ Unpack a structure from a buffer. """
    size = _structure_size(structure)
    return _unpack_structure(buf[pos:pos + size], structure)
    
data_type={
    1:'dBT',
    2:'REF',#'dBZ',
    3:'VEL',#'V',
    4:'W',
    5:'SQI',
    6:'CPA',
    7:'ZDR',
    8:'LDR',
    9:'RHO',
    10:'PHI',
    11:'KDP',
    12:'CP',
    14:'HCL',
    15:'CF',
    16:'SNRH',
    17:'SNRV',
    19:'POTS',
    21:'COP',
    26:'VELSZ',
    27:'DR',
    32:'Zc',
    33:'Vc',
    34:'Wc',
    35:'ZDRc'
}

GENERIC_HEADER = (
    ("magic_number", "4s"),
    ("major_version", "H"),
    ("minor_version", "H"),
    ("generic_type", "I"),
    ("product_type", "I"),
    ("res1", "16s"),
)

SITE_CONFIG = (
    ("site_code", "8s"),
    ("site_name", "32s"),
    ("Latitude", "f"),
    ("Longitude", "f"),
    ("antenna_height", "I"),
    ("ground_height", "I"),
    ("frequency", "f"),
    ("beam_width_hori", "f"),
    ("beam_width_vert", "f"),
    ("RDA_version", "I"),
    ("radar_type", "H"),
    ("antenna_gain", "H"),
    ("trans_loss", "H"),
    ("recv_loss", "H"),
    ("other_loss", "H"),
    ("res2", "46s"),
)

TASK_CONFIG = (
    ("task_name", "32s"),
    ("task_dsc", "128s"),
    ("polar_type", "I"),
    ("scan_type", "I"),
    ("pulse_width", "I"),
    ("scan_start_time", "I"),
    ("cut_number", "I"),
    ("hori_noise", "f"),
    ("vert_noise", "f"),
    ("hori_cali", "f"),
    ("vert_cali", "f"),
    ("hori_tmp", "f"),
    ("vert_tmp", "f"),
    ("ZDR_cali", "f"),
    ("PHIDP_cali", "f"),
    ("LDR_cali", "f"),
    ("res3", "40s"),
)

CUT_CONFIG = (
    ("process_mode", "I"),
    ("wave_form", "I"),
    ("PRF1", "f"),
    ("PRF2", "f"),
    ("dealias_mode", "I"),
    ("azimuth", "f"),
    ("elev", "f"),
    ("start_angle", "f"),
    ("end_angle", "f"),
    ("angular_reso", "f"),
    ("scan_spd", "f"),
    ("log_reso", "I"),
    ("dop_reso", "I"),
    ("max_range1", "I"),
    ("max_range2", "I"),
    ("start_range", "I"),
    ("sample1", "I"),
    ("sample2", "I"),
    ("phase_mode", "I"),
    ("atmos_loss", "f"),
    ("nyquist_spd", "f"),
    ("moments_mask", "Q"), # i8 LONG Q->8 B 
    ("moments_size_mask", "Q"), # Don't use L, for 4 B in Win but 8 B in Mac, unsafe! ############## debug in 2026-Apr-6
    ("misc_filter_mask", "I"),
    ("SQI_thres", "f"),
    ("SIG_thres", "f"),
    ("CSR_thres", "f"),
    ("LOG_thres", "f"),
    ("CPA_thres", "f"),
    ("PMI_thres", "f"),
    ("DPLOG_thres", "f"),
    ("res_thres", "4s"),
    ("dBT_mask", "I"),
    ("dBZ_mask", "I"),
    ("vel_mask", "I"),
    ("sw_mask", "I"),
    ("DP_mask", "I"),
    ("res_mask", "12s"),
    ("scan_sync", "4s"),
    ("direction", "I"),
    ("ground_clutter_classifier_type", "H"),
    ("ground_clutter_filter_type", "H"),
    ("ground_clutter_filter_notch_width", "H"),
    ("ground_clutter_filter_window", "H"),
    ("res4", "72s"),
)

RADIAL_HEADER = (
    ("radial_state", "I"),
    ("spot_blank", "I"),
    ("seq_number", "I"),
    ("radial_number", "I"),
    ("elevation_number", "I"),
    ("azimuth", "f"),
    ("elevation", "f"),
    ("seconds", "I"),
    ("microseconds", "I"),
    ("data_length", "I"),
    ("moment_number", "I"),
    ("res5", "H"),
    ("hori_est_noise", "H"),
    ("vert_est_noise", "H"),
    ("zip_type", "s"),
    ("res6", "13s"),
)

MOMENT_HEADER = (
    ("data_type", "I"),
    ("scale", "I"),
    ("offset", "I"),
    ("bin_length", "H"),
    ("flags", "H"),
    ("block_length", "I"),
    ("res", "12s"),
)


