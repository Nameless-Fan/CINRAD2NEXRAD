# Script: Radar data format transfer: CINRAD format -> NEXRAD level II ar2v format  
# Author: bo_fan@qq.com 
# Version: 20230912 v1.0 D:\Seafile\Share\notebooks\20230601ReadFile\CINRAD2NEXRAD.py
#          2026-Apr-6 improved     
# Script: Radar data format transfer: CINRAD format -> NEXRAD level II ar2v format  
# Author: bofan@smail.nju.edu.cn (Thanks for any bug report)  
# Version: 20230912 v1.0      
# Path: D:\Seafile\Share\notebooks\20230601ReadFile\CINRAD2NEXRAD.py

def cinrad2nexrad(inpath,outpath,site_name_4s=None):
    from ..read.cinrad import cinrad,data_type
    import pyart.io.nexrad_level2 as ar2v
    from pyart.io.nexrad_level2 import _structure_size

    # 1. Read CINRAD radar product
    radar=cinrad(inpath)

    # 2. Write NEXRAD level II product
    ## 2.1 Write volume header 
    volume_header=[] # ar2v.volume_header

    def convert_to_dict(tuple_list):
        # Create a dictionary using the dict() constructor and a list comprehension
        dictionary = dict((key, value) for key, value in tuple_list)
        # Return the completed dictionary
        return dictionary

    volume_header=convert_to_dict(ar2v.VOLUME_HEADER)
    volume_header['tape']='AR2V0006.'
    volume_header['extension']='001'
    volume_header['date']=int(radar.info['task']['scan_start_time']/24/3600)
    volume_header['time']=int((radar.info['task']['scan_start_time']/24/3600-volume_header['date'])*24*3600*1000)
    if site_name_4s is None:
        site_name_4s=str(radar.info['site']['site_name'][0:4],encoding='ASCII')
    volume_header['icao']=site_name_4s #str(radar.info['site']['site_name'][0:4],encoding='ASCII')

    rec=[] # various dict of 'header','msg_header','VOL','ELV','RAD' and other fields
    list_index=[2,3,7,10,9] # REF, VEL, ZDR, PHI, RHO

    for i in range(len(radar.rec['data'])):
        rec.append({})

        ## 2.2.2 Write data 
        for j in ['VOL','ELV','RAD']:
            field=j
            if field=='VOL':
                block=convert_to_dict(ar2v.VOLUME_DATA_BLOCK)
                block['block_type']='R'
                block['data_name']=field
                block['lrtup']=44
                block['version_major']=1
                block['version_minor']=0
                block['lat']=radar.info['site']['Latitude']
                block['lon']=radar.info['site']['Longitude']
                block['height']=radar.info['site']['antenna_height']
                block['feedhorn_height']=radar.info['site']['antenna_height'] # ?
                block['refl_calib']=radar.info['task']['hori_cali'] # ? e.g. KCLX: -44.6 
                block['power_h']=269.2 # ? e.g. KCLX
                block['power_v']=266.7 # ? e.g. KCLX
                block['diff_refl_calib']=radar.info['task']['ZDR_cali'] # e.g. KCLX: 0.35
                block['init_phase']=radar.info['task']['PHIDP_cali']+180 # e.g. KCLX: 60
                block['vcp']=212 
                block['spare']='' # ?
            elif field=='ELV':
                block=convert_to_dict(ar2v.ELEVATION_DATA_BLOCK)
                block['block_type']='R'
                block['data_name']=field
                block['lrtup']=12
                block['atmos']=int(radar.info['cut'][radar.rec['header'][i]['elevation_number']-1]['atmos_loss'])
                block['refl_calib']=radar.info['task']['hori_cali'] # ? e.g. KCLX: -43.8 
            elif field=='RAD':
                block=convert_to_dict(ar2v.RADIAL_DATA_BLOCK)
                block['block_type']='R'
                block['data_name']=field
                block['lrtup']=28
                block['unambig_range']=int(1.9e8/radar.info['cut'][radar.rec['header'][i]['elevation_number']-1]['PRF1']/1000)
                block['noise_h']=radar.info['task']['hori_noise']
                block['noise_v']=radar.info['task']['vert_noise']
                block['nyquist_vel']=int(radar.info['cut'][radar.rec['header'][i]['elevation_number']-1]['nyquist_spd'])
                block['spare']='' # ?
            rec[i][field]=block

        for j in range(len(radar.rec['data'][i]['moment_header'])):
            if radar.rec['data'][i]['moment_header'][j]['data_type'] not in list_index:
                continue
            field=data_type[radar.rec['data'][i]['moment_header'][j]['data_type']]
            block=convert_to_dict(ar2v.GENERIC_DATA_BLOCK)
            block['block_type']='D'
            block['data_name']=field
            block['reserved']=0 # ?
            block['ngates']=int(radar.rec['data'][i]['moment_header'][j]['block_length']/
                                radar.rec['data'][i]['moment_header'][j]['bin_length'])
            block['first_gate']=radar.info['cut'][radar.rec['header'][i]['elevation_number']-1]['start_range']
            if field=='VEL':
                block['gate_spacing']=radar.info['cut'][radar.rec['header'][i]['elevation_number']-1]['dop_reso']
            else:
                block['gate_spacing']=radar.info['cut'][radar.rec['header'][i]['elevation_number']-1]['log_reso']
            block['thresh']=100 # ?
            block['snr_thres']=16 # ?
            block['flags']=0 # ?
            block['word_size']=8*radar.rec['data'][i]['moment_header'][j]['bin_length']
            block['scale']=radar.rec['data'][i]['moment_header'][j]['scale']
            block['offset']=radar.rec['data'][i]['moment_header'][j]['offset']

            rec[i][field]=block
            rec[i][field]['data']=radar.rec['data'][i]['moment_data'][j]

        ## 2.2.1 Write data header 
        msg_header=convert_to_dict(ar2v.MSG_31)
        msg_header['id']=volume_header['icao']
        msg_header['collect_ms']=int((radar.rec['header'][i]['seconds']/24/3600-int(radar.rec['header'][i]['seconds']/24/3600))*24*3600*1000+radar.rec['header'][i]['microseconds'])
        msg_header['collect_date']=int(radar.rec['header'][i]['seconds']/24/3600)
        msg_header['azimuth_number']=radar.rec['header'][i]['radial_number']
        msg_header['azimuth_angle']=radar.rec['header'][i]['azimuth']
        msg_header['compress_flag']=0
        msg_header['spare_0']=0 # ?
        msg_header['radial_length']=[]
        msg_header['azimuth_resolution']=0 # ?
        msg_header['radial_spacing']=radar.rec['header'][i]['radial_state'] # 0-4 the same, >4 wrong
        msg_header['elevation_number']=radar.rec['header'][i]['elevation_number']
        msg_header['cut_sector']=1
        msg_header['elevation_angle']=radar.rec['header'][i]['elevation']
        msg_header['radial_blanking']=0
        msg_header['azimuth_mode']=0
        msg_header['block_count']=0
        msg_header['block_pointer_1']=0
        msg_header['block_pointer_2']=0
        msg_header['block_pointer_3']=0
        msg_header['block_pointer_4']=0   
        msg_header['block_pointer_5']=0
        msg_header['block_pointer_6']=0
        msg_header['block_pointer_7']=0
        msg_header['block_pointer_8']=0
        msg_header['block_pointer_9']=0
        msg_header['block_pointer_10']=0
        k=1
        fields=sorted(list(rec[i].keys()),key=(['VOL', 'ELV', 'RAD']+[data_type[ii] for ii in list_index]).index)
        for field in fields:
            if k>9:
                print('Too many fields(>10) to record!')
                break
            if field=='VOL' and k==1:
                msg_header['block_pointer_1']=_structure_size(ar2v.MSG_31)
            elif field=='ELV' and k==2:
                msg_header['block_pointer_2']=msg_header['block_pointer_1']+_structure_size(ar2v.VOLUME_DATA_BLOCK)
            elif field=='RAD' and k==3:
                msg_header['block_pointer_3']=msg_header['block_pointer_2']+_structure_size(ar2v.ELEVATION_DATA_BLOCK)
            elif k==4:
                msg_header['block_pointer_4']=msg_header['block_pointer_3']+_structure_size(ar2v.RADIAL_DATA_BLOCK)
            else:
                msg_header['block_pointer_'+str(k)]=msg_header['block_pointer_'+str(k-1)]+_structure_size(ar2v.GENERIC_DATA_BLOCK)+int(rec[i][field]['word_size']/8*rec[i][field]['ngates'])
            k=k+1   
        msg_header['block_count']=k-1
        msg_header['radial_length']=msg_header['block_pointer_'+str(k-1)]+_structure_size(ar2v.GENERIC_DATA_BLOCK)+int(rec[i][fields[-1]]['word_size']/8*rec[i][fields[-1]]['ngates'])
        rec[i]['msg_header']=msg_header


        ## 2.2 Write message header 
        header=convert_to_dict(ar2v.MSG_HEADER)
        header['size']=int((msg_header['radial_length']+_structure_size(ar2v.MSG_HEADER))/2) # ?
        header['channels']=0 # ?
        header['type']=31
        header['seq_id']=0 # ?
        header['date']=0 # ?
        header['ms']=0 # ?
        header['segments']=0 # ?
        header['seg_num']=0 # ?
        rec[i]['header']=header

    print(volume_header['icao'],volume_header['icao'],radar.info['site']['Latitude'],radar.info['site']['Longitude'],
        radar.info['site']['antenna_height'],13,'CN',site_name_4s,'\n',)

    # 3.Write into NEXRAD ar2v binary file
    import pyart.io.nexrad_level2 as ar2v
    import struct 

    def dic2vars(dic):
        # print([dic[i] for i in dic.keys()])
        args=[]
        for i in dic.keys():
            if isinstance(dic[i],str):
                args.append(bytes(dic[i],encoding='ASCII'))
            else:
                args.append(dic[i])
        # print(args)
        return args
    def _pack_structure(structure,args):
        fmt='>'+''.join(i[1]for i in structure) # NEXRAD is big-endian
        lst=struct.pack(fmt,*args)
        return lst
    def lst_plus(lst,new_lst,pointer):
        if len(lst)<pointer:
            lst=lst+b'\00'*(pointer-len(lst))
            lst=lst+new_lst
        elif len(lst)>=pointer:
            lst=lst[:pointer]
            lst=lst+new_lst
        return lst
    with open(outpath,'wb') as f:

        # write Volume Header Record
        structure=ar2v.VOLUME_HEADER
        args=dic2vars(volume_header)
        lst=_pack_structure(structure,args)
        f.write(lst)

        # write Local Data Manager (LDM) Compressed Record
        ## Message Data
        for i in range(len(rec)):
            f.write(b'\x00'*12) 
            ### MSG Header
            structure=ar2v.MSG_HEADER
            args=dic2vars(rec[i]['header'])
            lst=_pack_structure(structure,args)
            f.write(lst)
            ### MSG 31 Header   
            structure=ar2v.MSG_31
            args=dic2vars(rec[i]['msg_header'])
            lst=_pack_structure(structure,args)
            # f.write(lst3[:-4]) # Block_Pointer_10 wrong to RVOL

            list1=['VOL','ELV','RAD','REF', 'VEL', 'SW', 'ZDR', 'PHI', 'RHO', 'CFP']
            set2=rec[i].keys()
            list3=list(list1 & set2)
            list3=sorted(list3, key=list1.index)

            for j in list3: # the sequence of fields in this MSG 31
                #### MSG 31 Data Block
                if j == "VOL":
                    structure=ar2v.VOLUME_DATA_BLOCK
                    args=dic2vars(rec[i][j])
                    new_lst=_pack_structure(structure,args)
                    lst=lst_plus(lst,new_lst,rec[i]['msg_header']['block_pointer_'+str(list3.index(j)+1)])
                elif j == "ELV":
                    structure=ar2v.ELEVATION_DATA_BLOCK
                    args=dic2vars(rec[i][j])
                    new_lst=_pack_structure(structure,args)
                    lst=lst_plus(lst,new_lst,rec[i]['msg_header']['block_pointer_'+str(list3.index(j)+1)])
                elif j == "RAD":
                    structure=ar2v.RADIAL_DATA_BLOCK
                    args=dic2vars(rec[i][j])
                    new_lst=_pack_structure(structure,args)
                    lst=lst_plus(lst,new_lst,rec[i]['msg_header']['block_pointer_'+str(list3.index(j)+1)])
                else:
                    structure=ar2v.GENERIC_DATA_BLOCK
                    tmp_rec=rec[i][j].copy()
                    tmp_rec.pop('data')
                    args=dic2vars(tmp_rec)
                    new_lst=_pack_structure(structure,args)
                    lst=lst_plus(lst,new_lst,rec[i]['msg_header']['block_pointer_'+str(list3.index(j)+1)])   
                    ##### MSG 31 Generic Data
                    tmp_data=rec[i][j]['data']
                    l,ngates=len(tmp_data),tmp_rec['ngates']
                    if tmp_rec['word_size'] == 16:
                        tmp_format='>'+'H'*ngates
                        lst=lst+struct.pack(tmp_format,*tmp_data)
                    elif tmp_rec['word_size'] == 8:
                        tmp_format='>'+'B'*ngates
                        lst=lst+struct.pack(tmp_format,*tmp_data)
                    else:
                        import warnings
                        warnings.warn(
                            'Unsupported bit size: %s. Returning array dtype "B"',
                            tmp_rec['word_size'])
            f.write(lst)
        f.closed
    return 

if __name__=='__main__':
    indir='D:/Seafile/Share/notebooks/20230523DisplayGUI/dual_pol_radar_plot/testdata/'
    fname='NJU.20220729.105032.AR2.bz2'
    inpath=indir+fname
    outdir='./'
    outpath=outdir+fname+'.ar2v'
    cinrad2nexrad(inpath,outpath,'NJUC')
# def cinrad2nexrad(inpath, outpath, site_name_4s=None):
#     from ..read.cinrad import cinrad, data_type
#     import pyart.io.nexrad_level2 as ar2v
#     from pyart.io.nexrad_level2 import _structure_size
#     import struct

#     radar = cinrad(inpath)

#     # =========================
#     # 工具函数（修复版）
#     # =========================

#     def convert_to_dict(structure):
#         # ✅ 用0初始化，而不是format字符串
#         return {key: 0 for key, _ in structure}

#     def dic2vars(dic, structure):
#         args = []
#         for key, fmt in structure:
#             val = dic[key]

#             # 字符串字段
#             if fmt.endswith('s'):
#                 if isinstance(val, str):
#                     val = val.encode('ASCII')
#                 elif not isinstance(val, (bytes, bytearray)):
#                     raise TypeError(f"{key} must be bytes")

#             # 整数类字段
#             elif fmt in ['B','H','I','L','Q','b','h','i','l','q']:
#                 val = int(val)

#             # float字段
#             elif fmt in ['f','d']:
#                 val = float(val)

#             args.append(val)

#         return args
#     def pack_structure(structure, dic):
#         fmt = '>' + ''.join(i[1] for i in structure)
#         args = dic2vars(dic, structure)
#         return struct.pack(fmt, *args)

#     # =========================
#     # Volume Header
#     # =========================

#     volume_header = convert_to_dict(ar2v.VOLUME_HEADER)

#     volume_header['tape'] = 'AR2V0006.'
#     volume_header['extension'] = '001'

#     scan_time = radar.info['task']['scan_start_time']
#     volume_header['date'] = int(scan_time // (24 * 3600))
#     volume_header['time'] = int((scan_time % (24 * 3600)) * 1000)

#     if site_name_4s is None:
#         site_name_4s = str(radar.info['site']['site_name'][0:4], encoding='ASCII')

#     volume_header['icao'] = site_name_4s

#     rec = []
#     list_index = [2, 3, 7, 10, 9]

#     # =========================
#     # 构建每个 radial
#     # =========================

#     for i in range(len(radar.rec['data'])):
#         rec.append({})

#         # -------- VOL / ELV / RAD --------
#         for field in ['VOL', 'ELV', 'RAD']:

#             if field == 'VOL':
#                 block = convert_to_dict(ar2v.VOLUME_DATA_BLOCK)
#                 block.update({
#                     'block_type': 'R',
#                     'data_name': field,
#                     'lrtup': 44,
#                     'version_major': 1,
#                     'version_minor': 0,
#                     'lat': radar.info['site']['Latitude'],
#                     'lon': radar.info['site']['Longitude'],
#                     'height': radar.info['site']['antenna_height'],
#                     'feedhorn_height': radar.info['site']['antenna_height'],
#                     'refl_calib': radar.info['task']['hori_cali'],
#                     'power_h': 269.2,
#                     'power_v': 266.7,
#                     'diff_refl_calib': radar.info['task']['ZDR_cali'],
#                     'init_phase': radar.info['task']['PHIDP_cali'] + 180,
#                     'vcp': 212,
#                     'spare': 0
#                 })

#             elif field == 'ELV':
#                 block = convert_to_dict(ar2v.ELEVATION_DATA_BLOCK)
#                 elev = radar.rec['header'][i]['elevation_number'] - 1
#                 block.update({
#                     'block_type': 'R',
#                     'data_name': field,
#                     'lrtup': 12,
#                     'atmos': int(radar.info['cut'][elev]['atmos_loss']),
#                     'refl_calib': radar.info['task']['hori_cali']
#                 })

#             elif field == 'RAD':
#                 block = convert_to_dict(ar2v.RADIAL_DATA_BLOCK)
#                 elev = radar.rec['header'][i]['elevation_number'] - 1
#                 block.update({
#                     'block_type': 'R',
#                     'data_name': field,
#                     'lrtup': 28,
#                     'unambig_range': int(299792458 / radar.info['cut'][elev]['PRF1'] / 2 / 1000),
#                     'noise_h': radar.info['task']['hori_noise'],
#                     'noise_v': radar.info['task']['vert_noise'],
#                     'nyquist_vel': int(radar.info['cut'][elev]['nyquist_spd']),
#                     'spare': 0
#                 })

#             rec[i][field] = block

#         # -------- moment blocks --------
#         for j in range(len(radar.rec['data'][i]['moment_header'])):

#             mh = radar.rec['data'][i]['moment_header'][j]
#             if mh['data_type'] not in list_index:
#                 continue

#             field = data_type[mh['data_type']]
#             block = convert_to_dict(ar2v.GENERIC_DATA_BLOCK)

#             elev = radar.rec['header'][i]['elevation_number'] - 1

#             block.update({
#                 'block_type': 'D',
#                 'data_name': field,
#                 'reserved': 0,
#                 'ngates': int(mh['block_length'] / mh['bin_length']),
#                 'first_gate': radar.info['cut'][elev]['start_range'],
#                 'gate_spacing': radar.info['cut'][elev]['dop_reso'] if field == 'VEL'
#                                 else radar.info['cut'][elev]['log_reso'],
#                 'thresh': 100,
#                 'snr_thres': 16,
#                 'flags': 0,
#                 'word_size': 8 * mh['bin_length'],
#                 'scale': mh['scale'],
#                 'offset': mh['offset']
#             })

#             block['data'] = radar.rec['data'][i]['moment_data'][j]
#             rec[i][field] = block

#         # =========================
#         # MSG_31 header（重写 pointer 逻辑）
#         # =========================

#         msg = convert_to_dict(ar2v.MSG_31)

#         sec = radar.rec['header'][i]['seconds']
#         usec = radar.rec['header'][i]['microseconds']

#         msg.update({
#             'id': volume_header['icao'],
#             'collect_date': sec // (24 * 3600),
#             'collect_ms': (sec % (24 * 3600)) * 1000 + usec // 1000,
#             'azimuth_number': radar.rec['header'][i]['radial_number'],
#             'azimuth_angle': radar.rec['header'][i]['azimuth'],
#             'compress_flag': 0,
#             'spare_0': 0,
#             'azimuth_resolution': radar.info['cut'][0]['angular_reso'],
#             'radial_spacing': radar.info['cut'][0]['angular_reso'],
#             'elevation_number': radar.rec['header'][i]['elevation_number'],
#             'cut_sector': 1,
#             'elevation_angle': radar.rec['header'][i]['elevation'],
#             'radial_blanking': 0,
#             'azimuth_mode': 0
#         })

#         # 顺序
#         order = ['VOL', 'ELV', 'RAD'] + [data_type[ii] for ii in list_index]
#         fields = [x for x in order if x in rec[i]]

#         pointer = _structure_size(ar2v.MSG_31)

#         for idx, field in enumerate(fields):
#             msg[f'block_pointer_{idx+1}'] = pointer

#             if field == 'VOL':
#                 size = _structure_size(ar2v.VOLUME_DATA_BLOCK)
#             elif field == 'ELV':
#                 size = _structure_size(ar2v.ELEVATION_DATA_BLOCK)
#             elif field == 'RAD':
#                 size = _structure_size(ar2v.RADIAL_DATA_BLOCK)
#             else:
#                 b = rec[i][field]
#                 size = _structure_size(ar2v.GENERIC_DATA_BLOCK) + int(b['ngates'] * b['word_size'] / 8)

#             pointer += size

#         msg['block_count'] = len(fields)
#         msg['radial_length'] = pointer

#         rec[i]['msg_header'] = msg

#         # =========================
#         # MSG_HEADER
#         # =========================

#         header = convert_to_dict(ar2v.MSG_HEADER)

#         header.update({
#             'size': int((pointer + _structure_size(ar2v.MSG_HEADER)) / 2),
#             'channels': 2,
#             'type': 31,
#             'seq_id': radar.rec['header'][i]['seq_number'],
#             'date': sec // (24 * 3600),
#             'ms': (sec % (24 * 3600)) * 1000 + usec // 1000,
#             'segments': 0,
#             'seg_num': 0
#         })

#         rec[i]['header'] = header

#     # =========================
#     # 写文件
#     # =========================

#     with open(outpath, 'wb') as f:

#         # volume header
#         f.write(pack_structure(ar2v.VOLUME_HEADER, volume_header))

#         for i in range(len(rec)):

#             f.write(b'\x00' * 12)

#             f.write(pack_structure(ar2v.MSG_HEADER, rec[i]['header']))
#             f.write(pack_structure(ar2v.MSG_31, rec[i]['msg_header']))

#             order = ['VOL', 'ELV', 'RAD'] + ['REF', 'VEL', 'SW', 'ZDR', 'PHI', 'RHO', 'CFP']

#             for field in order:
#                 if field not in rec[i]:
#                     continue

#                 block = rec[i][field]

#                 if field == 'VOL':
#                     size = _structure_size(ar2v.VOLUME_DATA_BLOCK)
#                 elif field == 'ELV':
#                     size = _structure_size(ar2v.ELEVATION_DATA_BLOCK)
#                 elif field == 'RAD':
#                     size = _structure_size(ar2v.RADIAL_DATA_BLOCK)
#                 else:
#                     tmp = block.copy()
#                     data = tmp.pop('data')

#                     f.write(pack_structure(ar2v.GENERIC_DATA_BLOCK, tmp))

#                     if tmp['word_size'] == 16:
#                         f.write(struct.pack('>' + 'H' * tmp['ngates'], *data))
#                     else:
#                         f.write(struct.pack('>' + 'B' * tmp['ngates'], *data))

#     return

# if __name__=='__main__':
#     indir='../../'
#     fname='Z_RADR_I_ZG170_20260329064800_O_DOR_YTD2_CAP_FMT.bin.bz2'
#     inpath=indir+fname
#     outdir=indir
#     outpath=outdir+fname+'.ar2v'
#     cinrad2nexrad(inpath,outpath,'ZG170')