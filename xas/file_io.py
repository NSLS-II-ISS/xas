from datetime import datetime
import os
from subprocess import call
import numpy as np
import pandas as pd
from . import xray



def load_dataset_from_files(db, uid):
    def load_adc_trace(filename=''):
        df=pd.DataFrame()
        keys = ['times', 'timens', 'counter', 'adc']
        if os.path.isfile(filename):
            df_raw = pd.read_table(filename, delim_whitespace=True, comment='#', names=keys, index_col=False)
            df['timestamp'] = df_raw['times'] + 1e-9 * df_raw['timens']
            df['adc'] = df_raw['adc'].apply(lambda x: (int(x, 16) >> 8) - 0x40000 if (int(x, 16) >> 8) > 0x1FFFF else int(x, 16) >> 8) * 7.62939453125e-05
            return df
        else:
            return -1

    def load_enc_trace(filename=''):
        df = pd.DataFrame()
        keys = ['times', 'timens', 'encoder', 'counter', 'di']
        if os.path.isfile(filename):
            df_raw = pd.read_table(filename, delim_whitespace=True, comment='#', names=keys, index_col=False)
            df['timestamp'] = df_raw['times'] + 1e-9 * df_raw['timens']
            df['encoder'] = df_raw['encoder'].apply(lambda x: int(x) if int(x) <= 0 else -(int(x) ^ 0xffffff - 1))
            return df
        else:
            return -1

    def load_trig_trace(filename=''):
        keys = ['times', 'timens', 'encoder', 'counter', 'di']
        if os.path.isfile(filename):
            df = pd.read_table(filename, delim_whitespace=True, comment='#', names=keys,
                               index_col=False)
            df['timestamp'] = df['times'] + 1e-9 * df['timens']
            df = df.iloc[::2]
            return df.iloc[:, [5, 3]]
        else:
            return -1

    arrays = {}
    record = db[uid]
    for stream in record['descriptors']:
        data = pd.DataFrame()
        stream_device = stream['name']
        stream_name = stream['data_keys'][stream['name']]['devname']
        stream_source = stream['data_keys'][stream['name']]['source']
        stream_file = stream['data_keys'][stream['name']]['filename']
        #print(stream_file)

        if stream_source == 'pizzabox-di-file':
            data = load_trig_trace(stream_file)
        if stream_source == 'pizzabox-adc-file':
            #print(stream_device)
            data = load_adc_trace(stream_file)
            stream_offset = f'{stream_device} offset'
            if stream_offset in db[uid]['start']:
                #print("subtracting offset")
                data.iloc[:, 1] = data.iloc[:, 1] - record['start'][stream_offset]
            stream_gain =  f'{stream_device} gain'
            if stream_gain in db[uid]['start']:
                #print("correcting for gain")
                data.iloc[:, 1] = data.iloc[:, 1]/(10**record['start'][stream_gain])


        if stream_source == 'pizzabox-enc-file':
            data = load_enc_trace(stream_file)
            #print(stream_name)
            if stream_name =='hhm_theta':
                data.iloc[:,1] = xray.encoder2energy(data['encoder'], 360000,
                                                       -float(record['start']['angle_offset']))
                stream_name = 'energy'
                #print(stream_name)

        arrays[stream_name] = data

    return arrays

def validate_file_exists(path_to_file,file_type = 'interp'):
    """The function checks if the file exists or not. If exists, it adds an index to the file name.
    """
    if file_type == 'interp':
        prefix = 'r'
    elif file_type == 'bin':
        prefix = 'b'
    if os.path.isfile(path_to_file):
        (path, extension) = os.path.splitext(path_to_file)
        iterator = 2

        while True:

            new_filename = '{}-{}{:04d}{}'.format(path, prefix, iterator,extension)
            if not os.path.isfile(new_filename):
                return new_filename
            iterator += 1
    else:
        return path_to_file

def make_user_dir(path):
    folder_exists = os.path.isdir(path)
    if not folder_exists:
        os.mkdir(path)
        call(['chmod', '777', path])
    return folder_exists

def validate_path_exists(db, uid):
    path_to_file = db[uid].start['interp_filename']
    (path, filename) = os.path.split(path_to_file)
    if not os.path.isdir(path):
        os.mkdir(path)
        call(['chmod', '777', path])


def _get_value_from_hdr_start(hdr, key):
    if key in hdr['start']:
        return hdr['start'][key]
    else:
        return ''

def create_file_header(db,uid):
    hdr = db[uid]
    facility = hdr['start']['Facility']
    beamline = hdr['start']['beamline_id']
    pi = hdr['start']['PI']
    proposal = hdr['start']['PROPOSAL']
    saf = hdr['start']['SAF']
    name = hdr['start']['name']

    comment = hdr['start']['comment']
    year = hdr['start']['year']
    cycle = hdr['start']['cycle']
    scan_id = hdr['start']['scan_id']
    real_uid = hdr['start']['uid']
    start_time = hdr['start']['time']
    stop_time = hdr['stop']['time']
    human_start_time = str(datetime.fromtimestamp(start_time).strftime('%m/%d/%Y  %H:%M:%S'))
    human_stop_time = str(datetime.fromtimestamp(stop_time).strftime('%m/%d/%Y  %H:%M:%S'))
    human_duration = str(datetime.fromtimestamp(stop_time - start_time).strftime('%M:%S'))

    trajectory_name = _get_value_from_hdr_start(hdr, 'trajectory_name')
    element = _get_value_from_hdr_start(hdr, 'element')
    edge = _get_value_from_hdr_start(hdr, 'edge')
    e0 = _get_value_from_hdr_start(hdr, 'e0')

    i0_par = _get_value_from_hdr_start(hdr, 'i0_par')
    it_par = _get_value_from_hdr_start(hdr, 'it_par')
    ir_par = _get_value_from_hdr_start(hdr, 'ir_par')
    iff_par = _get_value_from_hdr_start(hdr, 'iff_par')

    i0_volt = _get_value_from_hdr_start(hdr, 'i0_volt')
    it_volt = _get_value_from_hdr_start(hdr, 'it_volt')
    ir_volt = _get_value_from_hdr_start(hdr, 'ir_volt')

    i0_gain = _get_value_from_hdr_start(hdr, 'i0_gain')
    it_gain = _get_value_from_hdr_start(hdr, 'it_gain')
    ir_gain = _get_value_from_hdr_start(hdr, 'ir_gain')
    iff_gain = _get_value_from_hdr_start(hdr, 'iff_gain')

    aux_detector = _get_value_from_hdr_start(hdr, 'aux_detector')

    nslsii_status = _get_value_from_hdr_start(hdr, 'nslsii_status')
    nslsii_energy = _get_value_from_hdr_start(hdr, 'nslsii_energy')
    nslsii_current = _get_value_from_hdr_start(hdr, 'nslsii_current')

    harmonic_rejection = _get_value_from_hdr_start(hdr, 'harmonic_rejection')

    mono_offset = _get_value_from_hdr_start(hdr, 'mono_offset')
    mono_encoder_resolution = _get_value_from_hdr_start(hdr, 'mono_encoder_resolution')
    mono_scan_type = _get_value_from_hdr_start(hdr, 'mono_scan_type')
    mono_scan_mode = _get_value_from_hdr_start(hdr, 'mono_scan_mode')
    mono_direction = _get_value_from_hdr_start(hdr, 'mono_direction')

    sample_stage = _get_value_from_hdr_start(hdr, 'sample_stage')
    sample_x_position = _get_value_from_hdr_start(hdr, 'sample_x_position')
    sample_y_position = _get_value_from_hdr_start(hdr, 'sample_y_position')
    plot_hint = _get_value_from_hdr_start(hdr, 'plot_hint')

    # comments =f'# Facility: {facility}\n'\
    #           f'# Beamline: {beamline}\n'\
    #           f'# Year: {year}\n' \
    #           f'# Cycle: {cycle}\n' \
    #           f'# SAF: {saf}\n' \
    #           f'# PI: {pi}\n'\
    #           f'# Proposal: {proposal}\n'\
    #           f'# Scan ID: {scan_id}\n' \
    #           f'# UID: {real_uid}\n'\
    #           f'# Comment: {comment}\n'\
    #           f'# Trajectory name: {trajectory_name}\n'\
    #           f'# Element: {element}\n'\
    #           f'# Edge: {edge}\n'\
    #           f'# E0: {e0}\n'\
    #           f'# Start time: {human_start_time}\n'\
    #           f'# Stop time: {human_stop_time}\n' \
    #           f'# Total time: {human_duration}\n#\n# '

    comments = \
        f'# Beamline.name: {beamline} - (08ID) Inner Shell Spectroscopy\n' \
        f'# Beamline.x-ray_source: {facility} damping wiggler\n' \
        f'# Beamline.collimation: ISS dual mirror system\n' \
        f'# Beamline.focusing: ISS toroid mirror\n' \
        f'# Beamline.harmonic_rejection: {harmonic_rejection}\n' \
        f'# Detector.I0: {i0_par}\n' \
        f'# Detector.I1: {it_par}\n'\
        f'# Detector.I2: {ir_par}\n' \
        f'# Detector.IF: {iff_par}\n' \
        f'# Detector.I0_volt: {i0_volt}\n' \
        f'# Detector.I1_volt: {it_volt}\n' \
        f'# Detector.I2_volt: {ir_volt}\n' \
        f'# Detector.I0_amp_gain: {i0_gain}\n' \
        f'# Detector.I1_amp_gain: {it_gain}\n' \
        f'# Detector.I2_amp_gain: {ir_gain}\n' \
        f'# Detector.IF_amp_gain: {iff_gain}\n' \
        f'# Detector.aux: {aux_detector}\n' \
        f'# Element.symbol: {element}\n' \
        f'# Element.edge: {edge}\n' \
        f'# Facility.name: {facility}\n' \
        f'# Facility.mode: {nslsii_status}\n' \
        f'# Facility.energy: {nslsii_energy}\n'\
        f'# Facility.current: {nslsii_current}\n' \
        f'# Facility.GUP: {proposal}\n'\
        f'# Facility.SAF: {saf}\n' \
        f'# Facility.cycle: {year}-{cycle}\n' \
        f'# Mono.name: Si(111)\n' \
        f'# Mono.d_spacing: 3.1354951\n' \
        f'# Mono.encoder_resolution: {mono_encoder_resolution}\n' \
        f'# Mono.angle_offset: {mono_offset}\n' \
        f'# Mono.scan_mode: {mono_scan_mode}\n' \
        f'# Mono.scan_type: {mono_scan_type}\n' \
        f'# Mono.trajectory_name: {trajectory_name}\n' \
        f'# Mono.direction: {mono_direction}\n' \
        f'# Sample.name: {name}\n' \
        f'# Sample.comment: {comment}\n' \
        f'# Sample.stage: {sample_stage}\n' \
        f'# ISS.sample_x_position: {sample_x_position}\n' \
        f'# ISS.sample_y_position: {sample_y_position}\n' \
        f'# Scan.experimenters: {pi}\n' \
        f'# Scan.edge_energy: {e0}\n' \
        f'# Scan.start_time: {human_start_time}\n' \
        f'# Scan.end_time: {human_stop_time}\n' \
        f'# Scan.duration: {human_duration}\n' \
        f'# Scan.transient_id: {scan_id}\n' \
        f'# Scan.uid: {real_uid}\n' \
        f'# Scan.plot_hint: {plot_hint}\n' \
        f'# Column.1: energy eV\n' \
        f'# Column.2: i0\n' \
        f'# Column.3: it\n' \
        f'# Column.4: ir\n' \
        f'# Column.5: if\n'


    return  comments

def find_e0(db, uid):
    e0 = -1
    if 'e0' in db[uid]['start']:
        e0 = float(db[uid]['start']['e0'])
    return e0


def save_interpolated_df_as_file(path_to_file, df, comments):
    cols = df.columns.tolist()
    fmt = '%17.6f ' + '%12.6e ' + (' '.join(['%12.6e' for i in range(len(cols)-2)]))
    header = '# ' + ' '.join(cols)

    df = df[cols]
    np.savetxt(path_to_file,
               df.values,
               fmt=fmt,
               delimiter=" ",
               header=header,
               comments=comments)

    #print("changing permissions to 774")
    call(['chmod', '774', path_to_file])


# def dump_uid_bundle(base_name, uids, db):
#     path_to_interp_file = db[uids[0]].start()['interp_filename']
#     (path, extension) = os.path.splitext(path_to_interp_file)



def save_binned_df_as_file(path_to_file, df, comments):
    (path, extension) = os.path.splitext(path_to_file)
    path_to_file = path + '.dat'
    path_to_file = validate_file_exists(path_to_file,file_type = 'bin')
    #cols = df.columns.tolist()[::-1]
    cols = df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    fmt = '%12.6f ' + (' '.join(['%12.6e' for i in range(len(cols) - 1)]))
    header = '  '.join(cols)

    df = df[cols]
    np.savetxt(path_to_file,
               df.values,
               fmt=fmt,
               delimiter=" ",
               header=f'# {header}',
               comments=comments)

    #print("changing permissions to 774")
    call(['chmod', '774', path_to_file])




def load_interpolated_df_from_file(filename):
    ''' Load interp file and return'''

    if not os.path.exists(filename):
        raise IOError(f'The file {filename} does not exist.')
    header = read_header(filename)
    keys = header[header.rfind('#'):][1:-1].split()
    df = pd.read_csv(filename, delim_whitespace=True, comment='#', names=keys, index_col=False)
    if 'energy' in [i.lower() for i in df.columns]: # i am sorry /// Denis
        df = df.sort_values('energy'.lower())
    return df, header

def load_binned_df_from_file(filename):
    ''' Load interp file and return'''

    if not os.path.exists(filename):
        raise IOError(f'The file {filename} does not exist.')
    header = read_header(filename)
    keys = header[header.rfind('#'):][1:-1].split()
    df = pd.read_csv(filename, delim_whitespace=True, comment='#', names=keys, index_col=False).sort_values(keys[0])

    energy_key = None
    for col in df.columns:
        if 'energy' in col:
            energy_key = col
            break

    if energy_key:
        df = df.rename(columns={energy_key: 'energy'})

    if 'energy' not in df.columns:
        df['energy'] = np.arange(len(df.index))
    return df, header

def read_header(filename):
    header = ''
    line = '#'
    with open(filename) as myfile:
        while line[0] == '#':
            line = next(myfile)
            header += line
    return header[:-len(line)]

def convert_header_to_dict(header):
    lines = header.split('\n')
    lines = [line for line in lines if len(line) > 0]
    lines = [line.replace('# ', '') for line in lines]




stepscan_channel_dict = {
    'hhm_energy': 'energy',
    'motor_emission_energy' : 'energy',
    'apb_ave_ch1_mean': 'i0',
    'apb_ave_ch2_mean': 'it',
    'apb_ave_ch3_mean': 'ir',
    'apb_ave_ch4_mean': 'iff',
    'apb_ave_ch5_mean': 'aux1',
    'apb_ave_ch6_mean': 'aux2',
    'apb_ave_ch7_mean': 'aux3',
    'apb_ave_ch8_mean': 'aux4',
    'pil100k_stats1_total': 'pil100_ROI1',
    'pil100k_stats2_total': 'pil100_ROI2',
    'pil100k_stats3_total': 'pil100_ROI3',
    'pil100k_stats4_total': 'pil100_ROI4',
    'xs_channel1_rois_roi01_value' : 'xs_channel1_rois_roi01_value',
    'xs_channel2_rois_roi01_value' : 'xs_channel2_rois_roi01_value',
    'xs_channel3_rois_roi01_value' : 'xs_channel3_rois_roi01_value',
    'xs_channel4_rois_roi01_value' : 'xs_channel4_rois_roi01_value',
    'xs_ROI1': 'xs_ROI1',
    'xs_ROI2': 'xs_ROI2',
    'xs_ROI3': 'xs_ROI3',
    'xs_ROI4': 'xs_ROI4'}

xs_channel_list = [
    'xs_channel1_rois_roi01_value',
    'xs_channel1_rois_roi02_value',
    'xs_channel1_rois_roi03_value',
    'xs_channel1_rois_roi04_value',
    'xs_channel2_rois_roi01_value',
    'xs_channel2_rois_roi02_value',
    'xs_channel2_rois_roi03_value',
    'xs_channel2_rois_roi04_value',
    'xs_channel3_rois_roi01_value',
    'xs_channel3_rois_roi02_value',
    'xs_channel3_rois_roi03_value',
    'xs_channel3_rois_roi04_value',
    'xs_channel4_rois_roi01_value',
    'xs_channel4_rois_roi02_value',
    'xs_channel4_rois_roi03_value',
    'xs_channel4_rois_roi04_value']

xs_channel_comb_dict = {'xs_ROI1' : [#'xs_channel1_rois_roi01_value',
                                     'xs_channel2_rois_roi01_value',
                                     'xs_channel3_rois_roi01_value',]
                                     #'xs_channel4_rois_roi01_value'],
                        # 'xs_ROI2': [#'xs_channel1_rois_roi02_value',
                        #             'xs_channel2_rois_roi02_value',
                        #             'xs_channel3_rois_roi02_value',]
                        #             #'xs_channel4_rois_roi02_value'],
                        # 'xs_ROI3': [#'xs_channel1_rois_roi03_value',
                        #             'xs_channel2_rois_roi03_value',
                        #             'xs_channel3_rois_roi03_value',]
                        #             #'xs_channel4_rois_roi03_value'],
                        # 'xs_ROI4': [#'xs_channel1_rois_roi04_value',
                        #             'xs_channel2_rois_roi04_value',
                        #             'xs_channel3_rois_roi04_value',]
                                    #'xs_channel4_rois_roi04_value'],
                        }



def stepscan_remove_offsets(hdr, fill=True):
    try:
        df = hdr.table(fill=True)
    except:
        df = hdr.table()

    for channel_name, _ in stepscan_channel_dict.items():
        if channel_name.startswith('apb'):
            offset = hdr.descriptors[0]["configuration"]['apb_ave']['data'][channel_name.replace("_mean", "_offset")]
            df[channel_name] = df[channel_name] - offset
    return df


def stepscan_normalize_xs(df):
    for channel_name in xs_channel_list:
        if channel_name in df.columns:
            df[channel_name] = df[channel_name] / df['xs_settings_acquire_time']
    return df


def combine_xspress3_channels(df):

    # if xs_channel_list[0] in df.columns:
    #     df_out = df.copy()
    #     for k in xs_channel_comb_dict.keys():
    #         df_out[k] = 0
    #         for channel in xs_channel_comb_dict[k]:
    #             df_out[k] += df_out[channel]
    #     return df_out
    #
    # else:
    #     return df
    return df


def filter_df_by_valid_keys(df):
    d = {}
    for key, new_key in stepscan_channel_dict.items():
        if key in df.columns:
            d[new_key] = df[key].values
    return pd.DataFrame(d)


def write_df_to_file(path_to_file, df, comments):
    cols = df.columns.tolist()

    fmt = '%12.6f ' + (' '.join(['%12.6e' for i in range(len(cols) - 1)]))
    header = '  '.join(cols)

    # df = df[cols]
    np.savetxt(path_to_file,
               df.values,
               fmt=fmt,
               delimiter=" ",
               header=f'# {header}',
               comments=comments)
    call(['chmod', '774', path_to_file])


def save_stepscan_as_file(path_to_file, df, comments):
    # assuming stepscan_channel_dict keys and values are ordered as above
    # select all columns from df with names in stepscan_channel_dict.keys()
    # and rename

    # valid_keys = set(stepscan_channel_dict.keys()) & set(df.columns)
    # valid_keys = [key for key in stepscan_channel_dict.keys() if key in df.columns]
    # data = df[valid_keys]
    df_upd = filter_df_by_valid_keys(df)
    # data.columns = [stepscan_channel_dict[k] for k in data.columns]
    write_df_to_file(path_to_file, df_upd, comments)
    return df_upd


# # ########@
# fpath = '/nsls2/xf08id/users/2021/1/306253/'
# for i, file in enumerate(flist_actual):
#     this_df, _ = load_interpolated_df_from_file(fpath+file)
#     if i == 0:
#         master_energy = this_df['energy'].values
#         mu_all = np.zeros(master_energy.shape)
#     this_energy = this_df['energy'].values
#     this_i0 = this_df['i0']
#     this_sdd = this_df['CHAN1ROI1'] + this_df['CHAN2ROI1'] + this_df['CHAN3ROI1'] + this_df['CHAN4ROI1']
#     this_mu = - (this_sdd / this_i0).values
#     mu_all += np.interp(master_energy, this_energy, this_mu)
#
# data_dict = {'i0' : np.ones(mu_all.size), 'mu_all' : mu_all, 'energy' : master_energy}
# data_df = pd.DataFrame(data_dict)
# fname = '3_AfterReduction Ir on TiO2 data 0225.dat'
# save_binned_df_as_file(fpath+fname, data_df, '')
