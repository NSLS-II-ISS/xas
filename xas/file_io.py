from datetime import datetime
import os
from subprocess import call
import numpy as np
import pandas as pd
from . import xray



def load_dataset_from_files(db,uid):
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
        print(stream_file)

        if stream_source == 'pizzabox-di-file':
            data = load_trig_trace(stream_file)
        if stream_source == 'pizzabox-adc-file':
            print(stream_device)
            data = load_adc_trace(stream_file)
            stream_offset = f'{stream_device} offset'
            if stream_offset in db[uid]['start']:
                print("subtracting offset")
                data.iloc[:, 1] = data.iloc[:, 1] - record['start'][stream_offset]
            stream_gain =  f'{stream_device} gain'
            if stream_gain in db[uid]['start']:
                print("correcting for gain")
                data.iloc[:, 1] = data.iloc[:, 1]/(10**record['start'][stream_gain])


        if stream_source == 'pizzabox-enc-file':
            data = load_enc_trace(stream_file)
            print(stream_name)
            if stream_name =='hhm_theta':
                data.iloc[:,1] = xray.encoder2energy(data['encoder'], 360000,
                                                       -float(record['start']['angle_offset']))
                stream_name = 'energy'
                print(stream_name)

        arrays[stream_name] = data

    return arrays

def validate_file_exists(path_to_file,file_type = 'interp'):
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


def validate_path_exists(db, uid):
    path_to_file = db[uid].start['interp_filename']
    (path, filename) = os.path.split(path_to_file)
    if not os.path.isdir(path):
        os.mkdir(path)
        call(['chmod', '777', path])
    else:
        print('...........Path exists')

def create_file_header(db,uid):
    facility = db[uid]['start']['Facility']
    beamline = db[uid]['start']['beamline_id']
    pi = db[uid]['start']['PI']
    proposal = db[uid]['start']['PROPOSAL']
    saf = db[uid]['start']['SAF']

    comment = db[uid]['start']['comment']
    year = db[uid]['start']['year']
    cycle = db[uid]['start']['cycle']
    scan_id = db[uid]['start']['scan_id']
    real_uid = db[uid]['start']['uid']
    start_time = db[uid]['start']['time']
    stop_time = db[uid]['stop']['time']
    human_start_time = str(datetime.fromtimestamp(start_time).strftime('%m/%d/%Y  %H:%M:%S'))
    human_stop_time = str(datetime.fromtimestamp(stop_time).strftime('%m/%d/%Y  %H:%M:%S'))
    human_duration = str(datetime.fromtimestamp(stop_time - start_time).strftime('%M:%S'))

    if 'trajectory_name' in db[uid]['start']:
        trajectory_name = db[uid]['start']['trajectory_name']
    else:
        trajectory_name = ''

    if 'element' in db[uid]['start']:
        element = db[uid]['start']['element']
    else:
        element = ''

    if 'edge' in db[uid]['start']:
        edge = db[uid]['start']['edge']
    else:
        edge = ''

    if 'e0' in db[uid]['start']:
        e0 = db[uid]['start']['e0']
    else:
        e0 = ''

    comments ='# Facility: {}\n'\
                 '# Beamline: {}\n'\
                 '# Year: {}\n' \
                 '# Cycle: {}\n' \
                 '# SAF: {}\n' \
                 '# PI: {}\n'\
                 '# Proposal: {}\n'\
                 '# Scan ID: {}\n' \
                 '# UID: {}\n'\
                 '# Comment: {}\n'\
                 '# Trajectory name: {}\n'\
                 '# Element: {}\n'\
                 '# Edge: {}\n'\
                 '# E0: {}\n'\
                 '# Start time: {}\n'\
                 '# Stop time: {}\n' \
                 '# Total time: {}\n#\n# '.format(
                  facility,
                  beamline,
                  year,
                  cycle,
                  saf,
                  pi,
                  proposal,
                  scan_id,
                  real_uid,
                  comment,
                  trajectory_name,
                  element,
                  edge,
                  e0,
                  human_start_time,
                  human_stop_time,
                  human_duration)
    return  comments

def find_e0(db, uid):
    e0 = -1
    if 'e0' in db[uid]['start']:
        e0 = float(db[uid]['start']['e0'])
    return e0


def save_interpolated_df_as_file(path_to_file, df, comments):
    cols = df.columns.tolist()
    fmt = '%17.6f ' + '%12.6f ' + (' '.join(['%12.6e' for i in range(len(cols)-2)]))
    header = '  '.join(cols)

    df = df[cols]
    np.savetxt(path_to_file,
               df.values,
               fmt=fmt,
               delimiter=" ",
               header=header,
               comments= comments)

    #print("changing permissions to 774")
    call(['chmod', '774', path_to_file])


def save_binned_df_as_file(path_to_file, df, comments):
    (path, extension) = os.path.splitext(path_to_file)
    path_to_file = path + '.dat'
    path_to_file = validate_file_exists(path_to_file,file_type = 'bin')
    cols = df.columns.tolist()
    fmt = '%12.6f ' + (' '.join(['%12.6e' for i in range(len(cols) - 1)]))
    header = '  '.join(cols)
    df = df[cols]
    np.savetxt(path_to_file,
               df.values,
               fmt=fmt,
               delimiter=" ",
               header=header,
               comments=comments)

    #print("changing permissions to 774")
    call(['chmod', '774', path_to_file])




def load_interpolated_df_from_file(filename):
    ''' Load interp file and return'''

    if not os.path.exists(filename):
        raise IOError(f'The file {filename} does not exist.')
    header = read_header(filename)
    keys = header[header.rfind('#'):][1:-1].split()
    df = pd.read_table(filename, delim_whitespace=True, comment='#', names=keys, index_col=False).sort_values(keys[1])
    return df, header


def read_header(filename):
    header = ''
    line = '#'
    with open(filename) as myfile:
        while line[0] == '#':
            line = next(myfile)
            header += line
    return header[:-len(line)]


