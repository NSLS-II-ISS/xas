import pandas as pd
import numpy as np
from . import xray

def load_apb_dataset_from_db(db, uid):
    hdr = db[uid]
    apb_dataset = list(hdr.data(stream_name='apb_stream', field='apb_stream'))[0]
    energy_dataset =  list(hdr.data(stream_name='pb9_enc1',field='pb9_enc1'))[0]
    angle_offset = -float(hdr['start']['angle_offset'])

    # ch_offset_keys = [key for key in hdr.start.keys() if key.startswith('ch') and key.endswith('_offset')]
    # ch_offsets = np.array([hdr.start[key] for key in ch_offset_keys])

    ch_offsets = get_ch_properties(hdr.start, 'ch', '_offset')
    ch_gains = get_ch_properties(hdr.start, 'ch', '_amp_gain')

    apb_dataset.iloc[:, 1:] -= ch_offsets
    apb_dataset.iloc[:, 1:] /= 1e6
    apb_dataset.iloc[:, 1:] /= (10**ch_gains)

    return apb_dataset, energy_dataset, angle_offset



def get_ch_properties(hdr_start, start, end):
    ch_keys = [key for key in hdr_start.keys() if key.startswith(start) and key.endswith(end)]
    return np.array([hdr_start[key] for key in ch_keys])



def translate_apb_dataset(apb_dataset, energy_dataset, angle_offset):
    data_dict= {}
    for column in apb_dataset.columns:
        if column != 'timestamp':
            adc = pd.DataFrame()
            adc['timestamp'] = apb_dataset['timestamp']
            adc['adc'] = apb_dataset[column]

            data_dict[column]=adc

    energy = pd.DataFrame()
    energy['timestamp'] = energy_dataset['ts_s'] + 1e-9 * energy_dataset['ts_ns']
    enc  = energy_dataset['encoder'].apply(lambda x: int(x) if int(x) <= 0 else -(int(x) ^ 0xffffff - 1))


    energy['encoder'] = xray.encoder2energy(enc, 360000, angle_offset)

    data_dict['energy'] = energy
    return data_dict







