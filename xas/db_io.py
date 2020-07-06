import pandas as pd
from . import xray

def load_apb_dataset_from_db(db, uid):
    hdr = db[uid]
    apb_dataset = list(hdr.data(stream_name='apb_stream', field='apb_stream'))[0]
    energy_dataset =  list(hdr.data(stream_name='pb9_enc1',field='pb9_enc1'))[0]
    angle_offset = -float(hdr['start']['angle_offset'])
    return apb_dataset, energy_dataset, angle_offset

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







