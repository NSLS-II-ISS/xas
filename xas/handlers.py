from databroker.assets.handlers_base import HandlerBase
from databroker.assets.handlers import Xspress3HDF5Handler
import pandas as pd
import os
import numpy as np
import time as ttime

class APBBinFileHandler(HandlerBase):
    "Read electrometer *.bin files"
    def __init__(self, fpath):
        # It's a text config file, which we don't store in the resources yet, parsing for now
        fpath_txt = f'{os.path.splitext(fpath)[0]}.txt'

        with open(fpath_txt, 'r') as fp:
            content = fp.readlines()
            content = [x.strip() for x in content]

        _ = int(content[0].split(':')[1])
        Gains = [int(x) for x in content[1].split(':')[1].split(',')]
        Offsets = [int(x) for x in content[2].split(':')[1].split(',')]
        FAdiv = float(content[3].split(':')[1])
        FArate = float(content[4].split(':')[1])
        trigger_timestamp = float(content[5].split(':')[1].strip().replace(',', '.'))

        raw_data = np.fromfile(fpath, dtype=np.int32)

        columns = ['timestamp', 'i0', 'it', 'ir', 'iff', 'aux1', 'aux2', 'aux3', 'aux4']
        num_columns = len(columns) + 1  # TODO: figure out why 1
        raw_data = raw_data.reshape((raw_data.size // num_columns, num_columns))

        derived_data = np.zeros((raw_data.shape[0], raw_data.shape[1] - 1))
        derived_data[:, 0] = raw_data[:, -2] + raw_data[:, -1]  * 8.0051232 * 1e-9  # Unix timestamp with nanoseconds
        for i in range(num_columns - 2):
            derived_data[:, i+1] = raw_data[:, i] #((raw_data[:, i] ) - Offsets[i]) / Gains[i]

        self.df = pd.DataFrame(data=derived_data, columns=columns)
        self.raw_data = raw_data

    def __call__(self):
        #print(f'Returning {self.df}')
        return self.df.values



class PizzaBoxEncHandlerTxtPD(HandlerBase):
    "Read PizzaBox text files using info from filestore."
    def __init__(self, fpath):
        # self.df = pd.read_table(fpath, names=['ts_s', 'ts_ns', 'encoder', 'index', 'state'], sep=' ')
        self.data = np.genfromtxt(fpath)

    def __call__(self):
        # return self.df
        return self.data



class ISSXspress3HDF5Handler(Xspress3HDF5Handler):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._roi_data = None
        # self._num_channels = None

    def _get_dataset(self):
        super()._get_dataset()

        if self._roi_data is not None:
            return
        print(f'{ttime.ctime()} Determining number of channels...')
        shape = self.dataset.shape
        if len(shape) != 3:
            raise RuntimeError(f'The ndim of the dataset is not 3, but {len(shape)}')
        _num_channels = shape[1]

        self._roi_data = {}
        all_keys = self._file['/entry/instrument/detector/NDAttributes'].keys()
        for chan in range(1, _num_channels + 1):
            base = f'CHAN{chan}ROI'
            keys = [k for k in all_keys if k.startswith(base) and not k.endswith('LM')]
            for key in keys:
                roi_num = int(key.replace(base, ''))
                self._roi_data[(chan, roi_num)] = self._file['/entry/instrument/detector/NDAttributes'][key][()]

    def __call__(self, data_type:str='spectrum', channel:int=1, roi_num:int=1):
        # print(f'{ttime.ctime()} XS dataset retrieving starting...')
        self._get_dataset()

        if data_type=='spectrum':
            # actually outputs spectrum:
            # return self._dataset[:, channel - 1, :].squeeze()
            return np.array([[0]])

        elif data_type=='roi':
            return self._roi_data[(channel, roi_num)].squeeze()

        else:
            raise KeyError(f'data_type={data_type} not supported')



PIL100k_HDF_DATA_KEY = 'entry/instrument/NDAttributes'
class ISSPilatusHDF5Handler(Xspress3HDF5Handler): # Denis: I used Xspress3HDF5Handler as basis since it has all the basic functionality and I more or less understand how it works
    specs = {'PIL100k_HDF5'} | HandlerBase.specs
    HANDLER_NAME = 'PIL100k_HDF5'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, key=PIL100k_HDF_DATA_KEY, **kwargs)
        self._roi_data = None
        self.hdfrois = [f'ROI{i + 1}' for i in range(4)]
        self.chanrois = [f'pil100k_ROI{i + 1}' for i in range(4)]


    def _get_dataset(self):
        if self._dataset is not None:
            return

        _data_columns = [self._file[self._key + f'/_{chanroi}Total'][()] for chanroi in self.hdfrois]
        self._roi_data = np.vstack(_data_columns).T
        self._image_data = self._file['entry/data/data'][()]
        # self._roi_data = pd.DataFrame(data_columns, columns=self.chanrois)
        # self._dataset = data_columns

    def __call__(self, data_type:str='image', roi_num:int=0):
        # print(f'{ttime.ctime()} XS dataset retrieving starting...')
        self._get_dataset()

        if data_type=='image':
            # print(output.shape, output.squeeze().shape)
            return self._image_data

        elif data_type=='roi':
            return self._roi_data[:, roi_num - 1].squeeze()

        else:
            raise KeyError(f'data_type={data_type} not supported')


