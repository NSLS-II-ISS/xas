from databroker.assets.handlers_base import HandlerBase
from databroker.assets.handlers import Xspress3HDF5Handler
import pandas as pd
import os
import numpy as np
import time as ttime
from xas.file_io import _shift_root
from xas.xia import decode_xmap_buffers

class ISSHandlerBaseShiftRoot(HandlerBase):
    def shift_root(self, fpath):
        return _shift_root(fpath)


class APBBinFileHandler(ISSHandlerBaseShiftRoot):
    "Read electrometer *.bin files"
    def __init__(self, fpath):
        fpath = self.shift_root(fpath)
        # It's a text config file, which we don't store in the resources yet, parsing for now
        # fpath_txt = f'{os.path.splitext(fpath)[0]}.txt'
        #
        # with open(fpath_txt, 'r') as fp:
        #     content = fp.readlines()
        #     content = [x.strip() for x in content]
        #
        # _ = int(content[0].split(':')[1])
        # Gains = [int(x) for x in content[1].split(':')[1].split(',')]
        # Offsets = [int(x) for x in content[2].split(':')[1].split(',')]
        # FAdiv = float(content[3].split(':')[1])
        # FArate = float(content[4].split(':')[1])
        # trigger_timestamp = float(content[5].split(':')[1].strip().replace(',', '.'))

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



class PizzaBoxEncHandlerTxtPD(ISSHandlerBaseShiftRoot):
    "Read PizzaBox text files using info from filestore."
    def __init__(self, fpath):
        fpath = self.shift_root(fpath)
        # self.df = pd.read_table(fpath, names=['ts_s', 'ts_ns', 'encoder', 'index', 'state'], sep=' ')
        self.data = np.genfromtxt(fpath)

    def __call__(self):
        # return self.df
        return self.data



class APBTriggerFileHandler(ISSHandlerBaseShiftRoot):
    "Read APB trigger *.bin files"
    def __init__(self, fpath):
        fpath = self.shift_root(fpath)
        raw_data = np.fromfile(fpath, dtype=np.int32)
        raw_data = raw_data.reshape((raw_data.size // 3, 3))
        columns = ['timestamp', 'transition']
        derived_data = np.zeros((raw_data.shape[0], 2))
        derived_data[:, 0] = raw_data[:, 1] + raw_data[:, 2]  * 8.0051232 * 1e-9  # Unix timestamp with nanoseconds
        derived_data[:, 1] = raw_data[:, 0]

        self.df = pd.DataFrame(data=derived_data, columns=columns)
        self.raw_data = raw_data

    def __call__(self):
        return self.df


class ISSXspress3HDF5Handler(Xspress3HDF5Handler, ISSHandlerBaseShiftRoot):

    def __init__(self, filename, **kwargs):
        filename = self.shift_root(filename)
        super().__init__(filename, **kwargs)
        self._roi_data = None
        # self._num_channels = None

    def _get_dataset(self):
        super()._get_dataset()
        # if self._dataset is not None:
        #     return
        # self._dataset = np.array([[0]])

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


    def close(self):
        super().close()
        self._roi_data = None

    def __call__(self, data_type:str='spectrum', channel:int=1, roi_num:int=1):
        # print(f'{ttime.ctime()} XS dataset retrieving starting...')
        self._get_dataset()

        if data_type=='spectrum':
            # actually outputs spectrum:
            # return self._dataset[:, channel - 1, :].squeeze()
            # return self._dataset
            return np.array([[0]])

        elif data_type=='roi':
            return self._roi_data[(channel, roi_num)].squeeze()

        else:
            raise KeyError(f'data_type={data_type} not supported')



PIL100k_HDF_DATA_KEY = 'entry/instrument/NDAttributes'
class ISSPilatusHDF5Handler(Xspress3HDF5Handler, ISSHandlerBaseShiftRoot): # Denis: I used Xspress3HDF5Handler as basis since it has all the basic functionality and I more or less understand how it works
    specs = {'PIL100k_HDF5'} | HandlerBase.specs
    HANDLER_NAME = 'PIL100k_HDF5'

    def __init__(self, filename, **kwargs):
        filename = self.shift_root(filename)
        super().__init__(filename, key=PIL100k_HDF_DATA_KEY, **kwargs)
        self._roi_data = None
        self.hdfrois = [f'ROI{i + 1}' for i in range(4)]
        self.chanrois = [f'pil100k_ROI{i + 1}' for i in range(4)]

    def close(self):
        # print('Closing pilatus h5 file')
        super().close()
        self._image_data = None
        self._roi_data = None

    def _get_dataset(self):
        if self._dataset is not None:
            return

        _data_columns = [self._file[self._key + f'/_{chanroi}Total'][()] for chanroi in self.hdfrois]
        self._roi_data = np.vstack(_data_columns).T
        # self._image_data = None
        self._image_data = self._file['entry/data/data'][()]
        # self._roi_data = pd.DataFrame(data_columns, columns=self.chanrois)
        # self._dataset = data_columns

    def __call__(self, data_type:str='image', roi_num:int=0):
        # print(f'{ttime.ctime()} XS dataset retrieving starting...')
        self._get_dataset()

        # print('In Pilatus handler')


        if data_type=='image':
            # print(output.shape, output.squeeze().shape)
            return self._image_data

        elif data_type=='roi':
            return self._roi_data[:, roi_num - 1].squeeze()

        else:
            raise KeyError(f'data_type={data_type} not supported')


class ISSPilatusHDF5HandlerLegacy(Xspress3HDF5Handler, ISSHandlerBaseShiftRoot): # Denis: I used Xspress3HDF5Handler as basis since it has all the basic functionality and I more or less understand how it works
    specs = {'PIL100k_HDF5'} | HandlerBase.specs
    HANDLER_NAME = 'PIL100k_HDF5'

    def __init__(self, filename, **kwargs):
        filename = self.shift_root(filename)
        super().__init__(filename, key=PIL100k_HDF_DATA_KEY, **kwargs)
        self._roi_data = None
        self.hdfrois = [f'ROI{i + 1}' for i in range(4)]
        self.chanrois = [f'pil100k_ROI{i + 1}' for i in range(4)]

    def close(self):
        # print('Closing pilatus h5 file')
        super().close()
        self._image_data = None
        self._roi_data = None

    def _get_dataset(self):
        if self._dataset is not None:
            return

        _data_columns = [self._file[self._key + f'/_{chanroi}Total'][()] for chanroi in self.hdfrois]
        data_columns = np.vstack(_data_columns).T
        self._roi_data = pd.DataFrame(data_columns, columns=self.chanrois)
        self._image_data = self._file['entry/data/data'][()]
        self._dataset = 'bla'
        # self._dataset = data_columns

    def __call__(self, *args, frame=None,  **kwargs):
        self._get_dataset()
        return_dict = {chanroi: self._roi_data[chanroi][frame] for chanroi in self.chanrois}
        return_dict['image'] = self._image_data[frame, :, :].squeeze()
        return return_dict


class ISSXIAHDF5HandlerLegacy(Xspress3HDF5Handler,
                              ISSHandlerBaseShiftRoot):  # Denis: I used Xspress3HDF5Handler as basis since it has all the basic functionality and I more or less understand how it works
    #    specs = {'PIL100k_HDF5'} | HandlerBase.specs
    HANDLER_NAME = 'XIA_XMAP_HDF5'

    def __init__(self, filename, **kwargs):
        filename = self.shift_root(filename)
        super().__init__(filename, **kwargs)
        self._roi_data = None
        # self._num_channels = None

    def _get_dataset(self):
        super()._get_dataset()
        # if self._dataset is not None:
        #     return
        # self._dataset = np.array([[0]])

        data, mapmode = decode_xmap_buffers(self._file['entry']['data']['data'])
        if mapmode != 2:
            raise RuntimeError('Data must be saved in ROI mapping mode!')
        self._dataset = data.counts
        #        if self._roi_data is not None:
        #            return
        #        print(f'{ttime.ctime()} Determining number of channels...')
        shape = self.dataset.shape
        if len(shape) != 3:
            raise RuntimeError(f'The ndim of the dataset is not 3, but {len(shape)}')
        _num_channels = shape[1]
        _maxroi = shape[2]
        self._roi_data = {}
        #        all_keys = self._file['/entry/instrument/detector/NDAttributes'].keys()
        for chan in range(_num_channels):

            #            base = f'CHAN{chan}ROI'
            #            keys = [k for k in all_keys if k.startswith(base) and not k.endswith('LM')]
            for roi_num in range(_maxroi):
                #                roi_num = int(key.replace(base, ''))
                self._roi_data[(chan + 1, roi_num)] = self._dataset[:, chan, roi_num]

    #                self._roi_data[(chan, roi_num)] = self._file['/entry/instrument/detector/NDAttributes'][key][()]

    def close(self):
        super().close()
        self._roi_data = None

    def __call__(self, data_type: str = 'spectrum', channel: int = 1, roi_num: int = 1):
        # print(f'{ttime.ctime()} XS dataset retrieving starting...')
        self._get_dataset()

        if data_type == 'spectrum':
            # actually outputs spectrum:
            # return self._dataset[:, channel - 1, :].squeeze()
            # return self._dataset
            return np.array([[0]])

        elif data_type == 'roi':
            return self._roi_data[(channel, roi_num)].squeeze()

        else:
            raise KeyError(f'data_type={data_type} not supported')

def register_all_handlers(db):
    db.reg.register_handler('PIZZABOX_ENC_FILE_TXT_PD',
                            PizzaBoxEncHandlerTxtPD, overwrite=True)
    db.reg.register_handler('APB',
                            APBBinFileHandler, overwrite=True)
    db.reg.register_handler('APB_TRIGGER',
                            APBTriggerFileHandler, overwrite=True)
    db.reg.register_handler('PIL100k_HDF5',
                            ISSPilatusHDF5Handler, overwrite=True)
    db.reg.register_handler(ISSXspress3HDF5Handler.HANDLER_NAME,
                            ISSXspress3HDF5Handler, overwrite=True)
    db.reg.register_handler(ISSXIAHDF5HandlerLegacy.HANDLER_NAME,
                            ISSXIAHDF5HandlerLegacy, overwrite=True)


def register_all_handlers_legacy(db):
    db.reg.register_handler('PIZZABOX_ENC_FILE_TXT_PD',
                            PizzaBoxEncHandlerTxtPD, overwrite=True)
    db.reg.register_handler('APB',
                            APBBinFileHandler, overwrite=True)
    db.reg.register_handler('APB_TRIGGER',
                            APBTriggerFileHandler, overwrite=True)
    db.reg.register_handler('PIL100k_HDF5',
                            ISSPilatusHDF5HandlerLegacy, overwrite=True)
    db.reg.register_handler(ISSXspress3HDF5Handler.HANDLER_NAME,
                            ISSXspress3HDF5Handler, overwrite=True)
    db.reg.register_handler(ISSXIAHDF5HandlerLegacy.HANDLER_NAME,
                            ISSXIAHDF5HandlerLegacy, overwrite=True)


def get_iss_db():
    from databroker import Broker
    db = Broker.named('iss-local')
    register_all_handlers(db)
    return db