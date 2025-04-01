import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from xas.xas_logger import get_logger
import time as ttime



def interpolate(dataset, key_base = None, sort=True):
    logger = get_logger()

    interpolated_dataset = {}
    min_timestamp = max([dataset.get(key).iloc[0, 0] for key in dataset])
    max_timestamp = min([dataset.get(key).iloc[len(dataset.get(key)) - 1, 0] for key in
                         dataset if len(dataset.get(key).iloc[:, 0]) > 5])
    if key_base is None:
        all_keys = []
        time_step = []
        for key in dataset.keys():
            all_keys.append(key)
            # time_step.append(np.mean(np.diff(dataset[key].timestamp)))
            time_step.append(np.median(np.diff(dataset[key].timestamp)))
        key_base = all_keys[np.argmax(time_step)]
    timestamps = dataset[key_base].iloc[:,0]

    condition = timestamps < min_timestamp
    timestamps = timestamps[np.sum(condition):]

    condition = timestamps > max_timestamp
    timestamps = timestamps[: (len(timestamps) - np.sum(condition) - 1)]

    interpolated_dataset['timestamp'] = timestamps.values

    for key in dataset.keys():
        print(f'Interpolating stream {key}...')
        logger.info(f'({ttime.ctime()}) Interpolating stream {key}...')

        time = dataset.get(key).iloc[:, 0].values
        val = dataset.get(key).iloc[:, 1].values
        if len(dataset.get(key).iloc[:, 0]) > 5 * len(timestamps):
            time = [time[0]] + [np.mean(array) for array in np.array_split(time[1:-1], len(timestamps))] + [time[-1]]
            val = [val[0]] + [np.mean(array) for array in np.array_split(val[1:-1], len(timestamps))] + [val[-1]]
            # interpolated_dataset[key] = np.array([timestamps, np.interp(timestamps, time, val)]).transpose()

        # interpolated_dataset[key] = np.array([timestamps, np.interp(timestamps, time, val)]).transpose()
        interpolator_func = interp1d(time, np.array([v for v in val]), axis=0)
        val_interp = interpolator_func(timestamps)
        if len(val_interp.shape) == 1:
            interpolated_dataset[key] = val_interp
        else:
            interpolated_dataset[key] = [v for v in val_interp]
        print(f'Interpolation of stream {key} is complete')
        logger.info(f'({ttime.ctime()}) Interpolation of stream {key} is complete')

    intepolated_dataframe = pd.DataFrame(interpolated_dataset)
    if sort:
        return intepolated_dataframe.sort_values('energy')
    else:
        return intepolated_dataframe



def interpolate_with_interp(dataset, key_base = None, sort=True):
    # logger = get_logger()

    interpolated_dataset = {}
    min_timestamp = max([dataset.get(key).iloc[0, 0] for key in dataset])
    max_timestamp = min([dataset.get(key).iloc[len(dataset.get(key)) - 1, 0] for key in
                         dataset if len(dataset.get(key).iloc[:, 0]) > 5])
    if key_base is None:
        all_keys = []
        time_step = []
        for key in dataset.keys():
            all_keys.append(key)
            # time_step.append(np.mean(np.diff(dataset[key].timestamp)))
            time_step.append(np.median(np.diff(dataset[key].timestamp)))
        key_base = all_keys[np.argmax(time_step)]
    timestamps = dataset[key_base].iloc[:,0]

    condition = timestamps < min_timestamp
    timestamps = timestamps[np.sum(condition):]

    condition = timestamps > max_timestamp
    timestamps = timestamps[: (len(timestamps) - np.sum(condition) - 1)]

    interpolated_dataset['timestamp'] = timestamps.values

    for key in dataset.keys():
        print(f'Interpolating stream {key}...')
        # logger.info(f'({ttime.ctime()}) Interpolating stream {key}...')
        if key == 'pil100k2_image':
            print(f'---------------------------{key}----------------------------')
            time = dataset.get(key).iloc[:, 0].values.astype(np.float64)
            val = np.stack(dataset.get(key).iloc[:, 1].values)

            shape_length = val.shape[0]

            val_flat = val.reshape(shape_length, -1)
            interpolated_flat = np.empty((len(timestamps), val_flat.shape[1]), dtype=val.dtype)

            for i in range(val_flat.shape[1]):
                interpolated_flat[:, i] = np.interp(timestamps, time, val_flat[:, i])

            interpolated_reshaped = interpolated_flat.reshape(len(timestamps), 195, 487).astype('object')

            interpolated_dataset[key] = [v for v in interpolated_reshaped]
        else:
            time = dataset.get(key).iloc[:, 0].values
            val = dataset.get(key).iloc[:, 1].values
            if len(dataset.get(key).iloc[:, 0]) > 5 * len(timestamps):
                time = [time[0]] + [np.mean(array) for array in np.array_split(time[1:-1], len(timestamps))] + [time[-1]]
                val = [val[0]] + [np.mean(array) for array in np.array_split(val[1:-1], len(timestamps))] + [val[-1]]

            interpolator_func = interp1d(time, np.array([v for v in val]), axis=0)
            val_interp = interpolator_func(timestamps)
            interpolated_dataset[key] = val_interp


        # interpolator_func = interp1d(time, np.array([v for v in val]), axis=0)
        # val_interp = interpolator_func(timestamps)
        # if len(val_interp.shape) == 1:
        #     interpolated_dataset[key] = val_interp
        # else:
        #     interpolated_dataset[key] = [v for v in val_interp]
        # print(f'Interpolation of stream {key} is complete')
        # logger.info(f'({ttime.ctime()}) Interpolation of stream {key} is complete')

    intepolated_dataframe = pd.DataFrame(interpolated_dataset)
    if sort:
        return intepolated_dataframe.sort_values('energy')
    else:
        return intepolated_dataframe