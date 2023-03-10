import numpy as np
import pandas as pd
from scipy.interpolate import interp1d


def interpolate(dataset, key_base = None, sort=True):
    interpolated_dataset = {}
    min_timestamp = max([dataset.get(key).iloc[0, 0] for key in dataset])
    max_timestamp = min([dataset.get(key).iloc[len(dataset.get(key)) - 1, 0] for key in
                         dataset if len(dataset.get(key).iloc[:, 0]) > 5])
    if key_base is None:
        all_keys = []
        time_step = []
        for key in dataset.keys():
            all_keys.append(key)
            time_step.append(np.mean(np.diff(dataset[key].timestamp)))
        key_base = all_keys[np.argmax(time_step)]
    timestamps = dataset[key_base].iloc[:,0]

    condition = timestamps < min_timestamp
    timestamps = timestamps[np.sum(condition):]

    condition = timestamps > max_timestamp
    timestamps = timestamps[: (len(timestamps) - np.sum(condition) - 1)]

    interpolated_dataset['timestamp'] = timestamps.values

    for key in dataset.keys():
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

    intepolated_dataframe = pd.DataFrame(interpolated_dataset)
    # intepolated_dataframe = pd.DataFrame(np.vstack((timestamps, np.array([interpolated_dataset[key][:, 1] for
    #                                                                         key in interpolated_dataset.keys()]))).transpose())
    # keys = ['timestamp']
    # keys.extend(interpolated_dataset.keys())
    # intepolated_dataframe.columns = keys

    # intepolated_dataframe['mu_t'] = np.log( intepolated_dataframe['i0'] / intepolated_dataframe['it'] )
    # intepolated_dataframe['mu_f'] = intepolated_dataframe['iff'] / intepolated_dataframe['i0']
    # intepolated_dataframe['mu_r'] = np.log( intepolated_dataframe['it'] / intepolated_dataframe['ir'] )

    if sort:
        return intepolated_dataframe.sort_values('energy')
    else:
        return intepolated_dataframe

