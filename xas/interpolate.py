import numpy as np
import pandas as pd



def interpolate(dataset,key_base = 'i0', sort=True):
    interpolated_dataset = {}
    min_timestamp = max([dataset.get(key).iloc[0, 0] for key in dataset])
    max_timestamp = min([dataset.get(key).iloc[len(dataset.get(key)) - 1, 0] for key in
                         dataset if len(dataset.get(key).iloc[:, 0]) > 5])

    try:
        if key_base not in dataset.keys():
            raise ValueError('Could not find "{}" in the loaded scan. Pick another key_base'
                             ' for the interpolation.'.format(key_base))
    except ValueError as err:
        print(err.args[0], '\nAborted...')
        return

    timestamps = dataset[key_base].iloc[:,0]

    condition = timestamps < min_timestamp
    timestamps = timestamps[np.sum(condition):]

    condition = timestamps > max_timestamp
    timestamps = timestamps[: len(timestamps) - np.sum(condition)]

    for key in dataset.keys():
        if len(dataset.get(key).iloc[:, 0]) > 5 * len(timestamps):
            time = [np.mean(array) for array in np.array_split(dataset.get(key).iloc[:, 0].values, len(timestamps))]
            val = [np.mean(array) for array in np.array_split(dataset.get(key).iloc[:, 1].values, len(timestamps))]
            interpolated_dataset[key] = np.array([timestamps, np.interp(timestamps, time, val)]).transpose()
        else:
            interpolated_dataset[key] = np.array([timestamps, np.interp(timestamps, dataset.get(key).iloc[: ,0].values,
                                                                        dataset.get(key).iloc[:,1])]).transpose()
    intepolated_dataframe = pd.DataFrame(np.vstack((timestamps, np.array([interpolated_dataset[array][:, 1] for
                                                                            array in interpolated_dataset]))).transpose())
    keys = ['timestamp']
    keys.extend(interpolated_dataset.keys())
    intepolated_dataframe.columns = keys
    if sort:
        return intepolated_dataframe.sort_values('energy')
    else:
        return intepolated_dataframe

