import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import databroker
from xas.file_io import load_binned_df_from_file
from xas.analysis import check_scan_from_file, average_scangroup_from_files

import matplotlib
from functools import reduce
# matplotlib.use("TkAgg")


def sort_scangroups(df_uid):
    """
    returns dictionary with scan group number as keys and uids as values
    """
    scan_groups = df_uid["scan_group"]
    all_uids = df_uid["uid"]

    sorted_scans = {}

    for sg, uid in zip(scan_groups, all_uids):
        if sg not in sorted_scans:
            sorted_scans[sg] = []
        sorted_scans[sg].append(uid)

    return sorted_scans


class ScanGroup:
    def __init__(self, uids, data, metadata) -> None:
        self.uids = uids  # list of uid strings
        self.data = data  # list of dataframes
        self.metadata = metadata  # list of pandas series

        # md_df keys may be different for older scans
        self.i0_mV = [
            np.array(df["i0"] * 10 ** md_df["ch1_amp_gain"] * 1000)
            for df, md_df in zip(self.data, self.metadata)
        ]
        self.it_mV = [
            np.array(df["it"] * 10 ** md_df["ch2_amp_gain"] * 1000)
            for df, md_df in zip(self.data, self.metadata)
        ]
        self.ir_mV = [
            np.array(df["ir"] * 10 ** md_df["ch3_amp_gain"] * 1000)
            for df, md_df in zip(self.data, self.metadata)
        ]
        self.iff_mV = [
            np.array(df["iff"] * 10 ** md_df["ch4_amp_gain"] * 1000)
            for df, md_df in zip(self.data, self.metadata)
        ]

        self.saturated_scan_inds = None
        self.zero_amp_inds = None

    @property
    def num_scans(self):
        return len(self.uids)

    @property
    def energy_array_list(self):
        _energy_array_list = [np.array(df["energy"]) for df in self.data]
        return _energy_array_list

    @property
    def current_array_list(self):
        """
        list of arrays, each containing i0, it, ir, and iff for each scan in ScanGroup
        """
        _current_array_list = [
            np.array(currents)
            for currents in zip(self.i0_mV, self.it_mV, self.ir_mV, self.iff_mV)
        ]
        return _current_array_list

    def check_saturation(self):
        _saturated_scans = []
        for ind, i_array in enumerate(self.current_array_list):
            if np.any(i_array < -3200):
                print(
                    f"saturated - need to check i0, it, iff for scan: \n {self.uids[ind]}"
                )
                _saturated_scans.append(ind)
        self.saturated_scan_inds = _saturated_scans

    def check_amplitude(self):
        _zero_amp = []
        for ind, i_array in enumerate(self.current_array_list):
            if np.all(i_array > 0):
                print(
                    f"no useful signal - need to check i0 for scan: \n {self.uids[ind]}"
                )
                _zero_amp.append(ind)
        self.zero_amp_inds = _zero_amp

    def plot_scans(self, ax=None):
        if ax is None:
            ax = plt.axes()
        for e_array, i_array in zip(self.energy_array_list, self.current_array_list):
            print(i_array)
            ax.plot(e_array, i_array.T)
            plt.show()


def find_all_scans_for_element(element, db):
    return db.search({'element' : element})

def search_db_for_scans(db, df_uid, elements_to_search, number_of_scans):
    for element in elements_to_search:
        element_uids = list(find_all_scans_for_element(element, db))[::-1][:number_of_scans]
        n_uids = len(element_uids)
        for i, uid in enumerate(element_uids):
            if i % 100 == 0: print(f'Progress: {i} / {n_uids}')
            # start = db[uid].start
            start = db[uid].metadata['start']
            stop = db[uid].metadata['stop']
            if stop['exit_status'] == 'fail':
                continue

            for key in df_uid.keys():
                try:
                    if key == 'filename':
                        val = start['interp_filename']
                        val = val[:-3] + 'dat'
                    else:
                        val = start[key]
                except KeyError:
                    val = None
                df_uid[key].append(val)

def search_db_for_entries(elements_to_search: list[str], number_of_scans: int):
    """
    Search database for scans up to `number_of_scans` for each element in `elements_to_search`.
    
    Returns `df_uid`, a DataFrame containing the metadata for all found scans.  
    """
    db = databroker.catalog['iss-local']
    # db_old = databroker.catalog['iss']
    df_uid = {'element': [], 'edge': [], 'uid': [], 'year': [],
              'cycle': [], 'PROPOSAL': [], 'time': [], 'name': [],
              'filename': [], 'scan_group_uid': [], 'sample_uid': []}
    search_db_for_scans(db, df_uid, elements_to_search, number_of_scans)
    df_uid = pd.DataFrame(df_uid)
    df_uid['scan_group'] = None
    return df_uid


def filter_df_uid_by_strings(df_uid: pd.DataFrame, strings_to_drop=['test', 'bla', 'calibration']):
    """ Remove scans containing any of `strings_to_drop` from `df_uid` and return result """
    return df_uid[~df_uid['name'].str.contains('|'.join(strings_to_drop))]


def reduce_name(name):
    if 'pos' in name:
        if ' (pos ' in name:
            idx = name.find(' (pos')
        else:
            idx = name.find(' pos')
        reduced_name = name[:idx]
    else:
        reduced_name = name[:-5]
    return reduced_name


def get_relevant_scans_for_row(df_uid, row, time_window=600, ls_dist_thresh=4):
    # filter by element/proposal
    filter_keys = ['element', 'edge', 'year', 'cycle', 'PROPOSAL']
    filter_list = [df_uid[k] == row[k] for k in filter_keys]

    # filter by time
    filter_list.append((row['time'] - df_uid['time']).abs() < time_window)

    # filter by ungrouped
    filter_list.append([(i is None) for i in df_uid['scan_group']])

    # filter by name
    filter_list.append((~(row['reduced_name'] != df_uid['reduced_name'])).tolist())

    relevant_scans = reduce(lambda x, y: x & y, filter_list)
    return relevant_scans


def group_scans(df_uid: pd.DataFrame, time_window=600):
    """
    Sort scans in `df_uid` into scangroups based on metadata, reduced scan name,
    and scan time proximity.
    """
    df_uid['scan_group'] = None
    scan_group_id = 0
    for i in df_uid.index:
        row = df_uid.loc[i]

        if row['scan_group'] is None:
            relevant_scans = get_relevant_scans_for_row(df_uid, row, time_window=time_window)

            if any(relevant_scans):
                df_uid.loc[relevant_scans, ('scan_group',)] = scan_group_id
            else:
                df_uid.loc[i, ('scan_group',)] = scan_group_id
            scan_group_id += 1


def check_scans_in_df_uid(df_uid: pd.DataFrame):
    """
    Check data quality for each scan in `df_uid`. New columns are
    added to `df_uid` for "mut", "mur", and "muf" with bools indicating
    good data quality (`True`) or poor data quality (`False`).
    """
    df_uid['mut_good'] = None
    df_uid['mur_good'] = None
    df_uid['muf_good'] = None
    db = databroker.catalog['iss-local']
    for i in df_uid.index:
        uid = df_uid.loc[i]['uid']
        md = db[uid].metadata['start']
        filename = df_uid.loc[i]['filename']
        print(i, end=' ')
        try:
            mu_good = check_scan_from_file(filename, md)
            print(mu_good)
            df_uid.loc[i, ('mut_good', )] = mu_good['mut']
            df_uid.loc[i, ('mur_good', )] = mu_good['mur']
            df_uid.loc[i, ('muf_good', )] = mu_good['muf']
        except Exception as e:
            print(e)
            pass

# df_uid that contains data
def populate_data_from_files(df_uid: pd.DataFrame) -> None:
    """
    Attempts to load data from `"filename"` column for each index of df_uid.
    If data is loaded it will be stored in "data" column at the same index
    as `pd.DataFrame` converted to `dict` (`df.to_dict()`).
    """
    df_uid['data'] = None
    for i in df_uid.index:
        filename = df_uid.loc[i]['filename']
        try:
            df, _ = load_binned_df_from_file(filename)
            # df = pd.DataFrame({"test": np.random.rand(10)})
            df_uid.at[i, "data"] = df.to_dict()
        except Exception as e:
            print(e)
            pass

# populate_df_uid_with_data(df_uid) # optional - populate with data and transfer for further local assessment of outlier rejection

def average_scangroups_in_df_uid(df_uid: pd.DataFrame):
    for scangroup in df_uid['scan_group'].unique():
        sub_df_uid = df_uid[df_uid['scan_group'] == scangroup]
        if np.all(sub_df_uid[['mut_good', 'muf_good', 'mur_good']].values):
            filenames = sub_df_uid['filename'].tolist()
            avg_mu_df, outliers_dict, dev_from_mean_dict = average_scangroup_from_files(filenames)
            print(scangroup, outliers_dict)
            # plt.subplot(311)
            plt.plot(avg_mu_df['energy'], avg_mu_df['mut'])

            # plt.subplot(312)
            plt.plot(avg_mu_df['energy'], avg_mu_df['mur'])

            # plt.subplot(313)
            plt.plot(avg_mu_df['energy'], avg_mu_df['muf'])



def test_df_uid():
    df_uid = search_db_for_entries(["Fe", "Cu"], 5000)
    df_uid = filter_df_uid_by_strings(df_uid)
    df_uid['reduced_name'] = df_uid['name'].apply(reduce_name)
    group_scans(df_uid)
    check_scans_in_df_uid(df_uid)
    populate_data_from_files(df_uid)
    df_uid.to_json("test_df_uid_5000.json")
    print(df_uid["data"])
    plt.figure(1, clear=True)
    average_scangroups_in_df_uid(df_uid)
    plt.show()


if __name__ == "__main__":
    test_df_uid()



######################
# test ScanGroup class
######################
# if __name__ == "__main__":
#
#     def test():
#         PATH = "/home/charles/Desktop/search_and_merge"
#
#         df_uid = pd.read_json(f"{PATH}/test_df_uid.json")
#         df_data = pd.read_json(f"{PATH}/test_data.json")  # use uids for keys
#         df_metadata = pd.read_json(f"{PATH}/test_metadata.json")  # use uids for keys
#
#         sorted_sgs = sort_scangroups(df_uid)
#
#         for sg in sorted_sgs:
#             uid_list = sorted_sgs[sg]
#             data_list = [pd.DataFrame(df_data[uid][0]) for uid in sorted_sgs[sg]]
#             md_list = [df_metadata[uid] for uid in sorted_sgs[sg]]
#
#             print(f"Scan group #{sg}")
#             scan_group = ScanGroup(uid_list, data_list, md_list)
#             scan_group.check_saturation()
#             scan_group.check_amplitude()
#             scan_group.plot_scans()
#
#     test()
