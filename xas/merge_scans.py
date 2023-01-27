import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
from dataclasses import dataclass
import databroker
from xas.file_io import load_binned_df_from_file
from xas.analysis import average_scangroup, average_scangroup_from_files, check_scan_from_file
from xas.analysis import prenormalize_data, standardize_energy_grid

import matplotlib
from functools import reduce
matplotlib.use("TkAgg")



def find_all_scans_for_element(element, db):
    return db.search({'element' : element})

def save_scans_into_dictionary(elements_to_search, number_of_scans=5000, scan_shift=10, outpath=''):
    output = {}
    db = databroker.catalog['iss-local']
    for element in elements_to_search:
        idx1 = number_of_scans + scan_shift
        idx2 = scan_shift
        print(idx1, idx2)
        element_uids = list(find_all_scans_for_element(element, db))[::-1][-idx1 : -idx2]
        n_uids = len(element_uids)
        for i, uid in enumerate(element_uids):
            if i % 100 == 0: print(f'Progress for {element}: {i} / {n_uids}')
            # start = db[uid].start
            start = db[uid].metadata['start']
            stop = db[uid].metadata['stop']
            if (stop is None) or (stop['exit_status'] == 'fail'):
                continue
            filename = start['interp_filename'][:-3] + 'dat'
            filename.replace('/nsls2/xf08id/users', '/nsls2/data/iss/legacy/processed')
            # print(filename)
            try:
                df, _ = load_binned_df_from_file(filename)
                # df = pd.DataFrame({"test": np.random.rand(10)})
                output[uid] = {'metadata' : {**start},
                               'data' : df.to_dict()}
            except Exception as e:
                print(e)
                pass
            # if filename.startswith('/nsls2/xf08id/users')
    with open(outpath, 'w') as f:
        json.dump(output, f)

    return output

# save_scans_into_dictionary(['Fe', 'Co', 'Ni'], number_of_scans=2000,
#                            outpath='/nsls2/data/iss/legacy/Sandbox/Charles/test_data_2023_01_27.json')


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
    Sort scans in `df_uid` into scan_groups based on metadata, reduced scan name,
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


def check_scans_quality(df_uid: pd.DataFrame):
    """
    Check data quality for each scan in `df_uid`. New columns are
    added to `df_uid` for "mut_good", "mur_good", and "muf_good" with bools indicating
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
def populate_from_files(df_uid: pd.DataFrame) -> None:
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

def average_scan_groups_from_files(df_uid: pd.DataFrame):
    for scan_group in df_uid['scan_group'].unique():
        sub_df_uid = df_uid[df_uid['scan_group'] == scan_group]
        if np.all(sub_df_uid[['mut_good', 'muf_good', 'mur_good']].values):
            filenames = sub_df_uid['filename'].tolist()
            avg_mu_df, outliers_dict, dev_from_mean_dict = average_scangroup_from_files(filenames)
            print(scan_group, outliers_dict)
            # plt.subplot(311)
            plt.plot(avg_mu_df['energy'], avg_mu_df['mut'])

            # plt.subplot(312)
            plt.plot(avg_mu_df['energy'], avg_mu_df['mur'])

            # plt.subplot(313)
            plt.plot(avg_mu_df['energy'], avg_mu_df['muf'])


def get_df_uid_from_files():
    df_uid = search_db_for_entries(["Fe", "Cu"], 5000)
    df_uid = filter_df_uid_by_strings(df_uid)
    df_uid['reduced_name'] = df_uid['name'].apply(reduce_name)
    group_scans(df_uid)
    check_scans_quality(df_uid)
    populate_from_files(df_uid)
    return df_uid   
    

def get_test_df_uid():
    df_uid = search_db_for_entries(["Fe", "Cu"], 5000)
    df_uid = filter_df_uid_by_strings(df_uid)
    df_uid['reduced_name'] = df_uid['name'].apply(reduce_name)
    group_scans(df_uid)
    check_scans_quality(df_uid)
    populate_from_files(df_uid)
    df_uid.to_json("test_df_uid_5000.json")
    print(df_uid["data"])
    plt.figure(1, clear=True)
    average_scan_groups_from_files(df_uid)
    plt.show()


def drop_empty_dfs(df_uid: pd.DataFrame):
    for i, df in zip(df_uid.index, df_uid["data"]):
        if df.empty:
            df_uid.drop(i, inplace=True)


def calculate_mus(data_col: pd.Series) -> None:
    """Try calculating mut, muf, and mur for each dataframe in 
    df_uid["data"] using `it`, `iff`, `ir`, and `i0` for current keys."""
    for i, df in zip(data_col.index, data_col):
        try:
            df["mut"] = -np.log(df["it"] / df["i0"])
            df["muf"] = df["iff"] / df["i0"]
            df["mur"] = -np.log(df["ir"] / df["it"])
        except KeyError:
            # no key errors after dropping empty dfs
            df["mut"] = KeyError
            df["muf"] = KeyError
            df["mur"] = KeyError
            print(f"key error at {i}")


def redo_mu_good(df_uid: pd.DataFrame, columns=("mut", "muf", "mur")):
    """recheck criteria that mu is always finite (no `NaN` or `inf`) for each channel"""
    for i, scan in df_uid.iterrows():
        for col in columns:
            if not np.all(np.isfinite(scan["data"][col])):
                df_uid.loc[i, col + "_good"] = False


def calc_outliers(df_uid: pd.DataFrame):
    """
    Attempt outlier rejection for each scan group. If successful cols
    are added to df_uid for outlier tagging and deviation from mean for each mu.
    """
    df_uid[["mut_outlier", "muf_outlier", "mur_outlier"]] = None
    df_uid[["mut_devmean", "muf_devmean", "mur_devmean"]] = None

    n_scans = len(df_uid)
    n_sgs = len(df_uid["scan_group"].unique())
    print(f"Number of scans: {n_scans}, number of groups: {n_sgs}")

    n_failed = 0
    n_all_bad_mu = 0

    for sg in df_uid["scan_group"].unique():
        sg_df_uid = df_uid.loc[df_uid["scan_group"] == sg]

        # filter "bad scans" before averaging
        sg_df_uid = sg_df_uid[(sg_df_uid["mut_good"] == True) & (sg_df_uid["muf_good"] == True) & (sg_df_uid["mur_good"] == True)]
        sg_datalist = sg_df_uid["data"].to_list()
        if sg_datalist:
            try:
                _, outlier_dict, dev_dict = average_scangroup(sg_datalist)
                outlier_df = pd.DataFrame(outlier_dict, index=sg_df_uid.index)
                outlier_df.rename(columns={"mut": "mut_outlier", "muf": "muf_outlier", "mur": "mur_outlier"}, inplace=True)
                dev_df = pd.DataFrame(dev_dict, index=sg_df_uid.index)
                dev_df.rename(columns={"mut": "mut_devmean", "muf": "muf_devmean", "mur": "mur_devmean"}, inplace=True)
                df_uid.loc[sg_df_uid.index, ["mut_outlier", "muf_outlier", "mur_outlier"]] = outlier_df
                df_uid.loc[sg_df_uid.index, ["mut_devmean", "muf_devmean", "mur_devmean"]] = dev_df
            except ValueError:
                print(f"averaging failed for group {sg} due to NaN or inf being inverted")
                n_failed += 1
        else:
            n_all_bad_mu += 1


def calc_prenorm(scan_group_df: pd.DataFrame, columns=["mut", "muf", "mur"]):
    """
    Use `prenormalize_data` function on scan_group and add prenorm spectra to each
    scan's "data" dataframe
    """
    for col in columns:
        all_col_data = np.array([df[col] for df in scan_group_df["data"]])
        all_col_prenorm = prenormalize_data(all_col_data)
        for df, prenorm_data in zip(scan_group_df["data"], all_col_prenorm):
            df[col + "_prenorm"] = prenorm_data


def plot_prenormed_scangroups(df_uid: pd.DataFrame):
    for sg in df_uid["scan_group"].unique():
        sg_df_uid = df_uid.loc[df_uid["scan_group"] == sg]
        try:
            calc_prenorm(sg_df_uid)
            fig, (ax1, ax2) = plt.subplots(1, 2)
            fig.suptitle(f"Scan group #{sg}")
            ax1.set_title("original data")
            ax2.set_title("after prenormalization")
            for _, scan in sg_df_uid.iterrows():
                df = scan["data"]    
                ax1.plot(df["energy"], df["mut"])
                ax1.plot(df["energy"], df["muf"])
                ax1.plot(df["energy"], df["mur"])
                
            for _, scan in sg_df_uid.iterrows():
                df = scan["data"]    
                ax2.plot(df["energy"], df["mut_prenorm"])
                ax2.plot(df["energy"], df["muf_prenorm"])
                ax2.plot(df["energy"], df["mur_prenorm"])
        except Exception:
            print(f"error plotting prenorm for {sg}")
        plt.show()


def plot_scangroup_outliers(scan_group_df: pd.DataFrame, ax_og: plt.Axes, ax_pre: plt.Axes, columns=["mut", "muf", "mur"]):
    """plot original and prenormalized data for a single scan group, with outliers indicated in red"""
    calc_prenorm(scan_group_df)
    ax_og.set_title("original data")
    ax_pre.set_title("after prenormalization")
    for _, scan in scan_group_df.iterrows():
        for col in columns:
            if scan[col + "_outlier"] and not np.isnan(scan[col + "_outlier"]):
                ax_og.plot(scan["data"]["energy"], scan["data"][col], c="r")
                ax_pre.plot(scan["data"]["energy"], scan["data"][col + "_prenorm"], c="r")
                outlier_devmean = scan[col + "_devmean"]
                print(col + f"outlier devmean = {outlier_devmean}")
            else:
                ax_og.plot(scan["data"]["energy"], scan["data"][col], c="k")
                ax_pre.plot(scan["data"]["energy"], scan["data"][col + "_prenorm"], c="k")


def plot_each_scangroups_outliers(df_uid: pd.DataFrame):
    """
    Plot each scan group that has any identified outliers individually.
    """
    for sg in df_uid["scan_group"].unique():
        sg_df_uid = df_uid.loc[df_uid["scan_group"] == sg]
        fig, (ax1, ax2) = plt.subplots(1, 2)
        fig.suptitle(f"scan group #{sg}")
        outlier_labels = sg_df_uid[["mut_outlier", "muf_outlier", "mur_outlier"]].values
        if np.any(outlier_labels):
            try:
                plot_scangroup_outliers(sg_df_uid, ax1, ax2)
                plt.show()
            except ValueError as ve:
                print(ve)
        plt.close(fig)


def plot_all_scangroups_outliers(df_uid: pd.DataFrame):
    """
    Plot all scan groups that have outliers on one figure. Outliers are labeled
    in red. 
    """
    fig, (ax1, ax2) = plt.subplots(1, 2)
    for sg in df_uid["scan_group"].unique():
        sg_df_uid = df_uid.loc[df_uid["scan_group"] == sg]
        outlier_labels = sg_df_uid[["mut_outlier", "muf_outlier", "mur_outlier"]].values
        devmean_values = sg_df_uid[["mut_devmean", "muf_devmean", "mur_devmean"]].values
        if np.any(outlier_labels):
            print(outlier_labels)
            print(devmean_values)
            try:
                plot_scangroup_outliers(sg_df_uid, ax1, ax2)
            except ValueError as ve:
                print(ve)
    plt.show()


def trimmed_zscores(all_data: np.ndarray, trim_fraction=0.2):
    trim_data = stats.trimboth(all_data, trim_fraction, axis=0)
    trim_mean = np.mean(trim_data, axis=0)
    trim_std = np.std(trim_data, ddof=1, axis=0)
    trim_zscores = ((all_data - trim_mean) / trim_std) ** 2
    return trim_zscores


def zscore_outlier_rejection(sg_df_uid: pd.DataFrame, columns=["mut", "muf", "mur"], threshold=25):
    good_scans_df = sg_df_uid[(sg_df_uid["mut_good"] == True) & (sg_df_uid["muf_good"] == True) & (sg_df_uid["mur_good"] == True)]
    dfs = [df for df in good_scans_df["data"]]
    if not dfs:  # empty list / no good data
        print("no good data")
        return
    dfs = standardize_energy_grid(dfs)
    if len(dfs) < 5:
        print("too few scans")
        return
    for col in columns:
        # put all scan data for col into array
        col_data = np.array([df[col] for df in dfs])
        col_data = prenormalize_data(col_data)
        col_zscores = trimmed_zscores(col_data)
        avg_zscores = np.mean(col_zscores, axis=1)

        print(col, avg_zscores)
        sg_df_uid.loc[good_scans_df.index, col + "_devmean"] = avg_zscores
        sg_df_uid.loc[good_scans_df.index, col + "_outlier"] = avg_zscores > threshold

    return sg_df_uid



def test_df_uid(path):
    matplotlib.use("TkAgg")
    df_uid = pd.read_json(path)
    df_uid["data"] = df_uid["data"].apply(pd.DataFrame)  # convert dicts
    df_uid = df_uid.reset_index(drop=True)
    drop_empty_dfs(df_uid)
    calculate_mus(df_uid["data"])
    redo_mu_good(df_uid)

    df_uid["mut_outlier"] = None
    df_uid["muf_outlier"] = None
    df_uid["mur_outlier"] = None
    df_uid["mut_devmean"] = None
    df_uid["muf_devmean"] = None
    df_uid["mur_devmean"] = None
    for sg in df_uid["scan_group"].unique():
        sg_df_uid = df_uid.loc[df_uid["scan_group"] == sg]
        # print(sg_df_uid["data"][0])
        print(f"calculating outlier rejection for {sg}")
        sg_df_uid = zscore_outlier_rejection(sg_df_uid)
        print(sg_df_uid)
        df_uid.loc[df_uid["scan_group"] == sg] = sg_df_uid
    plot_all_scangroups_outliers(df_uid)

    # calc_outliers(df_uid)
    # plot_all_scangroups_outliers(df_uid)
    return df_uid


if __name__ == "__main__":
    pass



#####################
# old ScanGroup class
#####################
def sort_scan_groups(df_uid):
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
        list of arrays, each containing i0, it, ir, and iff for each scan in Scan_Group
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



######################
# test Scan_Group class
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
#         sorted_sgs = sort_scan_groups(df_uid)
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
