import matplotlib

matplotlib.use("TkAgg")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from xas.file_io import load_binned_df_from_file, save_binned_df_as_file
import os


GAIN_DICT = {
    "i0": ["ch1_amp_gain"],
    "it": ["ch2_amp_gain"],
    "ir": ["ch3_amp_gain"],
    "iff": ["ch4_amp_gain"],
}


def degain(
    df: pd.DataFrame,
    md: dict,
    df_keys=["i0", "it", "ir", "iff"],
):
    """
    Use channel gain to calculate current in mV.
    Returns otherwise identical dataframe with current columns in mV.
    """
    df_mV = df.copy(deep=True)
    for k in df_keys:
        # use first key shared by GAIN_DICT[k] and md
        md_gain_key = list(set(GAIN_DICT[k]) & set(md.keys()))[0]
        df_mV[k] = df[k] * 10 ** md[md_gain_key] * 1000
    return df_mV


def check_saturation(df: pd.DataFrame, keys=["i0", "it", "ir", "iff"], threshold=3200):
    """
    Check currents (in mV) are not saturated relative to threshold.

    Returns `dict` with current channel keys (default `["i0", "it", "ir", "iff"]`) and `bool`
    values. `True` corresponds to unsaturated/good current. `False` corresponds to saturated current.
    """
    sat_dict = dict()
    for k in keys:
        if np.any(np.abs(df[k]) > threshold):
            sat_dict[k] = False
        else:
            sat_dict[k] = True
    return sat_dict


def check_amplitude(df: pd.DataFrame, keys=["i0", "it", "ir", "iff"], threshold=20):
    """
    Check currents (in mV) have sufficient non-zero amplitude relative to threshold.

    Returns `dict` with current channel keys (default `["i0", "it", "ir", "iff"]`) and `bool`
    values. `True` corresponds to good amplitude current. `False` corresponds to poor amplitude.
    """
    low_amp_dict = dict()
    for k in keys:
        if np.all(np.abs(df[k]) < threshold):
            low_amp_dict[k] = False
        else:
            low_amp_dict[k] = True
    return low_amp_dict


def check_mu_values(df: pd.DataFrame):
    """
    Calculate mu values for each channel (transmission, reference, fluorescence) and
    check all mu values are valid (no `NaN`, `inf`, etc.).

    Returns `dict` with mu channel keys (`"mut"`, `"mur"`, `"muf"`) and `bool` values.
    `True` corresponds to all valid values. `False` indicates at least one value is invalid
    for that channel.
    """
    mut = -np.log(df["it"] / df["i0"])
    mur = -np.log(df["ir"] / df["it"])
    muf = df["iff"] / df["i0"]
    valid_values = {
        "mut": np.all(np.isfinite(mut)),
        "mur": np.all(np.isfinite(mur)),
        "muf": np.all(np.isfinite(muf)),
    }
    return valid_values


def check_scan(df: pd.DataFrame, md: dict):
    """
    Combine `degain`, `check_saturation`, `check_amplitude`, and `check_mu_values`
    to determine data (mu) quality for each channel (transmission, fluorescence, reference).

    Returns `dict` with mu channel keys (`"mut"`, `"mur"`, `"muf"`) and `bool` values.
    `True` indicates good data quality. `False` indicates poor data quality.
    """

    df_mV = degain(df, md)
    unsaturated_currents = check_saturation(df_mV)
    good_amp_currents = check_amplitude(df_mV)
    valid_mu = check_mu_values(df_mV)
    mu_good = {mu: False for mu in ["mut", "muf", "mur"]}  # default all False
    if good_amp_currents["i0"] and unsaturated_currents["i0"]:
        if unsaturated_currents["it"]:
            mu_good["mut"] = valid_mu["mut"]
        if unsaturated_currents["iff"]:
            mu_good["muf"] = valid_mu["muf"]
    if (
        good_amp_currents["it"]
        and unsaturated_currents["it"]
        and unsaturated_currents["ir"]
    ):
        mu_good["mur"] = valid_mu["mur"]
    return mu_good


def check_scan_from_file(filename, md):

    # stable function for loading data from beamline file
    df, _ = load_binned_df_from_file(filename)

    return check_scan(df, md)

def standardize_energy_grid(dfs: list[pd.DataFrame], energy_key="energy", master_idx=0):
    """
    Interpolate data in each df to match master energy grid. Master grid can be
    selected via `master_idx` (default=0 or first DataFrame).
    """
    energy_master = dfs[master_idx][energy_key]
    dfs_out = [dfs[master_idx]]
    dfs_copy = dfs[:]
    dfs_copy.pop(master_idx)
    for df in dfs_copy:
        # if df.equals(dfs[master_idx]):
        #     continue
        _df = {energy_key: energy_master}
        for column in df.columns:
            if column != energy_key:
                _df[column] = np.interp(energy_master, df[energy_key], df[column])
        dfs_out.append(pd.DataFrame(_df))
    return dfs_out


def prenormalize_data(data: np.ndarray, energy: np.ndarray):
    """
    Prenormalize spectra to adjust for offset and scaling.
    Uses least square regression on median data to determine
    scale and offset for each spectrum.

    Returns array of prenormalized data.
    """
    n_curves, n_pts = data.shape
    data_out = np.zeros(data.shape)
    data_median = np.median(data, axis=0)

    energy_scaled = (energy - energy.min()) / (energy.max() - energy.min())

    for i in range(n_curves):
        basis = np.vstack(
            (
                np.ones(n_pts),
                data[i, :],
                energy_scaled,
                energy_scaled**2,
                energy_scaled**3,
            )
        ).T
        c, _, _, _ = np.linalg.lstsq(basis, data_median, rcond=None)
        data_out[i, :] = basis @ c
    return data_out


#
# def check_for_outliers(all_data: np.ndarray, trim_fraction=0.2, threshold=25):
#     """
#     Check set of scans for outliers. Scan data should be in rows in `all_data` array.
#
#     Uses deviation from trimmed mean to determine outliers. Returns boolean array
#     with `True` marking a scan as outlier, and array of deviation from trimmed mean values.
#     """
#     # all_data = np.array([df["mut"] for df in dfs])
#     n_pts = all_data.shape[1]
#     trim_data = stats.trimboth(all_data, trim_fraction, axis=0)
#     trim_mean = np.mean(trim_data, axis=0)
#     trim_std = np.std(trim_data, axis=0)
#     deviation_from_mean = (
#         np.sum(((all_data - trim_mean) / trim_std) ** 2, axis=1) / n_pts
#     )
#     return deviation_from_mean > threshold, deviation_from_mean
#
#
# def average_scangroup(
#     dfs: list[pd.DataFrame],
#     columns=["mut", "muf", "mur"],
#     energy_key="energy",
#     master_idx=0,
# ):
#     """
#     Takes in scangroup (list of DataFrames), standardizes energy grids,
#     checks each channel for outliers, and averages all non-outlier data.
#
#     ### Returns
#     pd.DataFrame(avg_data)
#         DataFrame with averaged data
#     outliers_dict
#         dict of boolean arrays indicating outliers for each channel
#     dev_from_mean_dict
#         dict of arrays with deviation from mean values for each channel
#
#     If there are too few scans in group (< 5) outlier rejection is not performed
#     so `outlier_dict` and `dev_from_mean_dict` are returned empty.
#     """
#     dfs = standardize_energy_grid(dfs, master_idx=master_idx)
#     avg_data = {energy_key: dfs[master_idx][energy_key]}
#     outliers_dict = {}
#     dev_from_mean_dict = {}
#     for col in columns:
#         all_data = np.array([df[col] for df in dfs])
#         all_data_prenorm = prenormalize_data(all_data)
#         # plt.plot(dfs[0][energy_key], all_data.T, "k-", alpha=0.25)
#         # plt.plot(dfs[0][energy_key], all_data_prenorm.T, "k-", alpha=1)
#
#         n_scans = all_data.shape[0]
#         if n_scans >= 5:  # only check for outliers on datasets with at least 5 scans
#             outliers, dev_from_mean = check_for_outliers(all_data_prenorm)
#             outliers_dict[col] = outliers
#             dev_from_mean_dict[col] = dev_from_mean
#         else:
#             # if less than 5 scans set outliers to all False
#             outliers = np.full(n_scans, False)
#             outliers_dict[col] = outliers
#         avg_data[col] = np.mean(all_data[~outliers], axis=0)
#
#     return pd.DataFrame(avg_data), outliers_dict, dev_from_mean_dict
#
#
# def average_scangroup_from_files(filenames):
#     """
#     Load data into DataFrames from files, calculate mu for each channel, and return
#     result of `average_scangroup` on DataFrames.
#     """
#     dfs = []
#     for filename in filenames:
#         df, _ = load_binned_df_from_file(filename)
#         df_mu = {
#             "energy": df["energy"],
#             "mut": -np.log(df["it"] / df["i0"]),
#             "muf": df["iff"] / df["i0"],
#             "mur": -np.log(df["ir"] / df["it"]),
#         }
#         dfs.append(pd.DataFrame(df_mu))
#     return average_scangroup(dfs)
#
#
# def test():
#     PATH = "/home/charles/Desktop/search_and_merge"
#     df_data = pd.read_json(f"{PATH}/test_data.json")  # use uids for keys
#     df_metadata = pd.read_json(f"{PATH}/test_metadata.json")  # use uids for keys
#     test_df = pd.DataFrame(df_data.iloc[:, 1][0])
#     test_md = dict(df_metadata.iloc[:, 1])
#
#     test_df_mV = degain(test_df, test_md)
#     print(test_df_mV, "\n")
#     print("unsaturated scans:  ", check_saturation(test_df_mV))
#     print("good amplitude:  ", check_amplitude(test_df_mV))
#     print("good data quality:  ", check_scan(test_df, test_md))
#
#     plt.plot(test_df["energy"], -np.log(test_df["it"] / test_df["i0"]))
#     plt.plot(test_df["energy"], -np.log(test_df["ir"] / test_df["it"]))
#     plt.plot(test_df["energy"], test_df["iff"] / test_df["i0"])
#
#     plt.show()
#
#
# if __name__ == "__main__":
#     test()
#


# def filenames_from_dir(path, base="", ext=""):
#     filenames = []
#     if type(base) == str:
#         (_, _, all_filenames) = next(os.walk(path))
#         for f in all_filenames:
#             if f.startswith(base) and f.endswith(ext) and ("merge" not in f):
#                 filenames.append(os.path.join(path, f))
#     elif type(base) == list:
#         for f in base:
#             filenames.append(os.path.join(path, f))
#     else:
#         print("Invalid file type. Return None")
#         return None
#     return np.sort(filenames)


def average_scans(path, base, ext=""):
    filenames = filenames_from_dir(path, base=base, ext=ext)
    header_av = "# merge\n"
    path_av = os.path.join(path, base + " merged" + ext)
    # define energy grid and read the data to merge
    n = len(filenames)
    energy_lo = np.zeros(n)
    energy_hi = np.zeros(n)
    df_list = []
    for i, f in enumerate(filenames):
        df, header = load_binned_df_from_file(f)
        df_list.append(df)
        energy_lo[i] = df.iloc[:, 0].min()
        energy_hi[i] = df.iloc[:, 0].max()
        _, new_line = os.path.split(f)
        header_av += "# " + new_line + "\n"

    energy_lo = energy_lo.max()
    energy_hi = energy_hi.min()
    E = np.array(df_list[0].iloc[:, 0])
    E = E[(E >= energy_lo) & (E <= energy_hi)]

    # create new df
    idx = pd.Index(np.arange(E.size))
    columns = df_list[0].columns
    df_av = pd.DataFrame(index=idx, columns=columns)
    df_av.iloc[:, 0] = E
    df_av.iloc[:, 1:] = 0

    # average
    for each_df in df_list:
        data = np.array(each_df.iloc[:, 1:])
        n_cols = data[0, :].size
        E_cur = np.array(each_df.iloc[:, 0])
        data_int = np.array([np.interp(E, E_cur, data[:, i]) for i in range(n_cols)]).T
        df_av.iloc[:, 1:] += data_int
    df_av.iloc[:, 1:] /= n

    columns_str = "  ".join(columns)
    fmt = "%12.6f " + (" ".join(["%12.6e" for i in range(len(columns) - 1)]))
    np.savetxt(
        path_av,
        df_av.values,
        delimiter=" ",
        fmt=fmt,
        header=f"# {columns_str}",
        comments=header_av,
    )
    # dfgd
    # save_binned_df_as_file(path_av, df_av, header_av)

    # plt.figure(1)
    # plt.clf()
    # for each_df in df_list:
    #     plt.plot(each_df['energy'], each_df['iff'])
    # plt.plot(df_av['energy'], df_av['iff'], 'r-')
    # dfgef


# def test_interpolation(fpath):


# def read_offsets_from_folder(folder):
#     files = filenames_from_dir(folder,base='', ext='.dat')
#     gains = np.array([3, 4, 5, 6, 7])
#     data_dict = {'apb_ave_ch1_mean' : [[], [], [], [], []],
#                  'apb_ave_ch2_mean' : [[], [], [], [], []],
#                  'apb_ave_ch3_mean' : [[], [], [], [], []],
#                  'apb_ave_ch4_mean' : [[], [], [], [], []],
#                  'apb_ave_ch5_mean' : [[], [], [], [], []],
#                  'apb_ave_ch6_mean' : [[], [], [], [], []],
#                  'apb_ave_ch7_mean' : [[], [], [], [], []],
#                  'apb_ave_ch8_mean' : [[], [], [], [], []]}
#
#     for file in files:
#         df = pd.read_csv(file)
#         this_gain = int(file[-5])
#         idx_gain = np.where(this_gain == gains)[0][0]
#         for key in data_dict.keys():
#             offset_value = df[key].values[0]
#             # print(key, idx_gain)
#             data_dict[key][idx_gain].append(offset_value)
#
#     return gains, data_dict
#
#
#
# ggg = dd['apb_ave_ch3_mean'][1]
# plt.figure(1)
# plt.close()
# plt.plot(ggg, 'k.')
# plt.plot([0, len(ggg)], [np.median(ggg), np.median(ggg)], 'r-')
#
# out_dict = {}
# for key in dd.keys():
#     bla = {}
#     for i, gain in enumerate(gains):
#         bla[int(gain)] = float(np.median(dd[key][i]))
#     out_dict[key] = bla
