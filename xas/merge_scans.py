import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import matplotlib

matplotlib.use("TkAgg")


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


if __name__ == "__main__":

    def test():
        PATH = "/home/charles/Desktop/search_and_merge"

        df_uid = pd.read_json(f"{PATH}/test_df_uid.json")
        df_data = pd.read_json(f"{PATH}/test_data.json")  # use uids for keys
        df_metadata = pd.read_json(f"{PATH}/test_metadata.json")  # use uids for keys

        sorted_sgs = sort_scangroups(df_uid)

        for sg in sorted_sgs:
            uid_list = sorted_sgs[sg]
            data_list = [pd.DataFrame(df_data[uid][0]) for uid in sorted_sgs[sg]]
            md_list = [df_metadata[uid] for uid in sorted_sgs[sg]]

            print(f"Scan group #{sg}")
            scan_group = ScanGroup(uid_list, data_list, md_list)
            scan_group.check_saturation()
            scan_group.check_amplitude()
            scan_group.plot_scans()

    test()
