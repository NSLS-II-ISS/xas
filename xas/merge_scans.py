import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import matplotlib

matplotlib.use("TkAgg")


def sort_uids(df_uid):
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
        self.uids = uids
        self.data = data
        self.metadata = metadata

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

    def check_saturation(self):
        # for ind, currents in enumerate(
        #     zip(self.i0_mV, self.it_mV, self.ir_mV, self.iff_mV)
        # ):
        #     current_array = np.array(currents)
        #     if np.any(current_array < -3200):
        #         print("saturated - need to check i0, it, iff")
        for df, md_df in zip(self.data, self.metadata):
            # md_df keys may be different in older scans
            i0_mV = df["i0"] * 10 ** md_df["ch1_amp_gain"] * 1000
            it_mV = df["it"] * 10 ** md_df["ch2_amp_gain"] * 1000
            ir_mV = df["ir"] * 10 ** md_df["ch3_amp_gain"] * 1000
            iff_mV = df["iff"] * 10 ** md_df["ch4_amp_gain"] * 1000
            i_array = np.array([i0_mV, it_mV, ir_mV, iff_mV])
            if np.any(i_array < -3200):
                print("saturated - need to check i0, it, iff")
            if np.all(i_array > 0):
                print("no useful signal - need to check i0 only")


if __name__ == "__main__":

    def test():
        PATH = "/home/charles/Desktop/search_and_merge"

        df_uid = pd.read_json(f"{PATH}/test_df_uid.json")
        df_data = pd.read_json(f"{PATH}/test_data.json")  # use uids for keys
        df_metadata = pd.read_json(f"{PATH}/test_metadata.json")  # use uids for keys

        sorted_uids = sort_uids(df_uid)

        for sg in sorted_uids:
            uid_list = sorted_uids[sg]
            data_list = [pd.DataFrame(df_data[uid][0]) for uid in sorted_uids[sg]]
            md_list = [df_metadata[uid] for uid in sorted_uids[sg]]

            print(f"Scan group #{sg}")
            scan_group = ScanGroup(uid_list, data_list, md_list)
            scan_group.check_saturation()

    test()
