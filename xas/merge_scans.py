import numpy as np
import pandas as pd


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


if __name__ == "__main__":

    def test():
        PATH = "/home/charles/Desktop/search_and_merge"

        df_uid = pd.read_json(f"{PATH}/test_df_uid.json")
        df_data = pd.read_json(f"{PATH}/test_data.json")  # use uids for keys
        df_metadata = pd.read_json(f"{PATH}/test_metadata.json")  # use uids for keys

        sorted_uids = sort_uids(df_uid)

        for sg in sorted_uids:
            print(sg)
            print(sorted_uids[sg])

    test()
