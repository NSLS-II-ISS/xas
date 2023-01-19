"""
Functions for testing outlier detection estimators on XAS data
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats
mpl.use("TkAgg")
from sklearn.covariance import EllipticEnvelope
from sklearn.neighbors import LocalOutlierFactor
from sklearn.ensemble import IsolationForest
from sklearn.svm import OneClassSVM

TEST_PATH = "/home/charles/Desktop/test_data/outlier_muf.pkl"


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
        basis = np.vstack((np.ones(n_pts), data[i, :], energy_scaled, energy_scaled**2, energy_scaled**3)).T
        c, _, _, _ = np.linalg.lstsq(basis, data_median, rcond=None)
        data_out[i, :] = basis @ c
    return data_out


def outlier_plot(x_vals, data: np.ndarray, outlier_labels: np.ndarray, ax=None):
    """Return pyplot axis with outliers plotted in red and inliers plotted in black"""
    if ax is None:
        ax = plt.gca()
    line_colors = ["r" if ol < 0 else "k" for ol in outlier_labels]
    ax.set_prop_cycle(color=line_colors)
    lines = ax.plot(x_vals, data.T)
    return lines


def add_toy_outlier(data: np.ndarray) -> np.ndarray:
    """
    Add toy "outlier" to row-based data array.
    Toy outlier is made by adding scaled Gaussian noise
    to the data average.
    """
    avg_data = np.mean(data, axis=0)
    noise = np.random.randn(avg_data.size) / 1000
    toy_outlier = avg_data + noise
    new_data = np.vstack((data, toy_outlier))
    return new_data


def load_muf_data(path) -> tuple[np.ndarray, np.ndarray]:
    """
    Load muf data from pickle file. 
    Return tuple of energy array and muf spectra (in array rows).
    """
    test_dfs: pd.Series = pd.read_pickle(path)
    test_dfs = test_dfs.to_list()
    muf_data = np.array([scan["muf"] for scan in test_dfs])
    energy = test_dfs[0]["energy"]
    return energy, muf_data


def trim_by_row(data: np.ndarray, trim_fraction: float) -> np.ndarray:
    """
    Trim rows from array based on their relative distance to the 
    median row data. 
    Distance for each row is defined by the average of the sqaure of the deviation 
    from the median for each point. 
    """
    med = np.median(data, axis=0)
    dist_from_med = (data - med) ** 2
    avg_row_dist = np.mean(dist_from_med, axis=1)
    return data[avg_row_dist <= np.quantile(avg_row_dist, 1 - trim_fraction)]


def fit_trimmed(estimator, data: np.ndarray, trim_fraction: float) -> None:
    """Shorthand function to fit estimator with trimmed data"""
    trim_data = stats.trimboth(data, trim_fraction / 2)
    estimator.fit(trim_data) 


def fit_row_trimmed(estimator, data: np.ndarray, trim_fraction: float) -> None:
    trim_data = trim_by_row(data, trim_fraction)
    estimator.fit(trim_data)     


def calc_mod_chisq(data: np.ndarray) -> np.ndarray:
    mad = np.median(np.abs(data - np.median(data, axis=0)).ravel()) / 0.67449
    ksi = (data - np.median(data, axis=0)) / mad
    n_pts = data.shape[1]
    mod_chisq = 1/n_pts * np.sum(ksi**2, axis=1)
    return mod_chisq


def modified_chisq_rejection(data: np.ndarray, threshold=25) -> np.ndarray:
    mod_chisq = calc_mod_chisq(data)
    results = mod_chisq > threshold
    return np.array([-1 if res else 1 for res in results])


def compare_LOF_modchisq(df_path):
    df_uid = pd.read_pickle(df_path)
    for sg in df_uid["scan_group"].unique():
        print(f"scan group {sg}")
        sg_df_uid = df_uid[df_uid["scan_group"] == sg]
        test_dfs = sg_df_uid["data"]
        test_dfs.reset_index(drop=True, inplace=True)
        mu_data = np.array([scan["mut"] for scan in test_dfs])
        energy = test_dfs[0]["energy"]
        if mu_data.shape[0] <= 4:
            continue
        mu_prenorm = prenormalize_data(mu_data, energy)
        
        est = LocalOutlierFactor(novelty=True)
        fit_trimmed(est, mu_prenorm, 0.4)
        LOF_labels = est.predict(mu_prenorm)
        print(f"LOF labels: {LOF_labels} \n")

        modchisq_labels = modified_chisq_rejection(mu_prenorm)
        print(f"Modified chi-sq values: {calc_mod_chisq(mu_prenorm)}")
        print(f"Modified chi-sq labels: {modchisq_labels} \n")

        if np.any(np.concatenate((LOF_labels, modchisq_labels)) < 0):
            fig, (ax1, ax2) = plt.subplots(1, 2)
            outlier_plot(energy, mu_prenorm, LOF_labels, ax=ax1)
            outlier_plot(energy, mu_prenorm, modchisq_labels, ax=ax2)
            ax1.set_title("trimmed LocalOutlierFactor")
            ax2.set_title("Modified Chi-Sq")
            plt.show()


def main(estimator, add_outlier=False):
    test_dfs: pd.Series = pd.read_pickle(TEST_PATH)
    test_dfs = test_dfs.to_list()
    muf_data = np.array([scan["muf"] for scan in test_dfs])
    energy = test_dfs[0]["energy"]

    muf_prenorm = prenormalize_data(muf_data)
    
    # reduce amount of scans
    # muf_prenorm = muf_prenorm[0:5]

    if add_outlier:
        muf_prenorm = add_toy_outlier(muf_prenorm)

    if estimator is None:
        outlier_labels = np.ones(muf_prenorm.shape[0])
    else:
        outlier_labels = estimator.fit_predict(muf_prenorm)
    
    print(outlier_labels)
    print(estimator.decision_function(muf_prenorm))
    print(estimator.score_samples(muf_prenorm))
    ax = outlier_plot(energy, muf_prenorm, outlier_labels)
    
    plt.show()


if __name__ == "__main__":
    main(IsolationForest())


###
# for sg in df_uid["scan_group"].unique():
#     print(f"scan group {sg}")
#     sg_df_uid = df_uid[df_uid["scan_group"] == sg]
#     test_dfs = sg_df_uid["data"]
#     test_dfs.reset_index(drop=True, inplace=True)
#     muf_data = np.array([scan["muf"] for scan in test_dfs])
#     energy = test_dfs[0]["energy"]
#     if muf_data.shape[0] <= 1:
#         continue
#     muf_prenorm = prenormalize_data(muf_data)
#     outlier_labels = EllipticEnvelope().fit_predict(muf_prenorm)
#     print(outlier_labels)
#     if np.any(outlier_labels < 0):
#         ax = outlier_plot(energy, muf_prenorm, outlier_labels)
#         plt.show()

# def plot_scangroup(sg_df, channels=("mut", "muf", "mur")):
#     for i, scan in sg_df.iterrows():
#         for ch in channels:
#             plt.plot(scan["data"]["energy"], scan["data"][ch])
#     plt.show()
