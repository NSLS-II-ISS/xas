"""
Functions for testing sklearn outlier detection estimators on XAS data
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


def prenormalize_data(data: np.ndarray):
    """
    Prenormalize spectra to adjust for offset and scaling.
    Uses least square regression on median data to determine 
    scale and offset for each spectrum.

    Returns array of prenormalized data. 
    """
    n_curves, n_pts = data.shape
    data_out = np.zeros(data.shape)
    data_median = np.median(data, axis=0)

    for i in range(n_curves):
        basis = np.vstack((np.ones(n_pts), data[i, :])).T
        c, _, _, _ = np.linalg.lstsq(basis, data_median, rcond=None)
        data_out[i, :] = basis @ c
    return data_out


def outlier_plot(x_vals, data: np.ndarray, outlier_labels: np.ndarray):
    ax = plt.axes()
    line_colors = ["r" if ol < 0 else "k" for ol in outlier_labels]
    ax.set_prop_cycle(color=line_colors)
    ax.plot(x_vals, data.T)
    return ax


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


def fit_trimmed(estimator, data: np.ndarray, trim_fraction: float):
    """Shorthand function to fit estimator with trimmed data"""
    trim_data = stats.trimboth(data, trim_fraction)
    estimator.fit(trim_data) 


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
