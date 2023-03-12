from code import interact
import json
import numpy as np
import pandas as pd

import matplotlib

# matplotlib.use('TkAgg')  # something wrong with my (Charlie's) system Qt
import matplotlib.pyplot as plt

from scipy.ndimage import center_of_mass, rotate
from scipy import linalg
from sklearn.covariance import MinCovDet

from .fitting import fit_gaussian_with_estimation


def percentile_threshold_filter(im2d, pmin=5, pmax=99.5):
    """set values below percentile equal to 0"""
    filt_im = im2d.copy()
    mask = (filt_im > np.percentile(filt_im, pmin)) & (
        filt_im < np.percentile(filt_im, pmax)
    )
    filt_im[mask] = 0
    return filt_im


def fit_plane(X, Y, Z):
    """fit plane to 3d data using least square regression"""
    X, Y, Z = np.array(X), np.array(Y), np.array(Z)
    assert X.size == Y.size == Z.size
    A = np.column_stack([X, Y, np.ones(X.size)])
    C, _, _, _ = linalg.lstsq(A, Z)

    # equation for plane: z = a*x + b*y + c
    a, b, c = C
    return a, b, c


def project_pt2line(x_pt, y_pt, p_xy):
    slope, intercept = p_xy
    x_proj = (x_pt + slope * y_pt - slope * intercept) / (slope**2 + 1)
    y_proj = slope * x_proj + intercept
    return x_proj, y_proj


def get_roi(metadata, roi="roi1"):
    rois = metadata["detectors"]["Pilatus 100k"]["config"]["roi"]
    roi_ = rois[roi]
    return roi_


def get_image_array(data):
    try:
        det_image = data["pil100k_image"]
        image_array = np.array(list(det_image)).squeeze()
    except:
        det_image = data["data_vars"]["pil100k_image"]["data"]
        image_array = np.array(det_image).squeeze()
    return image_array


def crop_roi(image_stack, roi):
    """crop 3D image stack to 2D region of interest"""
    min_x = roi["x"]
    dx = roi["dx"]
    max_x = roi["x"] + dx

    min_y = roi["y"]
    dy = roi["dy"]
    max_y = roi["y"] + dy

    crop_image = image_stack[:, min_y:max_y, min_x:max_x]
    return crop_image


def get_calib_energies(data):
    try:
        energies = list(data["hhm_energy"])
    except:
        energies = data["data_vars"]["hhm_energy"]["data"]
    return energies


def get_p_xy(image_stack):
    image_total = np.sum(image_stack, axis=0)
    ys, xs = image_total.shape
    _x = np.arange(xs)
    _y = np.arange(ys)
    x, y = np.meshgrid(_x, _y)
    x = x.ravel()
    y = y.ravel()
    intensity = image_total.ravel()
    p_xy = np.polyfit(x, y, 1, w=intensity)
    return p_xy


def reduce_image_alt(image2d_crop, p_xy):
    ys, xs = image2d_crop.shape
    x = np.arange(xs)
    y = np.arange(ys)
    x_array, y_array = np.meshgrid(x, y)
    x_array, y_array = x_array.ravel(), y_array.ravel()
    x_array_p, _ = project_pt2line(x_array, y_array, p_xy)
    image_array = image2d_crop.ravel()

    intensity = np.zeros(x.size)

    for i in range(x.size):
        weights = 1 - np.abs(x_array_p - x[i])
        weights[weights < 0] = 0
        intensity[i] = weights @ image_array
    return x, intensity


def run_calibration(image_stack_roi, energies, n_poly=2, output_diagnostics=False):
    assert image_stack_roi.shape[0] == len(
        energies
    ), "number of calibration images must match number of calibration energies"

    p_xy = get_p_xy(image_stack_roi)
    x_pix = None
    intensity_total = None
    intensity_total_fit = None
    x_pix_centers = []
    fwhms = []

    for image in image_stack_roi:
        # apply percentile filter
        # filtered_image = percentile_threshold_filter(image, 99)
        x_pix, intensity = reduce_image_alt(image, p_xy)
        x_pix_center, fwhm, _, _, intensity_fit = fit_gaussian_with_estimation(
            x_pix, intensity
        )

        if intensity_total is None:
            intensity_total = np.zeros(intensity.size)
            intensity_total_fit = np.zeros(intensity.size)

        intensity_total += intensity
        intensity_total_fit += intensity_fit

        x_pix_centers.append(x_pix_center)
        fwhms.append(fwhm)

    p_xe = np.polyfit(x_pix_centers, energies, n_poly)

    if output_diagnostics:
        return p_xy, p_xe, x_pix, intensity_total, intensity_total_fit, x_pix_centers
    else:
        return p_xy, p_xe


def plot_calibration_diagnostics(
    image_total, x_pix, intensity_total, intensity_total_fit, x_pix_centers, p_xy, p_xe
):
    y_pix_centers = np.polyval(p_xy, x_pix_centers)
    ax1 = plt.subplot(221)
    ax1.imshow(image_total)
    ax1.plot(x_pix_centers, y_pix_centers, "o", c="r")

    # generate energy map with same dimensions as roi image
    energy_map = np.zeros(image_total.shape)
    for y, x in np.ndindex(energy_map.shape):
        energy_map[y, x] = pixel2energy(x, y, p_xy, p_xe)
    ax2 = plt.subplot(222)
    ax2.imshow(energy_map, cmap="gray")

    _x = np.arange(0, energy_map.shape[1], 1)
    _y = np.polyval(p_xy, _x)
    ax2.plot(_x, _y, "-", c="r")

    plt.subplot(223)
    plt.plot(np.polyval(p_xe, x_pix), intensity_total, "k.-")
    plt.plot(np.polyval(p_xe, x_pix), intensity_total_fit, "r-")

    energy_hi = np.polyval(p_xe, np.array(x_pix_centers) - np.array(fwhms) / 2)
    energy_lo = np.polyval(p_xe, np.array(x_pix_centers) + np.array(fwhms) / 2)
    fwhms_energy = energy_hi - energy_lo

    plt.subplot(224)
    plt.plot(energies, fwhms_energy, "k.-")


def pixel2energy(x, y, p_xy, p_xe):
    x_p, y_p = project_pt2line(x, y, p_xy)
    energy = np.polyval(p_xe, x_p)
    return energy


def process_von_hamos_calibration(uid, db, output_diagnostics=False):
    hdr = db[uid]
    md = hdr.start
    t = hdr.table(fill=True)
    roi = get_roi(md)
    image_stack = get_image_array(t)
    image_stack = crop_roi(image_stack, roi)
    energies = get_calib_energies(t)
    return run_calibration(image_stack, energies, output_diagnostics=output_diagnostics)


def apply_von_hamos_calibration_to_image_stack(image_stack, uid_calibration, db):
    p_xy, p_xe = process_von_hamos_calibration(uid_calibration, db)
    intensity = []
    for image in image_stack:
        x_pix, _intensity = reduce_image_alt(image, p_xy)
        intensity.append(_intensity)
    energy = np.polyval(p_xe, x_pix)
    return energy, intensity


# def process_von_hamos_scan(uid, uid_calibration, db):


# process_von_hamos_calibration(uid, db)


def plot_calibration(image_stack, roi, calib_energies, polynom_xy, polynom_xe, ax=None):
    if ax is None:
        ax = plt.axes()

    im_sum = np.sum(image_stack, axis=0)
    ax.imshow(im_sum, vmin=0, vmax=100)

    x, dx, y, dy = roi.values()
    # plot box around roi
    ax.vlines([x, x + dx], y, y + dy, colors="y")
    ax.hlines([y, y + dy], x, x + dx, colors="y")

    # plot fitted x-y line across roi
    ax.plot(
        [x, x + dx],
        [np.polyval(polynom_xy, 0) + y, np.polyval(polynom_xy, dx) + y],
        "r-",
    )

    energy_map = np.zeros(image_stack.shape[1:])
    for _y, _x in np.ndindex(energy_map.shape):
        energy_map[_y, _x] = pixel2energy(_x - x, _y - y, polynom_xy, polynom_xe)

    ax.contour(energy_map, np.flip(calib_energies), colors="orange")

    return ax

    # return energy_centers, intensity


def test_calibration(calib_image_stack_roi, calib_energies, polynom_xy, polynom_xe):
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

    fwhm_array = np.zeros(calib_image_stack_roi.shape[0])

    for i, im2d in enumerate(calib_image_stack_roi):
        energy_centers, intensity = reduce_image(
            im2d, polynom_xy, polynom_xe, plot_spectrum=False
        )
        ax1.plot(energy_centers, intensity, c="k")

        # Ecen = energy_centers[np.argmax(intensity)]
        Ecen, fwhm, I_cor, I_fit, I_fit_raw = fit_gaussian_with_estimation(
            energy_centers, intensity
        )
        fwhm_array[i] = fwhm

        ax1.axvline(calib_energies[i], c="b")
        ax1.axvline(Ecen, c="g", ls="--")
        ax1.plot(energy_centers, I_fit_raw, "r--")

    ax2.set_xlabel("Energy (eV)")
    ax2.set_ylabel("Resolution (fwhm)")
    ax2.plot(calib_energies, fwhm_array, "ko")

    plt.show()

    return ax1, ax2


### testing/attempts

# def reduce_image(image2d_crop, p_xy, p_xe, plot_spectrum=False):
#     x_grid, y_grid = np.meshgrid(np.arange(image2d_crop.shape[1]), np.arange(image2d_crop.shape[0]))
#     x_array, y_array = x_grid.ravel(), y_grid.ravel()
#     image_array = image2d_crop.ravel()
#     energy_array = pixel2energy(x_array, y_array, p_xy, p_xe)
#
#     energy_edges, energy_step = np.linspace(energy_array.min() - 1e-6, energy_array.max() + 1e-6, image2d_crop.shape[1] + 1, retstep=True)
#     energy_centers = energy_edges[:-1] + energy_step/2
#     inds = np.digitize(energy_array, energy_edges)
#     inds -= 1  # For Denis
#
#     intensity = np.zeros(energy_centers.size)
#
#     for i in range(energy_centers.size):
#         intensity[i] = np.sum(image_array[inds == i])
#
#     if plot_spectrum:
#         plt.plot(energy_centers, intensity)
#         plt.show()
#
#     return energy_centers, intensity

# def gen_count_array(image):
#     count_array = np.empty((0, 2), int)
#     x_grid, y_grid = np.meshgrid(np.arange(image.shape[1]), np.arange(image.shape[0]))
#     for x, y, val in zip(x_grid.ravel(), y_grid.ravel(), image.ravel()):
#         for _ in range(val):
#             count_array = np.vstack([count_array, np.array([x, y])])
#     return count_array


# def mahalanobis_dist_filter(image, percentile):
#     """
#     sets points with mahalanobis distance above percentile threshold to zero
#     """
#     filt_image = image.copy()
#     norm_image = image.copy() ** 3
#     norm_image = np.round(norm_image * 100 / norm_image.max()).astype(np.int64)
#     data_pts = gen_count_array(norm_image)
#     robust_cov = MinCovDet().fit(data_pts)
#     mahalanobis_dist = robust_cov.mahalanobis(data_pts)
#     mahalanobis_dist_threshold = np.percentile(mahalanobis_dist, percentile)
#     image_mean = robust_cov.location_
#
#     unique_data_pts, unique_ind = np.unique(data_pts, axis=0, return_index=True)
#     unique_dist = mahalanobis_dist[unique_ind]
#     # print(data_pts.shape, unique_data_pts.shape)
#     # print(mahalanobis_dist.shape, unique_dist.shape)
#     unique_data_pts_filt = unique_data_pts[unique_dist > mahalanobis_dist_threshold, :]
#     filt_image[unique_data_pts_filt[:, 1], unique_data_pts_filt[:, 0]] = 0
#
#     # plt.imshow(image)
#     # plt.show()
#     # plt.plot(image_mean[0], image_mean[1], 'ro')
#     # plt.imshow(filt_image)
#     # plt.show()
#
#     # for pt, dist in zip(unique_data_pts, unique_dist):
#     #     if dist > mahalanobis_dist_threshold:
#     #         filt_image[pt[1], pt[0]] = 0
#
#     return filt_image, image_mean


if __name__ == "__main__":

    def file_io(data_file, md_file):
        """load data and metadata from local files"""
        data = pd.read_json(data_file)
        try:
            with open(md_file) as metadata:
                md = json.load(metadata)[1]
        except:
            with open(md_file) as metadata:
                md = json.load(metadata)
        return data, md

    def test():
        PATH = "/home/charles/Desktop/XES_calib"
        DATA_FILE = "Cu_calibration.json"
        MD_FILE = "Cu_calibration_md.json"

        data, md = file_io(f"{PATH}/{DATA_FILE}", f"{PATH}/{MD_FILE}")

        try:
            roi = get_roi(md)
        except:
            roi = {"x": 100, "dx": 250, "y": 80, "dy": 20}

        pix_array = get_image_array(data)

        pix_roi = crop_roi(pix_array, roi)
        # # rotate to test x-y fitting
        # for i, image in enumerate(pix_roi):
        #     pix_roi[i] = rotate(image, 3, reshape=False)

        energies = get_calib_energies(data)

        p_xy, p_xe = run_calibration(pix_roi, energies, plot_calibration=False)

        ax_test = plt.axes()
        plot_calibration(pix_array, roi, energies, p_xy, p_xe, ax=ax_test)
        plt.show()

        test_calibration(pix_roi, energies, p_xy, p_xe)

    test()
