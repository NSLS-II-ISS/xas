import matplotlib.pyplot as plt
import numpy as np

from scipy.optimize import curve_fit
from xas.math import gauss
import math

def get_mus(db, uid):
    data = db[uid].table()
    x = data['giantxy_x']
    y = data['giantxy_y']
    mut = np.log(data['apb_ave_ch1_mean']/data['apb_ave_ch2_mean'])
    muf = data['apb_ave_ch4_mean']/data['apb_ave_ch1_mean']
    return x,y, mut, muf


def analyze_spiral_scan(db, uid, conc, ax, save_fig_path=None):
    x, y, mut, muf = get_mus(db, uid)

    if conc>5:
        x_max, y_max = _analyze_measurement(x, y, mut)
        plot_xyz(x, y, mut, x_max, y_max, ax, save_fig_path=save_fig_path, label='mut')
    else:
        x_max, y_max = _analyze_measurement(x, y, muf)
        plot_xyz(x, y, muf, x_max, y_max, ax, save_fig_path=save_fig_path, label='muf')

    return x_max, y_max



def _analyze_measurement(x, y, z, r1=5, r2=(13.4/2-1)):
    x_im_center = x.iloc[0]
    y_im_center = y.iloc[0]

    xy_mask = (np.sqrt(np.abs(x - x_im_center) ** 2 +
                       np.abs(y - y_im_center) ** 2) > r1)

    # find golder center
    x_ho_com = com(x, z.max() - z, xy_mask)
    y_ho_com = com(y, z.max() - z, xy_mask)

    # find the area corresponding to the sample position
    xy_mask_recen = (np.sqrt(np.abs(x - x_ho_com) ** 2 +
                             np.abs(y - y_ho_com) ** 2) < r2)

    # x_max = x[xy_mask_recen][np.argmax(z[xy_mask_recen])]
    # y_max = y[xy_mask_recen][np.argmax(z[xy_mask_recen])]

    x_max = com(x, (z - z.min()) ** 3, xy_mask_recen) # square it to bring the center closer to the maximum
    y_max = com(y, (z - z.min()) ** 3, xy_mask_recen) # square it to bring the center closer to the maximum

    return x_max, y_max



def com(a_orig, w_orig, mask=None):
    a = a_orig.copy()
    w = w_orig.copy()
    if mask is not None:
        a = a[mask]
        w = w[mask]
    return np.sum(a * w)/np.sum(w)




def plot_xyz(x, y, z, x_max, y_max, ax, save_fig_path=None, label=None):

    if save_fig_path:
        plt.ioff()
        fig, ax = plt.subplots(1, 1, figsize=(3, 3))

    if ax:
        ax.tricontourf(x, y, z, 50)
        ax.plot(x_max, y_max, 'mx', ms=12, markeredgewidth=2)
        if label:
            plt.title(label)
            plt.xlabel('giantxy_x')
            plt.ylabel('giantxy_y')
    if save_fig_path:
        plt.savefig(save_fig_path, dpi=200)
        plt.ion()
        plt.close(fig)


def show_spiral_result(db,uid):

    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(6,3))
    fig.set_tight_layout(True)
    analyze_spiral_scan(db, uid, 10, ax1)
    analyze_spiral_scan(db, uid, 1, ax2)

    ax1.set_xlabel('giant_x')
    ax1.set_ylabel('giant_y')
    ax1.set_title('Transmission')

    ax2.set_xlabel('giant_x')
    ax2.set_ylabel('giant_y')
    ax2.set_title('Fluorescence')

    plt.savefig(r'/home/xf08id/Desktop/spiral_scan_1.png', dpi=600)


def determine_beam_position_from_fb_image(image, line = 420, center_point = 655, n_lines = 1, truncate_data=True, should_print_diagnostics=True):
    # print('>>>>> analyzing image')
    image = image.astype(np.int16)

    sum_lines = sum(image[:, [i for i in range(int(line - math.floor(n_lines/2)),
                                               int(line + math.ceil(n_lines/2)))]].transpose())
    if len(sum_lines) > 0:
        sum_lines = sum_lines - np.mean(sum_lines[:200]) # empirically we determined that first 200 pixels are BKG
    index_max = sum_lines.argmax()
    max_value = sum_lines.max()
    min_value = sum_lines.min()
    idx_to_fit = np.where(sum_lines > max_value / 2)
    x = np.arange(960)
    # sdfsd
    if max_value >= 10 and max_value <= n_lines * 100 and ((max_value - min_value) / n_lines) > 5:
        try:
            if truncate_data:
                coeff, var_matrix = curve_fit(gauss, x[idx_to_fit], sum_lines[idx_to_fit], p0=[1, index_max, 5])
            else:
                # print('actually got to fitting')
                coeff, var_matrix = curve_fit(gauss, x, sum_lines, p0=[1, index_max, 5])
            return coeff[1]
        except:
            if should_print_diagnostics:
                print('>>>> FEEDBACK - failed - Fitting failure')
            return None
    else:
        if should_print_diagnostics:
            print('>>>> FEEDBACK - failed - image is either empty or saturated')
        return None


import cv2

def generate_sample_camera_calibration(stage_xs, stage_ys, imgs_orig, dist_tol=10, dist_target=90, assume_ring_scan=True):
# def generate_sample_camera_calibration(imgs_orig, dist_tol=30, dist_target=75,
#                                            assume_ring_scan=True):
    sift = cv2.xfeatures2d.SIFT_create()

    imgs = []
    img_keypoints = []
    img_descriptors = []
    for _img in imgs_orig:
        img = cv2.normalize(cv2.cvtColor(_img.astype(np.float32), cv2.COLOR_GRAY2BGR), None, 0, 255, cv2.NORM_MINMAX).astype('uint8')
        imgs.append(img)
        keypoints, descriptors = sift.detectAndCompute(img, None)
        img_keypoints.append(keypoints)
        img_descriptors.append(descriptors)

    pix_xy1 = []
    pix_xy2 = []
    stage_xy = []

    n_imgs = len(imgs_orig)
    idx_pairs = [(i, i+1) for i in range(n_imgs-1)]
    if assume_ring_scan: idx_pairs.append((n_imgs-1, 0))

    # plt.figure(1, clear=True)
    pix_dists = []
    bf = cv2.BFMatcher(cv2.NORM_L1, crossCheck=False)
    for idx1, idx2 in idx_pairs:
        descriptors_1 = img_descriptors[idx1]
        descriptors_2 = img_descriptors[idx2]
        keypoints_1 = img_keypoints[idx1]
        keypoints_2 = img_keypoints[idx2]

        stage_dx = stage_xs[idx2] - stage_xs[idx1]
        stage_dy = stage_ys[idx2] - stage_ys[idx1]

        matches = bf.match(descriptors_1, descriptors_2)
        pix_xy1_all = np.array([keypoints_1[v.queryIdx].pt for v in matches])
        pix_xy2_all = np.array([keypoints_2[v.trainIdx].pt for v in matches])
        pix_dxy_all = pix_xy2_all - pix_xy1_all
        pix_dist = np.sqrt(pix_dxy_all[:, 0]**2 + pix_dxy_all[:, 1]**2)
        pix_dists.append(pix_dist)
        # plt.hist(pix_dist, int(pix_dist.max()))
        # pix_idx_filtered, = np.where((pix_dist >= dist_target - dist_tol) & (pix_dist <= dist_target + dist_tol))
        pix_idx_filtered, = np.where((np.abs(pix_dist - dist_target) <=  dist_tol))

        _pix_xy1 = [keypoints_1[matches[pix_idx].queryIdx].pt for pix_idx in pix_idx_filtered]
        _pix_xy2 = [keypoints_2[matches[pix_idx].trainIdx].pt for pix_idx in pix_idx_filtered]

        pix_xy1 += _pix_xy1
        pix_xy2 += _pix_xy2

        stage_xy += [(stage_dx, stage_dy)] * pix_idx_filtered.size

        plt.figure()

        plt.subplot(131)
        plt.hist(pix_dist, 200);

        plt.subplot(132)
        plt.imshow(imgs[idx1], cmap='gray')
        for i in range(len(_pix_xy1)):
            plt.plot(_pix_xy1[i][0], _pix_xy1[i][1], 'o')

        plt.subplot(133)
        plt.imshow(imgs[idx2], cmap = 'gray')
        for i in range(len(_pix_xy2)):
            plt.plot(_pix_xy2[i][0], _pix_xy2[i][1], 'o')

    pix_xy1 = np.array(pix_xy1)
    pix_xy2 = np.array(pix_xy2)
    stage_xy = np.array(stage_xy)

    basis_pix_motor = np.hstack((pix_xy1, stage_xy))
    target_pix = pix_xy2

    basis_pix_pix = np.hstack((pix_xy1, stage_xy))
    target_motor = stage_xy

    A_pix_motor_2_pixT, _, _, _ = np.linalg.lstsq(basis_pix_motor, target_pix, rcond=-1)
    A_pix_pix_2_motorT, _, _, _ = np.linalg.lstsq(basis_pix_pix, target_motor, rcond=-1)

    A_pix_motor_2_pix = A_pix_motor_2_pixT.T
    A_pix_pix_2_motor = A_pix_pix_2_motorT.T

    plt.figure(25, clear=True)
    plt.subplot(221)
    plt.plot((basis_pix_motor @ A_pix_motor_2_pixT)[:, 0], pix_xy2[:, 0], 'k.')
    plt.axis('square')

    return A_pix_motor_2_pix, A_pix_pix_2_motor


A_pix_motor_2_pix, A_pix_pix_2_motor = generate_sample_camera_calibration(dxs, dys, OUTPUT['camera_sp1'])





# -458 rel_spiral_square a467cce4-30a3-4796-9b87-f176982b15d1
# -459 rel_spiral_square efcb65d4-6852-4941-9a2d-125d17ae1a31
# -460 rel_spiral_square d23e1168-f97e-4f58-9a40-2160a926afd8
# -461 rel_spiral_square 283a1052-f0e6-4b3f-8a0d-64551375d03c