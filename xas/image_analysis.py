import matplotlib.pyplot as plt
import numpy as np

from scipy.optimize import curve_fit
from xas.math import gauss
import math

import matplotlib.path as mpltPath
from .fitting import Nominal2ActualConverter

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
from sklearn.covariance import MinCovDet
import time as ttime

def find_points_moving_with_stage(stage_xs, stage_ys, imgs_orig, pix_per_mm=15, pix_per_mm_tol=7, assume_ring_scan=True, plotting=False):

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

    pix_xy1 = np.array([[], []]).T
    pix_xy2 = np.array([[], []]).T
    stage_xy = np.array([[], []]).T

    n_imgs = len(imgs_orig)
    idx_pairs = [(i, i+1) for i in range(n_imgs-1)]
    if assume_ring_scan: idx_pairs.append((n_imgs-1, 0))

    # plt.figure(1, clear=True)
    # pix_dists = []
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
        pix_dxy_all_norm = pix_dxy_all / np.sqrt(stage_dx**2 + stage_dy**2)
        pix_dist_norm = np.sqrt(pix_dxy_all_norm[:, 0]**2 + pix_dxy_all_norm[:, 1]**2)
        pix_dist_mask = np.abs(pix_dist_norm - pix_per_mm) <= pix_per_mm_tol

        robust_cov = MinCovDet().fit(pix_dxy_all_norm[pix_dist_mask])
        mahal_dist = robust_cov.mahalanobis(pix_dxy_all_norm[pix_dist_mask])
        inlier_mask = mahal_dist <= np.percentile(mahal_dist, 50)

        _pix_xy1 = pix_xy1_all[pix_dist_mask][inlier_mask]
        _pix_xy2 = pix_xy2_all[pix_dist_mask][inlier_mask]

        pix_xy1 = np.vstack((pix_xy1, _pix_xy1))
        pix_xy2 = np.vstack((pix_xy2, _pix_xy2))
        stage_xy = np.vstack((stage_xy, np.tile(np.array([stage_dx, stage_dy]), (np.sum(inlier_mask), 1))))

    if plotting:
        plt.figure()
        plt.subplot(121)
        plt.imshow(imgs[idx1], cmap='gray')
        for _px, _py in _pix_xy1:
            plt.plot(_px, _py, '.')

        plt.subplot(122)
        plt.imshow(imgs[idx2], cmap='gray')
        for _px, _py in _pix_xy2:
            plt.plot(_px, _py, '.')

    return {'pix_xy1':      pix_xy1,
            'pix_xy2':      pix_xy2,
            'stage_xy':     stage_xy,
            'timestamp':    ttime.time()}


class CameraCalibration:

    def __init__(self, pix_xy1, pix_xy2, stage_xy, npoly=2):
        self.pix_xy1 = np.array(pix_xy1)
        self.pix_xy2 = np.array(pix_xy2)
        self.stage_xy = np.array(stage_xy)
        self.npoly = npoly
        self.generate_calibration()

    @property
    def info(self):
        return f'points: {self.pix_xy1.shape[0]}, deg: {self.npoly}'

    def update_npoly(self, npoly):
        self.npoly = npoly
        self.generate_calibration()

    def form_pix_subbasis(self, pix_xy, include_offset=False):
        if pix_xy.ndim == 1:
            m = 1
            x, y = pix_xy
        else:
            m = pix_xy.shape[0]
            x, y = pix_xy.T
        basis_list = []
        if include_offset:
            basis_list += [np.ones((m, 1))]
        basis_list += [np.vstack([x ** i * y ** (j - i) for i in range(j + 1)]).T for j in range(1, self.npoly + 1)]
        return np.hstack(basis_list)

    def form_basis_pix_pix(self, pix_xy1, pix_xy2):
        basis1 = self.form_pix_subbasis(pix_xy1)
        basis2 = self.form_pix_subbasis(pix_xy2, include_offset=False)
        return np.hstack((basis1, basis2))

    def form_basis_pix_motor(self, pix_xy, stage_xy):
        if stage_xy.ndim == 1:
            stage_xy = stage_xy[None, :]
        basis = self.form_pix_subbasis(pix_xy)
        return np.hstack((basis, stage_xy))

    def generate_calibration(self):
        basis_pix_motor = self.form_basis_pix_motor(self.pix_xy1, self.stage_xy)
        target_pix = self.pix_xy2

        basis_pix_pix = self.form_basis_pix_pix(self.pix_xy1, self.pix_xy2)
        target_motor = self.stage_xy

        A_pix_motor_2_pixT, _, _, _ = np.linalg.lstsq(basis_pix_motor, target_pix, rcond=-1)
        A_pix_pix_2_motorT, _, _, _ = np.linalg.lstsq(basis_pix_pix, target_motor, rcond=-1)

        self.A_pix_motor_2_pix = A_pix_motor_2_pixT.T
        self.A_pix_pix_2_motor = A_pix_pix_2_motorT.T

    def compute_new_pixel(self, old_xy, stage_xy):
        old_xy = np.array(old_xy)
        stage_xy = np.array(stage_xy)
        basis_pix_motor = self.form_basis_pix_motor(old_xy, stage_xy)
        return basis_pix_motor @ self.A_pix_motor_2_pix.T

    def compute_stage_motion(self, old_xy, new_xy):
        old_xy = np.array(old_xy)
        new_xy = np.array(new_xy)
        basis_pix_pix = self.form_basis_pix_pix(old_xy, new_xy)
        return  basis_pix_pix @ self.A_pix_pix_2_motor.T

    def check_quality(self, num=None):
        pix_xy2_pred = self.compute_new_pixel(self.pix_xy1, self.stage_xy)
        stage_xy_pred = self.compute_stage_motion(self.pix_xy1, self.pix_xy2)

        plt.figure(num=num, clear=True)
        plt.subplot(321)
        plt.plot(self.pix_xy2[:, 0], self.pix_xy2[:, 1], 'k.')
        plt.plot(pix_xy2_pred[:, 0], pix_xy2_pred[:, 1], 'r.')
        plt.axis('square')

        plt.subplot(322)
        plt.plot(self.stage_xy[:, 0], self.stage_xy[:, 1], 'k.')
        plt.axis('square')
        plt.plot(stage_xy_pred[:, 0], stage_xy_pred[:, 1], 'r.')

        plt.subplot(323)
        plt.plot(pix_xy2_pred[:, 0] - self.pix_xy2[:, 0], pix_xy2_pred[:, 1] - self.pix_xy2[:, 1], 'k.')
        plt.axis('square')

        plt.subplot(324)
        plt.plot(stage_xy_pred[:, 0] - self.stage_xy[:, 0], stage_xy_pred[:, 1] - self.stage_xy[:, 1], 'k.')
        plt.axis('square')

        pix_dist_err = np.sqrt(np.sum((pix_xy2_pred - self.pix_xy2) ** 2, axis=1))
        pix_dist_err_95 = np.percentile(pix_dist_err, 95)
        plt.subplot(325)
        cts, _, _ = plt.hist(pix_dist_err)
        plt.vlines(pix_dist_err_95, cts.min(), cts.max(), colors='r')
        plt.title(f'95% error level = {pix_dist_err_95:.2f}')

        stage_dist_err = np.sqrt(np.sum((stage_xy_pred - self.stage_xy) ** 2, axis=1))
        stage_dist_err_95 = np.percentile(stage_dist_err, 95)
        plt.subplot(326)
        cts, _, _ = plt.hist(stage_dist_err)
        plt.vlines(stage_dist_err_95, cts.min(), cts.max(), colors='r')
        plt.title(f'95% error level = {stage_dist_err_95:.2f}')



class CameraCalibrationFF:

    def __init__(self, pix_xy1, pix_xy2, stage_xy, npoly=2, plotting=False, fignum=1):
        self.pix_xy1 = np.array(pix_xy1)
        self.pix_xy2 = np.array(pix_xy2)
        self.stage_xy = np.array(stage_xy)
        self.npoly = npoly
        self.generate_calibration(plotting=plotting, fignum=fignum)

    @property
    def info(self):
        return f'points: {self.pix_xy1.shape[0]}, deg: {self.npoly}'

    @property
    def calibration_data_dict(self):
        return {'pix_xy1':  self.pix_xy1,
                'pix_xy2':  self.pix_xy2,
                'stage_xy': self.stage_xy,
                'npoly' :   self.npoly}

    def update_npoly(self, npoly):
        self.npoly = npoly
        self.generate_calibration()

    def form_pix_subbasis(self, pix_xy, include_offset=False):
        if pix_xy.ndim == 1:
            m = 1
            x, y = pix_xy
        else:
            m = pix_xy.shape[0]
            x, y = pix_xy.T
        basis_list = []
        if include_offset:
            basis_list += [np.ones((m, 1))]
        basis_list += [np.vstack([x ** i * y ** (j - i) for i in range(j + 1)]).T for j in
                       range(1, self.npoly + 1)]
        return np.hstack(basis_list)

    # def form_basis_pix_pix(self, pix_xy1, pix_xy2):
    #     basis1 = self.form_pix_subbasis(pix_xy1)
    #     basis2 = self.form_pix_subbasis(pix_xy2, include_offset=False)
    #     return np.hstack((basis1, basis2))
    #
    # def form_basis_pix_motor(self, pix_xy, stage_xy):
    #     if stage_xy.ndim == 1:
    #         stage_xy = stage_xy[None, :]
    #     basis = self.form_pix_subbasis(pix_xy)
    #     return np.hstack((basis, stage_xy))

    def generate_calibration(self, plotting=False, fignum=1):

        dxy = self.pix_xy2 - self.pix_xy1
        dxy = dxy / np.sum(self.stage_xy, axis=1)[:, None]
        xy = (self.pix_xy2 + self.pix_xy1) / 2

        basis_mot2dpix = self.form_pix_subbasis(xy, include_offset=True)
        dxy_pred = np.zeros(dxy.shape)

        move_mask = self.stage_xy[:, 1] == 0  # X motion
        A_dpixdmotx, _, _, _ = np.linalg.lstsq(basis_mot2dpix[move_mask, :], dxy[move_mask, :], rcond=-1)
        dxy_pred[move_mask, :] = (basis_mot2dpix @ A_dpixdmotx)[move_mask, :]

        move_mask = self.stage_xy[:, 0] == 0  # Y motion
        A_dpixdmoty, _, _, _ = np.linalg.lstsq(basis_mot2dpix[move_mask, :], dxy[move_mask, :], rcond=-1)
        dxy_pred[move_mask, :] = (basis_mot2dpix @ A_dpixdmoty)[move_mask, :]

        self.A_dpixdmotx = A_dpixdmotx
        self.A_dpixdmoty = A_dpixdmoty

        if plotting:
            plt.figure(fignum, clear=True)

            plt.subplot(211)
            move_mask = self.stage_xy[:, 1] == 0 # X motion
            plt.plot(xy[move_mask, 0], xy[move_mask, 1], 'k.')

            for i in np.where(move_mask)[0]:
                plt.arrow(xy[i, 0] - dxy[i, 0] / 2,
                          xy[i, 1] - dxy[i, 1] / 2,
                          dxy[i, 0], dxy[i, 1], color='k', head_width=5)
                plt.arrow(xy[i, 0] - dxy_pred[i, 0] / 2,
                          xy[i, 1] - dxy_pred[i, 1] / 2,
                          dxy_pred[i, 0], dxy_pred[i, 1], color='r', head_width=5)



            plt.subplot(212)
            move_mask = self.stage_xy[:, 0] == 0 # Y motion
            plt.plot(xy[move_mask, 0], xy[move_mask, 1], 'k.')

            for i in np.where(move_mask)[0]:
                plt.arrow(xy[i, 0] - dxy[i, 0] / 2,
                          xy[i, 1] - dxy[i, 1] / 2,
                          dxy[i, 0], dxy[i, 1], color='k', head_width=5)
                plt.arrow(xy[i, 0] - dxy_pred[i, 0] / 2,
                          xy[i, 1] - dxy_pred[i, 1] / 2,
                          dxy_pred[i, 0], dxy_pred[i, 1], color='r', head_width=5)

            plt.figure(fignum + 10, clear=True)

            move_mask = self.stage_xy[:, 1] == 0 # X motion

            plt.subplot(221)
            plt.plot(xy[move_mask, 0], dxy[move_mask, 0], 'g.')
            plt.plot(xy[move_mask, 0], dxy_pred[move_mask, 0], 'k.')

            plt.plot(xy[move_mask, 0], dxy[move_mask, 1], 'b.')
            plt.plot(xy[move_mask, 0], dxy_pred[move_mask, 1], 'k.')

            plt.subplot(222)
            plt.plot(xy[move_mask, 1], dxy[move_mask, 0], 'g.')
            plt.plot(xy[move_mask, 1], dxy_pred[move_mask, 0], 'k.')

            plt.plot(xy[move_mask, 1], dxy[move_mask, 1], 'b.')
            plt.plot(xy[move_mask, 1], dxy_pred[move_mask, 1], 'k.')

            move_mask = self.stage_xy[:, 0] == 0 # Y motion

            plt.subplot(223)
            plt.plot(xy[move_mask, 0], dxy[move_mask, 0], 'm.')
            plt.plot(xy[move_mask, 0], dxy_pred[move_mask, 0], 'k.')

            plt.plot(xy[move_mask, 0], dxy[move_mask, 1], 'c.')
            plt.plot(xy[move_mask, 0], dxy_pred[move_mask, 1], 'k.')

            plt.subplot(224)
            plt.plot(xy[move_mask, 1], dxy[move_mask, 0], 'm.')
            plt.plot(xy[move_mask, 1], dxy_pred[move_mask, 0], 'k.')

            plt.plot(xy[move_mask, 1], dxy[move_mask, 1], 'c.')
            plt.plot(xy[move_mask, 1], dxy_pred[move_mask, 1], 'k.')

    def C_dpixdmot(self, xy):
        basis = self.form_pix_subbasis(xy, include_offset=True)
        c_motx = basis @ self.A_dpixdmotx
        c_moty = basis @ self.A_dpixdmoty
        return np.vstack((c_motx, c_moty)).T

    def C_dmotdpix(self, xy):
        return np.linalg.inv(self.C_dpixdmot(xy))

    def _compute_new_pixel(self, old_xy, stage_xy, step_size=1):
        new_xy = old_xy.copy()
        n_steps = int(np.ceil(np.linalg.norm(stage_xy) / step_size))
        dmot = stage_xy / n_steps
        for i in range(n_steps):
            new_xy += self.C_dpixdmot(new_xy) @ dmot
        return new_xy

    def compute_new_pixel(self, old_xy, stage_xy, step_size=1):
        old_xy = np.array(old_xy, dtype=np.float64)
        stage_xy = np.array(stage_xy, dtype=np.float64)
        if old_xy.ndim == 1:
            old_xy = old_xy[None, :]
        npts = old_xy.shape[0]

        if stage_xy.ndim == 1:
            stage_xy = stage_xy[None, :]

        if stage_xy.shape[0] == 1:
            stage_xy = np.tile(stage_xy, (npts, 1))
        else:
            assert stage_xy.shape[0] == npts, 'The number of stage shifts should be either 1 or match the number of pixel points'

        new_pixels = []
        for i in range(npts):
            new_pixel = self._compute_new_pixel(old_xy[i, :], stage_xy[i, :], step_size=step_size)
            new_pixels.append(new_pixel)
        return np.array(new_pixels)

    def _compute_stage_motion(self, old_xy, new_xy, step_size=10):
        cur_xy = old_xy.copy()
        n_steps = int(np.ceil(np.linalg.norm(new_xy - old_xy) / step_size))
        dxy = (new_xy - old_xy) / n_steps
        stage_motion = np.array([0, 0], dtype=np.float64)
        for i in range(n_steps):
            stage_motion += self.C_dmotdpix(cur_xy) @ dxy
            cur_xy += dxy
        return stage_motion

    def compute_stage_motion(self, old_xy, new_xy, step_size=10):
        old_xy = np.array(old_xy, dtype=np.float64)
        new_xy = np.array(new_xy, dtype=np.float64)
        if old_xy.ndim == 1:
            old_xy = old_xy[None, :]
        npts = old_xy.shape[0]
        if new_xy.ndim == 1:
            new_xy = new_xy[None, :]

        if new_xy.shape[0] == 1:
            new_xy = np.tile(new_xy, (npts, 1))
        else:
            assert new_xy.shape[0] == npts, 'The number of new and old points should match!'

        stage_motions = []
        for i in range(npts):
            stage_motion = self._compute_stage_motion(old_xy[i, :], new_xy[i, :], step_size=step_size)
            stage_motions.append(stage_motion)
        return np.array(stage_motions)


    def check_quality(self, fignum=None):
        pix_xy2_pred = self.compute_new_pixel(self.pix_xy1, self.stage_xy)
        stage_xy_pred = self.compute_stage_motion(self.pix_xy1, self.pix_xy2)

        plt.figure(num=fignum, clear=True)
        plt.subplot(321)
        plt.plot(self.pix_xy2[:, 0], self.pix_xy2[:, 1], 'k.')
        plt.plot(pix_xy2_pred[:, 0], pix_xy2_pred[:, 1], 'r.')
        plt.axis('square')

        plt.subplot(322)
        plt.plot(self.stage_xy[:, 0], self.stage_xy[:, 1], 'k.')
        plt.axis('square')
        plt.plot(stage_xy_pred[:, 0], stage_xy_pred[:, 1], 'r.')

        plt.subplot(323)
        plt.plot(pix_xy2_pred[:, 0] - self.pix_xy2[:, 0], pix_xy2_pred[:, 1] - self.pix_xy2[:, 1], 'k.')
        plt.axis('square')

        plt.subplot(324)
        plt.plot(stage_xy_pred[:, 0] - self.stage_xy[:, 0], stage_xy_pred[:, 1] - self.stage_xy[:, 1], 'k.')
        plt.axis('square')

        pix_dist_err = np.sqrt(np.sum((pix_xy2_pred - self.pix_xy2) ** 2, axis=1))
        pix_dist_err_95 = np.percentile(pix_dist_err, 95)
        plt.subplot(325)
        cts, _, _ = plt.hist(pix_dist_err)
        plt.vlines(pix_dist_err_95, cts.min(), cts.max(), colors='r')
        plt.title(f'95% error level = {pix_dist_err_95:.2f}')

        stage_dist_err = np.sqrt(np.sum((stage_xy_pred - self.stage_xy) ** 2, axis=1))
        stage_dist_err_95 = np.percentile(stage_dist_err, 95)
        plt.subplot(326)
        cts, _, _ = plt.hist(stage_dist_err)
        plt.vlines(stage_dist_err_95, cts.min(), cts.max(), colors='r')
        plt.title(f'95% error level = {stage_dist_err_95:.2f}')

    def compute_grid_lines(self, nlines, xmax=1292, ymax=964, stage_step=-5):
        grid_lines = []

        x = np.hstack([np.linspace(1, xmax - 1, nlines), np.ones(nlines)])
        y = np.hstack([np.ones(nlines), np.linspace(1, ymax - 1, nlines)])
        s = np.vstack([np.array([[0, stage_step]] * nlines), np.array([[stage_step, 0]] * nlines)])

        for _x, _y, _s in zip(x, y, s):
            __x = _x
            __y = _y
            grid_line = [[__x, __y]]
            while (__x > 0) and (__x < xmax) and (__y > 0) and (__y < ymax):
                __x, __y = self.compute_new_pixel([__x, __y], _s).squeeze()
                grid_line.append([__x, __y])
            grid_lines.append(grid_line)
        return grid_lines


def create_mask_for_roi(image, roi):
    xmax, ymax = image.shape

    xgrid, ygrid = np.meshgrid(np.arange(0, xmax), np.arange(0, ymax))
    xgrid = xgrid.ravel()
    ygrid = ygrid.ravel()

    polygon = mpltPath.Path(np.array(roi)[:, [1, 0]])
    mask = polygon.contains_points([(x, y) for (x, y) in zip(xgrid, ygrid)])
    mask = mask.reshape((ymax, xmax)).T
    return mask

def create_mask_for_active_rois(image, roi_dict, active_rois=None):
    if active_rois is None: active_rois = list(roi_dict.keys())
    mask = np.zeros(image.shape, dtype=bool)
    for roi in active_rois:
        mask |= create_mask_for_roi(image, roi_dict[roi])
    return mask

def make_edge_roi_dict(roi_dict, dr: float=8.0):
    edge_roi_dict = {}
    for key, _polygon in roi_dict.items():
        polygon = np.array(_polygon, dtype=np.float64)
        edge_polygon = []
        for i, (x, y) in enumerate(polygon):
            i_prev = i - 1
            i_next = i + 1
            if i_next == polygon.shape[0]:
                i_next = 0

            v_prev = polygon[i_prev, :] - polygon[i, :]
            v_next = polygon[i_next, :] - polygon[i, :]

            v_prev /= np.linalg.norm(v_prev)
            v_next /= np.linalg.norm(v_next)

            v_out = -(v_prev + v_next)  # outward vector that points away from other two vertices
            v_out /= np.linalg.norm(v_out)

            dx, dy = np.round(dr * v_out)
            edge_polygon.append((x + dx, y + dy))
        edge_roi_dict[key] = edge_polygon
    return edge_roi_dict

def create_mask_roi_dict(image, roi_dict: dict=None):
    mask_roi_dict = {}
    for key, roi_polygon in roi_dict.items():
        mask_roi_dict[key] = create_mask_for_roi(image, roi_polygon)
    return mask_roi_dict

def create_mask_edge_dict(image, roi_dict: dict=None, mask_roi_dict: dict=None, dr: float=8.0, allow_roi_crosstalk=False):
    if mask_roi_dict is None:
        mask_roi_dict = create_mask_roi_dict(image, roi_dict)

    mask_edge_dict = {}
    edge_roi_dict = make_edge_roi_dict(roi_dict, dr=dr)

    if allow_roi_crosstalk:
        mask_all_rois = np.zeros(image.shape, dtype=bool)
    else:
        mask_all_rois = np.sum(list(mask_roi_dict.values()), axis=0).astype(bool)

    for key, roi_polygon in edge_roi_dict.items():
        mask_edge_dict[key] = create_mask_for_roi(image, roi_polygon) & ~mask_all_rois

    return mask_edge_dict



def estimate_image_background(image, mask_roi, mask_bkg, gauss_fwhm=20, plotting=False):

    a, b = image.shape
    x, y = np.meshgrid(np.arange(b), np.arange(a))
    x, y = x.ravel(), y.ravel()
    mask_roi_arr = mask_roi.ravel()
    mask_bkg_arr = mask_bkg.ravel()
    image_arr = image.ravel()

    tresh_lo = np.percentile(image_arr[mask_bkg_arr], 2)
    tresh_hi = np.percentile(image_arr[mask_bkg_arr], 98)
    zinger_mask = (image_arr >= tresh_lo) & (image_arr <= tresh_hi)

    mask_bkg_zinger_arr = mask_bkg_arr & zinger_mask

    dx = x[mask_roi_arr, None] - x[None, mask_bkg_zinger_arr]
    dy = y[mask_roi_arr, None] - y[None, mask_bkg_zinger_arr]
    # dx = x[:, None] - x[None, :]
    # dy = y[:, None] - y[None, :]
    dr = np.sqrt(dx ** 2 + dy ** 2)
    sigma = gauss_fwhm / 2.355
    A = np.exp(-0.5 * (dr / sigma) ** 2)
    A /= np.sum(A, axis=1)[:, None]
    if plotting:
        bkg_arr = np.zeros(image_arr.shape)
        bkg_arr[mask_roi_arr] = A @ image_arr[mask_bkg_zinger_arr]
        bkg = bkg_arr.reshape(a, b)
        cor = image.copy()
        cor[mask_roi] = bkg[mask_roi]
        plt.figure()
        plt.imshow(cor, vmin=np.percentile(cor, 5), vmax=np.percentile(cor, 95))

    return np.sum(A @ image_arr[mask_bkg_zinger_arr])

def get_total_counts_from_roi_mask(image, roi_mask):
    return np.sum(image[roi_mask])



def reduce_johann_images(interpolated_df, hdr, detector_key='pil100k2', image_key='image'):
    _roi_dict = hdr.start['detectors']['Pilatus 100k New']['config']['roi_polygon']
    trajectory_dict = hdr.start['spectrometer_relative_trajectory']

    crystals = list(trajectory_dict.keys())
    roi_dict = {}
    for crystal in crystals:
        roi_dict[crystal] = _roi_dict[crystal]

    det_image_key = f'{detector_key}_{image_key}'
    roi_mask_dict = create_mask_roi_dict(interpolated_df[det_image_key][0], roi_dict)

    for crystal in crystals:
        roi_key = f'{detector_key}_{crystal}'
        interpolated_df[roi_key] = interpolated_df[det_image_key].apply(lambda x: get_total_counts_from_roi_mask(x, roi_mask_dict[crystal]))

    return interpolated_df



def plot_polygon_roi(roi_dict, key='main', ax=None, **plot_kwargs):
    roi = roi_dict[key]
    if ax is None:
        ax = plt.gca()
    for i in range(len(roi)):
        if i < (len(roi) - 1):
            j = i + 1
        else:
            j = 0
        point_i = roi[i]
        point_j = roi[j]
        ax.plot([point_i[0], point_j[0]], [point_i[1], point_j[1]], **plot_kwargs)

def show_polygon_roi_for_pil100k_image(image, roi_dict): #, detector_key='pil100k2', image_key='image'):
    # det_image_key = f'{detector_key}_{image_key}'
    # total_image = np.sum(np.array([v for v in df[det_image_key].values]), axis=0)

    # trajectory_dict = hdr.start['spectrometer_relative_trajectory']
    # crystals = list(trajectory_dict.keys())

    # _roi_dict = hdr.start['detectors']['Pilatus 100k New']['config']['roi_polygon']

    plt.figure()
    plt.imshow(image, vmin=np.percentile(image, 5), vmax=np.percentile(image, 95))

    for crystal in roi_dict:
        plot_polygon_roi(roi_dict, key=crystal, color=_crystal_info_dict[crystal]['color'])








