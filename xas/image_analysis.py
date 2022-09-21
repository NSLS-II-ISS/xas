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
from sklearn.covariance import MinCovDet

def generate_sample_camera_calibration(stage_xs, stage_ys, imgs_orig, pix_per_mm=15, pix_per_mm_tol=5, assume_ring_scan=True):
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

    pix_xy1 = np.array([[], []]).T
    pix_xy2 = np.array([[], []]).T
    stage_xy = np.array([[], []]).T

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


        # plt.figure()
        #
        # plt.subplot(131)
        # plt.hist(pix_dist_norm, 200);
        #
        # plt.subplot(132)
        # plt.imshow(imgs[idx1], cmap='gray')
        # for i in range(_pix_xy1.shape[0]):
        #     if inlier_mask[i]:
        #         plt.plot(_pix_xy1[i, 0], _pix_xy1[i, 1], 'o')
        #
        # plt.subplot(133)
        # plt.imshow(imgs[idx2], cmap = 'gray')
        # for i in range(_pix_xy2.shape[0]):
        #     if inlier_mask[i]:
        #         plt.plot(_pix_xy2[i, 0], _pix_xy2[i, 1], 'o')

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

    return A_pix_motor_2_pix, A_pix_pix_2_motor, pix_xy1, pix_xy2, stage_xy, pix_dists


# A_pix_motor_2_pix, A_pix_pix_2_motor, pix_xy1, pix_xy2, pix_dists = generate_sample_camera_calibration(dxs, dys, OUTPUT['camera_sp1'], dist_target_x=65, dist_target_y=90)

# A_pix_motor_2_pix, A_pix_pix_2_motor, pix_xy1, pix_xy2, stage_xy, pix_dists = generate_sample_camera_calibration(dxs, dys, OUTPUT['camera_sp2'], dist_target_x=78, dist_target_y=90, dist_tol=40)
# A_pix_motor_2_pix, A_pix_pix_2_motor, pix_xy1, pix_xy2, stage_xy, pix_dists = generate_sample_camera_calibration(dxs, dys, OUTPUT['camera_sp2'], pix_per_mm=15, pix_per_mm_tol=7)



class CameraCalibration:

    def __init__(self, pix_xy1, pix_xy2, stage_xy, npoly=2):
        self.pix_xy1 = pix_xy1
        self.pix_xy2 = pix_xy2
        self.stage_xy = stage_xy
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


# npoly=1
# CC = CameraCalibration(pix_xy1, pix_xy2, stage_xy, npoly=npoly)
# CC.check_quality(25)
#
# npoly=2
# CC = CameraCalibration(pix_xy1, pix_xy2, stage_xy, npoly=npoly)
# CC.check_quality(26)
#
# npoly=3
# CC = CameraCalibration(pix_xy1, pix_xy2, stage_xy, npoly=npoly)
# CC.check_quality(27)

# CC.compute_new_pixel((100, 100), (5, 0))



# idx_sel = (stage_xy[:, 0] == 5) | (stage_xy[:, 0] == -5)
# dxy_norm = ((pix_xy2 - pix_xy1) / stage_xy[:, 0][:, None])[idx_sel, :]
#
#
# plt.figure(12, clear=True)
# plt.subplot(211)
# plt.plot(dxy_norm[:, 0], dxy_norm[:, 1], 'k.')
# # plt.plot(dxy_norm[inlier_mask, 0], dxy_norm[inlier_mask, 1], 'r.')
#
#
# idx_sel = (stage_xy[:, 1] == 5) | (stage_xy[:, 1] == -5)
# dxy_norm = ((pix_xy2 - pix_xy1) / stage_xy[:, 1][:, None])[idx_sel, :]
#
# # robust_cov = MinCovDet().fit(dxy_norm)
# # mahal_dist = robust_cov.mahalanobis(dxy_norm)
# # inlier_mask = mahal_dist <= np.percentile(mahal_dist, 50)
#
# # plt.figure(12, clear=True)
# plt.subplot(212)
# plt.plot(dxy_norm[:, 0], dxy_norm[:, 1], 'k.')
# # plt.plot(dxy_norm[inlier_mask, 0], dxy_norm[inlier_mask, 1], 'r.')
#
#
# # plt.hist(np.sqrt(dxy_norm[:, 0]**2 + dxy_norm[:, 1]**2))
#
#
# plt.figure(13, clear=True)
# bla = (pix_xy2 - pix_xy1)/5
# blabla = np.sqrt(bla[:, 0]**2 + bla[:, 1]**2)
# plt.hist(blabla)
#
# pix_xy_test = np.vstack((np.linspace(100, 1200 ,101), 500*np.ones(101)))
# stage_xy_test = np.vstack((-5*np.ones(101), 0*np.ones(101)))
# basis_test = np.vstack((pix_xy_test, stage_xy_test))
#
#
# plt.figure(3, clear=True)
# plt.subplot(221)
# plt.plot(pix_xy1[:, 0] , pix_xy2[:, 0], 'k.')
# plt.plot(pix_xy_test.T[:, 0], (A_pix_motor_2_pix @ basis_test).T[:, 0], 'r.')
#
# plt.subplot(222)
# plt.plot(pix_xy1[:, 1], pix_xy2[:, 1], 'k.')
# plt.plot(pix_xy_test.T[:, 1], (A_pix_motor_2_pix @ basis_test).T[:, 1], 'r.')
#
# plt.subplot(223)
# plt.plot(stage_xy[:, 0], pix_xy2[:, 0] - pix_xy1[:, 0], 'k.')
#
# plt.subplot(224)
# plt.plot(stage_xy[:, 0], pix_xy2[:, 1] - pix_xy1[:, 1], 'k.')
#
#


# -458 rel_spiral_square a467cce4-30a3-4796-9b87-f176982b15d1
# -459 rel_spiral_square efcb65d4-6852-4941-9a2d-125d17ae1a31
# -460 rel_spiral_square d23e1168-f97e-4f58-9a40-2160a926afd8
# -461 rel_spiral_square 283a1052-f0e6-4b3f-8a0d-64551375d03c