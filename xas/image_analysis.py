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


def determine_beam_position_from_fb_image(image, line = 420, center_point = 655, n_lines = 1, truncate_data=True):

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
            return None
    else:
        return None




# -458 rel_spiral_square a467cce4-30a3-4796-9b87-f176982b15d1
# -459 rel_spiral_square efcb65d4-6852-4941-9a2d-125d17ae1a31
# -460 rel_spiral_square d23e1168-f97e-4f58-9a40-2160a926afd8
# -461 rel_spiral_square 283a1052-f0e6-4b3f-8a0d-64551375d03c