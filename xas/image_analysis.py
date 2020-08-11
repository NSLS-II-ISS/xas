def get_mus(uid):
    data = db[uid].table()
    x = data['giantxy_x']
    y = data['giantxy_y']
    mut = np.log(data['apb_ave_ch1_mean']/data['apb_ave_ch2_mean'])
    muf = data['apb_ave_ch4_mean']/data['apb_ave_ch1_mean']
    return x,y, mut, muf


def analyze_spiral_scan(uid, conc, ax_t=None, ax_f=None):
    x, y, mut, muf = get_mus(uid)

    if conc>5:
        x_max, y_max = _analyze_measurement(x, y, mut)
    else:
        x_max, y_max = _analyze_measurement(x, y, muf)

    if ax_t:
        plot_xyz(x, y, mut, x_max, y_max, ax)

    if ax_f:
        plot_xyz(x, y, muf, x_max, y_max, ax)

    return x_max, y_max



def _analyze_measurement(x, y, z):
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

    x_max = com(x, (z - z.min()) ** 2, xy_mask_recen) # square it to bring the center closer to the maximum
    y_max = com(y, (z - z.min()) ** 2, xy_mask_recen) # square it to bring the center closer to the maximum
    return x_max, y_max



def com(a_orig, w_orig, mask=None):
    a = a_orig.copy()
    w = w_orig.copy()
    if mask is not None:
        a = a[mask]
        w = w[mask]
    return np.sum(a * w)/np.sum(w)




def plot_xyz(x, y, z, x_max, y_max, ax, r1=5, r2=(13.4/2-1)):

    ax.tricontourf(x, y, z, 50)
    ax.plot(x_max, y_max, 'mx', ms=25, markeredgewidth=5)
