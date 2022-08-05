from lmfit import Model
import numpy as np

def gaussian(x, amp, cen, sigma, bkg):
    """1-d gaussian: gaussian(x, amp, cen, sigma, bkg)"""
    return (bkg + amp * np.exp(-(x-cen)**2 / (2*sigma**2)))



def get_normalized_gaussian_scan(db, uid, return_norm_param=False):
    t = db[uid].table()
    E = t['hhm_energy'].values
    I = t['pil100k_stats1_total'].values
#     I_smooth = hampel(I)
#     I_smooth = I
    offset = np.mean(np.hstack((I[:2], I[-2:])))
    I -= offset
    scale = I.max()
    I /= scale
    if return_norm_param:
        return E, np.array(I), scale, offset
    return E, np.array(I)


def estimate_center_and_width_of_peak(E, I):
    E_cen = E[np.argmax(np.abs(I))]
    x = np.abs(I - 0.5)
    e_low = E < E_cen
    e_high = E > E_cen
    x1 = E[e_low][np.argmin(x[e_low])]
    x2 = E[e_high][np.argmin(x[e_high])]
    fwhm = np.abs(x1 - x2)
    return E_cen, fwhm



def fit_gaussian(E, I, Ecen0, fwhm0):
    gmodel = Model(gaussian)
    result = gmodel.fit(I, x=E, amp=1, cen=Ecen0, sigma=fwhm0/2.355, bkg=0)
    Ecen = result.params['cen'].value
    fwhm = np.abs(result.params['sigma'].value * 2.355)
    I_fit = (result.best_fit - result.params['bkg'].value) / result.params['amp']
    I_cor = (I - result.params['bkg'].value) / result.params['amp']
    I_fit_raw = result.best_fit
    return Ecen, fwhm, I_cor, I_fit, I_fit_raw

def fit_gaussian_with_estimation(E, I):
    Ecen0, fwhm0 = estimate_center_and_width_of_peak(E, I)
    return fit_gaussian(E, I, Ecen0, fwhm0)


def fit_linear_surf(x, y, z, plotting=False):
    A = np.hstack((x[:, None], y[:, None], np.ones((x.size, 1))))
    c, _, _, _ = np.linalg.lstsq(A, z, rcond=-1)
    # if plotting:
    #     try:
    #         mplot3d
    #     except NameError:
    #         from mpl_toolkits import mplot3d
    #     fig = plt.figure()
    #     ax = plt.axes(projection='3d')
    #     ax.scatter3D(x, y, z, marker='.', color='k')
    #     ax.scatter3D(x, y, A @ c, marker='.', color='r')
    return c

class Nominal2ActualConverter:

    def __init__(self, x_nominal, x_actual, n_poly=2):
        self.p_n2a = np.polyfit(x_nominal, x_actual, n_poly)
        self.p_a2n = np.polyfit(x_actual, x_nominal, n_poly)
        self.x_nominal = x_nominal
        self.x_actual = x_actual
        self.n_poly = n_poly

    def nom2act(self, e_nom):
        return np.polyval(self.p_n2a, e_nom)

    def act2nom(self, e_act):
        return np.polyval(self.p_a2n, e_act)




