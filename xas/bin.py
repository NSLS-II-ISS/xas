import numpy as np
from . import xray
import numexpr as ne
import pandas as pd
import time as ttime
# def get_transition_grid(E_step, E_range, n, ascend=True):
#     dE = (E_range*2/n - 2*E_step)/(n-1)
#     steps = E_step + np.arange(n)*dE
#     if not ascend:
#         steps = steps[::-1]
#     return np.cumsum(steps)

def get_transition_grid(dE_start, dE_end, E_range, round_up=True):
    if round_up:
        n = np.ceil(2 * E_range / (dE_start + dE_end))
    else:
        n = np.floor(2 * E_range / (dE_start + dE_end))
    delta = (E_range*2/n - 2*dE_start)/(n-1)
    steps = dE_start + np.arange(n)*delta
    # if not ascend:
    #     steps = steps[::-1]
    return np.cumsum(steps)


def xas_energy_grid(energy_range, e0, edge_start, edge_end, preedge_spacing, xanes_spacing, exafs_k_spacing,
                                 E_range_before=15, E_range_after = 20, n_before = 10, n_after = 20):


    # energy_range_lo= np.min(energy_range)
    # energy_range_hi = np.max(energy_range)
    energy_range_lo = np.min([e0 - 300, np.min(energy_range)])
    energy_range_hi = np.max([e0 + 2500, np.max(energy_range)])

    # preedge = np.arange(energy_range_lo, e0 + edge_start-1, preedge_spacing)
    preedge = np.arange(energy_range_lo, e0 + edge_start, preedge_spacing)

    # before_edge = np.arange(e0+edge_start,e0 + edge_start+7, 1)
    before_edge = preedge[-1] + get_transition_grid(preedge_spacing, xanes_spacing, E_range_before, round_up=False)

    edge = np.arange(before_edge[-1], e0+edge_end-E_range_after, xanes_spacing)

    # after_edge = np.arange(e0 + edge_end - 7, e0 + edge_end, 0.7)


    eenergy = xray.k2e(xray.e2k(e0+edge_end, e0), e0)
    post_edge = np.array([])

    while (eenergy < energy_range_hi):
        kenergy = xray.e2k(eenergy, e0)
        kenergy += exafs_k_spacing
        eenergy = xray.k2e(kenergy, e0)
        post_edge = np.append(post_edge, eenergy)

    after_edge = edge[-1] + get_transition_grid(xanes_spacing, post_edge[1] - post_edge[0], post_edge[0] - edge[-1], round_up=True)
    energy_grid = np.unique(np.concatenate((preedge, before_edge, edge, after_edge, post_edge)))
    energy_grid = energy_grid[(energy_grid >= np.min(energy_range)) & (energy_grid <= np.max(energy_range))]
    return  energy_grid


def _generate_convolution_bin_matrix(sample_points, data_x):
    fwhm = _compute_window_width(sample_points)
    delta_en = _compute_window_width(data_x)

    mat = _generate_sampled_gauss_window(data_x.reshape(1, -1),
                                        fwhm.reshape(-1, 1),
                                        sample_points.reshape(-1, 1))
    mat *= delta_en.reshape(1, -1)
    mat /= np.sum(mat, axis=1)[:, None]
    return mat


_GAUSS_SIGMA_FACTOR = 1 / (2 * (2 * np.log(2)) ** .5)

def _generate_sampled_gauss_window(x, fwhm, x0):
    sigma = fwhm * _GAUSS_SIGMA_FACTOR
    a = 1 / (sigma * (2 * np.pi) ** .5)
    data_y = ne.evaluate('a * exp(-.5 * ((x - x0) / sigma) ** 2)')
    # data_y = np.exp(-.5 * ((x - x0) / sigma) ** 2)
    # data_y /= np.sum(data_y)
    return data_y

def _compute_window_width(sample_points):
    '''Given smaple points compute windows via approx 1D voronoi

    Parameters
    ----------
    sample_points : array
        Assumed to be monotonic

    Returns
    -------
    windows : array
        Average of distances to neighbors
    '''
    d = np.diff(sample_points)
    fw = (d[1:] + d[:-1]) / 2
    return np.concatenate((fw[0:1], fw, fw[-1:]))




def bin(interpolated_dataset, e0, edge_start=-30, edge_end=50, preedge_spacing=5,
                        xanes_spacing= -1, exafs_k_spacing = 0.04, skip_binning=False ):
    if skip_binning:
        binned_df = interpolated_dataset
        col = binned_df.pop("energy")
        n = len(binned_df.columns)
        binned_df.insert(n, col.name, col)
        binned_df = binned_df.sort_values('energy')
    else:
        print(f'({ttime.ctime()}) Binning the data: BEGIN')
        if  xanes_spacing==-1:
            if e0 < 14000:
                xanes_spacing = 0.2
            elif e0 >= 14000 and e0 < 21000:
                xanes_spacing = 0.3
            elif e0 >= 21000:
                xanes_spacing = 0.4
            else:
                xanes_spacing = 0.3

        interpolated_energy_grid = interpolated_dataset['energy'].values
        binned_energy_grid = xas_energy_grid(interpolated_energy_grid, e0, edge_start, edge_end,
                              preedge_spacing, xanes_spacing, exafs_k_spacing)

        convo_mat = _generate_convolution_bin_matrix(binned_energy_grid, interpolated_energy_grid)
        ret = {'energy' : binned_energy_grid}
        for k, v in interpolated_dataset.items():
            if k != 'energy':
                data_array = v.values
                if len(data_array[0].shape) == 0:
                   ret[k] =  convo_mat @ data_array
                else:
                   data_ndarray = np.array([i for i in data_array], dtype = np.float64)
                   data_conv = np.tensordot(convo_mat, data_ndarray, axes=(1, 0))
                   ret[k] =  [i for i in data_conv]

        binned_df = pd.DataFrame(ret)
    binned_df = binned_df.drop('timestamp', 1)
    print(f'({ttime.ctime()}) Binning the data: DONE')
    return binned_df


# def bin_pilatus_images(interpolated_dataset, db, uid):
#     pass
# # make handler return images
# work from there



#
# fname = '/nsls2/xf08id/users/2020/3/300001/Cu_sine_1s_ud_10x 0001.raw'
# N = 10
# output_fname = '/nsls2/xf08id/Sandbox/ISS_beamline_paper/Cu_sine_1s_ud_x10.dat'
# T = 1
# T_offset = -0.1
#
# columns = ['timestamp',  'i0',  'it',  'ir',  'iff',  'aux1',  'aux2',  'aux3',  'aux4',  'energy',  'mu_t',  'mu_f',  'mu_r']
#
# data = np.genfromtxt(fname)
# df = pd.DataFrame(data, columns=columns)
# df.sort_values('timestamp', inplace=True)
# df['timestamp'] -= df['timestamp'].min()
# # T = df['timestamp'].max() / N / 2
#
#
# mus = []
#
# plt.figure(1)
# plt.clf()
#
# plt.figure(2)
# plt.clf()
#
# for i in range(N*2):
#     print(i)
#     df_loc = df[((df['timestamp']-T_offset) >= i*T) &
#                 ((df['timestamp']-T_offset) < (i+1)*T)]
#
#     plt.figure(1)
#     plt.plot(df_loc['timestamp'], df_loc['energy'])
#
#     df_int = bin(df_loc, 8979)
#     _e = df_int['energy'].values
#     _mu = np.log(df_int['it'] / df_int['ir']).values
#     plt.figure(2)
#     plt.plot(_e, _mu)
#     if i == 0:
#         energy = _e.copy()
#     _mu = np.interp(energy, _e, _mu)
#     _mu -= np.mean(_mu[energy<8900])
#     _mu /= np.mean(_mu[energy>9800])
#     mus.append(_mu)
#
# mus = np.array(mus).T
#
# ddd = np.hstack((energy[:, None], mus))
# np.savetxt(output_fname, ddd, header=('energy ' + 'd u '*N))
#
# #
# #
# # df_down_1 =
# # df_up_1 =
# #
# # df_down_5 = df[(df['timestamp'] >= 8*T) & (df['timestamp'] < 9*T)]
# # df_up_5 = df[(df['timestamp'] >= 9*T) & (df['timestamp'] < 10*T)]
# #
####


# SCRATCH

# fnames = ['/nsls2/xf08id/data/2021/03/09/en_a9f590ab',
#           '/nsls2/xf08id/data/2021/03/09/en_7377b91b',
#           '/nsls2/xf08id/data/2021/03/09/en_c6436896',
#           '/nsls2/xf08id/data/2021/03/04/en_ae93dac5',
#           '/nsls2/xf08id/data/2021/03/04/en_ae2b1fb6',
#           '/nsls2/xf08id/data/2021/03/04/en_79d333e4']

# fnames = ['/mnt/xf08ida-ioc1/test_5000',
#           '/mnt/xf08ida-ioc1/test_50000',
#           '/mnt/xf08ida-ioc1/test_89500',
#           '/mnt/xf08ida-ioc1/test_150000',
#           '/mnt/xf08ida-ioc1/test_200000']
#
# for f in fnames:
#     _d = np.genfromtxt(f)
#     print(1/np.median(np.diff(_d[:, 1]*1e-9)))



###


