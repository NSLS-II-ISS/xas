import numpy as np
from . import xray
import numexpr as ne
import pandas as pd

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


    energy_range_lo= np.min(energy_range)
    energy_range_hi = np.max(energy_range)

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
    return  np.unique(np.concatenate((preedge, before_edge, edge, after_edge, post_edge)))


def _generate_convolution_bin_matrix(sample_points, data_x):
    fwhm = _compute_window_width(sample_points)
    delta_en = _compute_window_width(data_x)

    mat = _generate_sampled_gauss_window(data_x.reshape(1, -1),
                                        fwhm.reshape(-1, 1),
                                        sample_points.reshape(-1, 1))
    mat *= delta_en.reshape(1, -1)
    return mat


_GAUSS_SIGMA_FACTOR = 1 / (2 * (2 * np.log(2)) ** .5)

def _generate_sampled_gauss_window(x, fwhm, x0):
    sigma = fwhm * _GAUSS_SIGMA_FACTOR
    a = 1 / (sigma * (2 * np.pi) ** .5)
    data_y = ne.evaluate('a * exp(-.5 * ((x - x0) / sigma) ** 2)')
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
                        xanes_spacing= -1, exafs_k_spacing = 0.04 ):

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
    ret = {k: convo_mat @ v.values for k, v in interpolated_dataset.items() if k != 'energy'}
    ret['energy'] = binned_energy_grid
    binned_df = pd.DataFrame(ret)
    binned_df = binned_df.drop('timestamp', 1)

    return binned_df