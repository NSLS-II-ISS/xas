import numpy as np
from . import xray
import numexpr as ne
import pandas as pd



def bin(interpolated_dataset, e0, edge_start=-30, edge_end=40, preedge_spacing=5,
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

    # Constants for converting from hwhm -> gaussian parameters
    GAUSS_SIGMA_FACTOR = 1 / (2 * (2 * np.log(2)) ** .5)

    def _generate_sampled_gauss_window(x, fwhm, x0):
        sigma = fwhm * GAUSS_SIGMA_FACTOR
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

    def _generate_convolution_bin_matrix(sample_points, data_x):
        fwhm = _compute_window_width(sample_points)
        delta_en = _compute_window_width(data_x)

        mat = _generate_sampled_gauss_window(data_x.reshape(1, -1),
                                            fwhm.reshape(-1, 1),
                                            sample_points.reshape(-1, 1))
        mat *= delta_en.reshape(1, -1)
        return mat

    def xas_energy_grid(energy_range, e0, edge_start, edge_end, preedge_spacing, xanes_spacing, exafs_k_spacing):
        energy_range_lo= np.min(energy_range)
        energy_range_hi = np.max(energy_range)

        preedge = np.arange(energy_range_lo, e0 + edge_start-1, preedge_spacing)

        before_edge = np.arange(e0+edge_start,e0 + edge_start+7, 1)

        edge = np.arange(e0+edge_start+7, e0+edge_end-7, xanes_spacing)

        after_edge = np.arange(e0 + edge_end - 7, e0 + edge_end, 0.7)


        eenergy = xray.k2e(xray.e2k(e0+edge_end, e0), e0)
        post_edge = np.array([])

        while (eenergy < energy_range_hi):
            kenergy = xray.e2k(eenergy, e0)
            kenergy += exafs_k_spacing
            eenergy = xray.k2e(kenergy, e0)
            post_edge = np.append(post_edge, eenergy)
        return  np.concatenate((preedge, before_edge, edge, after_edge, post_edge))

    interpolated_energy_grid = interpolated_dataset['energy'].values
    binned_energy_grid = xas_energy_grid(interpolated_energy_grid, e0, edge_start, edge_end,
                          preedge_spacing, xanes_spacing, exafs_k_spacing)


    convo_mat = _generate_convolution_bin_matrix(binned_energy_grid, interpolated_energy_grid)
    ret = {k: convo_mat @ v.values for k, v in interpolated_dataset.items() if k != 'energy'}
    ret['energy'] = binned_energy_grid
    binned_df = pd.DataFrame(ret)
    binned_df = binned_df.drop('timestamp', 1)

    return binned_df