import json
import numpy as np
import pandas as pd

import matplotlib
from matplotlib import projections
matplotlib.use('TkAgg')  # something wrong with my (Charlie's) system Qt 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from scipy.ndimage import center_of_mass, rotate
from scipy import linalg


PATH = '/home/charles/Desktop/XES_calib'
DATA_FILE = 'Cu_calibration.json'
MD_FILE = 'Cu_calibration_md.json'


def percentile_threshold(im2d, p):
   """ set values below percentile equal to 0 """
   im2d[im2d < np.percentile(im2d, p)] = 0
   return im2d


def fit_plane(X, Y, Z):
    """ fit plane to 3d data using least square regression """
    X, Y, Z = np.array(X), np.array(Y), np.array(Z)
    assert(X.size == Y.size == Z.size)
    A = np.column_stack([X, Y, np.ones(X.size)])
    C, _, _, _ = linalg.lstsq(A, Z)
    
    # equation for plane: z = a*x + b*y + c
    a, b, c = C[0], C[1], C[2]
    return a, b, c


def fit_xy(X, Y):
    X,Y = np.array(X), np.array(Y)
    slope, intercept = np.polyfit(X, Y, 1)
    return slope, intercept


def project_pt2line(x_pt, y_pt, slope, intercept):
    x_proj = (x_pt + slope*y_pt - slope*intercept) / (slope**2 + 1)
    y_proj = slope * x_proj + intercept
    return x_proj, y_proj


def file_io(data_file, md_file):
    """ load data and metadata from local files """
    data = pd.read_json(data_file)
    with open(md_file) as metadata:
        md = json.load(metadata)[1]
    return data, md


def get_roi(metadata, roi='roi1'):
    # Load ROI
    rois = metadata['detectors']['Pilatus 100k']['config']['roi']
    # print(rois)
    roi_ = rois[roi]
    return roi_


def get_image_array(data):
    det_image = data['pil100k_image']
    image_array = np.array(list(det_image)).squeeze()
    return image_array


def crop_roi(image_stack, roi):
    """ crop 3D image stack to 2D region of interest """
    min_x = roi['x']
    dx = roi['dx']
    max_x = roi['x'] + dx

    min_y = roi['y']
    dy = roi['dy']
    max_y = roi['y'] + dy

    crop_image = image_stack[:, min_y:max_y, min_x:max_x]
    return crop_image


def get_calib_energies(data):
    energies = list(data['hhm_energy'])
    return energies


def run_calibration(pix_roi, energies, plot_calibration=False):

    # assert(pix_roi.shape[0] == len(energies))

    filtered_roi = np.zeros(pix_roi.shape)

    COM = {'x': [], 'y': []}  # centers of mass
    for i, image in enumerate(pix_roi):
        # apply percentile filter
        filtered_image = percentile_threshold(image, 99)
        filtered_roi[i] = filtered_image

        y_com, x_com = center_of_mass(filtered_image)
        COM['x'].append(x_com)
        COM['y'].append(y_com)


    slope_xy, int_xy = fit_xy(COM['x'], COM['y'])

    energy_fit = np.poly1d(np.polyfit(COM['x'], energies, 2))
    # a, b, c = fit_plane(COM['x'], COM['y'], energies)
    def energy_function(x, y):
        x_p, y_p = project_pt2line(x, y, slope_xy, int_xy)
        energy = energy_fit(x_p)
        return energy

    if plot_calibration:

        # generate energy map with same dimensions as roi image
        energy_map = np.zeros(pix_roi[0].shape)
        for y, x in np.ndindex(energy_map.shape):
            # apply energy function to each point on map
            energy_map[y, x] = energy_function(x, y)

        ax1 = plt.subplot(211)
        ax1.imshow(np.sum(pix_roi, axis=0))
        for x, y in zip(COM['x'], COM['y']):
            ax1.plot(x, y, 'o', c='r')
        for x, energy in zip(COM['x'], energies):
            plt.axvline(x, c='orange')
            plt.text(x+1, 3, f'{energy}', c='orange')

        ax2 = plt.subplot(212)
        # ax2.contour(energy_map, np.flip(energies), colors='orange')
        ax2.imshow(energy_map, cmap='gray')

        # plot best fit line for x-y
        _x = np.arange(0, energy_map.shape[1], 1)
        _y = slope_xy * _x + int_xy
        ax2.plot(_x, _y, '-', c='r')
        plt.show()

        # 3D plot
        plt.clf()
        ax = plt.axes(projection='3d')
        ax.scatter(COM['x'], COM['y'], energies, c='darkblue')
        x_grid, y_grid = np.meshgrid(np.arange(0, energy_map.shape[1]), np.arange(0, energy_map.shape[0]))
        ax.plot_surface(x_grid, y_grid, energy_map, cmap='viridis', alpha=0.5)
        plt.show()

    return energy_function


def reduce_image(image2d, energy_func, plot_spectrum=False):

    calc_energies = {'point': [], 'energy': []}
    for y, x in np.ndindex(image2d.shape):
        energy = energy_func(x, y)
        calc_energies['point'].append((x,y))
        calc_energies['energy'].append(energy)

    energy_array = np.array(calc_energies['energy'])
    # print(energy_array.min(), energy_array.max())

    bin_edges, bin_size = np.linspace(energy_array.min(), energy_array.max(), image2d.shape[1], retstep=True)
    inds = np.digitize(energy_array, bin_edges)

    energy_intensities = np.zeros(bin_edges.size)

    for n, pt in enumerate(calc_energies['point']):
        xpt, ypt = pt[0], pt[1]
        intensity = image2d[ypt, xpt]
        energy_intensities[inds[n] - 1] += intensity

    if plot_spectrum:
        plt.plot(bin_edges + bin_size/2, energy_intensities)

        etest = np.arange(7995, 8085+1, 5)
        for e in etest:
            plt.axvline(e, c='r')

        plt.show()


if __name__ == '__main__':
    
    def test():
        
        data, md = file_io(f'{PATH}/{DATA_FILE}', f'{PATH}/{MD_FILE}')
        
        roi = get_roi(md)
        pix_array = get_image_array(data)

        pix_roi = crop_roi(pix_array, roi)
        # rotate to test x-y fitting
        for i, image in enumerate(pix_roi):
            pix_roi[i] = rotate(image, 3, reshape=False)
        
        energies = get_calib_energies(data)

        e_func = run_calibration(pix_roi, energies, plot_calibration=True)
        
        im2d = np.sum(pix_roi, axis=0)
        reduce_image(im2d, e_func, plot_spectrum=True)

    test()
