from code import interact
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
from sklearn.covariance import MinCovDet

from fitting import fit_gaussian_with_estimation


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
    a, b, c = C
    return a, b, c


def fit_xy(X, Y):
    slope, intercept = np.polyfit(X, Y, 1)
    return slope, intercept


def project_pt2line(x_pt, y_pt, p_xy):
    slope, intercept = p_xy
    x_proj = (x_pt + slope*y_pt - slope*intercept) / (slope**2 + 1)
    y_proj = slope * x_proj + intercept
    return x_proj, y_proj


def file_io(data_file, md_file):
    """ load data and metadata from local files """
    data = pd.read_json(data_file)
    try:
        with open(md_file) as metadata:
            md = json.load(metadata)[1]
    except:
        with open(md_file) as metadata:
            md = json.load(metadata)
    return data, md


def get_roi(metadata, roi='roi1'):
    # Load ROI
    rois = metadata['detectors']['Pilatus 100k']['config']['roi']
    # print(rois)
    roi_ = rois[roi]
    return roi_


def get_image_array(data):
    try:
        det_image = data['pil100k_image']
        image_array = np.array(list(det_image)).squeeze()
    except:
        det_image = data['data_vars']['pil100k_image']['data']
        image_array = np.array(det_image).squeeze()
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
    try:
        energies = list(data['hhm_energy'])
    except:
        energies = data['data_vars']['hhm_energy']['data']
    return energies


def gen_count_array(image):
    count_array = np.empty((0, 2), int)
    x_grid, y_grid = np.meshgrid(np.arange(image.shape[1]), np.arange(image.shape[0]))
    for x, y, val in zip(x_grid.ravel(), y_grid.ravel(), image.ravel()):
        for _ in range(val):
            count_array = np.vstack([count_array, np.array([x, y])])
    return count_array
        

def mahalanobis_dist_filter(image, p):
    data_array = gen_count_array(image)
    robust_cov = MinCovDet().fit(data_array)
    mahalanobis_dist = robust_cov.mahalanobis(data_array)
    mahalanobis_dist_threshold = np.percentile(mahalanobis_dist, p)
    
    mean_x = np.mean(data_array[:, 0])
    mean_y = np.mean(data_array[:, 1])
    # plt.imshow(image)
    # plt.plot(mean_x, mean_y, 'ro')
    # plt.show()
    

    for pt, dist in zip(data_array, mahalanobis_dist):
        if dist > mahalanobis_dist_threshold:
            image[pt[1], pt[0]] = 0

    # plt.imshow(image)
    # plt.plot(mean_x, mean_y, 'ro')
    # plt.show()
    
    return image

def run_calibration(image_stack_roi, energies, n_poly=2, plot_calibration=False):

    assert image_stack_roi.shape[0] == len(energies), "number of calibration images must match number of calibration energies"

    filtered_roi_stack = np.zeros(image_stack_roi.shape)

    COM = {'x': [], 'y': []}  # centers of mass
    for i, image in enumerate(image_stack_roi):
        # apply percentile filter
        filtered_image = percentile_threshold(image, 99)
        filtered_roi_stack[i] = filtered_image

        y_com, x_com = center_of_mass(filtered_image)
        COM['x'].append(x_com)
        COM['y'].append(y_com)

    # slope_xy, int_xy = fit_xy(COM['x'], COM['y'])
    polynom_xy = np.polyfit(COM['x'], COM['y'], 1)

    # energy_fit = np.poly1d(np.polyfit(COM['x'], energies, n_poly))
    polynom_xe = np.polyfit(COM['x'], energies, n_poly)

    # a, b, c = fit_plane(COM['x'], COM['y'], energies)
    # def energy_function(x, y):
    #     x_p, y_p = project_pt2line(x, y, polynom_xy)
    #     energy = np.polyval(polynom_xe, x_p)
    #     return energy

    if plot_calibration:

        # generate energy map with same dimensions as roi image
        energy_map = np.zeros(image_stack_roi[0].shape)
        for y, x in np.ndindex(energy_map.shape):
            # apply energy function to each point on map
            energy_map[y, x] = pixel2energy(x, y, polynom_xy, polynom_xe)

        ax1 = plt.subplot(211)
        ax1.imshow(np.sum(image_stack_roi, axis=0))

        ax1.plot(COM['x'], COM['y'], 'o', c='r')

        for x, energy in zip(COM['x'], energies):
            plt.axvline(x, c='orange')
            plt.text(x+1, 3, f'{energy}', c='orange')

        ax2 = plt.subplot(212)
        # ax2.contour(energy_map, np.flip(energies), colors='orange')
        ax2.imshow(energy_map, cmap='gray')

        # plot best fit line for x-y
        _x = np.arange(0, energy_map.shape[1], 1)
        _y = np.polyval(polynom_xy, _x)
        ax2.plot(_x, _y, '-', c='r')
        plt.show()

        # # 3D plot
        # plt.clf()
        # ax = plt.axes(projection='3d')
        # ax.scatter(COM['x'], COM['y'], energies, c='darkblue')
        # x_grid, y_grid = np.meshgrid(np.arange(0, energy_map.shape[1]), np.arange(0, energy_map.shape[0]))
        # ax.plot_surface(x_grid, y_grid, energy_map, cmap='viridis', alpha=0.5)
        # plt.show()

    return polynom_xy, polynom_xe


def plot_calibration(image_stack, roi, calib_energies, polynom_xy, polynom_xe, ax=None):
    
    if ax is None:
        ax = plt.axes()

    im_sum = np.sum(image_stack, axis=0)
    ax.imshow(im_sum, vmin=0, vmax=100)

    x, dx, y, dy = roi.values()
    # plot box around roi
    ax.vlines([x, x + dx], y, y + dy, colors='y')
    ax.hlines([y, y + dy], x, x + dx, colors='y')

    # plot fitted x-y line across roi
    ax.plot([x, x+dx], [np.polyval(polynom_xy, 0) + y, np.polyval(polynom_xy, dx) + y], 'r-')

    energy_map = np.zeros(image_stack.shape[1:])
    for _y, _x in np.ndindex(energy_map.shape):
        energy_map[_y, _x] = pixel2energy(_x - x, _y - y, polynom_xy, polynom_xe)
    
    ax.contour(energy_map, np.flip(calib_energies), colors='orange')

    return ax


def pixel2energy(x, y, p_xy, p_xe):
    x_p, y_p = project_pt2line(x, y, p_xy)
    energy = np.polyval(p_xe, x_p)
    return energy


def reduce_image(image2d_crop, p_xy, p_xe, plot_spectrum=False):

    x_grid, y_grid = np.meshgrid(np.arange(image2d_crop.shape[1]), np.arange(image2d_crop.shape[0]))
    x_array, y_array = x_grid.ravel(), y_grid.ravel()
    image_array = image2d_crop.ravel()
    energy_array = pixel2energy(x_array, y_array, p_xy, p_xe)

    # energy_map = np.zeros(image2d_crop.shape)
    # for y, x in np.ndindex(energy_map.shape):
    #     # apply energy function to each point on map
    #     energy_map[y, x] = pixel2energy(x, y, p_xy, p_xe)

    # calc_energies = {'point': [], 'energy': []}
    # for y, x in np.ndindex(image2d.shape):
    #     energy = energy_func(x, y)
    #     calc_energies['point'].append((x,y))
    #     calc_energies['energy'].append(energy)

    # energy_array = np.array(calc_energies['energy'])
    # # print(energy_array.min(), energy_array.max())

    energy_edges, energy_step = np.linspace(energy_array.min() - 1e-6, energy_array.max() + 1e-6, image2d_crop.shape[1] + 1, retstep=True)
    energy_centers = energy_edges[:-1] + energy_step/2
    inds = np.digitize(energy_array, energy_edges)
    inds -= 1  # For Denis 

    intensity = np.zeros(energy_centers.size)

    for i in range(energy_centers.size):
        intensity[i] = np.sum(image_array[inds == i])

    # for n, pt in enumerate(calc_energies['point']):
    #     xpt, ypt = pt[0], pt[1]
    #     intensity = image2d[ypt, xpt]
    #     energy_bin_intensities[inds[n] - 1] += intensity

    if plot_spectrum:
        plt.plot(energy_centers, intensity)
        plt.show()

    return energy_centers, intensity


def test_calibration(calib_image_stack_roi, calib_energies, polynom_xy, polynom_xe):
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

    fwhm_array = np.zeros(calib_image_stack_roi.shape[0])

    # ax1 = plt.subplot(211)
    for i, im2d in enumerate(calib_image_stack_roi):
        energy_centers, intensity = reduce_image(im2d, polynom_xy, polynom_xe, plot_spectrum=False)
        ax1.plot(energy_centers, intensity, c='k')

        Ecen, fwhm, I_cor, I_fit, I_fit_raw = fit_gaussian_with_estimation(energy_centers, intensity)
        fwhm_array[i] = fwhm

        ax1.axvline(calib_energies[i], c='b')
        ax1.axvline(Ecen, c='g', ls='--')
        ax1.plot(energy_centers, I_fit_raw, 'r--')
    
    # ax2 = plt.subplot(212)
    ax2.set_xlabel('Energy (eV)')
    ax2.set_ylabel('Resolution (fwhm)')
    ax2.plot(calib_energies, fwhm_array, 'ko')

    plt.show()

    return ax1, ax2



if __name__ == '__main__':
    
    def test():
        
        PATH = '/home/charles/Desktop/XES_calib'
        DATA_FILE = 'Cu_calibration.json'
        MD_FILE = 'Cu_calibration_md.json'


        data, md = file_io(f'{PATH}/{DATA_FILE}', f'{PATH}/{MD_FILE}')
        
        roi = get_roi(md)
        # roi = {'x': 100, 'dx': 250, 'y': 80, 'dy': 20}
        print(roi)
        pix_array = get_image_array(data)

        pix_roi = crop_roi(pix_array, roi)
        # rotate to test x-y fitting
        # for i, image in enumerate(pix_roi):
        #     pix_roi[i] = rotate(image, 3, reshape=False)
        
        energies = get_calib_energies(data)

        p_xy, p_xe = run_calibration(pix_roi, energies, plot_calibration=True)

        ax_test = plt.axes()
        plot_calibration(pix_array, roi, energies, p_xy, p_xe, ax=ax_test)
        plt.show()

        test_calibration(pix_roi, energies, p_xy, p_xe)
        
    test()
