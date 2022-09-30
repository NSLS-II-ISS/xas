import json
import matplotlib
from matplotlib import projections
matplotlib.use('TkAgg')  # something wrong with my (Charlie's) system Qt 
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.ndimage import center_of_mass, rotate
from scipy import stats, linalg

from mpl_toolkits.mplot3d import Axes3D

PATH = '/home/charles/Desktop/XES_calib'
DATA_FILE = 'Cu_calibration.json'
MD_FILE = 'Cu_calibration_md.json'


def percentile_threshold(im2d, p):
   """set values below percentile equal to 0"""
   im2d[im2d < np.percentile(im2d, p)] = 0
   return im2d


def fit_plane(X, Y, Z):
    """ fit plane to 3d data using least square regression"""
    X, Y, Z = np.array(X), np.array(Y), np.array(Z)
    assert(X.size == Y.size == Z.size)
    A = np.column_stack([X, Y, np.ones(X.size)])
    C, _, _, _ = linalg.lstsq(A, Z)
    
    # equation for plane: z = a*x + b*y + c
    a, b, c = C[0], C[1], C[2]
    return a, b, c


def main(data_file=f'{PATH}/{DATA_FILE}', md_file=f'{PATH}/{MD_FILE}'):
    
    # Load data and metadata
    data = pd.read_json(data_file)
    with open(md_file) as metadata:
        md = json.load(metadata)[1]
    
    det_image = data['pil100k_image']
    pix_array = np.array(det_image.to_list()).squeeze()

    # Load ROI
    rois = md['detectors']['Pilatus 100k']['config']['roi']
    # print(rois)
    roi1 = rois['roi1']

    min_x = roi1['x']
    dx = roi1['dx']
    max_x = roi1['x'] + dx

    min_y = roi1['y']
    dy = roi1['dy']
    max_y = roi1['y'] + dy

    pix_roi = pix_array[:, min_y:max_y, min_x:max_x]
    for i, image in enumerate(pix_roi):
        pix_roi[i] = rotate(image, 3, reshape=False)

    filtered_roi = np.zeros(pix_roi.shape)

    COM = {'x': [], 'y': []}  # centers of mass
    for i, image in enumerate(pix_roi):
        # apply percentile filter
        filtered_image = percentile_threshold(image, 99)
        filtered_roi[i] = filtered_image

        y_com, x_com = center_of_mass(filtered_image)
        COM['x'].append(x_com)
        COM['y'].append(y_com)

    energies = data['hhm_energy']
    assert(pix_roi.shape[0] == len(energies))

    # slope, intercept, _, _, _ = stats.linregress(COM['x'], energies)
    a, b, c = fit_plane(COM['x'], COM['y'], energies)
    def energy_function(x, y):
        # return intercept + slope * x
        return a*x + b*y + c

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
    ax2.imshow(energy_map, cmap='OrRd')

    plt.show()

    # 3D plot
    plt.clf()
    ax = plt.axes(projection='3d')
    ax.scatter(COM['x'], COM['y'], energies, c='darkblue')
    x_grid, y_grid = np.meshgrid(np.arange(0, energy_map.shape[1]), np.arange(0, energy_map.shape[0]))
    ax.plot_surface(x_grid, y_grid, energy_map, cmap='OrRd', alpha=0.5)
    plt.show()

if __name__ == '__main__':
    main()
