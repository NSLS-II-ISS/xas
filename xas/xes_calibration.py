import json
import matplotlib
matplotlib.use('TkAgg')  # something wrong with my (Charlie's) system Qt 
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.ndimage import center_of_mass
from scipy.interpolate import interp1d
from scipy import stats

PATH = '/home/charles/Desktop/XES_calib'
DATA_FILE = 'Cu_calibration.json'
MD_FILE = 'Cu_calibration_md.json'


def percentile_threshold(im2d, p):
   """set values below percentile equal to 0"""
   im2d[im2d < np.percentile(im2d, p)] = 0
   return im2d

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
    print(pix_roi[0])

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

    slope, intercept, _, _, _ = stats.linregress(COM['x'], energies)
    def energy_function(x):
        return intercept + slope * x

    # generate energy map with same dimensions as roi image
    energy_map = np.zeros(pix_roi[0].shape)
    for y, x in np.ndindex(energy_map.shape):
        # apply energy function to each point on map
        energy_map[y, x] = energy_function(x)

    ax1 = plt.subplot(211)
    ax1.imshow(np.sum(pix_roi, axis=0))
    for x, y in zip(COM['x'], COM['y']):
        ax1.plot(x, y, 'o', c='r')

    ax2 = plt.subplot(212)
    ax2.imshow(energy_map, cmap='OrRd')

    plt.show()

if __name__ == '__main__':
    main()
