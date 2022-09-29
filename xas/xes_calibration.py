import json
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')  # something wrong with my (Charlie's) system Qt 
import pandas as pd
from scipy.ndimage import center_of_mass

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
    max_x = roi1['x'] + roi1['dx']

    min_y = roi1['y']
    max_y = roi1['y'] + roi1['dy']

    pix_roi = pix_array[:, min_y:max_y, min_x:max_x]

    filtered_roi = np.zeros(pix_roi.shape)

    COM = {'x': [], 'y': []}  # centers of mass
    for i, image in enumerate(pix_roi):
        # apply percentile filter
        filtered_image = percentile_threshold(image, 99)
        filtered_roi[i] = filtered_image

        y_com, x_com = center_of_mass(filtered_roi[i])
        COM['x'].append(x_com)
        COM['y'].append(y_com)

    plt.imshow(np.sum(pix_roi, axis=0))
    for x, y in zip(COM['x'], COM['y']):
        plt.plot(x, y, marker='o', c='r')

    plt.show()

    energies = data['hhm_energy']
    
    assert(pix_roi.shape[0] == len(energies))
    n_energies = len(energies)




if __name__ == '__main__':
    main()
