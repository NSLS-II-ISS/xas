

import numpy as np
import pyFAI as pyFAI





def get_ai(dist=40, poni1=480, poni2=100,
           pixel1=0.172, pixel2=0.172,
           rot1=0, rot2=0, rot3=np.pi/2,
           energy=11300):
    wavelength = 12.3984 / energy * 1e-3 * 1e-10
    ai = pyFAI.AzimuthalIntegrator(dist=dist*1e-3,
                                   poni1=poni1 * pixel1 * 1e-3, poni2=poni2 * pixel2 * 1e-3,
                                   pixel1=pixel1 * 1e-3, pixel2=pixel2 * 1e-3,
                                   rot1=rot1, rot2=rot2, rot3=rot3,
                                   wavelength=wavelength)
    return ai



def integrate_pil100k_image_stack(df_images, energies, dist, center_x, center_y):
    for image
