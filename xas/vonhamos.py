import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from xas.file_io import load_binned_df_from_file, write_df_to_file, load_extended_data_from_file
import time as ttime
from xas.fitting import fit_gaussian_with_estimation, fit_gaussian, Nominal2ActualConverter
import pandas as pd
import os

from scipy.ndimage import center_of_mass, rotate
from scipy import linalg
from sklearn.covariance import MinCovDet

PILATUS_KEY = 'Pilatus 100k New'

DET2KEY = {'Pilatus 100k' : 'pil100k',
           'Pilatus 100k New': 'pil100k2'}


def percentile_threshold_filter(im2d, pmin=5, pmax=99.5):
    """ set values below percentile equal to 0 """
    filt_im = im2d.copy()
    mask = (filt_im > np.percentile(filt_im, pmin)) & (filt_im < np.percentile(filt_im, pmax))
    filt_im[mask] = 0
    return filt_im


def fit_plane(X, Y, Z):
    """ fit plane to 3d data using least square regression """
    X, Y, Z = np.array(X), np.array(Y), np.array(Z)
    assert (X.size == Y.size == Z.size)
    A = np.column_stack([X, Y, np.ones(X.size)])
    C, _, _, _ = linalg.lstsq(A, Z)

    # equation for plane: z = a*x + b*y + c
    a, b, c = C
    return a, b, c


def project_pt2line(x_pt, y_pt, p_xy):
    slope, intercept = p_xy
    x_proj = (x_pt + slope * y_pt - slope * intercept) / (slope ** 2 + 1)
    y_proj = slope * x_proj + intercept
    return x_proj, y_proj


def get_roi(metadata, roi='roi1', detector=PILATUS_KEY):
    rois = metadata['detectors'][detector]['config']['roi']
    roi_ = rois[roi]
    return roi_


def get_image_array(extended_data):
    try:
        det_image = extended_data['pil100k2_image']
        image_array = np.array(list(det_image)).squeeze()
        if image_array.ndim == 2: # this is to handle "over-squeezing"
            image_array = image_array[None, :, :]
    except:
        det_image = extended_data['data_vars']['pil100k2_image']['data']
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


def get_p_xy(image_stack):
    image_total = np.sum(image_stack, axis=0)
    ys, xs = image_total.shape
    _x = np.arange(xs)
    _y = np.arange(ys)
    x, y = np.meshgrid(_x, _y)
    x = x.ravel()
    y = y.ravel()
    intensity = image_total.ravel()
    p_xy = np.polyfit(x, y, 1, w=intensity)
    return p_xy


def reduce_image_alt(image2d_crop, p_xy):
    ys, xs = image2d_crop.shape
    x = np.arange(xs)
    y = np.arange(ys)
    x_array, y_array = np.meshgrid(x, y)
    x_array, y_array = x_array.ravel(), y_array.ravel()
    x_array_p, _ = project_pt2line(x_array, y_array, p_xy)
    image_array = image2d_crop.ravel()

    intensity = np.zeros(x.size)

    for i in range(x.size):
        weights = 1 - np.abs(x_array_p - x[i])
        weights[weights < 0] = 0
        intensity[i] = weights @ image_array
    return x, intensity


def run_calibration(image_stack_roi, energies, n_poly=2, output_diagnostics=False):
    assert image_stack_roi.shape[0] == len(
        energies), "number of calibration images must match number of calibration energies"

    p_xy = get_p_xy(image_stack_roi)
    x_pix = None
    intensity_total = None
    intensity_total_fit = None
    x_pix_centers = []
    fwhms = []

    for image in image_stack_roi:
        # apply percentile filter
        # filtered_image = percentile_threshold_filter(image, 99)
        x_pix, intensity = reduce_image_alt(image, p_xy)
        x_pix_center, fwhm, _, _, intensity_fit = fit_gaussian_with_estimation(x_pix, intensity)

        if intensity_total is None:
            intensity_total = np.zeros(intensity.size)
            intensity_total_fit = np.zeros(intensity.size)

        intensity_total += intensity
        intensity_total_fit += intensity_fit

        x_pix_centers.append(x_pix_center)
        fwhms.append(fwhm)

    p_xe = np.polyfit(x_pix_centers, energies, n_poly)

    if output_diagnostics:
        plot_calibration_diagnostics(np.sum(image_stack_roi, axis=0), energies,
                                     x_pix, intensity_total, intensity_total_fit, x_pix_centers, fwhms, p_xy, p_xe)
        return p_xy, p_xe, x_pix, intensity_total, intensity_total_fit, x_pix_centers
    else:
        return p_xy, p_xe


def plot_calibration_diagnostics(image_total, energies,
                                 x_pix, intensity_total, intensity_total_fit,
                                 x_pix_centers, fwhms, p_xy, p_xe):
    y_pix_centers = np.polyval(p_xy, x_pix_centers)
    ax1 = plt.subplot(221)
    ax1.imshow(image_total)
    ax1.plot(x_pix_centers, y_pix_centers, 'o', c='r')

    # generate energy map with same dimensions as roi image
    energy_map = np.zeros(image_total.shape)
    for y, x in np.ndindex(energy_map.shape):
        energy_map[y, x] = pixel2energy(x, y, p_xy, p_xe)
    ax2 = plt.subplot(222)
    ax2.imshow(energy_map, cmap='gray')

    _x = np.arange(0, energy_map.shape[1], 1)
    _y = np.polyval(p_xy, _x)
    ax2.plot(_x, _y, '-', c='r')

    plt.subplot(223)
    plt.plot(np.polyval(p_xe, x_pix), intensity_total, 'k.-')
    plt.plot(np.polyval(p_xe, x_pix), intensity_total_fit, 'r-')

    energy_hi = np.polyval(p_xe, np.array(x_pix_centers) - np.array(fwhms) / 2)
    energy_lo = np.polyval(p_xe, np.array(x_pix_centers) + np.array(fwhms) / 2)
    fwhms_energy = energy_hi - energy_lo

    plt.subplot(224)
    plt.plot(energies, fwhms_energy, 'k.-')


def pixel2energy(x, y, p_xy, p_xe):
    x_p, y_p = project_pt2line(x, y, p_xy)
    energy = np.polyval(p_xe, x_p)
    return energy


def get_cropped_image_stack(data, roi_coords):
    image_stack = get_image_array(data)
    image_stack = crop_roi(image_stack, roi_coords)
    return image_stack

def process_calibration_for_roi(df, md, roi='roi1', roi_dict=None, detector=PILATUS_KEY, output_diagnostics=False):
    if roi_dict is None:
        roi_dict = md['detectors'][detector]['config']['roi']
    image_stack = get_cropped_image_stack(df, roi_dict[roi])
    energies = get_calib_energies(df)
    return run_calibration(image_stack, energies, output_diagnostics=output_diagnostics)

def trivial_calibration():
    p_xy, p_xe = [1, 0], [1, 0]
    return p_xy, p_xe

def process_calibration_for_roi_uid(uid, db, **kwargs):
    hdr = db[uid]
    md = hdr.start
    df = hdr.table(fill=True)
    return process_calibration_for_roi(df, md, **kwargs)

def scan_and_calibration_roi_match(md, uid_calibration, db, roi='roi1', detector=PILATUS_KEY):
    md_calibration = db[uid_calibration].start
    roi_scan = get_roi(md, roi=roi, detector=detector)
    roi_calibration = get_roi(md_calibration, roi=roi, detector=detector)
    return (roi_scan == roi_calibration)

def apply_calibration_for_roi(df, extended_data, md, uid_calibration, db, roi='roi1', roi_dict=None, detector=PILATUS_KEY, droi=5):


    if roi_dict is None:
        roi_dict = md['detectors'][detector]['config']['roi']
        roi_dict_calibration = None
        enforce_roi = False
    else:
        roi_dict_calibration = roi_dict
        enforce_roi = True

    do_trivial_calibration = False

    reason = ''
    if db is None:
        do_trivial_calibration = True
        reason += '- databroker is not defined in the von hamos processing pipeline\n'

    if (not do_trivial_calibration) and (not uid_calibration):
        do_trivial_calibration = True
        reason += '- calibration uid is not defined\n'

    if (not do_trivial_calibration) and (not scan_and_calibration_roi_match(md, uid_calibration, db, roi=roi, detector=detector)) and (not enforce_roi):
        do_trivial_calibration = True
        reason += f'- ROI coordinates mismatch between the data scan and the calibration scan\n'

    if do_trivial_calibration:
        print(f'Could not apply energy calibration to {roi}. Reason(s):\n{reason}')
        p_xy, p_xe = trivial_calibration()
    else:
        p_xy, p_xe = process_calibration_for_roi_uid(uid_calibration, db, roi=roi, roi_dict=roi_dict_calibration, detector=detector)

    image_stack = get_cropped_image_stack(extended_data, roi_dict[roi])
    # need to do something about the bkg intensity using droi
    intensity = []
    for image in image_stack:
        x_pix, _intensity = reduce_image_alt(image, p_xy)
        intensity.append(_intensity)

    intensity = np.array(intensity)
    energy = np.polyval(p_xe, x_pix)

    vh_roi_data = {f'energy' : energy, f'{DET2KEY[detector]}' : intensity}
    return vh_roi_data, roi_dict


def process_von_hamos_scan(df, extended_data, comments, hdr, path_to_file, db=None,
                           detector=PILATUS_KEY, roi_keys=None, roi_dict=None,
                           droi=5, save_dat=True):
    comments += f'# Spectrometer.type: von Hamos\n' \
                f'# Spectrometer.detector: {detector}'

    md = hdr.start

    if roi_keys is None:
        roi_keys = ['roi1']

    vh_data_dict = {}

    for roi_key in roi_keys:

        if md['scan_for_calibration_purpose']:
            uid_calibration = md['uid']
        else:
            uid_calibration = md['spectrometer_config']['energy_calibration_uid']
        vh_data_dict[roi_key], roi_dict = apply_calibration_for_roi(df, extended_data, md, uid_calibration, db, roi=roi_key, roi_dict=roi_dict, detector=detector, droi=droi)

    extended_data = {**extended_data, 'von_hamos_data' : vh_data_dict}

    for k, roi in roi_dict.items():
        for c, v in roi.items():
            comments += f'# Spectrometer.detector.{k}.{c}: {v}\n'

    if save_dat:
        file_paths = save_vh_data_to_file(path_to_file, df, vh_data_dict, comments)
        # if vh_scan.kind == 'xes':
        #     file_paths = save_vh_scan_to_file(path_to_file, vh_scan, comments)
        # else:
        #     file_paths = []
    else:
        file_paths = []

    return extended_data, comments, file_paths

def make_vh_dfs(df, vh_data_dict):
    try:
        hhm_energy = df['energy'].values
    except:
        hhm_energy = None
    i0 = df['i0'].values
    vh_dfs = []
    suffixes = []
    for roi_key, roi_data in vh_data_dict.items():
        _df = {}
        for data_key, data_array in roi_data.items():
            if data_key == 'energy':
                _df['energy'] = data_array
            else:
                for k, intensity in enumerate(data_array):
                    if hhm_energy is not None:
                        energy_key = f'_{hhm_energy[k]:.2f}'.replace('.', '_') # note the underscore
                    else:
                        energy_key = ''
                    _df_key = f'{data_key}_{roi_key}{energy_key}'
                    _df[_df_key] = intensity
                    _df_key = f'i0{energy_key}'
                    _df[_df_key] = i0[k] * np.ones(intensity.size)
        _df = pd.DataFrame(_df)
        vh_dfs.append(_df)
        suffixes.append(f'vh_{roi_key}')
    return vh_dfs, suffixes



def save_vh_data_to_file(path_to_file, df, vh_data_dict, comments):
    dfs, suffixes = make_vh_dfs(df, vh_data_dict)
    (path, extension) = os.path.splitext(path_to_file)
    paths = []
    for df, suffix in zip(dfs, suffixes):
        path_frame = f'{path} {suffix}{extension}'
        print(f'VON HAMOS PROCESSING: data will be saved in {path_frame}')
        write_df_to_file(path_frame, df, comments)
        paths.append(path_frame)

        # for visualization and debugging
        path_frame = f'{path} {suffix}.rixs'
        write_df_to_file(path_frame, df, comments)
    return paths
# legacy
# def save_vh_scan_to_file(path_to_file, vh_scan, comments):
#     dfs, suffixes = vh_scan.make_dfs()
#     (path, extension) = os.path.splitext(path_to_file)
#     paths = []
#     for df, suffix in zip(dfs, suffixes):
#         path_frame = f'{path} {suffix}{extension}'
#         print(f'VON HAMOS PROCESSING: data will be saved in {path_frame}')
#         write_df_to_file(path_frame, df, comments)
#         paths.append(path_frame)
#     return paths



def process_von_hamos_scan_legacy(df, extended_data, comments, hdr, path_to_file, detector=PILATUS_KEY, roi='auto', droi=5, save_dat=True):

    # if ('spectrometer_scan_kind' in hdr.start.keys()) and (hdr.start['spectrometer_scan_kind'] == 'calibration'):
    #     vh_scan = VonHamosCalibration(df, extended_data)
    # else:
    vh_scan = VonHamosScan(df, extended_data)

    if roi == 'auto':
        roi_dict = hdr.start['detectors'][detector]['config']['roi']
    else:
        # TODO: assert roi to be a dict
        roi_dict = roi

    vh_scan.set_roi(roi_dict, droi=droi)
    vh_scan.integrate_images()

    # if ('spectrometer_scan_kind' in hdr.start.keys()) and (hdr.start['spectrometer_scan_kind'] == 'calibration'):
    #     try:
    #         vh_scan.calibrate(['roi1'])
    #     except Exception as e:
    #         print(e)
    #
    # if ('spectrometer_calibration_uid' in hdr.start.keys()):
    #     spectrometer_calibration_uid = hdr.start['spectrometer_calibration_uid']
    #     vh_calibration = vh_scan.from_db(hdr.start[spectrometer_calibration_uid])

    extended_data = vh_scan.augment_extended_data(extended_data)

    comments += f'# Spectrometer.type: von Hamos\n' \
                f'# Spectrometer.detector: {detector}'
    for k, roi in roi_dict.items():
        for c, v in roi.items():
            comments += f'# Spectrometer.detector.{k}.{c}: {v}\n'

    if save_dat:
        if vh_scan.kind == 'xes':
            file_paths = save_vh_scan_to_file(path_to_file, vh_scan, comments)
        else:
            file_paths = []
    else:
        file_paths = []

    return extended_data, comments, file_paths


















# class VonHamosScan:
#
#     @classmethod
#     def from_db(cls, db, uid, ext_data_path='extended_data'):
#         hdr = db[uid]
#         fpath = hdr.start['interp_file']
#         fpath, _ = os.path.splitext(fpath)
#         fpath += '.dat'
#         return cls.from_h5(fpath, ext_data_path=ext_data_path)
#     #
#     @classmethod
#     def from_h5(cls, path_to_file, ext_data_path='extended_data'):
#
#         df, _ = load_binned_df_from_file(path_to_file)
#         folder, file = os.path.split(path_to_file)
#         folder = os.path.join(folder, ext_data_path)
#         filename, _ = os.path.splitext(file)
#         filename += '.h5'
#         path_to_ext_file = os.path.join(folder, filename)
#         extended_data = load_extended_data_from_file(path_to_ext_file)
#         return cls(df, extended_data)
#
#     def __init__(self, df, extended_data, image_key='pil100k_image'):
#
#         if 'energy' in df.columns:
#             self.energy = df.energy.values
#             self.kind = 'rixs'
#         else:
#             self.energy = None
#             self.kind = 'xes'
#
#         self.i0 = df.i0.values
#         # self.iff = df.iff.values
#         self.images = extended_data[image_key]
#
#         if self.kind == 'rixs':
#             if not np.isclose(df.energy.values[0], self.energy[0], 1e-4):
#                 self.images = self.images[::-1, :, :]
#
#         self.images = self.images#/np.abs(self.i0)[:, None, None]
#         # self.muf = self.iff#/self.i0
#
#         self.total_image = np.sum(self.images, axis=0)
#         self.energy_converter = None
#
#
#     def set_roi(self, roi_dict, droi=5):
#         self.roi_dict = roi_dict
#         self.droi = droi
#         self.remove_bkg_from_images()
#
#     def show_roi(self, fignum=1, vmin=None, vmax=None):
#         if vmin is None: vmin = self.total_image.min()
#         if vmax is None: vmax = np.percentile(np.unique(self.total_image), 50)
#         ysize, xsize = self.total_image.shape
#
#         fig, ax = plt.subplots(num=fignum, clear=True)
#         ax.imshow(self.total_image, vmin=vmin, vmax=vmax)
#         for k, roi in self.roi_dict.items():
#             rect = patches.Rectangle((roi['x'], roi['y']), roi['dx'], roi['dy'],
#                                      linewidth=1, edgecolor='r', facecolor='none')
#             ax.add_patch(rect)
#
#         plt.xlim(0, xsize)
#         plt.ylim(0, ysize)
#
#     def _convert_roi_to_indexes(self, roi, droi=False):
#         ysize, xsize = self.total_image.shape
#         x, y, dx, dy = roi['x'], roi['y'], roi['dx'], roi['dy']
#         if droi:
#             x -= self.droi
#             y -= self.droi
#             dx += 2 * self.droi
#             dy += 2 * self.droi
#         x1 = x
#         x2 = np.min([x + dx + 1, xsize])
#         y1 = y
#         y2 = np.min([y + dy + 1, ysize])
#         return x1, x2, y1, y2
#
#     def remove_bkg_from_images(self):
#         for k, roi in self.roi_dict.items():
#             x1, x2, y1, y2 = self._convert_roi_to_indexes(roi)
#             xw1, xw2, yw1, yw2 = self._convert_roi_to_indexes(roi, droi=True)
#
#
#             bkg_mask = np.zeros(self.total_image.shape, dtype=bool)
#             bkg_mask[yw1:yw2, xw1:xw2] = True
#             bkg_mask[y1:y2, x1:x2] = False
#
#             # image_bkg = image.copy()
#             # image_bkg[~bkg_mask] = -100
#
#             ysize, xsize = self.total_image.shape
#             y_mesh, x_mesh = np.meshgrid(np.arange(ysize), np.arange(xsize), indexing='ij')
#
#             y_bkg = y_mesh[bkg_mask].ravel()
#             x_bkg = x_mesh[bkg_mask].ravel()
#
#             data_bkg = np.array([im[bkg_mask].ravel() for im in self.images])
#             # dfg
#             # mask = mask_by_percentiles(i_bkg)
#             mask = np.all(data_bkg >= 0, axis=0)  # & (i_bkg < thresh)
#             y_bkg = y_bkg[mask]
#             x_bkg = x_bkg[mask]
#             data_bkg = data_bkg[:, mask]
#
#             c = fit_linear_surf(y_bkg, x_bkg, data_bkg.T)
#
#             A_fit = np.hstack((y_mesh.ravel()[:, None], x_mesh.ravel()[:, None], np.ones((y_mesh.ravel().size, 1))))
#             self.images_bkg = (A_fit @ c).reshape(self.images.shape)
#             self.images_no_bkg = self.images - self.images_bkg
#
#     def integrate_images(self):
#         self.xes = {}
#         for k, roi in self.roi_dict.items():
#             x1, x2, y1, y2 = self._convert_roi_to_indexes(roi)
#             pixel = np.arange(x1, x2)
#             intensity = np.sum(self.images[:, y1: y2, x1: x2], axis=1)
#             intensity_bkg = np.sum(self.images_bkg[:, y1: y2, x1: x2], axis=1)
#             intensity_no_bkg = np.sum(self.images_no_bkg[:, y1: y2, x1: x2], axis=1)
#             self.xes[k] = {'pixel' : pixel,
#                            'intensity' : intensity,
#                            'intensity_bkg' : intensity_bkg,
#                            'intensity_no_bkg' : intensity_no_bkg}
#
#     # def append_calibration(self, calibration):
#     #     self.energy_converter = calibration.energy_converter
#
#     def augment_extended_data(self, extended_data):
#         aug_data = {}
#         for k, v in self.xes.items():
#             key = f'pil100k_{k}_vh'
#             aug_data[key] = v
#         return {**extended_data, **aug_data}
#
#     # @property
#     # def emission_energy(self):
#     #     try:
#     #         return self.energy_converter.nom2act(self.pixel)
#     #     except:
#     #         raise Exception('No energy converter')
#     #
#     def make_dfs(self, rois=['roi1', 'roi2']):
#
#         # if self.kind == 'xes':
#         dfs = []
#         suffixes = []
#         for roi in rois:
#             for i in range(self.xes[roi]['intensity'].shape[0]):
#                 if self.energy is not None:
#                     energy_str = f' {str(int(self.energy[i]))} eV'
#                 else:
#                     energy_str = ''
#                 suffixes.append(f'{roi} frame {(i+1):04d}{energy_str}')
#
#                 data_dict = {'energy' : self.xes[roi]['pixel'],
#                      'i0': np.ones(self.xes[roi]['pixel'].shape) * self.i0[i],
#                      f'pil100k_VH__intensity' : self.xes[roi]['intensity'][i, :],
#                      f'pil100k_VH__bkg' : self.xes[roi]['intensity_bkg'][i, :],
#                      f'pil100k_VH__intensity_no_bkg' : self.xes[roi]['intensity_no_bkg'][i, :]}
#                 dfs.append(pd.DataFrame(data_dict))
#         return dfs, suffixes
#         # elif self.kind == 'rixs':
#         #     pass
#
#
# class VonHamosCalibration(VonHamosScan):
#
#     def __init__(self, *args):
#         super().__init__(*args)
#
#     def calibrate(self, rois=None, n_poly=2, plotting=False):
#         if rois is None:
#             rois = self.roi_dict.keys()
#         for roi in rois:
#             self.calibrate_roi(roi, n_poly=2, plotting=False)
#
#
#     def calibrate_roi(self, roi, n_poly=2, plotting=False):
#         xes = self.xes[roi]
#         pixel = xes['pixel']
#         intensity = xes['intensity'].copy()
#
#         intensity /= np.max(intensity, axis=0)[None, :]
#
#         n_scans = intensity.shape[1]
#
#         pixel_cen = np.zeros(n_scans)
#         pixel_fwhm = np.zeros(n_scans)
#
#         if plotting:
#             plt.figure()
#
#         for i in range(n_scans):
#             pixel_cen[i], pixel_fwhm[i], intensity_fit = self._fit_elastic_line(pixel, intensity)
#
#             if plotting:
#                 plt.plot(pixel, intensity[:, i] - i, 'k-')
#                 plt.plot(pixel, intensity_fit - i, 'r')
#
#         energy_converter = Nominal2ActualConverter(pixel_cen, self.energy, n_poly=n_poly)
#         self.xes[roi]['emission_energy'] = energy_converter.nom2act(pixel)
#
#
#
#
#     def _fit_elastic_line(self, x, y, threshold=5):
#         cen, fwhm, _, _, y_fit = fit_gaussian(x, y, x.min() + y.argmax(), 1)
#         if fwhm > threshold:
#             y_new = y - y_fit
#             y_new /= y_new.max()
#             # pix_cen, fwhm, _, _, I_fit_raw = fit_gaussian_with_estimation(self.pixel, new_spectrum)
#             cen, fwhm, _, _, y_fit = fit_gaussian(x, y_new, x.min() + y_new.argmax(), 1)
#         return cen, fwhm, y_fit









