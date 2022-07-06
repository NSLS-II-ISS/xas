import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from xas.file_io import load_binned_df_from_file, write_df_to_file
import time as ttime
from xas.fitting import fit_gaussian_with_estimation, fit_gaussian, Nominal2ActualConverter
import pandas as pd
import os
# from lmfit.models import GaussianModel, LinearModel
# _gmodel = GaussianModel() + LinearModel()



class VonHamosScan:

    # @classmethod
    # def from_db(cls, db, uid):
    #     df = db[uid].table(fill=True)
    #     return cls(df)
    #
    # @classmethod
    # def from_file(cls, db, path_to_file):
    #     _, header = load_binned_df_from_file(path_to_file)
    #     uid = [i for i in header.split('\n# ') if 'uid' in i][0].split(' ')[-1]
    #     return cls.from_db(db, uid)

    def __init__(self, df, extended_data, image_key='pil100k_image'):

        if 'energy' in df.columns:
            self.energy = df.energy.values
            self.kind = 'rixs'
        else:
            self.kind = 'xes'

        self.i0 = df.i0.values
        # self.iff = df.iff.values
        self.images = extended_data[image_key]

        if self.kind == 'rixs':
            if not np.isclose(df.energy.values[0], self.energy[0], 1e-4):
                self.images = self.images[::-1, :, :]

        self.images = self.images#/np.abs(self.i0)[:, None, None]
        # self.muf = self.iff#/self.i0

        self.total_image = np.sum(self.images, axis=0)
        self.energy_converter = None


    def set_roi(self, roi_dict, droi=5):
        self.roi_dict = roi_dict
        self.droi = droi

    def show_roi(self, fignum=1, vmin=None, vmax=None):
        if vmin is None: vmin = self.total_image.min()
        if vmax is None: vmax = np.percentile(np.unique(self.total_image), 50)
        ysize, xsize = self.total_image.shape

        fig, ax = plt.subplots(num=fignum, clear=True)
        ax.imshow(self.total_image, vmin=vmin, vmax=vmax)
        for k, roi in self.roi_dict.items():
            rect = patches.Rectangle((roi['x'], roi['y']), roi['dx'], roi['dy'],
                                     linewidth=1, edgecolor='r', facecolor='none')
            ax.add_patch(rect)

        plt.xlim(0, xsize)
        plt.ylim(0, ysize)

    def integrate_images(self):
        self.xes = {}
        ysize, xsize = self.total_image.shape
        for k, roi in self.roi_dict.items():
            x, y, dx, dy = roi['x'], roi['y'], roi['dx'], roi['dy']
            x1 = x
            x2 = np.min([x + dx + 1, xsize])
            y1 = y
            y2 = np.min([y + dy + 1, ysize])
            pixel = np.arange(x1, x2)
            intensity = np.sum(self.images[:, y1: y2, x1: x2], axis=1)
            self.xes[k] = {'pixel' : pixel, 'intensity' : intensity}

    def append_calibration(self, calibration):
        self.energy_converter = calibration.energy_converter

    def augment_extended_data(self, extended_data):
        aug_data = {}
        for k, v in self.xes.items():
            key = f'pil100k_{k}_vh'
            aug_data[key] = v
        return {**extended_data, **aug_data}

    @property
    def emission_energy(self):
        try:
            return self.energy_converter.nom2act(self.pixel)
        except:
            raise Exception('No energy converter')

    def make_dfs(self):

        if self.kind == 'xes':
            dfs = []
            for i in range(self.xes.shape[1]):
                d = {'energy' : self.pixel,
                     'i0': np.ones(self.pixel.shape) * self.i0[i],
                     'pil100k_VH_counts' : self.xes[:, i]}
                dfs.append(pd.DataFrame(d))
            return dfs
        elif self.kind == 'rixs':
            pass


class VonHamosCalibration(VonHamosScan):

    def __init__(self, *args):
        super().__init__(*args)

    def calibrate(self, n_poly=2, plotting=False):
        self.xes /= np.max(self.xes, axis=0)[None, :]

        n_scans = self.xes.shape[1]

        self.pixel_cen = np.zeros(n_scans)
        self.pixel_fwhm = np.zeros(n_scans)
        self.energy_converter = Nominal2ActualConverter(self.pixel_cen, self.energy, n_poly=n_poly)

        if plotting:
            plt.figure(clear=True)

            for i in range(n_scans):
                self.pixel_cen[i], self.pixel_fwhm[i], I_fit_raw = self._fit_elastic_line(self.pixel, self.xes[:, i])

                plt.plot(self.pixel, self.xes[:, i] - i, 'k-')
                plt.plot(self.pixel, I_fit_raw - i, 'r')



    def _fit_elastic_line(self, x, y, threshold=5):
        cen, fwhm, _, _, y_fit = fit_gaussian(x, y, x.min() + y.argmax(), 1)
        if fwhm > threshold:
            y_new = y - y_fit
            y_new /= y_new.max()
            # pix_cen, fwhm, _, _, I_fit_raw = fit_gaussian_with_estimation(self.pixel, new_spectrum)
            cen, fwhm, _, _, y_fit = fit_gaussian(x, y_new, x.min() + y_new.argmax(), 1)
        return cen, fwhm, y_fit


def process_von_hamos_scan(df, extended_data, comments, hdr, detector='Pilatus 100k', roi='auto', droi=5):
    vh_scan = VonHamosScan(df, extended_data)

    if roi == 'auto':
        roi_dict = hdr.start['detectors'][detector]['config']['roi']
    else:
        # TODO: assert roi to be a dict
        roi_dict = roi

    vh_scan.set_roi(roi_dict, droi=droi)
    vh_scan.integrate_images()

    df = vh_scan.augment_df(df)

    return df, comments


def save_vh_scan_to_file(path_to_file, vh_scan, comments):
    dfs = vh_scan.make_dfs()
    (path, extension) = os.path.splitext(path_to_file)
    for i, df in enumerate(dfs):
        path_frame = f'{path} frame {i + 1} {extension}'
        print(f'VON HAMOS PROCESSING: data will be saved in {path_frame}')
        write_df_to_file(path_frame, df, comments)







