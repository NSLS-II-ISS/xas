import numpy as np
import matplotlib.pyplot as plt
from xas.file_io import load_binned_df_from_file
import time as ttime
from xas.fitting import fit_gaussian_with_estimation, fit_gaussian, Nominal2ActualConverter

# from lmfit.models import GaussianModel, LinearModel
# _gmodel = GaussianModel() + LinearModel()

class VonHamosScan:

    def __init__(self, db, fname):

        df, header = load_binned_df_from_file(fname)
        uid = [i for i in header.split('\n# ') if 'uid' in i][0].split(' ')[-1]

        start = ttime.time()
        t = db[uid].table(fill=True)
        print(f'db read in {ttime.time() - start}')

        self.energy = df.energy.values
        self.i0 = df.i0.values
        self.iff = df.iff.values
        self.images = np.array([i.squeeze() for i in t.pil100k_image])

        if not np.isclose(t.hhm_energy.values[0], self.energy[0], 1e-4):
            self.images = self.images[::-1, :, :]

        self.images = self.images/np.abs(self.i0)[:, None, None]
        self.muf = self.iff/self.i0

        self.total_image = np.sum(self.images, axis=0)
        self.y = 0
        self.x = 0
        self.dy, self.dx = self.total_image.shape
        self.energy_converter = None


    def set_roi(self, y, dy, x, dx):
        self.y = y
        self.dy = dy
        self.x = x
        self.dx = dx

    def show_roi(self, fignum=1, vmin=None, vmax=None):
        if vmin is None: vmin = self.total_image.min()
        if vmax is None: vmax = np.percentile(np.unique(self.total_image), 50)
        ysize, xsize = self.total_image.shape

        plt.figure(fignum)
        plt.clf()
        plt.imshow(self.total_image, vmin=vmin, vmax=vmax)
        plt.vlines([self.x, self.x + self.dx], 0, ysize, colors='r')
        plt.hlines([self.y, self.y + self.dy], 0, xsize, colors='r')

        plt.xlim(self.x - 10, self.x + self.dx + 10)
        plt.ylim(self.y - 10, self.y + self.dy + 10)


    def integrate_images(self):
        self.pixel = np.arange(self.x, self.x +  self.dx + 1)
        self.rixs = np.mean(self.images[:, self.y : self.y + self.dy + 1, self.x : self.x + self.dx + 1], axis=1).T

    def append_calibration(self, calibration):
        self.energy_converter = calibration.energy_converter

    @property
    def emission_energy(self):
        try:
            return self.energy_converter.nom2act(self.pixel)
        except:
            raise Exception('No energy converter')




class VonHamosCalibration(VonHamosScan):

    def __init__(self, *args):
        super().__init__(*args)

    def calibrate(self):
        self.rixs /= np.max(self.rixs, axis=0)[None, :]

        n_scans = self.rixs.shape[1]

        self.pixel_cen = np.zeros(n_scans)
        self.pixel_fwhm = np.zeros(n_scans)

        plt.figure(2)
        plt.clf()

        for i in range(n_scans):
            self.pixel_cen[i], self.pixel_fwhm[i], I_fit_raw = self._fit_elastic_line(self.pixel, self.rixs[:, i])

            plt.plot(self.pixel, self.rixs[:, i] - i, 'k-')
            plt.plot(self.pixel, I_fit_raw - i, 'r')

        self.energy_converter = Nominal2ActualConverter(self.pixel_cen, self.energy, n_poly=2)

    def _fit_elastic_line(self, x, y, threshold=5):
        cen, fwhm, _, _, y_fit = fit_gaussian(x, y, x.min() + y.argmax(), 1)
        if fwhm > threshold:
            y_new = y - y_fit
            y_new /= y_new.max()
            # pix_cen, fwhm, _, _, I_fit_raw = fit_gaussian_with_estimation(self.pixel, new_spectrum)
            cen, fwhm, _, _, y_fit = fit_gaussian(x, y_new, x.min() + y_new.argmax(), 1)
        return cen, fwhm, y_fit


        #
        #

        # plt.vlines(self.pixel.min() + self.rixs[:, 10].argmax(), 0, 1)

        #
