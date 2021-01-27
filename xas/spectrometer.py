from xas.fitting import (get_normalized_gaussian_scan, estimate_center_and_width_of_peak, fit_gaussian,
                         Nominal2ActualConverter)
import numpy as np
from matplotlib import pyplot as plt

def analyze_elastic_scan(db, uid):
    E, I = get_normalized_gaussian_scan(db, uid)
    Ecen0, fwhm0 = estimate_center_and_width_of_peak(E, I)
    return (*fit_gaussian(E, I, Ecen0, fwhm0), E)



def analyze_many_elastic_scans(db, uids, E_nominal, plotting=False):
    E_actual = []
    resolution = []
    I_cors = []
    I_fits = []
    E_scans = []

    for uid in uids:
        Ecen, fwhm, I_cor, I_fit, E_scan = analyze_elastic_scan(db, uid)
        E_actual.append(Ecen)
        resolution.append(fwhm)
        I_cors.append(I_cor)
        I_fits.append(I_fit)
        E_scans.append(E_scan)

    E_actual = np.array(E_actual)
    resolution = np.sqrt(np.array(resolution)**2 - (1.3e-4*E_actual)**2)
    energy_converter = Nominal2ActualConverter(E_nominal, E_actual, n_poly=2)

    if plotting:
        plt.figure()
        plt.clf()

        plt.subplot(211)
        for i in range(len(I_cors)):
            plt.plot(E_scans[i], I_cors[i], 'k.-')
            plt.plot(E_scans[i], I_fits[i] , 'r-')

        plt.subplot(223)
        plt.plot(E_nominal, E_actual, 'k.', label='data')
        _e_nom_grid = np.linspace(E_nominal.min(), E_nominal.max(), E_nominal.size*10)
        _e_act_grid = energy_converter.nom2act(_e_nom_grid)
        plt.plot(_e_nom_grid, _e_act_grid, 'r-', label='fit')
        plt.xlabel('nominal energy')
        plt.ylabel('actual energy')
        plt.legend()

        plt.subplot(224)
        plt.plot(E_actual, resolution, 'k.-')
        plt.xlabel('actual energy')
        plt.ylabel('resolution')

    return energy_converter







class Crystal:

    def __init__(self, R, r, hkl, kind):
        self.R = R
        self.r = r
        self.hkl = hkl
        self.kind = kind
        self.d = self.lat_const / self.refl_order

    @property
    def lat_const(self):
        if self.kind == 'Si':
            return 5.431020511
        elif self.kind == 'Ge':
            return 5.658

    @property
    def refl_order(self):
        h, k, l = self.hkl
        return np.sqrt(h ** 2 + k ** 2 + l ** 2)

    def E2L(self, E):
        return 12398.4 / E

    def L2E(self, L):
        return 12398.4 / L

    def bragg_angle(self, E):
        L = self.E2L(E)
        return np.arcsin(L / (2 * self.d))

    def bragg_energy(self, ba):
        L = 2 * self.d * np.sin(ba)
        return self.L2E(L)

    def place_E(self, E):
        self.E = E
        self.ba = self.bragg_angle(E)
        self.ba_deg = np.rad2deg(self.ba)
        self._place()

    def place_ba(self, ba):
        self.ba = ba
        self.E = self.bragg_energy(ba)
        self._place()

    def _place(self):
        # self.Ox = self.R / 2 * np.cos(np.pi / 2 - self.ba)
        # self.Oy = self.R / 2 * np.sin(np.pi / 2 - self.ba)
        # self.Oz = 0
        #
        # self.OOx = 0
        # self.OOy = 2 * self.Oy
        # self.OOz = self.Oz
        #
        # self.x = self.Ox * 2
        # self.y = 0
        # self.z = 0
        #
        #
        # ksi = self.R / 2 * np.sqrt(2 - 2 * np.cos(2 * self.ba))
        # self.Dx = self.x - ksi * np.cos(np.pi - 2 * self.ba)
        # self.Dy = ksi * np.sin(np.pi - 2 * self.ba)
        # self.Dz = 0
        #
        # phi = np.linspace(0, np.pi * 2, 361)
        # self.rowland_circle_x = self.Ox + self.R / 2 * np.cos(phi)
        # self.rowland_circle_y = self.Oy + self.R / 2 * np.sin(phi)
        # self.rowland_circle_z = 0

        self.x = self.R/2 * (1 + np.cos(np.pi - 2*self.ba))
        self.y = self.R / 2 * np.sin(np.pi - 2 * self.ba)
        self.d_y = self.R  * np.sin(np.pi - 2 * self.ba)

