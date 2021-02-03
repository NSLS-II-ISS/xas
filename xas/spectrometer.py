from xas.fitting import (get_normalized_gaussian_scan, estimate_center_and_width_of_peak, fit_gaussian,
                         Nominal2ActualConverter)
from xas.file_io import  load_binned_df_from_file

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



def get_xes_data(db, uid, offset=-10):
    t = db[uid].table()
    xes = t['pil100k_stats1_total']/(np.abs(t['apb_ave_ch1_mean']) - offset)
    return xes


def parse_file_with_uids(file_with_uids):
    uids_herfd, uids_xes = [], []
    with open(file_with_uids, "r") as f:
        lines = f.readlines()
    for line in lines:
        words = line.split(' ')
        uids_herfd.append(words[-2])
        uids_xes.append(words[-1][:-2])
    return uids_herfd, uids_xes


def get_herfd_data(db, uid):
    filename = db[uid].start['interp_filename']
    df, _ = load_binned_df_from_file(filename)
    energy = df['energy']
    herfd = np.abs((df['pil100_ROI1'] - df['pil100_ROI2'])/ df['i0']).values
    return energy, herfd



def parse_rixs_scan(db, file_with_uids):
    uids_herfd, uids_xes = parse_file_with_uids(file_with_uids)
    xes_data = []
    herfd_data = []
    for uid_herfd, uid_xes in zip(uids_herfd, uids_xes):
        energy, herfd = get_herfd_data(db, uid_herfd)
        xes = get_xes_data(db, uid_xes)
        herfd_data.append(herfd)
        xes_data.append(xes)
    herfd_data = np.array(herfd_data).T
    xes_data = np.array(xes_data).T
    # linear algebra magic for scaling purposes
    c, _, _, _ = np.linalg.lstsq(xes_data[:,:1], xes_data, rcond=-1)
    herfd_data /= c
    return herfd_data, xes_data, energy


def convert_rixs_to_energy_transfer(Ein, Eout, herfd):
    dEmin = Ein.min() - Eout.max()
    dEmax = Ein.max() - Eout.min()

    transfer_step = np.min([np.min(np.abs(np.diff(Ein))), np.abs(Eout[1] - Eout[0])])*np.sqrt(2)
    dE = np.arange(dEmin, dEmax, transfer_step)
    # dE = np.linspace(dEmin, dEmax, 601)
    rixs = np.zeros((Ein.size, dE.size))

    for idx in range(Ein.size):
        each_Ein = Ein[idx]
        each_transfer = (each_Ein - Eout)
        idx_ord = np.argsort(each_transfer)
        rixs[idx, :] = np.interp(dE, each_transfer[idx_ord], herfd[idx_ord, idx], left=0, right=0)

    return dE, rixs

plt.figure()
plt.plot(energies_vtc_cubanes, xes_co3mno4_all, 'k-', lw=1, alpha=0.3)
plt.plot(energies_vtc_cubanes, np.mean(xes_co3mno4_all, axis=1), 'r-', lw=2)

# herfd_data, xes_data, energy_in = parse_rixs_scan(db, '/nsls2/xf08id/users/2021/1/308190/rixs_uids_ebce49.txt')
# energy_transfer, rixs = convert_rixs_to_energy_transfer(energy_in, energies_emission[:herfd_data.shape[1]], herfd_data.T)
# plt.figure()
# plt.figure(); plt.contourf(energy_in, energy_transfer, (rixs/rixs.max()).T, 251, vmin=0.0, vmax=0.15, cmap='jet')
# plt.axis('image')
# plt.xlim(7706, 7715); plt.ylim(56, 65)






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

