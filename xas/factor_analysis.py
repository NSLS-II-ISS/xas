import numpy as np
import matplotlib.pyplot as plt


def svd_analysis(energy, mus, emin=None, emax=None, n_cmp=None, ac_thresh=0.8, plotting=False):
    e_mask = np.ones(energy.size, dtype=bool)
    if emin is not None: e_mask &= (energy >= emin)
    if emax is not None: e_mask &= (energy <= emax)

    mus_masked = mus[e_mask, :]
    u, s, v = np.linalg.svd(mus_masked)
    ac_u = np.sum(u[1:, :] * u[:-1, :], axis=0)
    ac_v = np.sum(v.T[1:, :] * v.T[:-1, :], axis=0)

    if n_cmp is None:
        n_cmp = np.min([np.sum(ac_u > ac_thresh), np.sum(ac_v > ac_thresh)])

    if plotting:
        plt.figure()
        plt.subplot(321)
        plt.plot(energy[e_mask], mus_masked)
        plt.xlabel('energy, eV')
        plt.ylabel('mu norm')
        plt.title('Data')

        plt.subplot(323)
        plt.semilogy(s, 'ks-')
        plt.xlabel('component number')
        plt.ylabel('singular value')
        plt.title('Singular values')

        plt.subplot(325)
        plt.plot(ac_u, 'rs-', label='AC_U')
        plt.plot(ac_v, 'bs-', label='AC_V')
        plt.xlabel('component number')
        plt.ylabel('autocorrelaion value')
        plt.title('Autocorrelation')
        plt.legend()
        plt.xlim(0, n_cmp+5)

        plt.subplot(222)
        shift_y = np.tile(np.arange(n_cmp+3)[None, :], (np.sum(e_mask), 1))*0.3
        plt.plot(energy[e_mask], u[:, :(n_cmp+3)] - shift_y)
        plt.xlabel('energy, eV')
        plt.ylabel('LSV')
        plt.title('U -Left Singular Vectors')

        plt.subplot(224)
        shift_y = np.tile(np.arange(n_cmp + 3)[None, :], (mus.shape[1], 1))*0.3
        plt.plot(v.T[:, :(n_cmp + 3)] - shift_y)
        plt.xlabel('spectrum number')
        plt.ylabel('RSV')
        plt.title('V - Right Singular Vectors')


    return u, s, v, n_cmp


def evolving_svd_analysis(energy, mus, emin=None, emax=None, plotting=False, n_cmp=None):
    n_energies, n_curves = mus.shape
    all_cmp = np.min([n_energies, n_curves])
    ss_forward = np.zeros((all_cmp, n_curves))
    for i in range(n_curves):
        _, s, _, _n_cmp = svd_analysis(energy, mus[:, :i+1], emin=emin, emax=emax)
        ss_forward[:s.size, i] = s

    ss_backward = np.zeros((all_cmp, n_curves))
    for i in range(n_curves-1, -1, -1):
        _, s, _, _ = svd_analysis(energy, mus[:, i:], emin=emin, emax=emax)
        ss_backward[:s.size, i] = s

    if n_cmp is None:
        n_cmp = _n_cmp
    elif n_cmp == 'all':
        n_cmp = all_cmp

    if plotting:
        plt.figure()
        plt.subplot(211)
        plt.semilogy(ss_forward.T[:, :n_cmp])
        plt.semilogy(ss_forward.T[:, n_cmp:(n_cmp+3)], 'k--')
        plt.xlabel('number fof curves')
        plt.title('Forward')

        plt.subplot(212)
        plt.semilogy(ss_backward.T[:, :n_cmp])
        plt.semilogy(ss_backward.T[:, n_cmp:(n_cmp + 3)], 'k--')
        plt.xlabel('number fof curves')
        plt.title('Backward')

    return ss_forward, ss_backward



def script_for_elis_slide(energy, mus, emin=9600, emax=9780):
    fpath = '/nsls2/xf08id/Sandbox/2021_analysis_slides/'
    image_base_name = 'image_v1_'
    n_curves = mus.shape[1]
    e_mask = np.ones(energy.size, dtype=bool)
    if emin is not None: e_mask &= (energy >= emin)
    if emax is not None: e_mask &= (energy <= emax)

    font = {'size': 7}

    matplotlib.rc('font', **font)

    for i in range(4, n_curves):
    # for i in [30, 50, 80, n_curves]:
        mus_subset = mus[:, :i]
        u, s, v,_ = svd_analysis(energy, mus_subset, emin=emin, emax=emax, plotting=False)

        plt.figure(1, figsize=(18/2.54, 6/2.54))
        plt.clf()

        plt.subplot(141)
        plt.plot(energy[e_mask], mus_subset[e_mask, :])
        plt.xlim(emin, emax)
        plt.ylim(-0.1, 1.7)

        plt.subplot(142)
        plt.semilogy(s, 'k.-')
        plt.ylim(1e-3, 1e3)

        plt.subplot(143)
        shift_y = np.tile(np.arange(4)[None, :], (np.sum(e_mask), 1)) * 0.3
        plt.plot(energy[e_mask], u[:, :4] - shift_y)
        plt.xlim(emin, emax)
        plt.ylim(-1.2, 0.1)

        plt.subplot(144)
        shift_y = np.tile(np.arange(4)[None, :], (mus_subset.shape[1], 1)) * 0.3
        plt.plot(v.T[:, :4]-shift_y)
        plt.ylim(-1.2, 0.1)
        #plt.xlim(0, n_curves)
        plt.tight_layout()
        plt.savefig(fpath + image_base_name + str(i) + '.png', dpi=300)



#
# def process_rixs_von_hamos(db, uid, roi, subtract_bkg=False):
#    t = db[uid].table(fill=True)
#    try:
#        energies_in = t['hhm_energy'].values
#    except:
#        energies_in = None
#    images_pandas = t['pil100k_image']
#    images = np.array([im[0, :, :] for im in images_pandas])
#    spectra = []
#    for image in images:
#        spectrum = process_image_von_hamos(image, roi, subtract_bkg=subtract_bkg)
#        spectra.append(spectrum)
#    spectra = np.array(spectra)
#    return energies_in, spectra
#
#
# def process_image_von_hamos(image, roi, subtract_bkg=False):
#    roi_x1, roi_x2, roi_y1, roi_y2, int_axis = roi
#    image_roi = image[roi_x1: roi_x2, roi_y1: roi_y2]
#    if subtract_bkg:
#        if int_axis == 0:
#            bkg_roi = image[roi_x1: roi_x2, roi_y1 + roi_y2: roi_y2 + roi_y2]
#        if int_axis == 0:
#            bkg_roi = image[roi_x1 + roi_x2: roi_x2 + roi_x2, roi_y1: roi_y2]
#        bkg = np.sum(bkg_roi, axis=int_axis)
#        bkg = savgol_filter(bkg, 21, 3)
#    else:
#        bkg = None
#
#    spectrum = np.sum(image_roi, axis=int_axis)
#    if bkg is not None:
#        spectrum = spectrum - bkg
#
#    return spectrum
#
#
# # VTC
# roi_x1, roi_x2, roi_y1, roi_y2, int_axis = 10, 30, 0, 300, 0
# roi = (roi_x1, roi_x2, roi_y1, roi_y2, int_axis)
# uid_ti_foil = 'b9008469-e6e7-4345-8872-5b632a90e29a'
# _, spectra = process_rixs_von_hamos(db, uid_ti_foil, roi, subtract_bkg=True)
# spectrum_av_Ti = np.mean(spectra, axis=0)
# spectrum_av_Ti /= spectrum_av_Ti.max()
#
# uid_tin = 'ead8b8e4-3171-4bbc-ab34-8dfff3bcac8c'
# _, spectra = process_rixs_von_hamos(db, uid_tin, roi, subtract_bkg=True)
# spectrum_av_TiN = np.mean(spectra, axis=0)
# spectrum_av_TiN /= spectrum_av_TiN.max()
#
# uid_fetio3 = '2dd621a7-8690-4efc-82ab-0b18fbeb1fa8'
# _, spectra = process_rixs_von_hamos(db, uid_fetio3, roi, subtract_bkg=True)
# spectrum_av_FeTiO3 = np.mean(spectra, axis=0)
# spectrum_av_FeTiO3 /= spectrum_av_FeTiO3.max()
#
# uid_liti2o3 = 'b0aa6acb-a341-4bef-8b07-631009a610c3'
# _, spectra = process_rixs_von_hamos(db, uid_liti2o3, roi, subtract_bkg=True)
# spectrum_av_LiTi2O3 = np.mean(spectra, axis=0)
# spectrum_av_LiTi2O3 /= spectrum_av_LiTi2O3.max()
#
#
# fig, ax = plt.subplots(1,1)
# ax.plot(spectrum_av_Ti, label='Ti')
# ax.plot(spectrum_av_TiN, label='TiN')
# ax.plot(spectrum_av_FeTiO3, label='FeTiO3')
# ax.plot(spectrum_av_LiTi2O3, label='LiTi2O3')
# ax.legend()
#
#
#
# e = ((205 - energy)/(206-123)*(4960-4932) + 4932)
#
#
# fig, ax = plt.subplots(1,1, figsize=(8/2.55, 6/2.55))
# ax.plot(e, spectrum_av_Ti, 'k-', label='Ti foil')
# ax.plot(e, spectrum_av_TiN, 'b-', label='TiN')
# ax.plot(e, spectrum_av_FeTiO3, 'r-', label='FeTiO3')
#
# ax.plot(e[e>4940], spectrum_av_Ti[e>4940]*20, 'k--')
# ax.plot(e[e>4940], spectrum_av_TiN[e>4940]*20, 'b--')
# ax.plot(e[e>4940], spectrum_av_FeTiO3[e>4940]*20, 'r--')
#
# ax.set_xlim(4915, 4975)
# ax.set_ylim(-0.05, 1.05)
# ax.text(4965, 0.5, 'x20', va='center', ha='center')
#
# ax.vlines(4940, -0.05, 1.05, linestyles='--', colors='k')
# ax.legend(frameon=False)
# ax.set_xlabel('Emission energy')
# ax.set_ylabel('Intensity')
# plt.tight_layout()
#
