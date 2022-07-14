

import numpy as np
import pyFAI as pyFAI





def get_ai(dist=40, center_ver=93, center_hor=440,
           pixel1=0.172, pixel2=0.172,
           rot1=0, rot2=0, rot3=np.pi/2,
           energy=11300):
    wavelength = 12.3984 / energy * 1e-3 * 1e-10
    ai = pyFAI.AzimuthalIntegrator(dist=dist*1e-3,
                                   poni1=center_ver * pixel1 * 1e-3, poni2=center_hor * pixel2 * 1e-3,
                                   pixel1=pixel1 * 1e-3, pixel2=pixel2 * 1e-3,
                                   rot1=rot1, rot2=rot2, rot3=rot3,
                                   wavelength=wavelength)
    return ai

ai = get_ai()

image = np.array(data['pil100k_image'][0])
mask = image < np.percentile(image.ravel(), 5)
#
# # res = ai.integrate1d_ng(np.array(data_list[0]),
# #                         1000,
# #                         mask=mask,
# #                         unit="2th_deg")
# # jupyter.plot1d(res)
#
res2 = ai.integrate2d_ng(image,
                         300, 360,
                         mask=mask,
                         unit="r_mm")

jupyter.plot2d(res2)



# hdr = db['928a1184-5ea7-4eb3-aa44-3b5b725ce1b3']# Ir dimer sample
# hdr = db['c723a571-e3e2-4d15-9b36-cf573bf85834'] # neat MeCN
# pil100k_data = hdr.table(stream_name='pil100k_stream', fill=True)
#
#
# images_raw = pil100k_data['pil100k_image'][1]




from scipy.signal import medfilt2d
def process_image(image, nw=5):
    a, b = image.shape
    image_smooth = medfilt2d(image, kernel_size=nw)
    dw = int((nw - 1) / 2)
    thresh = np.percentile(image.ravel(), 5)
    I = np.zeros(image.shape)
    sigma = np.zeros(image.shape)
    gains = np.zeros(image.shape)
    for i in range(dw, a - dw):
        for j in range(dw, b - dw):
            subset = image[(i - dw) : (i + dw + 1), (j - dw) : (j + dw + 1)]
            subset_smooth = image_smooth[(i - dw): (i + dw + 1), (j - dw): (j + dw + 1)]
            # print(subset.shape)
            subset = subset.ravel()
            subset_smooth = subset_smooth.ravel()
            subset_smooth = subset_smooth[subset > thresh]
            subset = subset[subset > thresh]
            _I = np.mean(subset)
            _sigma = np.std(subset - subset_smooth)
            I[i, j] = _I
            sigma[i, j] = _sigma
            gains[i, j] = _sigma**2 / _I

    plt.figure(8)
    plt.clf()
    plt.subplot(221)
    plt.plot(I.ravel(), sigma.ravel()**2, 'k.')
    plt.plot([0, 20e3], [0, 20e3], 'r-')

    plt.subplot(222)
    plt.imshow(I)

    plt.subplot(223)
    plt.imshow(sigma**2)

    plt.subplot(224)
    plt.imshow(gains, vmin=0.9, vmax=1.1)

process_image(image)


def integrate_pil100k_image_stack(images_array, dist=40, center_ver=93, center_hor=440, deadtime_cor=False):
    ai = get_ai(dist=dist, center_ver=center_ver, center_hor=center_hor)
    s = []
    mask = None
    for image in images_array:
        if deadtime_cor:
            image *= np.exp(-image * 25 * 160e-9)
        if mask is None:
            mask = image < np.percentile(image.ravel(), 5)
        res = ai.integrate1d_ng(image,
                                350,
                                mask=mask,
                                unit="2th_deg")
        s.append(res[1])
        tth = res[0]
    return tth, np.array(s)


fname_base = 'Ir sample 2 scan AXS wide 600 um'
folder = '/nsls2/data/iss/legacy/processed/2022/2/300011/'

from xas.file_io import load_binned_df_from_file, load_extended_data_from_file
from xas.file_io import load_interpolated_df_and_extended_data_from_file


df = None
data = None
for i in range(1, 16):
    f = f'{folder}{fname_base} {i:04d}-r0002.dat'
    df_i, ext_data_i, _ = load_interpolated_df_and_extended_data_from_file(f)
    if df is None:
        df = df_i
        data = ext_data_i['pil100k_image']
    else:
        df += df_i
        data += ext_data_i['pil100k_image']

df = df/15
data = data/15



# df, _ = load_binned_df_from_file(f'{folder}{fname}.dat')
# data = load_extended_data_from_file(f'{folder}/extended_data/{fname}.h5')['pil100k_image']


# df, _ = load_binned_df_from_file('/nsls2/data/iss/legacy/processed/2022/2/300010/Neat MeCN AXS 0019-r0013.dat')
# data = pd.read_json('/nsls2/data/iss/legacy/processed/2022/2/300010/extended_data/Neat MeCN AXS 0019-r0013.json')

tth, s = integrate_pil100k_image_stack(data)
i0s = -df['i0']
i0s /= i0s[50]
s /= i0s[:, None]
qq = 4*np.pi / (12.3984 / 11.220) * np.sin(np.deg2rad(tth)/2)
iff = xview_gui.project[-1].flat
# iff = df['iff']/df['i0']
s = s.T
energies = df['energy'].values
# s /= s.max()

def preproc_s(s, n_lo=5, n_hi=5, n_fit=1):
    s_lo = s[:, :n_lo]
    u, _, _ = np.linalg.svd(s_lo)
    #
    basis = u[:, :n_fit]
    c, _, _, _ = np.linalg.lstsq(basis, s)

    s_b = basis @ c
    s_preproc = s - s_b

    nc=3

    plt.figure(5, clear=True)
    plt.subplot(221)
    # plt.plot(u[:, :n_fit])
    plt.imshow(s_preproc)

    plt.subplot(222)
    plt.plot(s_preproc[:, 161] * 25)

    plt.subplot(223)
    plt.plot(s_preproc[:, 300] * 25)

    plt.subplot(224)
    plt.plot(v[:, :nc])
    return s_preproc

s_preproc = preproc_s(s)

def rm_fluorescence(s, iff):
    c, _, _, _ = np.linalg.lstsq(iff[:, None], s.T)
    # c[c<0] = 0
    s_b = (iff[:, None] @ c).T
    plt.figure(6, clear=True)
    plt.subplot(221)
    plt.imshow(s)

    plt.subplot(222)
    plt.imshow(s - s_b)

    return s - s_b

s_preproc = rm_fluorescence(s_preproc, iff)


def process_stack_of_patterns(energy, s, iff, emin=11000, emax=11270, e0=11217):
    e_lo_mask = energy <= emin
    e_hi_mask = energy >= emax
    e_mask = e_lo_mask | e_hi_mask
    poly = np.polyfit(energy[e_lo_mask], s[e_lo_mask, :], 1)
    s_bkg = np.zeros(s.shape)
    s_proc = np.zeros(s.shape)
    # poly_iff = np.polyfit(energy[e_lo_mask], s[e_lo_mask, :], 1)

    basis = np.vstack((np.ones(energy.size), energy, energy**2, iff, iff*energy)).T
    print(basis.shape, s[:, 0].shape)
    for i in range(s.shape[1]):
        # preedge = np.polyval(poly[:, i], energy)
        # p_postedge = np.polyfit(energy[e_hi_mask], (s[:, i] - preedge)[, ], 1)
        # postedge = np.polyval(p_postedge, energy)
        # bkg_loc = preedge.copy()
        # bkg_loc[energy >=e0] = (preedge + postedge)[energy >=e0]
        # iff_bkg_loc = np.zeros(energy.size)
        # iff_bkg_loc[energy >=e0] = 1
        # basis = np.vstack((iff, energy)).T
        # c, _, _, _ = np.linalg.lstsq( iff[e_hi_mask, None], (s[:, i] - preedge)[e_hi_mask])
        # s_bkg[:, i] = preedge + c * iff
        c, _, _, _ = np.linalg.lstsq(basis[e_mask], s[e_mask, i] )
        s_bkg[:, i] = basis @ c
        # s_proc[:, i] = (s[:, i] - preedge)/postedge

        # s_bkg[:, i] = basis @ c + np.polyval(poly[:, i], energy)

    s_proc = s - s_bkg

    # s_proc /= np.mean(s_proc[e_hi_mask, :], axis=0)
    return s_proc, s_bkg



s_proc, s_bkg = process_stack_of_patterns(df['energy'].values, s.T, iff)

plt.figure(2, clear=True)

plt.subplot(221)
plt.imshow(np.log10(s.T))

plt.subplot(222)
plt.imshow(s_proc, vmin=-20, vmax=20)


plt.subplot(223)
# plt.plot(energies, s.T[:, 50])
# plt.plot(energies, s_bkg[:, 50])

plt.plot(energies, s_proc[:, 50])


plt.subplot(224)
nc=[130, 160]
# plt.plot(qq, s.T[160, :])
# plt.plot(qq, s_bkg[160, :])
# for i in nc:
#     plt.plot(qq, s_proc[i, :])

plt.plot(qq, np.mean(s_proc[125:150, :], axis=0))
plt.plot(qq, np.mean(s_proc[160:200, :], axis=0))

from scipy.optimize import nnls

def fit_stack_of_patterns(energy, s, iff, ff0, ff_real, ff_imag):

    s_fit = np.zeros(s.shape)
    s_bkg = np.zeros(s.shape)
    anomalous_scat = np.zeros((2, s.shape[1]))
    # c, _, _, _ = np.linalg.lstsq(basis, s)
    # s_fit = basis @ c
    # c_bkg = c.copy()
    # c_bkg[-2:, :] = 0
    # s_bkg = basis @ c_bkg
    # print(c.shape)
    for i in range(s.shape[1]):
        ff_lin = ff_real
        ff_quad = ff_real ** 2 + ff_imag ** 2
        ff_lin_use = ff_lin - ff_lin[0]
        ff_quad_use = ff_quad - ff_quad[0]

        basis = np.vstack((np.ones(energy.size), energy, energy ** 2,
                           ff_lin_use, ff_quad_use)).T

        c, _, _, _ = np.linalg.lstsq(basis, s[:, i])
        # c, _ = nnls(basis, s[:, i])
        s_fit[:, i] = basis @ c
        c_bkg = c.copy()
        c_bkg[-2:] = 0
        s_bkg[:, i] = basis @ c_bkg
        anomalous_scat[:, i] = c[-2:]
    return s_fit, s_bkg, anomalous_scat

import kkcalc


def kk_calculate_real_from_array(energy, mu_flat, element, merge_points=None, add_background=True, fix_distortions=False, curve_tolerance=None, curve_recursion=50):
    """Do all data loading and processing and then calculate the kramers-Kronig transform.
    Parameters
    ----------
    NearEdgeDataFile : string
    	Path to file containg near-edge data
    ChemicalFormula : string
    	A standard chemical formula string consisting of element symbols, numbers and parentheses.
    merge_points : list or tuple pair of `float` values, or None
    	The photon energy values (low, high) at which the near-edge and scattering factor data values
    	are set equal so as to ensure continuity of the merged data set.
    Returns
    -------
    This function returns a numpy array with columns consisting of the photon energy, the real and the imaginary parts of the scattering factors.
    """
    Stoichiometry = kkcalc.data.ParseChemicalFormula(element)
    Relativistic_Correction = kkcalc.kk.calc_relativistic_correction(Stoichiometry)
    Full_E, Imaginary_Spectrum = kkcalc.data.calculate_asf(Stoichiometry)
    _data = np.vstack((energy, mu_flat)).T
    NearEdge_Data = kkcalc.data.convert_data(_data, FromType='xanes', ToType='asf')
    Full_E, Imaginary_Spectrum = kkcalc.data.merge_spectra(NearEdge_Data, Full_E, Imaginary_Spectrum, merge_points=merge_points, add_background=add_background, fix_distortions=fix_distortions)
    Real_Spectrum = kkcalc.kk.KK_PP(Full_E, Full_E, Imaginary_Spectrum, Relativistic_Correction)
    if curve_tolerance is not None:
        output_data = kkcalc.kk.improve_accuracy(Full_E, Real_Spectrum, Imaginary_Spectrum, Relativistic_Correction, curve_tolerance, curve_recursion)
    else:
        Imaginary_Spectrum_Values = kkcalc.data.coeffs_to_ASF(Full_E, np.vstack((Imaginary_Spectrum, Imaginary_Spectrum[-1])))
        output_data = np.vstack((Full_E, Real_Spectrum, Imaginary_Spectrum_Values)).T
    return output_data


bla = kk_calculate_real_from_array(energies, iff, 'Ir')

import xraydb
ff_data = kkcalc.kk.kk_calculate_real('/nsls2/data/iss/legacy/processed/2022/2/300010/Ir sample 2 AXS try4 cont 0037-r0002 iff-i0.mu',
                               'Ir', input_data_type='xanes')
ff_real = np.interp(df['energy'], ff_data[:, 0], ff_data[:, 1])
ff_imag = np.interp(df['energy'], ff_data[:, 0], ff_data[:, 2])
ff0 = xraydb.f0('Ir', qq)




# s_proc, s_bkg = process_stack_of_patterns(df['energy'], s, iff)
s_fit, s_bkg, anomalous_scat = fit_stack_of_patterns(df['energy'].values, s_preproc, iff, ff0, ff_real, ff_imag)
s_proc = s_preproc - s_bkg

plt.figure(1, clear=True);

plt.subplot(221)
# plt.plot(s.T)
plt.imshow(s_preproc - s_bkg)
# plt.plot((s - s_bkg).T)

nc = [100, 200]
offset = -200
plt.subplot(222)
for i, n in enumerate(nc):
    plt.plot(df['energy'], s_preproc[:, n] - i * offset, 'k-')
    # plt.plot(df['energy'], s_bkg[:, n] - i * offset)
    plt.plot(df['energy'], s_fit[:, n] - i * offset, 'r-')
    plt.plot(df['energy'], s_bkg[:, n] - i * offset, 'm--')

# plt.plot(df['energy'], ff_real/70)
# plt.plot(df['energy'], iff/10 + 1)

# plt.plot(s[:, 600])
# plt.plot(s_bkg[:, 600])

plt.subplot(223)
plt.plot(df['energy'], s_proc[:, 200])
plt.plot(df['energy'], s_proc[:, 70])

plt.subplot(224)
plt.plot(qq, anomalous_scat[0, :])
plt.plot(qq, anomalous_scat[1, :])

# plt.plot(s_proc[130, :])
# plt.plot(df['iff']/df['i0'])
# plt.plot(df['energy'], df['iff']/df['i0'])








