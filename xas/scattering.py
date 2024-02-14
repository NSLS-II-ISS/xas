import matplotlib.pyplot as plt
import numpy as np
import pyFAI as pyFAI
import kkcalc, kkcalc.data, kkcalc.kk

from xas.xasproject import XASDataSet


def get_ai(dist=40, center_ver=37., center_hor=45.,
           pixel1=0.172, pixel2=0.172,
           rot1=0, rot2=0, rot3=np.pi/2,
           energy=24350):
    wavelength = 12.3984 / energy * 1e-3 * 1e-10
    ai = pyFAI.AzimuthalIntegrator(dist=dist*1e-3,
                                   poni1=center_ver * pixel1 * 1e-3, poni2=center_hor * pixel2 * 1e-3,
                                   pixel1=pixel1 * 1e-3, pixel2=pixel2 * 1e-3,
                                   rot1=rot1, rot2=rot2, rot3=rot3,
                                   wavelength=wavelength)
    return ai

ai = get_ai()
#
res2 = ai.integrate2d_ng(df.pil100k_image[0],
                         300, 360, mask=mask,
                         unit="r_mm")
fig, ax = plt.subplots(1, clear=True, num=1)
from pyFAI.gui import jupyter
jupyter.plot2d(res2, ax=ax)


def _generate_mask(img, log_thresh_lo=None, log_thresh_hi=None, mask_edges=True):
    mask = np.zeros(img.shape, dtype=bool)
    mask[np.isnan(np.log(img))] = 1
    if log_thresh_lo is not None: mask[np.log10(img) <= log_thresh_lo] = 1
    if log_thresh_hi is not None: mask[np.log10(img) >= log_thresh_hi] = 1
    if mask_edges:
        mask[0, :] = 1
        mask[-1, :] = 1
        mask[:, 0] = 1
        mask[:, -1] = 1
    return mask

def generate_mask(df, image_key='pil100k_image', **kwargs):
    img = df[image_key].mean()
    return _generate_mask(img, **kwargs)






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

# process_image(image)


def integrate_pil100k_image_stack(images_array, energy, dist=40, center_ver=93, center_hor=440, deadtime_cor=False, mask=None, npts=425,
                                  **kwargs):
    ai = get_ai(dist=dist, center_ver=center_ver, center_hor=center_hor, **kwargs)
    s = []
    # mask = None
    for image in images_array:
        if deadtime_cor:
            image *= np.exp(-image * 25 * 160e-9)
        if mask is None:
            mask = image < np.percentile(image.ravel(), 5)
        res = ai.integrate1d_ng(image,
                                npts,
                                mask=mask,
                                unit="2th_deg")
        s.append(res[1])
        tth = res[0]

    q_all = 4 * np.pi / (12398.4 / energy[np.newaxis, :]) * np.sin(np.deg2rad(tth[:, np.newaxis])/2)
    qmin = q_all.min(axis=0).max()
    qmax = q_all.max(axis=0).min()
    nq = np.sum((q_all >= qmin) & (q_all <= qmax), axis=0).min()
    q = np.linspace(qmin, qmax, nq)
    sq = np.zeros((q.size, energy.size))
    for i in range(energy.size):
        sq[:, i] = np.interp(q, q_all[:, i], s[i])

    return q, sq

# reduced
# df =  get_processed_df_from_uid('abaa4995-9bd1-4a5b-9b2d-4338ea30ab12', db, logger=None, draw_func_interp=None, draw_func_bin=None,
#                                 print_func=None, save_interpolated_file=False, return_processed_df=True, load_images=True)

from xas.file_io import load_binned_df_and_extended_data_from_file
# uids_red = ['d6b89a0d-1cf5-4b0b-8695-bc84cb5007e5',
#             '318cec96-5213-4773-baf2-da50cae02f2d',
#             '138e62df-c191-4831-9249-a5d252e8b6e8',
#             '4989c59b-76a0-46dd-8723-51fd6f1ff6e5',
#             'abaa4995-9bd1-4a5b-9b2d-4338ea30ab12',]
uids_red = ['b5d1704b-8c99-4f71-884f-14204166ab7b', # oxide
           '119523cd-0aae-4964-af3d-140458a51bcd',
           '4f94bb97-166e-4194-b145-bab39479bbf1',
           '23bfef24-7602-4a65-b761-cd9617152716',
           '8e71ee10-c9a3-4c49-9c59-11f169ca94d9' ]

dfs = []
for uid in uids_red:
    path_to_file = db[uid].start['interp_filename'][:-3] + 'dat'
    _df, _ext_data, _ = load_binned_df_and_extended_data_from_file(path_to_file)
    _df['pil100k_image'] = [i for i in _ext_data['pil100k_image']]
    dfs.append(_df)

total_image = np.mean(np.array([np.array([i for i in df['pil100k_image']]) for df in dfs]), axis=(0,1))
plt.figure(1, clear=True); plt.imshow(np.log10(total_image), vmin=3, vmax=4.9)

mask = _generate_mask(total_image, log_thresh_lo=3.05, log_thresh_hi=4.9)
mask[:70, :70] = 1
mask[:120, 340:] = 1
plt.figure(2, clear=True); plt.imshow(mask)

ai = get_ai(center_ver=39., center_hor=45.)
#
res2 = ai.integrate2d_ng(total_image,
                         300, 360, mask=mask,
                         unit="r_mm")
fig, ax = plt.subplots(1, clear=True, num=3)
from pyFAI.gui import jupyter
jupyter.plot2d(res2, ax=ax)

plt.figure(4, clear=True)
plt.plot(res2[0].T / (np.mean(res2[0].T[65:75, :], axis=0))[None, :])
plt.xlim(60, 90)
plt.ylim(0.6, 1.2)


energy_offset = 3.75
# dist=43
for df in dfs: #[:1]:
    df['i0_norm'] = df.i0 / df.i0.median()
    df['mu'] = df['iff'] / df['i0']
    q, sq = integrate_pil100k_image_stack(df.pil100k_image / (df['i0_norm']), df.energy.values + energy_offset,
                                          dist=43, center_ver=39, center_hor=45, mask=mask, npts=425,
                                          rot1=np.deg2rad(0.00), rot2=np.deg2rad(-0.1))
    df['sq'] = sq.T.tolist()

energy = df.energy.values + energy_offset

# qind1 = 350
# qind2 = 370

qind1 = 140
qind2 = 170

q_max = np.array([q[qind1:qind2][i] for i in np.argmax(sq[qind1:qind2, :], axis=0)])
peak_max = np.max(sq[qind1:qind2, :], axis=0) / np.mean(sq[qind1:qind2, :], axis=0)

plt.figure(5, clear=True)
plt.subplot(211)
plt.plot(sq)

plt.subplot(212)
plt.plot(q, sq / np.mean(sq[190:200, :], axis=0)[None, :])
# plt.plot(q_max, peak_max, 'k.-')
plt.xlim(q[120], q[370])
# plt.ylim(3000, 15000)
plt.ylim(0.5, 1.75)

plt.figure(6, clear=True)
plt.plot(energy, q_max)
plt.title(q_max.max() - q_max.min())

# fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, clear=True, num=5)
# ax.plot_surface(np.tile(q[:, None], (1, energy.size)),
#                 np.tile(energy[None, :], (q.size, 1)),
#                 sq,
#                 cmap='coolwarm',
#                 linewidth=0, antialiased=False)
#
# ax.set_xlim(2.5, 4)
# ax.set_zlim(10000, 15000)

plt.figure(7, clear=True)
# plt.contourf(energy, q, sq, 50, vmin=5e3, vmax=40e3)
plt.plot(energy, sq[[153, 220], :].T)

mu = np.mean(np.array([(df['iff'] / df['i0']).values for df in dfs]), axis=0)
sq = np.mean(np.array([np.array([i for i in df['sq']]) for df in dfs]), axis=0).T
imgs = np.mean(np.array([np.array([i for i in df['pil100k_image']]) for df in dfs]), axis=0)

def flatten_mu(energy, mu):
    ds = XASDataSet(mu=mu, energy=energy)
    ds.flatten()
    return ds.flat

mu_flat = flatten_mu(energy, mu)

# energy_mask = ((energy <=11200) & (energy <= 11500)) | ((energy >=11600) & (energy <= 12000))
energy_mask =  (energy <= 11200) | (energy >=11800)
energy_norm = energy/np.mean(energy)
A = np.vstack((energy_norm**0, energy_norm**1, energy_norm**2,
               mu_flat * energy_norm**0, mu_flat * energy_norm**1, mu_flat * energy_norm**2)).T
# A = np.vstack((np.ones(energy.size), energy/np.mean(energy), mu_flat)).T
A_masked = A[energy_mask, :]
p, _, _, _ = np.linalg.lstsq(A_masked, sq[:, energy_mask].T, rcond=-1)

sq_b = (A @ p).T

plt.figure(8, clear=True)
plt.subplot(211)
plt.plot(energy, sq[[155, 230], :].T)
plt.plot(energy, sq_b[[155, 230], :].T, 'k-')

plt.subplot(212)
plt.plot(q, p[-1, :])
plt.hlines([np.median(p[-1, q>=2.5])], 2.5, 5.7, colors='r')

q_sel = (q>=3.4) & (q<=4.44)
plt.figure(9, clear=True)
plt.subplot(211)
plt.plot(energy, sq[q_sel, :].T)
plt.plot(energy, sq_b[q_sel, :].T, 'k-')

plt.subplot(212)
plt.plot(energy, sq[q_sel, :].T - sq_b[q_sel, :].T)

f0 = 78
f1 = xraydb.f1_chantler('Pt', energy)
f2 = xraydb.f2_chantler('Pt', energy)
f_lin = f0 + f1
f_quad = f_lin**2 + f2**2

from scipy.optimize import nnls
basis = np.vstack((np.ones(energy.size), f_lin, f_quad)).T
sq_fit = np.zeros(sq.shape)

C = (1 - 0.88 * (energy / np.mean(energy)))

energy_mask = energy <= 11500
n_sel = [153, 220]
# p, _, _, _ = np.linalg.lstsq(basis[energy_mask, :], sq[n_sel, :][:, energy_mask].T, rcond=-1)
p = np.zeros((basis.shape[1], q.size))
for i in n_sel:
    p[:, i], _ = nnls(basis[energy_mask, :], sq[i, energy_mask].T * C[energy_mask])

sq_fit = (basis @ p).T


# C /= np.mean(C)

plt.figure(12, clear=True)
plt.plot(energy, sq[[153, 220], :].T * C[:, None])
plt.plot(energy, sq_fit[[153, 220], :].T, 'k-')

# plt.plot(energy, (f0 + f1))
# plt.plot(energy, (f0 + f1)**2 + f2**2)







# uids = list(range(12594, 12604))
uids = [ 'abaa4995-9bd1-4a5b-9b2d-4338ea30ab12',
         '4989c59b-76a0-46dd-8723-51fd6f1ff6e5',
         '138e62df-c191-4831-9249-a5d252e8b6e8',
         '318cec96-5213-4773-baf2-da50cae02f2d',
         'd6b89a0d-1cf5-4b0b-8695-bc84cb5007e5',
         'a257ea3f-0087-4fa2-9bdb-ec2b5d2008ce',
         'afcb00dc-551d-4fc1-b50c-3da84ca85c1c',
         'd07d25a0-add8-4b1b-acc0-eb26c31f972d',
         '72ef63d1-77f1-4730-88c0-eb17334ca0b4',
         '1df48476-96b6-4b66-a1ca-79704b6ba00e']

dfs = []
sq_list = []
for uid in uids:
     df =  get_processed_df_from_uid(uid, db, logger=None, draw_func_interp=None, draw_func_bin=None,
                                                                                print_func=None, save_interpolated_file=False, return_processed_df=True,
                                                                                load_images=True)
     q, _sq = integrate_pil100k_image_stack(df.pil100k_image / (df.i0 / df.i0.median()), df.energy.values, dist=50, center_ver=37, center_hor=45, mask=mask)
     sq_list.append(_sq)
     dfs.append(df)

sq_av = np.mean(np.array(sq_list), axis=0)

plt.figure(1, clear=True)
plt.imshow(sq_av)

plt.figure(2, clear=True)
plt.contourf(df.energy, q, sq_av, 50)

plt.figure(3, clear=True)
plt.plot(df.energy, sq_av[260:261, :].T)
plt.plot(df.energy, sq_av[400:401, :].T)
# plt.plot(q[250:275], sq_av[250:275, ::50])

# fname_base = 'Ir sample 2 scan AXS wide 600 um'

# reprocess the data
# process files for sample 2 with AXS wide trajectory
# from xas.process import process_interpolate_bin_from_uid
# # for i in range(255080, 255104+1):
# for i in range(255114, 255313 + 1):
#     print(i)
#     process_interpolate_bin_from_uid(i, db)

# from xas.file_io import load_binned_df_and_extended_data_from_file, save_extended_data_as_file
# folder = '/nsls2/data/iss/legacy/processed/2022/2/300011/'
# fname_base = 'Ir sample 1 scan AXS wide 600 um data'



# dist=40
# center_ver=105
# center_hor=438
# for i in range(2, 201):
#     f = f'{folder}{fname_base} {i:04d}-r0003.dat'
#     df_i, ext_data_i, _ = load_binned_df_and_extended_data_from_file(f)
#     q_i, sq_i = integrate_pil100k_image_stack(ext_data_i['pil100k_image'], df_i['energy'], dist=dist, center_ver=center_ver, center_hor=center_hor)
#     int_data_i = {'q' : q_i, 'sq' : sq_i, 'dist' : dist, 'center_ver' : center_ver, 'center_hor' : center_hor}
#     f_int = f'{f[:-4]}_int.dat'
#     save_extended_data_as_file(f_int, int_data_i, data_kind='default', ext_data_path='extended_data')
#
#     # if df is None:
#         df = df_i
#         data = ext_data_i['pil100k_image']

# df, ext_data, _ = load_binned_df_and_extended_data_from_file('/nsls2/data/iss/legacy/processed/2022/2/300011/Ir sample 2 scan AXS wide 600 um cont 0001-r0003.dat')
#
# from xas.file_io import load_binned_df_from_file, load_extended_data_from_file
# from xas.file_io import load_interpolated_df_and_extended_data_from_file
#
#
# df = None
# data = None
# for i in range(1, 16):
#     f = f'{folder}{fname_base} {i:04d}-r0002.dat'
#     df_i, ext_data_i, _ = load_interpolated_df_and_extended_data_from_file(f)
#     if df is None:
#         df = df_i
#         data = ext_data_i['pil100k_image']
#     else:
#         df += df_i
#         data += ext_data_i['pil100k_image']
#
# df = df/15
# data = data/15



# df, _ = load_binned_df_from_file(f'{folder}{fname}.dat')
# data = load_extended_data_from_file(f'{folder}/extended_data/{fname}.h5')['pil100k_image']


# df, _ = load_binned_df_from_file('/nsls2/data/iss/legacy/processed/2022/2/300010/Neat MeCN AXS 0019-r0013.dat')
# data = pd.read_json('/nsls2/data/iss/legacy/processed/2022/2/300010/extended_data/Neat MeCN AXS 0019-r0013.json')

# tth, s = integrate_pil100k_image_stack(data)
# i0s = -df['i0']
# i0s /= i0s[50]
# s /= i0s[:, None]
# qq = 4*np.pi / (12.3984 / 11.220) * np.sin(np.deg2rad(tth)/2)
# iff = xview_gui.project[-1].flat
# # iff = df['iff']/df['i0']
# s = s.T
# energies = df['energy'].values
# # s /= s.max()
#
# def preproc_s(s, n_lo=5, n_hi=5, n_fit=1):
#     s_lo = s[:, :n_lo]
#     u, _, _ = np.linalg.svd(s_lo)
#     #
#     basis = u[:, :n_fit]
#     c, _, _, _ = np.linalg.lstsq(basis, s)
#
#     s_b = basis @ c
#     s_preproc = s - s_b
#
#     nc=3
#
#     plt.figure(5, clear=True)
#     plt.subplot(221)
#     # plt.plot(u[:, :n_fit])
#     plt.imshow(s_preproc)
#
#     plt.subplot(222)
#     plt.plot(s_preproc[:, 161] * 25)
#
#     plt.subplot(223)
#     plt.plot(s_preproc[:, 300] * 25)
#
#     plt.subplot(224)
#     plt.plot(v[:, :nc])
#     return s_preproc
#
# s_preproc = preproc_s(s)
#
# def rm_fluorescence(s, iff):
#     c, _, _, _ = np.linalg.lstsq(iff[:, None], s.T)
#     # c[c<0] = 0
#     s_b = (iff[:, None] @ c).T
#     plt.figure(6, clear=True)
#     plt.subplot(221)
#     plt.imshow(s)
#
#     plt.subplot(222)
#     plt.imshow(s - s_b)
#
#     return s - s_b
#
# s_preproc = rm_fluorescence(s_preproc, iff)


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



# s_proc, s_bkg = process_stack_of_patterns(df['energy'].values, s.T, iff)
#
# plt.figure(2, clear=True)
#
# plt.subplot(221)
# plt.imshow(np.log10(s.T))
#
# plt.subplot(222)
# plt.imshow(s_proc, vmin=-20, vmax=20)
#
#
# plt.subplot(223)
# # plt.plot(energies, s.T[:, 50])
# # plt.plot(energies, s_bkg[:, 50])
#
# plt.plot(energies, s_proc[:, 50])
#
#
# plt.subplot(224)
# nc=[130, 160]
# # plt.plot(qq, s.T[160, :])
# # plt.plot(qq, s_bkg[160, :])
# # for i in nc:
# #     plt.plot(qq, s_proc[i, :])
#
# plt.plot(qq, np.mean(s_proc[125:150, :], axis=0))
# plt.plot(qq, np.mean(s_proc[160:200, :], axis=0))
#
# from scipy.optimize import nnls

def fit_stack_of_patterns(energy, sq, iff, f_real, f_imag):

    # s_fit = np.zeros(sq.shape)
    # s_bkg = np.zeros(sq.shape)
    # anomalous_scat = np.zeros((2, sq.shape[1]))

    f_lin = 2 * f_real
    f_quad = f_real ** 2 + f_imag ** 2
    f_lin_use = f_lin - f_lin[0]
    f_quad_use = f_quad - f_quad[0]

    basis = np.vstack((np.ones(energy.size), energy, energy**2, iff,
                       f_lin_use, f_quad_use)).T

    c, _, _, _ = np.linalg.lstsq(basis, sq.T)
    s_fit = (basis @ c).T
    c_bkg = c.copy()
    c_bkg[-2:, :] = 0
    anomalous_scat = c[-2:, :]
    s_bkg = (basis @ c_bkg).T
    # print(c.shape)
    # for i in range(sq.shape[1]):
    #
    #
    #     basis = np.vstack((np.ones(energy.size), energy, energy ** 2,
    #                        ff_lin_use, ff_quad_use)).T
    #
    #     c, _, _, _ = np.linalg.lstsq(basis, sq[:, i])
    #     # c, _ = nnls(basis, s[:, i])
    #     s_fit[:, i] = basis @ c
    #     c_bkg = c.copy()
    #     c_bkg[-2:] = 0
    #     s_bkg[:, i] = basis @ c_bkg
    #     anomalous_scat[:, i] = c[-2:]
    return s_fit, s_bkg, anomalous_scat




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


def kkcalc_data_calculate_asf(Stoichiometry):
    """Sum scattering factor data for a given chemical stoichiometry.

    Parameters
    ----------
    Stoichiometry : a list of elemental symbol,number pairs

    Returns
    -------
    total_E: 1D numpy array listing the starting photon energies of the segments that the spectrum is broken up into.
    total_Im_coeffs: nx5 numpy array in which each row lists the polynomial coefficients describing the shape of the spectrum in that segment.
    """
    ELEMENT_DATABASE = kkcalc.data.load_Element_Database()
    logger = kkcalc.data.logger
    numpy = np
    logger.info("Calculate material scattering factor data from the given stoichiometry")

    if len(Stoichiometry) is 0:
        logger.error("No elements described by input.")
        return None
    else:
        # get unique energy points
        total_E = numpy.array([])
        for element, n in Stoichiometry:
            total_E = numpy.concatenate((total_E, ELEMENT_DATABASE[str(element)]['E']))
        total_E = numpy.unique(total_E)
        # add weighted asf data sets for KK calculation
        total_Im_coeffs = numpy.zeros((len(total_E) - 1, 5))
        counters = numpy.zeros((len(Stoichiometry)), dtype=numpy.int64)
        for i, E in enumerate(total_E[1:]):
            sum_Im_coeffs = 0
            for j in range(len(counters)):
                sum_Im_coeffs += Stoichiometry[j][1] * ELEMENT_DATABASE[str(Stoichiometry[j][0])]['Im'][counters[j], :]
                counters[j] += ELEMENT_DATABASE[str(Stoichiometry[j][0])]['E'][counters[j] + 1] == E
            total_Im_coeffs[i, :] = sum_Im_coeffs
        return total_E, total_Im_coeffs

Full_E, Imaginary_Spectrum = kkcalc_data_calculate_asf(Stoichiometry)

Stoichiometry = kkcalc.data.ParseChemicalFormula('Pd')
Relativistic_Correction = kkcalc.kk.calc_relativistic_correction(Stoichiometry)
Full_E, Imaginary_Spectrum = kkcalc.data.calculate_asf(Stoichiometry)
_data = np.vstack((df.energy, df.flat)).T
NearEdge_Data = kkcalc.data.convert_data(_data, FromType='xanes', ToType='asf')
Full_E, Imaginary_Spectrum = kkcalc.data.merge_spectra(NearEdge_Data, Full_E, Imaginary_Spectrum, merge_points=None, add_background=True, fix_distortions=False)
Real_Spectrum = kkcalc.kk.KK_PP(Full_E, Full_E, Imaginary_Spectrum, Relativistic_Correction)
Imaginary_Spectrum_Values = kkcalc.data.coeffs_to_ASF(Full_E, np.vstack((Imaginary_Spectrum, Imaginary_Spectrum[-1])))
output_data = np.vstack((Full_E, Real_Spectrum, Imaginary_Spectrum_Values)).T

def compute_f_real_imag(energy, mu_flat, element):
    ff_data = kk_calculate_real_from_array(energy, mu_flat, element)
    ff_energy, ff_real, ff_imag = ff_data.T
    f_real = np.interp(energy, ff_energy, ff_real)
    f_imag = np.interp(energy, ff_energy, ff_imag)
    return f_real, f_imag


def integrate_scattering_dataset(ds, dist=40, center_ver=93, center_hor=440, plotting=True):
    q, sq = integrate_pil100k_image_stack(np.abs(ds.ext_data['pil100k_image']), ds.energy, dist=dist,
                                           center_ver=center_ver, center_hor=center_hor, deadtime_cor=False)
    ds.q = q
    ds.sq = sq

    if plotting:
        plt.figure()
        plt.plot(q, sq)

integrate_scattering_dataset(x[0])

def process_scattering_dataset(ds, dist=40, center_ver=93, center_hor=440, plotting=True):
    element = ds.md['element']
    f_real, f_imag = compute_f_real_imag(ds.energy, ds.flat, element)

    sq_fit, sq_bkg, anomalous_scat = fit_stack_of_patterns(ds.energy, ds.sq, ds.flat, f_real, f_imag)

    if plotting:

        fig, ax1 = plt.subplots(1)
        ax2 = ax1.twinx()
        ax1.plot(ds.energy, f_real, 'k-', label='real')
        ax2.plot(ds.energy, f_imag, 'r-', label='imag')
        ax1.legend()
        ax2.legend()

        plt.figure()
        plt.subplot(231)
        plt.imshow(ds.sq)

        plt.subplot(232)
        plt.imshow(sq_fit)

        plt.subplot(233)
        plt.imshow(ds.sq - sq_bkg)

        plt.subplot(234)
        plt.plot(ds.energy, ds.sq[200, :], 'k-')
        plt.plot(ds.energy, sq_fit[200, :], 'r-')
        plt.plot(ds.energy, sq_bkg[200, :], 'g-')

        plt.subplot(235)
        plt.plot(ds.q, anomalous_scat[0, :])
        plt.plot(ds.q, anomalous_scat[1, :])


# process_scattering_dataset(x[0])
# bla = kk_calculate_real_from_array(energies, iff, 'Ir')
#
# import xraydb
# ff_data = kkcalc.kk.kk_calculate_real('/nsls2/data/iss/legacy/processed/2022/2/300010/Ir sample 2 AXS try4 cont 0037-r0002 iff-i0.mu',
#                                'Ir', input_data_type='xanes')
# ff_real = np.interp(df['energy'], ff_data[:, 0], ff_data[:, 1])
# ff_imag = np.interp(df['energy'], ff_data[:, 0], ff_data[:, 2])
# ff0 = xraydb.f0('Ir', qq)
#
#
#
#
# # s_proc, s_bkg = process_stack_of_patterns(df['energy'], s, iff)
# s_fit, s_bkg, anomalous_scat = fit_stack_of_patterns(df['energy'].values, s_preproc, iff, ff0, ff_real, ff_imag)
# s_proc = s_preproc - s_bkg
#
# plt.figure(1, clear=True);
#
# plt.subplot(221)
# # plt.plot(s.T)
# plt.imshow(s_preproc - s_bkg)
# # plt.plot((s - s_bkg).T)
#
# nc = [100, 200]
# offset = -200
# plt.subplot(222)
# for i, n in enumerate(nc):
#     plt.plot(df['energy'], s_preproc[:, n] - i * offset, 'k-')
#     # plt.plot(df['energy'], s_bkg[:, n] - i * offset)
#     plt.plot(df['energy'], s_fit[:, n] - i * offset, 'r-')
#     plt.plot(df['energy'], s_bkg[:, n] - i * offset, 'm--')
#
# # plt.plot(df['energy'], ff_real/70)
# # plt.plot(df['energy'], iff/10 + 1)
#
# # plt.plot(s[:, 600])
# # plt.plot(s_bkg[:, 600])
#
# plt.subplot(223)
# plt.plot(df['energy'], s_proc[:, 200])
# plt.plot(df['energy'], s_proc[:, 70])
#
# plt.subplot(224)
# plt.plot(qq, anomalous_scat[0, :])
# plt.plot(qq, anomalous_scat[1, :])
#
# # plt.plot(s_proc[130, :])
# # plt.plot(df['iff']/df['i0'])
# # plt.plot(df['energy'], df['iff']/df['i0'])
#
#






