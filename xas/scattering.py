

import numpy as np
import pyFAI as pyFAI





def get_ai(dist=40, center_ver=100, center_hor=480,
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

# image = np.array(data_list[0])
# mask = image < np.percentile(image.ravel(), 5)
#
# # res = ai.integrate1d_ng(np.array(data_list[0]),
# #                         1000,
# #                         mask=mask,
# #                         unit="2th_deg")
# # jupyter.plot1d(res)
#
# res2 = ai.integrate2d_ng(image,
#                          300, 360,
#                          mask=mask,
#                          unit="r_mm")
#
# jupyter.plot2d(res2)



def integrate_pil100k_image_stack(df_images, energies=None, dist=40, center_ver=100, center_hor=480):
    ai = get_ai(dist=40, center_ver=100, center_hor=480)
    s = []
    mask = None
    for item in df_images.values:
        image = np.array(item[0])
        if mask is None:
            mask = image < np.percentile(image.ravel(), 5)
        res = ai.integrate1d_ng(image,
                                1000,
                                mask=mask,
                                unit="2th_deg")
        s.append(res[1])
        tth = res[0]
    return tth, np.array(s)


df, _ = load_binned_df_from_file('/nsls2/data/iss/legacy/processed/2022/2/300010/Ir sample 2 AXS try4 cont 0037-r0002.dat')
data = pd.read_json('/nsls2/data/iss/legacy/processed/2022/2/300010/extended_data/Ir sample 2 AXS try4 cont 0037-r0002.json')
tth, s = integrate_pil100k_image_stack(data)
s /= -df['i0'].values[:, None]
s /= s.max()

def process_stack_of_patterns(energy, s, iff, emin=11200, emax=11250):
    e_lo_mask = energy <= emin
    # e_hi_mask = energy >= emax
    poly = np.polyfit(energy[e_lo_mask], s[e_lo_mask, :], 1)
    s_bkg = np.zeros(s.shape)
    # for i in range(s.shape[1]):

    # s_proc = s - s_bkg

    basis = np.vstack((np.ones(energy.size), energy, energy**2, energy**3, iff)).T
    print(basis.shape, s[:, 0].shape)
    for i in range(s.shape[1]):
        np.polyval(poly[:, i], energy)
        c, _, _, _ = np.linalg.lstsq(basis, s[:, i])
        s_bkg[:, i] = basis @ c

    s_proc = s - s_bkg

    # s_proc /= np.mean(s_proc[e_hi_mask, :], axis=0)
    return s_proc, s_bkg

s_proc, s_bkg = process_stack_of_patterns(df['energy'], s, (df['iff']/df['i0']).values)

plt.figure(1, clear=True);

plt.subplot(221)
plt.plot(s_proc.T)

plt.subplot(222)
plt.plot(df['energy'], s[:, 700])
plt.plot(df['energy'], s_bkg[:, 700])

plt.subplot(223)
plt.plot(s_proc[:, 200])

plt.subplot(224)
plt.plot(s_proc[300, :])

plt.plot(df['energy'], df['iff']/df['i0'])