import h5py

from xas.fitting import (get_normalized_gaussian_scan, estimate_center_and_width_of_peak, fit_gaussian,
                         Nominal2ActualConverter)
from xas.file_io import  load_binned_df_from_file

import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import savgol_filter
import os
import time as ttime
from numpy.polynomial import Polynomial

from scipy.spatial.transform import Rotation
from scipy.optimize import fsolve


def compute_rowland_circle_geometry(x_src, y_src, R, bragg_deg, det_dR):
    bragg = np.deg2rad(bragg_deg)

    x_cr = -R * np.cos(np.pi / 2 - bragg)
    y_cr = 0

    x_det = +2 * x_cr * np.cos(bragg) * np.cos(bragg) - det_dR * np.cos(np.pi - 2 * bragg)
    y_det = -2 * x_cr * np.cos(bragg) * np.sin(bragg) - det_dR * np.sin(np.pi - 2 * bragg)

    return (x_cr + x_src), (y_cr + y_src), (x_det + x_src), (y_det + y_src)


def _rotate_xyz(omega, bragg, xyz):
    rmat = Rotation.from_euler('zyz', [-(90 - bragg), omega, (90 - bragg)], degrees=True).as_matrix()
    return rmat @ xyz

def _solve_omega_func(omega, bragg, xyz, dz):
    _, _, z_cr_rot = _rotate_xyz(omega, bragg, xyz)
    return z_cr_rot - dz

def _compute_omega(R, bragg, dz):
    x_cr, _, _, _ = compute_rowland_circle_geometry(0, 0, R, bragg, 0)
    omega = np.arcsin(dz / (np.abs(x_cr) * np.sin(np.deg2rad(bragg))))
    return np.rad2deg(omega)

def _compute_omega_from_x_cr(x_cr, bragg, dz):
    omega = np.arcsin(dz / (np.abs(x_cr) * np.sin(np.deg2rad(bragg))))
    return np.rad2deg(omega)

def _compute_rotated_rowland_circle_geometry(x_cr_main, y_cr_main, x_det, y_det, bragg_deg, dz, output_omega=False):

    # dz = -dz # quirk of this coordinate system
    xyz_main_cr = np.array([x_cr_main, y_cr_main, 0])
    # omega0 = dz / 10 # heuristic that seems to work OK
    # omega = fsolve(_solve_omega_func, omega0, args=(bragg_deg, xyz_main_cr, dz))
    omega = _compute_omega_from_x_cr(x_cr_main, bragg_deg, dz)

    x_cr_rot, y_cr_rot, z_cr_rot = _rotate_xyz(omega, bragg_deg, xyz_main_cr)
    roll_cr_rot = bragg_deg + np.rad2deg(np.arcsin(np.abs(y_cr_rot / np.abs(x_cr_main))))
    # roll_cr_rot = 90 - np.rad2deg(np.arctan(np.abs((y_det - y_cr_rot) / (x_det - x_cr_rot))))/2
    yaw_cr_rot = np.rad2deg(np.arctan(np.abs((z_cr_rot) / (x_det - x_cr_rot))))
    if output_omega:
        return x_cr_rot, y_cr_rot, roll_cr_rot, yaw_cr_rot, omega
    return x_cr_rot, y_cr_rot, roll_cr_rot, yaw_cr_rot

def compute_rotated_rowland_circle_geometry(x_src, y_src, R, bragg, dz, output_omega=False):
    x_cr_main, y_cr_main, x_det, y_det = compute_rowland_circle_geometry(x_src, y_src, R, bragg, 0)
    output = _compute_rotated_rowland_circle_geometry(x_cr_main, y_cr_main, x_det, y_det, bragg, dz, output_omega=output_omega)
    return output

# compute_rotated_rowland_circle_geometry(0, 0, 1000, 85, 139.5)
# bragg = 90
# R = 1000
# det_dR = 0
# # x_cr_main, y_cr_main, x_det, y_det = compute_rowland_circle_geometry(0, 0, R, bragg, det_dR)
# # _compute_rotated_rowland_circle_geometry(x_cr_main, y_cr_main, x_det, y_det, bragg, 139.5)
#
# compute_rotated_rowland_circle_geometry(0, 0, R, bragg, 139.5)



from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt

def plot_crystal_geometry(ax, xyz_src, xyz_cr, xyz_det, xyz_rcc, xyz_rc, xyz_cr_ana):
    ax.plot(*xyz_src, 'k.')
    # ax.plot(*xyz_cr, 'b.')
    ax.plot(*xyz_det, 'r.')
    # ax.plot(*xyz_rcc, 'k.', alpha=0.5)
    ax.plot(*xyz_rc, 'k-', alpha=0.5)
    plt.plot(*xyz_cr_ana, 'b-')

# def compute_omega(R, bragg_deg, dz):
#     bragg = np.deg2rad(bragg_deg)
#
#     SC1 = R * np.cos(np.pi / 2 - bragg)
#     omega_prime = np.arccos(dz / SC1)
#     C0C1 =




def visualize_rowland_circle_geometry(R, bragg, det_dR=0, fignum=1, dz1=139.5, dz2=129.5, r=50):
    x_src, y_src, z_src = 0, 0, 0
    x_cr_main, y_cr_main, x_det, y_det = compute_rowland_circle_geometry(x_src, y_src, R, bragg, det_dR)
    z_cr_main = 0
    z_det = 0

    phi_cr = np.arange(0, 360)
    x_cr_ana_main = np.zeros(phi_cr.size)
    y_cr_ana_main = r * np.sin(np.deg2rad(phi_cr))
    z_cr_ana_main = r * np.cos(np.deg2rad(phi_cr))
    rmat_cr_ana = Rotation.from_euler('z', [90-bragg], degrees=True).as_matrix().squeeze()
    # print(rmat_cr_ana)
    xyz_cr_ana_main = rmat_cr_ana @ np.array([x_cr_ana_main, y_cr_ana_main, z_cr_ana_main])
    xyz_cr_ana_main[0] += x_cr_main

    x_rcc_main = -R / 2 * np.cos(np.deg2rad(90 - bragg))
    y_rcc_main = R / 2 * np.sin(np.deg2rad(90 - bragg))
    z_rcc_main = 0

    phi_rc = np.arange(0, 360)
    x_rc_main = R / 2 * np.cos(np.deg2rad(phi_rc)) + x_rcc_main
    y_rc_main = R / 2 * np.sin(np.deg2rad(phi_rc)) + y_rcc_main
    z_rc_main = np.zeros(phi_rc.size)

    def compute_aux_geometry(_dz):
        x_cr_aux, y_cr_aux, roll_cr_aux, yaw_cr_aux, omega_aux = compute_rotated_rowland_circle_geometry(x_src, y_src, R, bragg, _dz, output_omega=True)
        z_cr_aux = _dz
        xyz_rcc_aux = _rotate_xyz(omega_aux, bragg, np.array([x_rcc_main, y_rcc_main, z_rcc_main]))
        xyz_rc_aux = _rotate_xyz(omega_aux, bragg, np.array([x_rc_main, y_rc_main, z_rc_main]))
        xyz_cr_ana_aux = _rotate_xyz(omega_aux, bragg, xyz_cr_ana_main)
        return x_cr_aux, y_cr_aux, z_cr_aux, roll_cr_aux, yaw_cr_aux, omega_aux, xyz_rcc_aux, xyz_rc_aux, xyz_cr_ana_aux

    x_cr_aux2, y_cr_aux2, z_cr_aux2, roll_cr_aux2, yaw_cr_aux2, omega_aux2, xyz_rcc_aux2, xyz_rc_aux2, xyz_cr_ana_aux2 = compute_aux_geometry(-dz1)
    x_cr_aux3, y_cr_aux3, z_cr_aux3, roll_cr_aux3, yaw_cr_aux3, omega_aux3, xyz_rcc_aux3, xyz_rc_aux3, xyz_cr_ana_aux3 = compute_aux_geometry(dz1)
    x_cr_aux4, y_cr_aux4, z_cr_aux4, roll_cr_aux4, yaw_cr_aux4, omega_aux4, xyz_rcc_aux4, xyz_rc_aux4, xyz_cr_ana_aux4 = compute_aux_geometry(-dz1 - dz2)
    x_cr_aux5, y_cr_aux5, z_cr_aux5, roll_cr_aux5, yaw_cr_aux5, omega_aux5, xyz_rcc_aux5, xyz_rc_aux5, xyz_cr_ana_aux5 = compute_aux_geometry(dz1 + dz2)

    # print(f'{omega_aux2=}, omega={_compute_omega(R, bragg, -dz)}, d={omega_aux2 - _compute_omega(R, bragg, -dz)}')
    # print('')

    print(f'MAIN CRYSTAL COORDS: x={x_cr_main : 0.3f}, y={y_cr_main : 0.3f}, z={z_cr_main : 0.3f}')
    print(f'AUX2 CRYSTAL COORDS: x={x_cr_aux2 : 0.3f}, y={y_cr_aux2 : 0.3f}, z={z_cr_aux2 : 0.3f}')
    print(f'AUX3 CRYSTAL COORDS: x={x_cr_aux3 : 0.3f}, y={y_cr_aux3 : 0.3f}, z={z_cr_aux3 : 0.3f}')
    print(f'AUX4 CRYSTAL COORDS: x={x_cr_aux4 : 0.3f}, y={y_cr_aux4 : 0.3f}, z={z_cr_aux4 : 0.3f}')
    print(f'AUX5 CRYSTAL COORDS: x={x_cr_aux5 : 0.3f}, y={y_cr_aux5 : 0.3f}, z={z_cr_aux5 : 0.3f}')

    fig = plt.figure(fignum, clear=True)
    ax = plt.axes(projection='3d')
    plot_crystal_geometry(ax,
                          (x_src, z_src, y_src),
                          (x_cr_main, z_cr_main, y_cr_main),
                          (x_det, z_det, y_det),
                          (x_rcc_main, z_rcc_main, y_rcc_main),
                          (x_rc_main, z_rc_main, y_rc_main),
                          (xyz_cr_ana_main[0], xyz_cr_ana_main[2], xyz_cr_ana_main[1]))


    plot_crystal_geometry(ax,
                          (x_src, z_src, y_src),
                          (x_cr_aux2, z_cr_aux2, y_cr_aux2),
                          (x_det, z_det, y_det),
                          xyz_rcc_aux2[[0, 2, 1]],
                          (xyz_rc_aux2[0], xyz_rc_aux2[2], xyz_rc_aux2[1]),
                          (xyz_cr_ana_aux2[0], xyz_cr_ana_aux2[2], xyz_cr_ana_aux2[1]))

    plot_crystal_geometry(ax,
                          (x_src, z_src, y_src),
                          (x_cr_aux3, z_cr_aux3, y_cr_aux3),
                          (x_det, z_det, y_det),
                          xyz_rcc_aux3[[0, 2, 1]],
                          (xyz_rc_aux3[0], xyz_rc_aux3[2], xyz_rc_aux3[1]),
                          (xyz_cr_ana_aux3[0], xyz_cr_ana_aux3[2], xyz_cr_ana_aux3[1]))

    plot_crystal_geometry(ax,
                          (x_src, z_src, y_src),
                          (x_cr_aux4, z_cr_aux4, y_cr_aux4),
                          (x_det, z_det, y_det),
                          xyz_rcc_aux4[[0, 2, 1]],
                          (xyz_rc_aux4[0], xyz_rc_aux4[2], xyz_rc_aux4[1]),
                          (xyz_cr_ana_aux4[0], xyz_cr_ana_aux4[2], xyz_cr_ana_aux4[1]))

    plot_crystal_geometry(ax,
                          (x_src, z_src, y_src),
                          (x_cr_aux5, z_cr_aux5, y_cr_aux5),
                          (x_det, z_det, y_det),
                          xyz_rcc_aux5[[0, 2, 1]],
                          (xyz_rc_aux5[0], xyz_rc_aux5[2], xyz_rc_aux5[1]),
                          (xyz_cr_ana_aux5[0], xyz_cr_ana_aux5[2], xyz_cr_ana_aux5[1]))

    ax.set_xlim(-1200, 200)
    ax.set_ylim(-700, 700)
    ax.set_zlim(-700, 700)
    plt.axis('equal')


def normalize_peak(y_orig, bkg1=2, bkg2=-2, nmax=1, return_scale=False):
    y = y_orig.copy()
    offset = np.mean(np.hstack((y[:bkg1], y[bkg2:])))
    y -= offset
    scale = np.mean(np.sort(y)[-nmax:])
    y /= scale
    if return_scale:
        return y, scale
    return y


_pilatus_roi_colors = {1: 'tab:blue',
                       2: 'tab:orange',
                       3: 'tab:green',
                       4: 'tab:purple'}

def analyze_elastic_fly_scan(db, uid, rois=None, plot_func=None):
    fname_bin = db[uid].start['interp_filename'][:-3] + 'dat'
    df, _ = load_binned_df_from_file(fname_bin)
    energy = df['energy'].values

    if rois is None: rois = [1]

    for i in rois:
        field = f'pil100k2_roi{i}'
        intensity = df[field].values
        intensity = normalize_peak(intensity)
        Ecen0, fwhm0 = estimate_center_and_width_of_peak(energy, intensity)
        Ecen, fwhm, intensity_cor, intensity_fit, intensity_fit_raw = fit_gaussian(energy, intensity, Ecen0, fwhm0)
        print(f'{field}: {Ecen=:0.3f}, {fwhm=:0.3f}')
        if plot_func is not None:
            roi_color = _pilatus_roi_colors[i]
            roi_label = f'roi{i}'
            plot_func(energy, intensity_cor, intensity_fit, Ecen, fwhm, roi_label=roi_label, roi_color=roi_color, )
        fwhm_return = fwhm
    return fwhm




# analyze_elastic_fly_scan(db, -1)

def analyze_elastic_scan(db, uid):
    E, I, scale, offset = get_normalized_gaussian_scan(db, uid, return_norm_param=True)
    Ecen0, fwhm0 = estimate_center_and_width_of_peak(E, I)
    Ecen, fwhm, I_cor, I_fit, I_fit_raw = fit_gaussian(E, I, Ecen0, fwhm0)
    return Ecen, fwhm, I_cor, I_fit, (I_fit_raw*scale + offset), E


# def analyze_linewidth_fly_scan(db, uid, x_key='johann_main_crystal_motor_cr_main_roll', rois=None, plot_func=None):
#     fname_bin = db[uid].start['interp_filename'][:-3] + 'dat'
#     df, _ = load_binned_df_from_file(fname_bin)
#     x = df[x_key].values
#     if rois is None: rois = [1]
#
#     for i in rois:
#         field = f'pil100k_roi{i}'
#         y = np.abs(df[field].values / df['i0'].values)
#         y_smooth = savgol_filter(y, 5, 3)
#         y_smooth = normalize_peak(y_smooth, bkg1=5, bkg2=-5, nmax=5)
#         cen, fwhm = estimate_center_and_width_of_peak(x, y_smooth)
#         print(f'{field}: {cen=:0.3f}, {fwhm=:0.3f}')
#         if plot_func is not None:
#             roi_color = _pilatus_roi_colors[i]
#             roi_label = f'roi{i}'
#             plot_func(energy, intensity_smooth, intensity_smooth, cen, fwhm, roi_label=roi_label, roi_color=roi_color, )
#         fwhm_return = fwhm
#     return fwhm

def analyze_linewidth_fly_scan(db, uid, x_key='johann_main_crystal_motor_cr_main_roll', rois=None,
                               plot_func=None, processing_kwargs=None, **plot_kwargs):
    if (processing_kwargs is not None) and ('interp_filename' in processing_kwargs):
        path_to_file = processing_kwargs['interp_filename']
    else:
        path_to_file = db[uid].start['interp_filename']
    fname_bin = path_to_file[:-3] + 'dat'
    df, _ = load_binned_df_from_file(fname_bin)
    x = df[x_key].values
    if rois is None: rois = [1]

    field_base = 'pil100k2' if 'pil100k2_roi1' in df.columns else 'pil100k'
    output = {}

    for i in rois:
        if type(i) == int:
            field = f'{field_base}_roi{i}'
        else:
            field = f'{field_base}_{i}'

        # print(field)
        y = np.abs(df[field].values / (df['i0'] / df['i0'].median()).values)

        y, y_scale = normalize_peak(y, bkg1=5, bkg2=-5, nmax=5, return_scale=True)

        y_max = np.mean(np.sort(y)[-5:])
        mask = y >= y_max * 0.5

        x_roi = x[mask]
        y_roi = y[mask]

        y_roi_log = np.log(y_roi)
        n_order = 7
        p = Polynomial.fit(x_roi, y_roi_log, n_order)

        y_roi_log_fit = p(x_roi)
        y_roi_fit = np.exp(y_roi_log_fit)

        y_log_fit = p(x)
        y_smooth = np.exp(y_log_fit)
        y_smooth_roi = np.exp(y_roi_log_fit)

        x_extrema = p.deriv().roots()
        x_maxima = x_extrema[p.deriv().deriv()(x_extrema) < 0]
        x_maxima = x_maxima[(x_maxima >= x_roi.min()) & (x_maxima <= x_roi.max())]

        x_maximum = np.real(x_maxima[np.argmin(np.abs(x_maxima - x_roi[np.argmax(y_roi)]))])
        y_maximum = np.exp(p(x_maximum))

        x1, x2 = estimate_hm_positions_of_peak(x_roi, y_roi_fit / y_maximum)
        cen = x_maximum
        fwhm = np.abs(x1 - x2)
        # y_smooth = savgol_filter(y, 5, 3)
        # y_smooth = normalize_peak(y_smooth, bkg1=5, bkg2=-5, nmax=5)
        # cen, fwhm = estimate_center_and_width_of_peak(x, y_smooth)
        print(f'{field}: {cen=:0.3f}, {fwhm=:0.3f}')
        x_com = np.sum(y_smooth_roi * x_roi) / np.sum(y_smooth_roi)

        # plt.plot(x, y/y_maximum, '.')
        # plt.plot(x_roi, y_roi_fit/y_maximum, '-')

        if plot_func is not None:
            # roi_color = _pilatus_roi_colors[i]
            # roi_label = f'roi{i}'
            plot_func(x, y / y_maximum,
                      x_fit=x_roi, y_fit=y_smooth_roi / y_maximum,
                      x_peak=cen, y_peak=1, fwhm=fwhm, x_label=x_key, y_label=field, **plot_kwargs)

        if len(rois) == 1:
            return fwhm, y_maximum * y_scale, x_maximum, x_com

        output[i] = {'fwhm': fwhm, 'max_value': y_maximum * y_scale, 'max_loc': x_maximum, 'com_loc': x_com}

    return output


def estimate_hm_positions_of_peak(x, y):
    x_cen = x[np.argmax(np.abs(y))]
    y_diff = np.abs(y - 0.5)
    x_low = x < x_cen
    x_high = x > x_cen
    x1 = np.interp(0.5, y[x_low], x[x_low])
    x2 = np.interp(0.5, y[x_high][::-1], x[x_high][::-1])
    return x1, x2



def fit_polynom_and_estimate_minimum(x, y, deg=3):
    x = np.array(x)
    y = np.array(y)

    p = Polynomial.fit(x, y, deg)
    x_extrema = p.deriv().roots()
    x_minima = x_extrema[p.deriv().deriv()(x_extrema) > 0]
    x_minima = x_minima[(x_minima>=x.min()) & (x_minima<=x.max())]
    if len(x_minima) == 0:
        x_minimum = x[np.argmin(y)]
    else:
        y_minima = p(x_minima)
        x_minimum = x_minima[np.argmin(y_minima)]

    y_minimum = p(x_minimum)
    x_grid = np.linspace(x.min(), x.max(), x.size * 10)
    y_grid = p(x_grid)
    # plt.figure(2, clear=True)
    # plt.plot(x, y, 'k.-')
    # plt.plot(x, p(x), 'r-')
    # plt.vlines([x_minimum, x_minimum], y.min(), y.max(), colors='m')

    return x_minimum, y_minimum, x_grid, y_grid


def get_optimal_crystal_alignment_position(x, y, optimum='minimum', plot_func=None, **plot_kwargs):
    x = np.array(x)
    y = np.array(y)
    if optimum == 'maximum':
        x_min, y_min, x_fit, y_fit = fit_polynom_and_estimate_minimum(x, -y)
        y_fit *= -1
        y_min *= -1
    else:
        x_min, y_min, x_fit, y_fit = fit_polynom_and_estimate_minimum(x, y)
    if plot_func is not None:
        plot_func(x, y, x_fit=x_fit, y_fit=y_fit, x_peak=x_min, y_peak=y_min, curve_index=0, **plot_kwargs)
    return x_min

_crystal_info_dict = { 'main': {'motor_name': 'johann_main_crystal_motor_cr_main_roll', 'color': 'tab:blue'},
                       'aux2': {'motor_name': 'johann_aux2_crystal_motor_cr_aux2_roll', 'color': 'tab:orange'},
                       'aux3': {'motor_name': 'johann_aux3_crystal_motor_cr_aux3_roll', 'color': 'tab:green'},
                       'aux4': {'motor_name': 'johann_aux4_crystal_motor_cr_aux4_roll', 'color': 'tab:red'},
                       'aux5': {'motor_name': 'johann_aux5_crystal_motor_cr_aux5_roll', 'color': 'tab:purple'}}

def convert_roll_to_energy_for_johann_fly_scan(df, hdr, LUT=None):
    energy_base = None
    if LUT is None:
        LUT = hdr.start['spectrometer_config']['fly_calibration_dict']['LUT']

    if LUT is not None:
        trajectory_dict = hdr.start['spectrometer_relative_trajectory']
        crystals = list(trajectory_dict.keys())

        for crystal in crystals:
            if energy_base is None:
                energy_base = f'energy_{crystal}'
            _roll = LUT[crystal]
            _energy = LUT['energy']
            e2roll_converter = Nominal2ActualConverter(_energy, _roll)
            _cr_roll_motor_key = _crystal_info_dict[crystal]['motor_name']
            _cr_roll = df[_cr_roll_motor_key]
            _cr_energy = e2roll_converter.act2nom(_cr_roll)
            df[f'energy_{crystal}'] = _cr_energy

    return df, energy_base

def interpolate_df_emission_energies_on_common_grid(df, hdr, energy_key, det_key='pil100k2',):
    trajectory_dict = hdr.start['spectrometer_relative_trajectory']
    crystals = list(trajectory_dict.keys())
    for crystal in crystals:
        if crystal != energy_key:
            this_xes_key = f'{det_key}_{crystal}'
            df[this_xes_key] = np.interp(df[energy_key], df[f'energy_{crystal}'], df[this_xes_key])

    df = df.rename(columns={energy_key: 'energy'})
    columns_to_drop = [c for c in df.columns if ('energy_' in c) or ('_crystal_motor_cr' in c)]
    df = df.drop(columns=columns_to_drop)
    rois = [f'pil100k2_{k}' for k in hdr.start.spectrometer_relative_trajectory.keys()]
    df['pil100k2_xes_total'] = np.sum(df[rois].values, axis=1)
    df['pil100k2_xes_mean'] = np.mean(df[rois].values, axis=1)
    return df


def filter_johann_image_kwargs(kwargs: dict=None):
    if kwargs is None:
        kwargs = {}
    return {k: v for k, v in kwargs.items() if k in ['poly_roi_dict']}

def filter_johann_calibration_kwargs(kwargs: dict=None):
    if kwargs is None:
        kwargs = {}
    return {k: v for k, v in kwargs.items() if k in ['LUT']}


pilatus_mask = np.ones((195, 487), dtype=bool)
pilatus_mask[15, 352] = False
pilatus_mask[158, 11] = False

def pilatus_image_com(image, roi):
    x, y, dx, dy = roi
    image_roi = image[x : x + dx, y : y + dy]
    image_roi_x = np.sum(image_roi, axis=1)
    image_roi_y = np.sum(image_roi, axis=0)
    roi_x = np.arange(dx)
    roi_y = np.arange(dy)
    com_x = np.sum(image_roi_x * roi_x) / np.sum(image_roi_x)
    com_y = np.sum(image_roi_y * roi_y) / np.sum(image_roi_y)
    return x + com_x, y + com_y

def pilatus_position_roi(image, dx, dy):
    image[pilatus_mask] = np.percentile(image, 10)
    com_x0, com_y0 = pilatus_image_com(image, (0, 0, 195, 487))
    print(com_x0, com_y0)
    com_x, com_y = pilatus_image_com(image, (int(com_x0 - dx/2*1.5),
                                             int(com_y0 - dy/2*1.5),
                                             int(dx*1.5), int(dy*1.5)))
    print(int(com_x0 - dx/2*1.5),
          int(com_y0 - dy/2*1.5),
          int(dx*1.5), int(dy*1.5))
    return com_x, com_y

# c_x, c_y = pilatus_position_roi(img, 30, 80)
# plt.figure(1)
# plt.clf()
# plt.imshow(img, vmin=0, vmax=100)
# plt.hlines(c_x, 0, 487)
# plt.vlines(c_y, 0, 195)



def analyze_many_elastic_scans(db, uids, E_nominal, plotting=False, short_output=True):
    E_actual = []
    resolution = []
    I_cors = []
    I_fits = []
    I_fit_raws = []
    E_scans = []

    for uid in uids:
        Ecen, fwhm, I_cor, I_fit, I_fit_raw, E_scan = analyze_elastic_scan(db, uid)
        E_actual.append(Ecen)
        resolution.append(fwhm)
        I_cors.append(I_cor)
        I_fits.append(I_fit)
        I_fit_raws.append(I_fit_raw)
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
    if short_output:
        return energy_converter
    else:
        return energy_converter, E_actual, resolution, I_fit_raws



def get_xes_data(db, uid, offset=-10):
    xes = 1
    if uid:
        t = db[uid].table()
        if 'apb_ave_ch1_mean' in t.columns:
            denom = ((np.abs(t['apb_ave_ch1_mean']) - offset)).values
        else:
            denom = 1
        xes = (t['pil100k_stats1_total']/denom).values.squeeze()
    return xes


def parse_file_with_uids(file_with_uids, xes_normalization=True):
    energies_out, uids_herfd, uids_xes = [], [], []
    with open(file_with_uids, "r") as f:
        lines = f.readlines()
    for line in lines:
        words = line.split(' ')
        if xes_normalization:
            energies_out.append(float(words[-3]))
            uids_herfd.append(words[-2])
            uids_xes.append(words[-1][:-2])
        else:
            energies_out.append(float(words[-2]))
            uids_herfd.append(words[-1][:-2])
            uids_xes.append(None)
    return np.array(energies_out), uids_herfd, uids_xes

def parse_rixslog(uids_bundle):
    energies_out, uids_norm = [], []
    f = h5py.File(uids_bundle, 'r')
    uids_herfd = list(f.keys())
    n_spectra = len(uids_herfd)
    energies_out = np.zeros(n_spectra)
    uids_norm = np.zeros(n_spectra, dtype='<U40')
    for i, uid_herfd in enumerate(uids_herfd):
        ds = f[uid_herfd]
        energies_out[i] = ds['emission_energy'][()]
        if 'uid_norm' in ds.keys():
            uids_norm[i] = ds['uid_norm'][()]
    f.close()
    energies_out = np.array(energies_out)
    idx_ord = np.argsort(energies_out)
    energies_out = energies_out[idx_ord]
    uids_herfd = np.array(uids_herfd)[idx_ord]
    uids_norm = np.array(uids_norm)[idx_ord]

    return energies_out, uids_herfd, uids_norm


    # for line in lines:
    #     words = line.split(' ')
    #     if xes_normalization:
    #         energies_out.append(float(words[-3]))
    #         uids_herfd.append(words[-2])
    #         uids_xes.append(words[-1][:-2])
    #     else:
    #         energies_out.append(float(words[-2]))
    #         uids_herfd.append(words[-1][:-2])
    #         uids_xes.append(None)
    # return np.array(energies_out), uids_herfd, uids_xes


def get_herfd_data(db, uid):
    filename = db[uid].start['interp_filename']
    path, extension = os.path.splitext(filename)
    if extension != '.dat':
        filename = path+'.dat'
    print(f' @@@@@@@ Reading {filename}')
    df, _ = load_binned_df_from_file(filename)
    energy = df['energy'].values
    # herfd = np.abs((df[ROI])/ df['i0']).values
    # if ROI_bkg is not None:
    #     herfd_bkg = np.abs((df[ROI_bkg]) / df['i0']).values
    #     return energy, (herfd-herfd_bkg)
    return energy, df



def parse_rixs_scan(db, file_with_uids, xes_normalization=True):
    energies_out, uids_herfd, uids_xes = parse_file_with_uids(file_with_uids, xes_normalization=xes_normalization)
    xes_data = []
    herfd_data = []
    for uid_herfd, uid_xes in zip(uids_herfd, uids_xes):
        energies_in, df = get_herfd_data(db, uid_herfd)
        herfd_data.append(herfd)
        if xes_normalization:
            xes = get_xes_data(db, uid_xes)
            xes_data.append(xes)

    herfd_data = np.array(herfd_data).T

    if xes_normalization:
        xes_data = np.array(xes_data).T
        # linear algebra magic for scaling purposes
        # c, _, _, _ = np.linalg.lstsq(xes_data[None, :], xes_data[0, None], rcond=-1)
        c = xes_data/xes_data[0]
        herfd_data /= c
    return herfd_data, xes_data, energies_in, energies_out


def parse_rixslog_scan(db, file_with_uids, xes_normalization=True):
    energies_out, uids_herfd, uids_norm = parse_rixslog(file_with_uids)
    xes_data = []
    herfd_dfs = []
    keys = []
    energies_in = None
    n_spectra = energies_out.size

    for i in range(n_spectra):
        energies_in, df = get_herfd_data(db, uids_herfd[i])
        keys = df.keys()
        herfd_dfs.append(df)
        xes = get_xes_data(db, uids_norm[i])
        xes_data.append(xes)

    result_dict = {}
    for key in keys:
        matrix = []
        for _df in herfd_dfs:
            arr = _df[key].values
            matrix.append(arr)
        result_dict[key] = np.array(matrix).T


    xes_data = np.array(xes_data).T

    if xes_normalization:
        # linear algebra magic for scaling purposes
        # c, _, _, _ = np.linalg.lstsq(xes_data[:,:1], xes_data, rcond=-1)
        c = xes_data/xes_data[0]
        for key in keys:
            if str.lower(key).startswith('pil'):
                result_dict[key] /= c
    result_dict['conecntration_norm'] = xes_data
    result_dict['energy_in'] = result_dict['energy'][:, 0]
    result_dict.pop('energy')
    result_dict['energy_out'] = energies_out

    return result_dict


def convert_rixs_to_energy_transfer(Ein, Eout, herfd):
    dEmin = Ein.min() - Eout.max()
    dEmax = Ein.max() - Eout.min()

    transfer_step = np.min([np.min(np.abs(np.diff(np.unique(Ein)))), np.abs(Eout[1] - Eout[0])])*np.sqrt(2)
    dE = np.arange(dEmin, dEmax, transfer_step)
    # dE = np.linspace(dEmin, dEmax, 601)
    rixs = np.zeros((Ein.size, dE.size))

    for idx in range(Ein.size):
        each_Ein = Ein[idx]
        each_transfer = (each_Ein - Eout)
        idx_ord = np.argsort(each_transfer)
        rixs[idx, :] = np.interp(dE, each_transfer[idx_ord], herfd[idx_ord, idx], left=0, right=0)

    return dE, rixs

# plt.figure()
# plt.plot(energies_vtc_cubanes, xes_co3mno4_all, 'k-', lw=1, alpha=0.3)
# plt.plot(energies_vtc_cubanes, np.mean(xes_co3mno4_all, axis=1), 'r-', lw=2)

# herfd_data, xes_data, energy_in = parse_rixs_scan(db, '/nsls2/xf08id/users/2021/1/308190/rixs_uids_ebce49.txt')
# energy_transfer, rixs = convert_rixs_to_energy_transfer(energy_in, energies_emission[:herfd_data.shape[1]], herfd_data.T)
# plt.figure()
# plt.figure(); plt.contourf(energy_in, energy_transfer, (rixs/rixs.max()).T, 251, vmin=0.0, vmax=0.15, cmap='jet')
# plt.axis('image')
# plt.xlim(7706, 7715); plt.ylim(56, 65)


# plt.figure(figsize=(8/2.54, 6/2.54))
# plt.contourf(energy_in, energy_transfer-2.5, (rixs/rixs.max()).T, 551, vmin=0.00, vmax=0.1, cmap='jet')
# plt.xlim(7706, 7715)
# plt.ylim(56, 65)
# plt.xlabel('Incident energy, eV')
# plt.ylabel('Energy transfer, eV')

def process_rixs_von_hamos(db, uid, roi, subtract_bkg=False):
    t = db[uid].table(fill=True)
    try:
        energies_in = t['hhm_energy'].values
    except:
        energies_in = None
    images_pandas = t['pil100k_image']
    images = np.array([im[0, :, :] for im in images_pandas])
    spectra = []
    for image in images:
        spectrum = process_image_von_hamos(image, roi, subtract_bkg=subtract_bkg)
        spectra.append(spectrum)
    spectra = np.array(spectra)
    return energies_in, spectra


def process_image_von_hamos(image, roi, subtract_bkg=False):
    roi_x1, roi_x2, roi_y1, roi_y2, int_axis = roi
    image_roi = image[roi_x1: roi_x2, roi_y1: roi_y2]
    if subtract_bkg:
        if int_axis == 0:
            bkg_roi = image[roi_x1: roi_x2, roi_y1 + roi_y2: roi_y2 + roi_y2]
        if int_axis == 1:
            bkg_roi = image[roi_x1 + roi_x2: roi_x2 + roi_x2, roi_y1: roi_y2]
        bkg = np.sum(bkg_roi, axis=int_axis)
        bkg = savgol_filter(bkg, 21, 3)
    else:
        bkg = None

    spectrum = np.sum(image_roi, axis=int_axis)
    if bkg is not None:
        spectrum = spectrum - bkg

    return spectrum

'''
# VTC
roi_x1, roi_x2, roi_y1, roi_y2, int_axis =  10, 30, 0, 300, 0
roi = (roi_x1, roi_x2, roi_y1, roi_y2, int_axis)
uid_ti_foil = 'b9008469-e6e7-4345-8872-5b632a90e29a'
_, spectra = process_rixs_von_hamos(db, uid_ti_foil, roi, subtract_bkg=True)
spectrum_av_Ti = np.mean(spectra, axis=0)
spectrum_av_Ti /= spectrum_av_Ti.max()

uid_tin = 'ead8b8e4-3171-4bbc-ab34-8dfff3bcac8c'
_, spectra = process_rixs_von_hamos(db, uid_tin, roi, subtract_bkg=True)
spectrum_av_TiN = np.mean(spectra, axis=0)
spectrum_av_TiN /= spectrum_av_TiN.max()

uid_fetio3 = '2dd621a7-8690-4efc-82ab-0b18fbeb1fa8'
_, spectra = process_rixs_von_hamos(db, uid_fetio3, roi, subtract_bkg=True)
spectrum_av_FeTiO3 = np.mean(spectra, axis=0)
spectrum_av_FeTiO3 /= spectrum_av_FeTiO3.max()

uid_liti2o3 = 'b0aa6acb-a341-4bef-8b07-631009a610c3'
_, spectra = process_rixs_von_hamos(db, uid_liti2o3, roi, subtract_bkg=True)
spectrum_av_LiTi2O3 = np.mean(spectra, axis=0)
spectrum_av_LiTi2O3 /= spectrum_av_LiTi2O3.max()


uid_TiO2 = '53c3ae46-6195-403e-9eca-93283cbf9c2b'
_, spectra = process_rixs_von_hamos(db, uid_TiO2, roi, subtract_bkg=True)
spectrum_av_TiO2 = np.mean(spectra, axis=0)
spectrum_av_TiO2 /= spectrum_av_TiO2.max()


uid_Ti2O3 = 'eefcc422-e82f-4279-a6fd-054a9902ca49'
_, spectra = process_rixs_von_hamos(db, uid_Ti2O3, roi, subtract_bkg=True)
spectrum_av_Ti2O3 = np.mean(spectra, axis=0)
spectrum_av_Ti2O3 /= spectrum_av_Ti2O3.max()

fig, ax = plt.subplots(1,1)
# ax.plot(spectrum_av_Ti[::-1], label='Ti')
# ax.plot(spectrum_av_TiN[::-1], label='TiN')
# ax.plot(spectrum_av_FeTiO3[::-1], label='FeTiO3')
# ax.plot(spectrum_av_LiTi2O3[::-1], label='LiTi2O3')
ax.plot(spectrum_av_TiO2[::-1], label='TiO2')
ax.plot(spectrum_av_Ti2O3[::-1], label='Ti2O3')
ax.legend()




# RIXS
roi_x1, roi_x2, roi_y1, roi_y2, int_axis = 0, 192, 200, 300, 1
roi = (roi_x1, roi_x2, roi_y1, roi_y2, int_axis)
# uid_ti_foil = 'fbba85fc-383b-42d0-9f51-2d7552da919c'
# uid_srtio3 = '6820fc1b-80ef-4c4b-8218-577c395f4162'

uid_fetio3 = 'd537688b-ca25-41fa-8f39-ae926db7860e'
energies_out = (np.arange(192) - 130) * 0.2 + 4513
# # energies_in, herfds = process_rixs_von_hamos(db, uid_ti_foil, roi)
# # energies_in, herfds = process_rixs_von_hamos(db, uid_srtio3, roi)
energies_in, herfds = process_rixs_von_hamos(db, uid_fetio3, roi)
#
# from scipy.signal import savgol_filter
#
herfds_smooth = savgol_filter(herfds.T, 1, 0).T

energies_transfer, rixs = convert_rixs_to_energy_transfer(energies_in, energies_out, herfds.T)
_, rixs_smooth = convert_rixs_to_energy_transfer(energies_in, energies_out, herfds_smooth.T)


plt.figure(1)
plt.clf()
# plt.contourf(energies_in, energies_transfer, np.log10((rixs/rixs.max())).T, 51, cmap='jet')
# plt.contour(energies_in, energies_transfer, np.log10((rixs/rixs.max())).T, 51, colors='k', linewidths=0.5)
plt.contourf(energies_in, energies_transfer, np.log10((rixs_smooth/rixs_smooth.max())).T, 31, cmap='jet')
plt.contour(energies_in, energies_transfer, np.log10((rixs_smooth/rixs_smooth.max())).T, 31, colors='k', linewidths=0.5)
plt.axis('image')
plt.ylim(454, 491)
plt.xlim(4964, 4986)
#
# plt.figure(2)
# plt.clf()
# plt.contourf(energies_in, energies_out, np.log10((herfds/herfds.max()).T), 999, cmap='jet')
# # plt.contour(energies_in, energies_out, herfds.T, 21, colors='k', linewidths=0.5, vmin=0, vmax=0.1)
# plt.axis('image')
# plt.ylim(4487, 4525.2)
# # plt.xlim(4964, 4976)
#
# herfd_kalpha_1 = np.sum(herfds_smooth[:, 65:82], axis=1)
# herfd_kalpha_2 = np.sum(herfds_smooth[:, 98:115], axis=1)
# herfd_kalpha_1 -= herfd_kalpha_1[0]
# herfd_kalpha_1 = herfd_kalpha_1 / herfd_kalpha_1.max()
# herfd_kalpha_2 -= herfd_kalpha_2[0]
# herfd_kalpha_2 = herfd_kalpha_2/ herfd_kalpha_2.max()
# plt.figure(3)
# plt.clf()
# plt.plot(energies_in, herfd_kalpha_1)
# plt.plot(energies_in, herfd_kalpha_2)



# plt.figure(1); plt.clf();
# plt.imshow(images[-1, roi_x1:roi_x2, roi_y1:roi_y2], vmin=0, vmax=1000)
'''

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

        self.x = self.R / 2 * (1 + np.cos(np.pi - 2*self.ba))
        self.y = self.R / 2 * np.sin(np.pi - 2 * self.ba)
        self.d_y = self.R  * np.sin(np.pi - 2 * self.ba)

    def compute_energy_from_positions(self, x, y, d_y):
        ba1 = - 0.5 * (np.arccos(x / (self.R / 2) - 1) - np.pi)
        ba2 = - 0.5 * (np.arcsin(y / (self.R / 2)) - np.pi)
        ba3 = - 0.5 * (np.arcsin(d_y / (self.R)) - np.pi)
        e1 = self.bragg_energy(ba1)
        e2 = self.bragg_energy(ba2)
        e3 = self.bragg_energy(ba3)
        return np.mean([e1, e2, e3])





files=('/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0002.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0003.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0004.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0005.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0006.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0007.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0008.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0009.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0010.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0011.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0012.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0013.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0014.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0015.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0016.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0017.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0018.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0019.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0020.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0021.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0022.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0023.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0024.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0025.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0026.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0027.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0028.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0029.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0030.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0031.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0032.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0033.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0034.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0035.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0036.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0037.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0038.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0039.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0040.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0041.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0042.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0043.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0044.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0045.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0046.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0047.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0048.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0049.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0050.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0051.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0052.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0053.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0054.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0055.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0056.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0057.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0058.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0059.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0060.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0061.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0062.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0063.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0064.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0065.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0066.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0067.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0068.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0069.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0070.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0071.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0072.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0073.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0074.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0075.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0076.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0077.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0078.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0079.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0080.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0081.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0082.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0083.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0084.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0085.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0086.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0087.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0088.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0089.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0090.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0091.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0092.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0093.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0094.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0095.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0096.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0097.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0098.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0099.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0100.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0101.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0102.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0103.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0104.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0105.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0106.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0107.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0108.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0109.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0110.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0111.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0112.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0113.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0114.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0115.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0116.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0117.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0118.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0119.dat',
    '/nsls2/xf08id/users/2021/1/300001/Co3O4 RIXS 0001-r0120.dat')


def read_files_for_rixs(files):

    herfds = []
    for file in files:
        data = np.genfromtxt(file)
        energy_in = data[:, 0]
        herfd = -data[:, 9] / data[:, 1]
        herfds.append(herfd)

    return energy_in, herfds
