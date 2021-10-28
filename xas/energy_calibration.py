import numpy as np
from .file_io import load_binned_df_from_file
# from isstools.xasproject.xasproject import XASDataSet
from xas.xasproject import XASDataSet
# from xas.trajectory import read_trajectory_limits
from lmfit import Parameters, minimize
import time as ttime
import xraydb
from xas import xray

def get_foil_spectrum(element, edge, db_proc):
    r = db_proc.search({'Sample_name' : f'{element} foil', 'Edge' : edge})
    if len(r) == 0:
        return None, None
    uid_proc = list(r)[-1]
    ds = db_proc[uid_proc].primary.read()
    energy = ds['Energy'].values
    mu = ds['mu_norm'].values
    return energy, mu


def compute_shift_between_spectra(energy, mu, energy_ref_roi, mu_ref_roi):

    def interpolated_spectrum(pars):
        e_shift = pars.valuesdict()['e_shift']
        x = np.interp(energy_ref_roi, energy - e_shift, mu)
        basis = np.vstack((x, np.ones(x.shape))).T
        c, _, _, _ = np.linalg.lstsq(basis, mu_ref_roi)
        return basis @ c

    def residuals(pars):
        return (interpolated_spectrum(pars) - mu_ref_roi)

    pars = Parameters()
    pars.add('e_shift', value=0)
    out = minimize(residuals, pars)
    e_shift = out.params['e_shift'].value
    mu_fit = interpolated_spectrum(out.params)
    return e_shift, mu_fit


def get_energy_offset(uid, db, db_proc, dE=25, plot_fun=None):
    start = db[uid].start
    fname_raw = start['interp_filename']
    if fname_raw.endswith('.raw'):
        fname_bin = fname_raw[:-4] + '.dat'

        for i in range(5):
            try:
                df, _ = load_binned_df_from_file(fname_bin)
            except:
                print(f'[Energy Calibration] Attempt to read data {i+1}')
                ttime.sleep(1)

        energy = df['energy'].values
        _mu = -np.log(df['ir'] / df['it']).values
        ds = XASDataSet(mu=_mu, energy=energy)
        mu = ds.flat

        element = start['element']
        edge = start['edge']
        e0 = float(start['e0'])
        energy_ref, mu_ref = get_foil_spectrum(element, edge, db_proc)
        mask = (energy_ref >= (e0 - dE)) & (energy_ref <= (e0 + dE))

        energy_ref_roi = energy_ref[mask]
        mu_ref_roi = mu_ref[mask]
        shift, mu_fit = compute_shift_between_spectra(energy, mu, energy_ref_roi, mu_ref_roi)

        if plot_fun is not None:
            # mu = np.interp(energy_ref_roi, energy, mu)
            plot_fun(energy_ref_roi, mu_ref_roi, mu_fit)

        return e0, e0+shift
        # return e0, shift, energy_ref_roi, mu_ref_roi, mu
        # return energy, mu_ref



# def process_calibration(element, edge, db, db_proc, hhm, trajectory_manager, dE=25, axis=None, canvas=None):
#     e_shift, en_ref, mu_ref, mu = get_energy_offset(-1, db, db_proc, dE=dE)
#     # energy_nominal = xraydb.xray_edge(element, edge).energy
#     # energy_actual = energy_nominal + e_shift
#     # offset_actual = xray.energy2encoder(energy_actual, hhm.pulses_per_deg) / hhm.pulses_per_deg
#     # offset_nominal = xray.energy2encoder(energy_nominal, hhm.pulses_per_deg) / hhm.pulses_per_deg
#     # angular_offset_shift = offset_actual - offset_nominal
#     # new_angular_offset = hhm.angle_offset.value - angular_offset_shift
#     # if hhm.set_new_angle_offset(new_angular_offset):
#     #     current_index =
#     #
#     #     return e_shift, en_ref, mu_ref, mu
#     #
#     #
#     #
#     #
#     # _offset_act = xray.energy2encoder(e0_act, hhm.pulses_per_deg)
#     # _offset_nom = xray.energy2encoder(e0_nom, hhm.pulses_per_deg)
#     # delta_offset = (_offset_act - _offset_nom) / hhm.pulses_per_deg
#     # new_offset = hhm.angle_offset.value - delta_offset
#     # yield from bps.mv(hhm.angle_offset, new_offset)
#     return e_shift, en_ref, mu_ref, mu



