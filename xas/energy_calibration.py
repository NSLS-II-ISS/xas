import numpy as np
from .file_io import load_binned_df_from_file
# from isstools.xasproject.xasproject import XASDataSet
from xas.xasproject import XASDataSet
# from xas.trajectory import read_trajectory_limits
from lmfit import Parameters, minimize
import time as ttime
import xraydb
from xas import xray
from .db_io import get_fly_uids_for_proposal

def get_foil_spectrum(element, edge, db_proc):
    r = db_proc.search({'Sample_name' : f'{element} foil', 'Edge' : edge})
    uid_proc = list(r)[0]
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


def get_energy_offset(uid, db, db_proc, dE=25, plot_fun=None, attempts=5, sleep_time=1, full_return=False):
    start = db[uid].start
    fname_raw = start['interp_filename']
    if fname_raw.endswith('.raw'):
        fname_bin = fname_raw[:-4] + '.dat'

        for i in range(attempts):
            try:
                df, _ = load_binned_df_from_file(fname_bin)
            except:
                print(f'[Energy Calibration] Attempt to read data {i+1}')
                ttime.sleep(sleep_time)
                df = None

        try:
            energy = df['energy'].values
            _mu = -np.log(df['ir'] / df['it']).values
            ds = XASDataSet(mu=_mu, energy=energy)
            mu = ds.flat

            element = start['element']
            edge = start['edge']
            e0 = float(start['e0'])
            # energy_ref, mu_ref = get_foil_spectrum(element, edge, db_proc)
            energy_ref, mu_ref = db_proc.foil_spectrum(element, edge)
            mask = (energy_ref >= (e0 - dE)) & (energy_ref <= (e0 + dE))

            energy_ref_roi = energy_ref[mask]
            mu_ref_roi = mu_ref[mask]
            shift, mu_fit = compute_shift_between_spectra(energy, mu, energy_ref_roi, mu_ref_roi)
            e_cor = e0 + shift
            if plot_fun is not None:
                # mu = np.interp(energy_ref_roi, energy, mu)
                plot_fun(energy_ref_roi, mu_ref_roi, mu_fit)

        except Exception as e:
            print(f'[Energy Calibration] Error: {e}')
            e0, e_cor, energy_ref_roi, mu_ref_roi, mu_fit = None, None, None, None, None

        if full_return:
            return e0, e_cor, energy_ref_roi, mu_ref_roi, mu_fit
        else:
            return e0, e_cor




        # return e0, shift, energy_ref_roi, mu_ref_roi, mu
        # return energy, mu_ref


def get_energy_offset_for_proposal(db, year, cycle, proposal, db_proc, dE=25):
    uids = get_fly_uids_for_proposal(db, year, cycle, proposal)
    n_uids = len(uids)
    print(f'Found {n_uids} scans')
    e0_list = []
    time_e0_list = []
    chisq_list = []
    for i, uid in enumerate(uids):
        print(f'Analyzing scan {i} of {n_uids} (UID = {uid}): ', end='')
        hdr = db[uid]
        if ('exit_status' in hdr.stop.keys()) and (hdr.stop['exit_status']):
            _e_ref, _e_obs, _, mu_ref, mu_fit = get_energy_offset(uid, db, db_proc, dE=dE, attempts=1, sleep_time=0, full_return=True)
            if _e_obs:
                e0_list.append(_e_obs)
                _chisq = np.sum((mu_ref - mu_fit) ** 2) / mu_ref.size
                chisq_list.append(_chisq)
                _time = hdr.start['time']
                time_e0_list.append(_time)
                print(f'_e_obs={_e_obs}, chisq={_chisq}')

    e0 = np.array(e0_list)
    time_e0 = np.array(time_e0_list)
    chisq = np.array(chisq_list)
    return time_e0, e0, chisq


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



