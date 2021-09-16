import numpy as np
from .file_io import load_binned_df_from_file
# from isstools.xasproject.xasproject import XASDataSet
from xas.xasproject import XASDataSet
from xas.trajectory import read_trajectory_limits
from lmfit import Parameters, minimize
import time as ttime
import xraydb


def get_foil_spectrum(element, db_proc):
    r = db_proc.search({'Sample_name': element + ' foil'})
    if len(r) == 0:
        return None, None
    uid_proc = list(r)[-1]
    ds = db_proc[uid_proc].primary.read()
    energy = ds['Energy'].values
    mu = ds['mu_norm'].values
    return energy, mu


def compute_shift_between_spectra(energy, mu, energy_ref_roi, mu_ref_roi, e0, mask):

    def residuals(pars):
        e_shift = pars.valuesdict()['e_shift']
        x = np.interp(energy_ref_roi, energy - e_shift, mu)
        basis = np.vstack((x, np.ones(x.shape))).T
        c, _, _, _ = np.linalg.lstsq(basis, mu_ref_roi)
        return (basis @ c - mu_ref_roi)

    pars = Parameters()
    pars.add('e_shift', value=0)
    out = minimize(residuals, pars)
    e_shift = out.params['e_shift'].value

    return e_shift


def get_energy_offset(uid, db, db_proc, dE=25):
    start = db[uid].start
    fname_raw = start['interp_filename']
    if fname_raw.endswith('.raw'):
        fname_bin = fname_raw[:-4] + '.dat'

        for i in range(5):
            try:
                df, _ = load_binned_df_from_file(fname_bin)
            except:
                ttime.sleep(1)

        energy = df['energy'].values
        _mu = -np.log(df['ir'] / df['it']).values
        ds = XASDataSet(mu=_mu, energy=energy)
        mu = ds.flat

        element = start['element']
        e0 = float(start['e0'])
        energy_ref, mu_ref = get_foil_spectrum(element, db_proc)
        mask = (energy_ref >= (e0 - dE)) & (energy_ref <= (e0 + dE))

        energy_ref_roi = energy_ref[mask]
        mu_ref_roi = mu_ref[mask]
        shift = compute_shift_between_spectra(energy, mu, energy_ref_roi, mu_ref_roi, e0, mask)

        mu = np.interp(energy_ref_roi, energy, mu)
        return shift, energy_ref_roi,mu_ref_roi, mu
        # return energy, mu_ref

def validate_calibration(element, edge,db_proc, hhm, ):
    # check if current trajectory is good for this calibration
    r = db_proc.search({'Sample_name': element + ' foil'})
    if len(r) == 0:
        return False, f'Error: No matching foil has been found'

    e_min, e_max = read_trajectory_limits(hhm)
    edge_energy = xraydb.xray_edge(element, edge).energy
    if not ((edge_energy > e_min) and (edge_energy < e_max)):
        return False, f'Error: invalid trajectory for this calibration'

    return True, ''

def process_calibration(db, db_proc, hhm, dE=25, axis=None, canvas=None):
    e_shift, en_ref,mu_ref,mu = get_energy_offset(-1, db, db_proc, dE=dE)








    # _offset_act = xray.energy2encoder(e0_act, hhm.pulses_per_deg)
    # _offset_nom = xray.energy2encoder(e0_nom, hhm.pulses_per_deg)
    # delta_offset = (_offset_act - _offset_nom) / hhm.pulses_per_deg
    # new_offset = hhm.angle_offset.value - delta_offset
    # yield from bps.mv(hhm.angle_offset, new_offset)
    # return e_shift,en_ref,mu_ref,mu