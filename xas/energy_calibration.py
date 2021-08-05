import numpy as np
from .file_io import load_binned_df_from_file
# from isstools.xasproject.xasproject import XASDataSet
from xas.xasproject import XASDataSet
from lmfit import Parameters, minimize


def get_foil_spectrum(element, db_proc):
    r = db_proc.search({'Sample_name': element + ' foil'})
    uid_proc = list(r)[-1]
    ds = db_proc[uid_proc].primary.read()
    energy = ds['Energy'].values
    mu = ds['mu_norm'].values
    return energy, mu


def compute_shift_between_spectra(energy, mu, energy_ref, mu_ref, e0, dE=25):
    mask = (energy_ref>=(e0-dE)) & (energy_ref<=(e0+dE))
    energy_ref_roi = energy_ref[mask]
    mu_ref_roi = mu_ref[mask]

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
        df, _ = load_binned_df_from_file(fname_bin)
        energy = df['energy'].values
        _mu = -np.log(df['ir'] / df['it']).values
        ds = XASDataSet(mu=_mu, energy=energy)
        mu = ds.flat

        element = start['element']
        e0 = float(start['e0'])
        energy_ref, mu_ref = get_foil_spectrum(element, db_proc)

        return compute_shift_between_spectra(energy, mu, energy_ref, mu_ref, e0, dE=dE)

        # return energy, mu_ref