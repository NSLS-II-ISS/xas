import numpy as np

def k2e(k, E0):
    """
    Convert from k-space to energy in eV

    Parameters
    ----------
    k : float
        k value
    E0 : float
        Edge energy in eV

    Returns
    -------
    out : float
        Energy value

    See Also
    --------
    :func:`isstools.conversions.xray.e2k`
    """
    return ((1000/(16.2009 ** 2)) * (k ** 2)) + E0

def e2k(E, E0):
    """
    Convert from energy in eV to k-space

    Parameters
    ----------
    E : float
        Current energy in eV
    E0 : float
        Edge energy in eV

    Returns
    -------
    out : float
        k-space value

    See Also
    --------
    :func:`isstools.conversions.xray.k2e`
    """
    return 16.2009 * (((E - E0)/1000) ** 0.5)

def encoder2energy(encoder, pulses_per_deg, offset = 0):
    """
    Convert from encoder counts to energy in eV

    Parameters
    ----------
    encoder : float or np.array()
        Encoder counts to convert
    pulses_per_deg: float
        Number of pulses per degree of the encoder
    offset : float
        Offset in degrees to adjust the conversion

    Returns
    -------
    out : float or np.array()
        Energy value or array of values

    See Also
    --------
    :func:`isstools.conversions.xray.energy2encoder`
    """
    return -12398.42 / (2 * 3.1356 * np.sin(np.deg2rad((encoder/pulses_per_deg) - float(offset))))

def energy2encoder(energy, pulses_per_deg, offset = 0):
    """
    Convert from energy in eV to encoder counts

    Parameters
    ----------
    energy : float or np.array()
        Energy in eV to convert
    offset : float
        Offset in degrees to adjust the conversion

    Returns
    -------
    out : float or np.array()
        Encoder counts value or array of values

    See Also
    --------
    :func:`isstools.conversions.xray.encoder2energy`
    """
    return pulses_per_deg * (np.degrees(np.arcsin(-12398.42 / (2 * 3.1356 * energy))) - float(offset))



def generate_energy_grid(e0, preedge_start, xanes_start, xanes_end, exafs_end, preedge_spacing,
                        xanes_spacing, exafsk_spacing, int_time_preedge = 1, int_time_xanes = 1, int_time_exafs = 1, k_power = 0):
    preedge = np.arange(e0 + preedge_start, e0 + xanes_start, preedge_spacing)
    preedge_int = np.ones(len(preedge)) * int_time_preedge

    edge = np.arange(e0 + xanes_start, e0 + xanes_end, xanes_spacing)
    edge_int = np.ones(len(edge)) * int_time_xanes

    iterator = exafsk_spacing
    kenergy = 0
    postedge = np.array([])

    energy_end = k2e(exafs_end, e0)
    exafs_int = []
    while(kenergy + e0 + xanes_end < energy_end):
        kenergy = k2e(iterator, e0) - e0

        postedge = np.append(postedge, e0 + xanes_end + kenergy)
        k_current = e2k(e0 + xanes_end + kenergy,e0)
        exafs_int.append(int_time_exafs * (k_current ** k_power))
        iterator += exafsk_spacing

    integration_times = np.append(np.append(preedge_int, edge_int), np.array(exafs_int))
    grid = np.append(np.append(preedge, edge), postedge)
    return grid, integration_times
    #return np.append(np.append(preedge, edge), postedge)







