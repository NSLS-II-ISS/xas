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

    # This is how it's defined in the IOC
    record(calcout, "XF:08IDA-OP{Mono:HHM-Ax:E}Mtr-SP") {
    field(INPB, "-1977.004107667")
    field(INPC, "3.141592653")
    field(INPE, "XF:08IDA-OP{Mono:HHM-Ax:E}Offset.VAL PP MS")

    "ASIN(B/A)*180/C - E")
    """
    return pulses_per_deg * (np.degrees(np.arcsin(-12398.42 / (2 * 3.1356 * energy))) - float(offset))

def energy2angle(energy,  offset = 0):
    return np.degrees(np.arcsin(-12398.42 / (2 * 3.1356 * energy))) - float(offset)


def generate_energy_grid(e0, preedge_start, xanes_start, xanes_end, exafs_end, preedge_spacing,
                         xanes_spacing, exafsk_spacing, dwell_time_preedge = 1, dwell_time_xanes = 1, dwell_time_exafs = 1, k_power = 0):
    preedge = np.arange(e0 + preedge_start, e0 + xanes_start, preedge_spacing)
    preedge_int = np.ones(len(preedge)) * dwell_time_preedge

    edge = np.arange(e0 + xanes_start, e0 + xanes_end, xanes_spacing)
    edge_int = np.ones(len(edge)) * dwell_time_xanes

    iterator = exafsk_spacing
    kenergy = 0
    postedge = np.array([])

    energy_end = k2e(exafs_end, e0)
    exafs_int = []
    while(kenergy + e0 + xanes_end < energy_end):
        kenergy = k2e(iterator, e0) - e0

        postedge = np.append(postedge, e0 + xanes_end + kenergy)
        k_current = e2k(e0 + xanes_end + kenergy,e0)
        exafs_int.append(dwell_time_exafs * (k_current ** k_power))
        iterator += exafsk_spacing

    integration_times = np.append(np.append(preedge_int, edge_int), np.array(exafs_int))
    grid = np.append(np.append(preedge, edge), postedge)
    return grid, integration_times
    #return np.append(np.append(preedge, edge), postedge)


def generate_energy_grid_from_dict(scan_parameters):
    energy, time_grid = generate_energy_grid(scan_parameters['e0'],
                                             scan_parameters['preedge_start'],
                                             scan_parameters['XANES_start'],
                                             scan_parameters['XANES_end'],
                                             scan_parameters['EXAFS_end'],
                                             scan_parameters['preedge_stepsize'],
                                             scan_parameters['XANES_stepsize'],
                                             scan_parameters['EXAFS_stepsize'],
                                             dwell_time_preedge= scan_parameters['preedge_dwelltime'],
                                             dwell_time_xanes= scan_parameters['XANES_dwelltime'],
                                             dwell_time_exafs= scan_parameters['EXAFS_dwelltime'],
                                             k_power = scan_parameters['k_power'])
    if scan_parameters['revert']:
        energy = energy[::-1]
        time_grid = time_grid[::-1]

    time = np.cumsum(time_grid)

    return energy, time_grid, time


def generate_emission_energy_grid(e0, preline_start, mainline_start, mainline_end, postline_end,
                                  preline_stepsize, mainline_stepsize, postline_stepsize,
                                  preline_dwelltime, mainline_dwelltime, postline_dwelltime, revert):

    energy_grid_preline = np.arange(e0 + preline_start, e0 + mainline_start, preline_stepsize)
    time_grid_preline = np.ones(energy_grid_preline.size) * preline_dwelltime

    energy_grid_mainline = np.arange(e0 + mainline_start, e0 + mainline_end, mainline_stepsize)
    time_grid_mainline = np.ones(energy_grid_mainline.size) * mainline_dwelltime

    energy_grid_postline = np.arange(e0 + mainline_end, e0 + postline_end + postline_stepsize, postline_stepsize)
    time_grid_postline = np.ones(energy_grid_postline.size) * postline_dwelltime

    energy_grid = np.hstack((energy_grid_preline, energy_grid_mainline, energy_grid_postline))
    time_grid = np.hstack((time_grid_preline, time_grid_mainline, time_grid_postline))

    if revert:
        energy_grid = energy_grid[::-1]
        time_grid = time_grid[::-1]



    return energy_grid, time_grid

def generate_emission_energy_grid_from_dict(scan_parameters):
    if 'preline_dwelltime' in scan_parameters.keys():
        preline_dwelltime = scan_parameters['preline_dwelltime']
    else:
        preline_dwelltime = 1

    if 'mainline_dwelltime' in scan_parameters.keys():
        mainline_dwelltime = scan_parameters['preline_dwelltime']
    else:
        mainline_dwelltime = 1

    if 'postline_dwelltime' in scan_parameters.keys():
        postline_dwelltime = scan_parameters['postline_dwelltime']
    else:
        postline_dwelltime = 1

    energy_grid, time_grid = generate_emission_energy_grid(scan_parameters['e0'],
                                                           scan_parameters['preline_start'],
                                                           scan_parameters['mainline_start'],
                                                           scan_parameters['mainline_end'],
                                                           scan_parameters['postline_end'],
                                                           scan_parameters['preline_stepsize'],
                                                           scan_parameters['mainline_stepsize'],
                                                           scan_parameters['postline_stepsize'],
                                                           preline_dwelltime,
                                                           mainline_dwelltime,
                                                           postline_dwelltime,
                                                           scan_parameters['revert'])
    return energy_grid, time_grid