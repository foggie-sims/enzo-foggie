"""
Read in pySB99 files and combine them into an hdf5 file of the format Enzo
expects for pre-SN feedback and local stellar radiation
(PreSNFeedback / UseLocalStellarRadiation, see src/enzo/ReadPreSNFeedbackTable.C).

Script version of combine_pySB_hdf5_with_RT_Fields.ipynb, without the plots.

Every output array has first index population metallicity and second index
population age. All quantities are smoothed in time with a moving window
(the raw models have spikes that can cause issues with interpolation) and
normalized from the 10^6 Msun SB99 cluster to 1 Msun.

This assumes Starburst99 has already been run for each metallicity in
INITIAL_METAL_FRACTION, with the results in the directories listed in
MODELDIRS. The table was originally generated with pySB99, the python wrapper
for Starburst99 (https://github.com/CalumHawcroft/Starburst).

Usage:
    python combine_pySB_hdf5_with_RT_Fields.py [--wind-dir DIR] [--sed-dir DIR] [--output FILE]
"""

import argparse

import h5py
import numpy as np
import scipy.integrate
import scipy.interpolate

# Make sure this lines up with how SB99 was run
POPULATION_AGE = np.arange(0.00e6, 50e6, 0.1e6)  # yr
# Make sure this lines up with the different-metallicity SB99 runs. Note that
# these are actual metal fractions so 0.02 is solar.
INITIAL_METAL_FRACTION = np.array([0.0, 0.0004, 0.002, 0.006, 0.014, 0.02])
MODELDIRS = ['pySB_Z0', 'pySB_Z0004', 'pySB_Z002', 'pySB_Z006', 'pySB_Z014', 'pySB_Z02']

CLUSTER_MASS = 1e6  # Msun; change this if SB99 tables were run with a different mass
WINDOW_SIZE = 10

# Constants
light_c = 2.998e10
SolarMass = 1.989e33
yr_s = 3.155e7
h_cgs = 6.62607e-27  # erg * s
c_cgs = 2.99792458e10  # cm /s
h_ev = 4.135667696e-15  # eV s
angstrom2cgs = 1e-8


def smooth(arr):
    """Moving-window average in time, edge-padded to keep the length."""
    window = np.ones(WINDOW_SIZE) / WINDOW_SIZE
    return np.convolve(np.pad(arr, WINDOW_SIZE // 2, mode='edge'), window, mode='valid')[:-1]


def event_rate(sed_with_nebular, sed_wavelength, mask, sigma):
    """Integrate photon counts * cross section over the masked band, for each age.

    sigma is either a scalar or an array over the masked wavelengths.
    Returns units of cm^2/s.
    """
    erg_axis = h_cgs * c_cgs / (sed_wavelength[mask] * angstrom2cgs)  # erg/photon
    counts = 10**sed_with_nebular[:, mask] / erg_axis  # erg s-1/A / (erg/photon)
    return scipy.integrate.simpson(counts * sigma, x=sed_wavelength[mask], axis=-1)


def get_average_draine_kappa(sed_with_nebular, sed_wavelength, mask, kappa,
                             kappa_wavelength, sigma=None):
    """Average Draine kappa over the masked band, for each age.

    sigma=None: luminosity weighted, for the Habing field.
    sigma given (scalar or array over the mask): weighted by sigma * photon
    counts, i.e. the same integrand as event_rate, for photo-rates.
    Returns cm^2/g of gas.
    """
    f = scipy.interpolate.interp1d(kappa_wavelength, kappa, kind='linear',
                                   bounds_error=False, fill_value='extrapolate')
    kappa_interp = f(sed_wavelength[mask])  # kappa at the SED wavelengths

    L_lambda = 10**sed_with_nebular[:, mask]  # erg/s/A
    if sigma is None:
        weight = L_lambda
    else:
        erg_axis = h_cgs * c_cgs / (sed_wavelength[mask] * angstrom2cgs)
        weight = L_lambda / erg_axis * sigma  # photons/s/A * cm^2

    num = scipy.integrate.simpson(kappa_interp * weight, x=sed_wavelength[mask], axis=-1)
    den = scipy.integrate.simpson(weight, x=sed_wavelength[mask], axis=-1)
    return num / den


def read_draine_kappa(path, kind="eff"):
    """Read a Draine (2003) kext_albedo_*.all dust table.

    Returns wavelength [micron], photon energy [eV] and dust opacity per
    gram of gas [cm^2/g], sorted by increasing wavelength.
    kind: "ext" (extinction), "abs" (absorption) or
          "eff" (absorption + (1 - <cos>) * scattering).
    """
    with open(path) as f:
        lines = f.read().splitlines()

    # Dust and gas mass per H nucleon, from the header.
    for line in lines:
        if "M_dust per H nucleon" in line:
            mdust_per_H = float(line.split()[0])
        if "M_gas/M_dust" in line:
            gas_to_dust = float(line.split()[0])
    mgas_per_H = mdust_per_H * gas_to_dust  # g of gas per H nucleon

    # Data rows start after the dashed line; columns 0-4 are numeric,
    # some rows carry a trailing comment (filter names) that is ignored.
    start = next(i for i, l in enumerate(lines) if l.startswith("-----------")) + 1
    rows = [l.split()[:5] for l in lines[start:] if l.strip()]
    lam, albedo, mean_cos, cext_per_H, _ = np.array(rows, dtype=float).T

    # Filter rows are inserted out of order: sort and drop repeats.
    lam, idx = np.unique(lam, return_index=True)
    albedo, mean_cos, cext_per_H = albedo[idx], mean_cos[idx], cext_per_H[idx]

    kappa_ext = cext_per_H / mgas_per_H                    # cm^2 per g of gas
    kappa_abs = kappa_ext * (1 - albedo)
    kappa_sca = kappa_ext * albedo
    kappa = {"ext": kappa_ext,
             "abs": kappa_abs,
             "eff": kappa_abs + (1 - mean_cos) * kappa_sca}[kind]

    return lam*1e4, kappa #wavelength in Angstrom, kappa


def process_model(wind_dir, sed_dir, draine_table):
    """Return a dict of smoothed, per-Msun time series for one metallicity."""
    out = {}

    ### Winds ###
    mass = np.loadtxt(wind_dir + '/windmass.txt')  # Wind mass is in units of Msun/yr
    metals = np.loadtxt(wind_dir + '/windmetals.txt')  # Wind metal mass is in units of Msun/yr
    momentum = 10**np.loadtxt(wind_dir + '/windmom.txt')  # Wind momentum is actually pdot, in units of dynes (g cm/s^2)
    Lbol = 10**np.loadtxt(wind_dir + '/bololum.txt')  # Bolometric luminosity is in units of ergs/s
    out['wind_mass'] = smooth(mass) / CLUSTER_MASS
    out['wind_metals'] = smooth(metals) / CLUSTER_MASS
    # Convert pdot to Msun/yr*km/s and normalize by star cluster size
    out['wind_mom'] = (smooth(momentum) / CLUSTER_MASS + smooth(Lbol) / light_c / CLUSTER_MASS) / (SolarMass / yr_s * 1e5)

    ### SEDs ###
    sed = np.load(sed_dir + "/pySB_SED_stellar.npy")
    sed_with_nebular = np.load(sed_dir + "/pySB_SED_stellar_and_nebular.npy")
    sed_wavelength = np.loadtxt(sed_dir + '/SED_wavelength.txt')
    print("SED wavelength ranges from", np.min(sed_wavelength), np.max(sed_wavelength))

    ### Opacities from Draine 2003 ###
    kappa_wavelength, kappa = read_draine_kappa(draine_table,'eff')

    ### H2 Photodissociation ###
    lw_mask = (sed_wavelength <= 1100.8) & (sed_wavelength >= 910.2)
    lw_sigma = 3.71e-18  # cm^2
    out['lw_events'] = smooth(event_rate(sed_with_nebular, sed_wavelength, lw_mask, lw_sigma)) / CLUSTER_MASS
    lw_kappa = get_average_draine_kappa(sed_with_nebular, sed_wavelength, lw_mask, kappa, kappa_wavelength,sigma=lw_sigma)
    print("lw_kappa is:",np.mean(lw_kappa))
    out['lw_kappa'] = smooth(lw_kappa)

    ### H- Photodetachment ###
    nu_th_hm = 0.755 / h_ev  # 0.755 eV threshold
    hm_max_wavelength = 3e8 / nu_th_hm * 10**10  # ~16400 Angstrom
    hm_mask = (sed_wavelength <= hm_max_wavelength) & (sed_wavelength >= 1250)
    sed_nu = 3e8 / (sed_wavelength[hm_mask] * 1e-10)  # Hz
    hm_sigma = 7.928e5 * np.power(sed_nu - nu_th_hm, 1.5) / np.power(sed_nu, 3.)
    out['hm_events'] = smooth(event_rate(sed_with_nebular, sed_wavelength, hm_mask, hm_sigma)) / CLUSTER_MASS
    hm_kappa = get_average_draine_kappa(sed_with_nebular, sed_wavelength, hm_mask, kappa, kappa_wavelength,sigma=hm_sigma)
    print("hm_kappa is:",np.mean(hm_kappa))
    out['hm_kappa'] = smooth(hm_kappa)

    ### CO Photodissociation ###
    co_mask = (sed_wavelength <= 1118) & (sed_wavelength >= 885)
    co_sigma = 1.7455e-17  # cm^2
    out['co_events'] = smooth(event_rate(sed_with_nebular, sed_wavelength, co_mask, co_sigma)) / CLUSTER_MASS
    out['co_kappa'] = smooth(get_average_draine_kappa(sed_with_nebular, sed_wavelength, co_mask, kappa, kappa_wavelength,sigma=co_sigma))

    # log10 of photon energy in eV, used for the CI and OI cross section fits
    x = np.log10(h_ev * c_cgs / (sed_wavelength * angstrom2cgs))

    ### CI - CII Photoionization ###
    ci_mask = sed_wavelength <= 1110  # Ionization potential
    ci_sigma = np.zeros_like(x)
    m0 = x < 1.05
    m1 = (x >= 1.05) & (x < 2.48)
    m2 = x >= 2.48
    ci_sigma[m0] = -20.
    ci_sigma[m1] = -16.332 + 0.320 * x[m1] - 0.636 * x[m1] * x[m1]
    ci_sigma[m2] = -14.549 - 0.391 * x[m2] - 0.403 * x[m2] * x[m2]
    ci_sigma = np.power(10., ci_sigma)  # cm^2
    out['ci_events'] = smooth(event_rate(sed_with_nebular, sed_wavelength, ci_mask, ci_sigma[ci_mask])) / CLUSTER_MASS
    out['ci_kappa'] = smooth(get_average_draine_kappa(sed_with_nebular, sed_wavelength, ci_mask, kappa, kappa_wavelength, sigma=ci_sigma[ci_mask]))

    ### OI - OII Photoionization ###
    oi_mask = sed_wavelength <= 910  # Ionization potential
    oi_sigma = np.zeros_like(x)
    m0 = x < 1.13
    m1 = (x >= 1.13) & (x < 1.22)
    m2 = (x >= 1.22) & (x < 1.26)
    m3 = (x >= 1.26) & (x < 1.43)
    m4 = (x >= 1.43) & (x < 1.73)
    m5 = x >= 1.73
    oi_sigma[m0] = -20.
    oi_sigma[m1] = -17.400
    oi_sigma[m2] = -32.040 + 12.000 * x[m2]
    oi_sigma[m3] = -16.968 + 0.038 * x[m3]
    oi_sigma[m4] = -15.341 - 1.100 * x[m4]
    oi_sigma[m5] = -12.054 - 3.000 * x[m5]
    oi_sigma = np.power(10., oi_sigma)  # cm^2
    out['oi_events'] = smooth(event_rate(sed_with_nebular, sed_wavelength, oi_mask, oi_sigma[oi_mask])) / CLUSTER_MASS
    out['oi_kappa'] = smooth(get_average_draine_kappa(sed_with_nebular, sed_wavelength, oi_mask, kappa, kappa_wavelength, sigma=oi_sigma[oi_mask]))

    ### Habing Field ###
    habing_mask = (sed_wavelength <= 2400.0) & (sed_wavelength >= 910.2)
    habing0 = 1.6e-3  # erg/s / cm^2 / G0
    habing_luminosity = scipy.integrate.simpson(10**sed_with_nebular[:, habing_mask],
                                                x=sed_wavelength[habing_mask], axis=-1)  # erg/s
    out['habing'] = smooth(habing_luminosity / habing0) / CLUSTER_MASS  # G0 cm^2
    habing_kappa = get_average_draine_kappa(sed_with_nebular, sed_wavelength, habing_mask, kappa, kappa_wavelength)
    print("habing_kappa=",np.mean(habing_kappa))
    out['habing_kappa'] = smooth(habing_kappa)

    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--wind-dir', default='/Users/ctrapp/Downloads/Starburst/',
                        help='Directory holding the per-metallicity wind/bololum text files')
    parser.add_argument('--sed-dir', default='/Users/ctrapp/Documents/GitHub/Starburst/',
                        help='Directory holding the per-metallicity pySB SED files')
    parser.add_argument('--output', default='preSN_feedback_SB99_RT.hdf5',
                        help='Output hdf5 table')
    parser.add_argument('--draine-table', default='/Users/ctrapp/Downloads/Starburst/kext_albedo_WD_MW_3.1_60_D03.all',
                        help='Directory holding drain dust absoprtions for birth cloud kappas')
    
    args = parser.parse_args()

    print(len(POPULATION_AGE), len(INITIAL_METAL_FRACTION))

    # Process each metallicity's SB99 run. Each entry of `models` is a dict
    # of 1D time series (one value per population age), one dict per metallicity.
    models = []
    for modeldir in MODELDIRS:
        wind_dir = args.wind_dir + modeldir
        sed_dir = args.sed_dir + modeldir
        models.append(process_model(wind_dir, sed_dir, args.draine_table))

    def stack_over_metallicity(key):
        """Combine one quantity from every metallicity into a 2D array
        of shape (n_metallicity, n_age), the layout Enzo expects."""
        rows = []
        for model in models:
            rows.append(model[key])
        return np.array(rows)

    # Save the arrays into an hdf5 file with the group and dataset structure Enzo expects.
    with h5py.File(args.output, 'w') as f:
        grp = f.create_group("indexer")
        grp.create_dataset("initial_metal_fraction", data=INITIAL_METAL_FRACTION)
        dset_age = grp.create_dataset("population_age", data=POPULATION_AGE)
        dset_age.attrs['time_unit'] = 'yr'

        grp = f.create_group("SB99_models")
        datasets = [
            ("wind_mass_rate",            'wind_mass',   'Msun/yr'),
            ("wind_metal_mass_rate",      'wind_metals', 'Msun/yr'),
            ("wind_and_Lbol_momentum",    'wind_mom',    'Msun/yr*km/s'),
            ("h2_photodissociation_rate", 'lw_events',   'cm**2/s'),
            ("h2_photodissociation_kappa",'lw_kappa' ,   'cm**2/g'),
            ("hm_photodetatchment_rate",  'hm_events',   'cm**2/s'),
            ("hm_photodetatchment_kappa", 'hm_kappa' ,   'cm**2/g'),
            ("co_photodissociation_rate", 'co_events',   'cm**2/s'),
            ("co_photodissociation_kappa",'co_kappa' ,   'cm**2/g'),
            ("ci_photoionization_rate",   'ci_events',   'cm**2/s'),
            ("ci_photoionization_kappa",  'ci_kappa' ,   'cm**2/g'),
            ("oi_photoionization_rate",   'oi_events',   'cm**2/s'),
            ("oi_photoionization_kappa",  'oi_kappa' ,   'cm**2/g'),
            ("habing_luminosity",         'habing',      'G0*cm**2'),
            ("habing_kappa",              'habing_kappa','cm**2/g')
        ]
        for name, key, unit in datasets:
            dset = grp.create_dataset(name, data=stack_over_metallicity(key))
            dset.attrs['unit'] = unit

    with h5py.File(args.output, 'r') as f:
        print(list(f['SB99_models'].keys()))


if __name__ == '__main__':
    main()
