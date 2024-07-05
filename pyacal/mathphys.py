"""Mathematical and Physical constants."""

import numpy as _np


# International System of Units
# =============================
# ref.: https://en.wikipedia.org/wiki/SI_base_unit


class BaseUnits:
    """."""
    # Base Units
    # ==========

    meter = 1.0
    kilogram = 1.0
    second = 1.0
    ampere = 1.0
    kelvin = 1.0
    mole = 1.0
    candela = 1.0


class Constants:
    """."""
    _u = BaseUnits
    # temporary auxiliary derived units
    _volt = (_u.kilogram * _u.meter**2) / (_u.ampere * _u.second**2)
    _coulomb = _u.second * _u.ampere
    _joule = _u.kilogram * _u.meter**2 / _u.second**2

    # physical constants
    # ==================

    # --- exact by definition --

    light_speed = 299792458 * (_u.meter / _u.second)
    gas_constant = 8.314462618 * (_joule / _u.mole / _u.kelvin)
    boltzmann_constant = 1.380649e-23 * (_joule / _u.kelvin)
    avogadro_constant = 6.02214076e23 * (1 / _u.mole)
    elementary_charge = 1.602176634e-19 * (_coulomb)
    reduced_planck_constant = 1.054571817e-34 * (_joule * _u.second)

    # --- measured ---

    # 2022-03-19 - https://physics.nist.gov/cgi-bin/cuu/Value?me|search_for=electron+mass
    electron_mass = 9.1093837015e-31 * (_u.kilogram)

    # 2022-03-19 - https://physics.nist.gov/cgi-bin/cuu/Value?mu0|search_for=vacuum+permeability
    vacuum_permeability = 1.25663706212e-6 * \
        (_volt * _u.second / _u.ampere / _u.meter)

    # --- derived ---

    # [Kg̣*m^2/s^2] - derived
    electron_rest_energy = electron_mass * _np.pow(light_speed, 2)

    # [V·s/(A.m)]  - derived
    vacuum_permitticity = 1.0/(vacuum_permeability * _np.pow(light_speed, 2))

    # [T·m^2/(A·s)] - derived
    vacuum_impedance = vacuum_permeability * light_speed

    # [m] - derived
    electron_radius = _np.pow(elementary_charge, 2) / \
        (4*_np.pi*vacuum_permitticity*electron_rest_energy)

    _joule_2_eV = _joule / elementary_charge

    # [m]/[GeV]^3 - derived
    rad_cgamma = 4*_np.pi*electron_radius / \
        _np.pow(electron_rest_energy/elementary_charge/1.0e9, 3) / 3

    # [m] - derived
    Cq = (55.0/(32*_np.sqrt(3.0))) * (reduced_planck_constant) * \
        light_speed / electron_rest_energy

    # [m^2/(s·GeV^3)] - derived
    Ca = electron_radius*light_speed / \
        (3*_np.pow(electron_rest_energy*_joule_2_eV/1.0e9, 3))


class DerivedUnits:
    """."""
    _u = BaseUnits
    _c = Constants
    newton = _u.kilogram * _u.meter / _u.second
    coulomb = _u.second * _u.ampere
    joule = newton * _u.meter
    watt = joule / _u.second
    volt = watt / _u.ampere
    weber = volt * _u.second
    tesla = weber / _u.meter**2
    pascal = _u.kilogram / (_u.meter * _u.second**2)

    radian = (_u.meter / _u.meter)
    (mA, uA) = (1e-3, 1e-6)
    (km, cm, mm, um, nm) = (1e3, 1e-2, 1e-3, 1e-6, 1e-9)
    (rad, mrad, urad, nrad) = (1e0, 1e-3, 1e-6, 1e-9)
    (minute, hour, day, year) = (60, 60*60, 24*60*60, 365.25*24*60*60)

    electron_volt = _c.elementary_charge * volt
    (eV, MeV, GeV) = (electron_volt, electron_volt*1e6, electron_volt*1e9)


class Conversions:
    """."""
    # conversions factors
    # ===================
    # conversion factors should be defined instead of conversion functions
    # whenever possible. The reason is that it is more general since a
    # conversion function over iterable would have to be defined, whereas some
    # iterables (from numpy, for example) defines multiplication by scalar.
    _u = BaseUnits
    _du = DerivedUnits

    radian_2_degree = (180.0/_np.pi)
    degree_2_radian = (_np.pi/180.0)

    rad_2_mrad = (_du.rad / _du.mrad)
    meter_2_mm = (_u.meter / _du.mm)
    joule_2_eV = (_du.joule / _du.electron_volt)
    pascal_2_bar = _du.pascal * 1.0e-5
    eV_2_GeV = (_du.eV / _du.GeV)
    joule_2_GeV = joule_2_eV * eV_2_GeV
    ev_2_joule = 1.0 / joule_2_eV
    GeV_2_eV = 1.0 / eV_2_GeV
    mrad_2_rad = 1.0 / rad_2_mrad
    mm_2_meter = 1.0 / meter_2_mm

    def beam_rigidity(energy):
    """Return beam rigidity, beta amd game, given its energy [GeV]."""
    brho, _, beta, gamma, _ = _beam.beam_rigidity(energy=energy)
    return brho, beta, gamma