# constants.py

"""
This module defines physical and cosmological constants used throughout the simulation.
"""

import numpy as np

# Physical constants
G = 6.67430e-11         # Gravitational constant (m^3 kg^-1 s^-2)
c = 2.99792458e8        # Speed of light in vacuum (m/s)
h_bar = 1.054571817e-34 # Reduced Planck constant (J·s)
ELECTRON_MASS = 9.10938356e-31  # Electron mass (kg)
PROTON_MASS = 1.6726219e-27     # Proton mass (kg)
mu_0 = 4 * np.pi * 1e-7         # Vacuum permeability (N/A^2)
epsilon_0 = 8.854187817e-12     # Vacuum permittivity (F/m)

# Cosmological constants
HUBBLE_CONSTANT = 2.3e-18       # Hubble constant (1/s)
DARK_MATTER_DENSITY = 0.3       # Fraction of critical density
DARK_ENERGY_DENSITY = 0.7       # Fraction of critical density
CRITICAL_DENSITY = 1.88e-26     # kg/m^3
L = 1.1056e-52             # Cosmological constant (1/m^2)

# Fundamental Physical Constants
h = h_bar * 2 * np.pi  # Planck constant (J⋅s)
k_B = 1.380649e-23  # Boltzmann constant (J/K)
e = 1.602176634e-19  # Elementary charge (C)
alpha = 1/137.035999084  # Fine structure constant

# Particle Physics Constants
NEUTRON_MASS = 1.67492749804e-27  # kg
ATOMIC_MASS_UNIT = 1.660539067e-27  # kg
HIGGS_MASS = 2.24e-25  # kg
TOP_QUARK_MASS = 3.0784e-25  # kg
Z_BOSON_MASS = 1.62557e-25  # kg
W_BOSON_MASS = 1.433e-25  # kg

# Astronomical Constants
SOLAR_MASS = 1.98847e30  # kg
EARTH_MASS = 5.97217e24  # kg
MOON_MASS = 7.342e22  # kg
MERCURY_MASS = 3.285e23  # kg
VENUS_MASS= 4.867e24  # kg
MARS_MASS = 6.39e23  # kg
JUPITER_MASS = 1.898e27  # kg
SATURN_MASS = 5.683e26  # kg
URANUS_MASS = 8.681e25  # kg
NEPTUNE_MASS = 1.024e26  # kg

AU = 1.495978707e11  # Astronomical Unit (m)
LIGHT_YEAR = 9.461e15  # m
PARSEC = 3.086e16  # m

# Additional Cosmological Constants
OMEGA_M = DARK_MATTER_DENSITY  # Matter density parameter
OMEGA_L = DARK_ENERGY_DENSITY  # Dark energy density parameter
OMEGA_K = 1.0 - OMEGA_M - OMEGA_L  # Curvature density parameter

# Derived Planck Units
PLANCK_LENGTH = np.sqrt(h_bar * G / c**3)  # m
PLANCK_MASS = np.sqrt(h_bar * c / G)  # kg
PLANCK_TIME = PLANCK_LENGTH / c  # s
PLANCK_TEMPERATURE = np.sqrt(h_bar * c**5 / (G * k_B**2))  # K
PLANCK_CHARGE = np.sqrt(4 * np.pi * epsilon_0 * h_bar * c)  # C

# Unit Conversions
def ev_to_joule(ev: float) -> float:
    """Convert electron volts to joules"""
    return ev * e

def kelvin_to_ev(T: float) -> float:
    """Convert temperature from Kelvin to electron volts"""
    return k_B * T / e

def parsec_to_meter(pc: float) -> float:
    """Convert parsecs to meters"""
    return pc * PARSEC

def solar_mass_to_kg(solar_masses: float) -> float:
    """Convert solar masses to kilograms"""
    return solar_masses * SOLAR_MASS

# Physical Calculations
def schwarzschild_radius(mass: float) -> float:
    """Calculate Schwarzschild radius for given mass"""
    return 2 * G * mass / c**2

def gravitational_binding_energy(mass: float, radius: float) -> float:
    """Calculate gravitational binding energy"""
    return 3 * G * mass**2 / (5 * radius)

def escape_velocity(mass: float, radius: float) -> float:
    """Calculate escape velocity"""
    return np.sqrt(2 * G * mass / radius)

# Orbital Parameters (m)
MERCURY_SEMI_MAJOR_AXIS = 5.790905e10
VENUS_SEMI_MAJOR_AXIS = 1.08208e11
EARTH_SEMI_MAJOR_AXIS = 1.495978707e11  # AU
MARS_SEMI_MAJOR_AXIS = 2.279391e11
JUPITER_SEMI_MAJOR_AXIS = 7.786235e11
SATURN_SEMI_MAJOR_AXIS = 1.43344e12
MOON_SEMI_MAJOR_AXIS = 3.844e8  # Relative to Earth

# Orbital Eccentricities
MERCURY_ECCENTRICITY = 0.20563069
VENUS_ECCENTRICITY = 0.00677323
EARTH_ECCENTRICITY = 0.01671022
MARS_ECCENTRICITY = 0.09341233
JUPITER_ECCENTRICITY = 0.04839266
SATURN_ECCENTRICITY = 0.05415060
MOON_ECCENTRICITY = 0.0549

# Orbital Periods (s)
MERCURY_PERIOD = 7.60052e6
VENUS_PERIOD = 1.941e7
EARTH_PERIOD = 3.15576e7  # 365.25 days
MARS_PERIOD = 5.93543e7
JUPITER_PERIOD = 3.74335e8
SATURN_PERIOD = 9.29424e8
MOON_PERIOD = 2.36068e6  # 27.322 days

# Mean Orbital Velocities (m/s)
def calculate_mean_orbital_velocity(semi_major_axis: float, mass_central: float) -> float:
    """Calculate mean orbital velocity using Kepler's laws"""
    return np.sqrt(G * mass_central / semi_major_axis)

MERCURY_VELOCITY = calculate_mean_orbital_velocity(MERCURY_SEMI_MAJOR_AXIS, SOLAR_MASS)
VENUS_VELOCITY = calculate_mean_orbital_velocity(VENUS_SEMI_MAJOR_AXIS, SOLAR_MASS)
EARTH_VELOCITY = calculate_mean_orbital_velocity(EARTH_SEMI_MAJOR_AXIS, SOLAR_MASS)
MARS_VELOCITY = calculate_mean_orbital_velocity(MARS_SEMI_MAJOR_AXIS, SOLAR_MASS)
JUPITER_VELOCITY = calculate_mean_orbital_velocity(JUPITER_SEMI_MAJOR_AXIS, SOLAR_MASS)
SATURN_VELOCITY = calculate_mean_orbital_velocity(SATURN_SEMI_MAJOR_AXIS, SOLAR_MASS)
MOON_VELOCITY = calculate_mean_orbital_velocity(MOON_SEMI_MAJOR_AXIS, EARTH_MASS)

# Escape Velocities (m/s)
def calculate_escape_velocity(mass: float, radius: float) -> float:
    """Calculate escape velocity from a body"""
    return np.sqrt(2 * G * mass / radius)

# Physical Radii (m)
SOLAR_RADIUS = 6.957e8
EARTH_RADIUS = 6.371e6
MOON_RADIUS = 1.737e6
MERCURY_RADIUS = 2.440e6
VENUS_RADIUS = 6.052e6
MARS_RADIUS = 3.390e6
JUPITER_RADIUS = 6.991e7
SATURN_RADIUS = 5.823e7

# Derived Quantities
def calculate_surface_gravity(mass: float, radius: float) -> float:
    """Calculate surface gravity of a body"""
    return G * mass / radius**2

def calculate_orbital_energy(mass: float, semi_major_axis: float, mass_central: float) -> float:
    """Calculate orbital energy"""
    return -G * mass * mass_central / (2 * semi_major_axis)

def calculate_angular_momentum(mass: float, semi_major_axis: float, mass_central: float, eccentricity: float) -> float:
    """Calculate specific angular momentum"""
def calculate_angular_momentum(semi_major_axis: float, mass_central: float, eccentricity: float) -> float:
    """Calculate specific angular momentum"""
    return np.sqrt(G * mass_central * semi_major_axis * (1 - eccentricity**2))
AU_TO_METERS = 1.495978707e11
LIGHT_YEAR_TO_METERS = 9.461e15
PARSEC_TO_METERS = 3.086e16

# Ensure no imports from 'galaxy.py' or 'particles.py'
