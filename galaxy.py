# galaxy.py

from constants import (G, PARSEC, SOLAR_MASS, MERCURY_MASS, VENUS_MASS, 
                      EARTH_MASS, MOON_MASS, MARS_MASS, JUPITER_MASS, CRITICAL_DENSITY)
from typing import Optional, Tuple
import numpy as np
from dataclasses import dataclass
from particles import Particle

@dataclass
class GalaxyProperties:
    """Class to store properties of a galaxy."""
    mass: float  # in kilograms
    radius: float  # in meters
    luminosity: float  # in watts
    dark_matter_density: float  # in kg/m^3
    dark_energy_density: float  # in kg/m^3
    scale_radius: float  # in meters

    @classmethod
    def default_properties(cls):
        """Initialize with default properties."""
        return cls(
            mass=1.0e42,  # Example mass in kg
            radius=5.0e20,  # Example radius in meters
            luminosity=3.828e26,  # Solar luminosity in watts
            dark_matter_density=1.0e-21,  # Example dark matter density in kg/m^3
            dark_energy_density=0.7 * CRITICAL_DENSITY,  # Assuming DARK_ENERGY_DENSITY is 0.7
            scale_radius=1.0e21  # Example scale radius in meters
        )

class Galaxy:
    def __init__(self, properties: GalaxyProperties, scale_radius: float, position: np.ndarray):
        self.properties = properties
        self.scale_radius = scale_radius
        self.position = np.array(position, dtype=np.float64)
        # Delay import to prevent circular import
        # from particles import Particle
        # Initialize particles as needed within methods

    def spiral_arm_phase(self, r: float) -> float:
        """Calculate the spiral arm phase based on radius."""
        if r <= 0:
            raise ValueError("Radius must be positive for logarithm calculation.")
        ln_r = np.log(r / self.scale_radius)
        theta = ln_r / np.tan(self.properties.pitch_angle)
        return theta * self.properties.spiral_arms

    def spiral_arm_density(self, r: float, theta: float) -> float:
        """Calculate spiral arm density enhancement."""
        if r <= 0:
            return 0.0  # No enhancement at the galaxy center
        phase = self.spiral_arm_phase(r)
        return np.cos(self.properties.spiral_arms * theta - phase)

    def _gas_density(self, r: float, theta: float, z: float) -> float:
        """Calculate gas density at a given cylindrical coordinate."""
        density = (self.properties.gas_fraction / (4 * np.pi * self.properties.radius**3)) * \
                  np.exp(-r / self.scale_radius) * \
                  (1 + self.spiral_arm_density(r, theta))
        # Include vertical (z) dependence
        density *= np.exp(-abs(z) / self.properties.scale_height)
        return density

    def _cylindrical_coordinates(self, position: np.ndarray) -> Tuple[float, float, float]:
        """Convert Cartesian coordinates to cylindrical coordinates."""
        x, y, z = position
        r = np.linalg.norm([x, y])  # Cylindrical radius
        theta = np.arctan2(y, x)
        # Constrain r and z to reasonable ranges
        r = np.clip(r, 0, 1e3 * self.scale_radius)
        z = np.clip(z, -1e3 * self.properties.scale_height, 1e3 * self.properties.scale_height)
        return r, theta, z

    def _field_components(self, B_strength: float, r: float, theta: float, z: float) -> np.ndarray:
        """
        Calculate the magnetic field components in cylindrical coordinates.
        Returns the field in Cartesian coordinates.
        """
        if r <= 0:
            raise ValueError("Radius must be positive for logarithm calculation.")
        ln_r = np.log(r / self.scale_radius)

        # Example calculations for B_r and B_theta
        B_r = B_strength * np.cos(ln_r)
        B_theta = B_strength * np.sin(ln_r)

        # Calculate B_z with vertical dependence
        B_z = B_strength * np.tanh(z / self.properties.scale_height) * 0.1

        # Convert cylindrical B-field components to Cartesian coordinates
        B_x = B_r * np.cos(theta) - B_theta * np.sin(theta)
        B_y = B_r * np.sin(theta) + B_theta * np.cos(theta)

        return np.array([B_x, B_y, B_z])

    def star_formation_density(self, position: np.ndarray) -> float:
        """Calculate star formation rate density at a given position."""
        r, theta, z = self._cylindrical_coordinates(position)
        base_rate = self._base_star_formation_rate(r, z)
        arm_enhancement = self.spiral_arm_density(r, theta)
        schmidt_law = self._schmidt_kennicutt_law(r)
        return base_rate * arm_enhancement * schmidt_law

    def _base_star_formation_rate(self, r: float, z: float) -> float:
        """Calculate the base star formation rate."""
        # Implement the base star formation rate logic here
        # Placeholder implementation:
        return 1.0

    def _schmidt_kennicutt_law(self, r: float) -> float:
        """Implement the Schmidt-Kennicutt law for star formation."""
        # Implement the Schmidt-Kennicutt law here
        # Placeholder implementation:
        return 1.0

def initialize_solar_system():
    """Initialize solar system bodies with slight z-components to avoid coplanar points"""
    return [
        Particle(SOLAR_MASS, [0, 0, 0], [0, 0, 0], "Sun"),
        Particle(MERCURY_MASS, [5.7e10, 0, 1e3], [0, 4.7e4, 0], "Mercury"),
        Particle(VENUS_MASS, [1.1e11, 0, -1e3], [0, 3.5e4, 0], "Venus"),
        Particle(EARTH_MASS, [1.496e11, 0, 2e3], [0, 2.978e4, 0], "Earth"),
        Particle(MOON_MASS, [1.496e11 + 3.844e8, 0, -2e3], [0, 3.0e4, 0], "Moon"),
        Particle(MARS_MASS, [2.279e11, 0, 1.5e3], [0, 2.41e4, 0], "Mars"),
        Particle(JUPITER_MASS, [7.786e11, 0, -1.5e3], [0, 1.307e4, 0], "Jupiter")
    ]
