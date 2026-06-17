# forces.py

import numpy as np
from typing import Optional, Dict
from constants import (G, L, c, h_bar, mu_0, epsilon_0)
from dataclasses import dataclass, asdict
import json

def hubble_expansion(a: float, H_0: float) -> float:
    """Calculate Hubble parameter with dark energy"""
    return H_0 * np.sqrt(a * (1 + L / 3))

def gravitational_force(m1: float, m2: float, r_vec: np.ndarray) -> np.ndarray:
    """Calculate gravitational force with full relativistic corrections"""
    r = np.linalg.norm(r_vec)
    if r == 0:
        return np.zeros(3)
    
    # Schwarzschild radius check
    r_s = 2 * G * (m1 + m2) / c**2
    if r < r_s:
        return np.zeros(3)  # Inside event horizon
        
    # Post-Newtonian corrections up to 2PN order
    pn1 = 1 + (3 * G * (m1 + m2)) / (r * c**2)  # 1PN
    pn2 = 1 + (41 * G**2 * (m1 + m2)**2) / (4 * r**2 * c**4)  # 2PN
    
    force_magnitude = G * m1 * m2 / r**2 * pn1 * pn2
    return force_magnitude * r_vec / r

def electromagnetic_force(q1: float, q2: float, r_vec: np.ndarray, 
                        v1: np.ndarray, v2: np.ndarray) -> np.ndarray:
    """Calculate electromagnetic force including magnetic effects"""
    r = np.linalg.norm(r_vec)
    if r == 0:
        return np.zeros(3)
    
    # Electric force (Coulomb)
    E_force = q1 * q2 / (4 * np.pi * epsilon_0 * r**2) * r_vec / r
    
    # Magnetic force (Biot-Savart)
    v_rel = v1 - v2
    B_force = mu_0 * q1 * q2 / (4 * np.pi * r**2) * np.cross(v_rel, r_vec)
    
    return E_force + B_force / c**2

def dark_matter_effect(position: np.ndarray, mass: float, 
                      dark_matter_density: float) -> np.ndarray:
    """Calculate dark matter force with NFW profile"""
    r = np.linalg.norm(position)
    if r == 0:
        return np.zeros(3)
    
    # NFW profile parameters
    rs = 20000 * 3.086e16  # Scale radius (20 kpc)
    rho_s = dark_matter_density  # Characteristic density
    
    # NFW density profile
    x = r / rs
    
    # Enclosed mass and force
    mass_enclosed = 4 * np.pi * rho_s * rs**3 * (np.log(1 + x) - x / (1 + x))
    force_magnitude = G * mass * mass_enclosed / r**2
    
    return -force_magnitude * position / r

def quantum_force(mass: float, position: np.ndarray, 
                 characteristic_length: float) -> np.ndarray:
    """Calculate quantum force from uncertainty principle"""
    r = np.linalg.norm(position)
    if r == 0:
        return np.zeros(3)
        
    # Heisenberg uncertainty principle based force
    force_magnitude = h_bar**2 / (mass * characteristic_length**3)
    return force_magnitude * position / r

def radiation_pressure(luminosity: float, r_vec: np.ndarray,
                      opacity: float = 1.0,
                      temperature: float = 5800.0) -> np.ndarray:
    """Calculate radiation pressure force"""
    r = np.linalg.norm(r_vec)
    if r == 0:
        return np.zeros(3)
    
    # Stefan-Boltzmann radiation
    sigma = 5.67e-8  # Stefan-Boltzmann constant
    pressure = (sigma * temperature**4) / c
    
    # Distance attenuation
    pressure *= luminosity / (4 * np.pi * r**2)
    
    # Opacity correction
    pressure *= opacity
    
    return pressure * r_vec / r

def total_force(particle1: Dict, particle2: Dict,
                dark_matter_density: float = 1e-21) -> np.ndarray:
    """Calculate total force between particles including all effects"""
    r_vec = particle2['position'] - particle1['position']
    
    # Gravitational force
    f_grav = gravitational_force(particle1['mass'], particle2['mass'], r_vec)
    
    # Add electromagnetic force if charges present
    f_em = np.zeros(3)
    if 'charge' in particle1 and 'charge' in particle2:
        f_em = electromagnetic_force(
            particle1['charge'], particle2['charge'],
            r_vec, particle1['velocity'], particle2['velocity']
        )
    
    # Add dark matter effect
    f_dm = dark_matter_effect(r_vec, particle1['mass'], dark_matter_density)
    
    # Add radiation pressure if luminosity present
    f_rad = np.zeros(3)
    if 'luminosity' in particle2:
        f_rad = radiation_pressure(particle2['luminosity'], r_vec)
    
    return f_grav + f_em + f_dm + f_rad

@dataclass
class SimulationConfig:
    # Required fields (no defaults)
    time_step: float
    total_steps: int
    dark_matter_density: float
    # Optional fields (with defaults)
    enable_relativity: bool = True
    enable_quantum_effects: bool = False
    enable_dark_energy: bool = True
    max_simulation_radius: float = 1e6  # Default value, can be adjusted
    spatial_resolution: float  # Define spatial_resolution
    current_time: float = 0.0
    hubble_constant: float = 70.0  # Default Hubble constant
    memory_limit: float = 16.0  # Default memory limit in GB
    
    visualization_config: Optional[Dict] = None
    output_config: Optional[Dict] = None

    def __post_init__(self):
        if self.visualization_config is None:
            self.visualization_config = {
                "enable": True,
                "trail_length": 50,
                "update_interval": 10
            }
        
        if self.output_config is None:
            self.output_config = {
                "save_interval": 100,
                "output_directory": "simulation_data",
                "format": "json"
            }

    def validate(self) -> bool:
        """Validate configuration parameters"""
        try:
            assert self.time_step > 0
            assert self.total_steps > 0
            assert self.dark_matter_density >= 0
            assert self.hubble_constant > 0
            assert self.spatial_resolution > 0
            assert self.max_simulation_radius > self.spatial_resolution
            assert self.memory_limit > 0
            return True
        except AssertionError as e:
            raise ValueError(f"Invalid configuration parameters: {e}")
    
    def estimate_memory_usage(self) -> float:
        """Estimate memory usage in GB"""
        particles_memory = (self.max_simulation_radius / self.spatial_resolution)**3 * 24  # bytes per particle
        return particles_memory * 1e-9
    
    def to_dict(self) -> dict:
        """Convert configuration to dictionary"""
        return asdict(self)
    
    @classmethod
    def from_dict(cls, config_dict: dict) -> 'SimulationConfig':
        """Create configuration from dictionary"""
        return cls(**config_dict)
    
    def save_config(self, filename: str):
        """Save configuration to file"""
        with open(filename, 'w') as f:
            json.dump(self.to_dict(), f)

    @classmethod
    def load_config(cls, filename: str) -> 'SimulationConfig':
        with open(filename, 'r') as f:
            config_dict = json.load(f)
        return cls.from_dict(config_dict)
@dataclass
class Forces:
    """Class representing combined forces acting on a particle.
    
    Attributes:
        gravitational: Gravitational force vector
        electromagnetic: Electromagnetic force vector
        dark_matter: Dark matter force vector
        radiation_pressure: Radiation pressure force vector
        name: Identifier name
        magnitude: Total force magnitude
        direction: Force direction unit vector
        memory_limit: Memory limit in bytes
        some_default_arg: Optional parameter
        another_default_arg: Optional boolean flag
    """
    gravitational: np.ndarray
    electromagnetic: np.ndarray  
    dark_matter: np.ndarray
    radiation_pressure: np.ndarray
    name: str
    magnitude: float
    direction: np.ndarray
    memory_limit: int
    some_default_arg: float = 0.0
    another_default_arg: bool = True
    
    def __post_init__(self):
        # Validate input arrays
        for field in ['gravitational', 'electromagnetic', 'dark_matter', 'radiation_pressure', 'direction']:
            value = getattr(self, field)
            if not isinstance(value, np.ndarray) or value.shape != (3,):
                raise ValueError(f"{field} must be a 3D numpy array")
        
        if not isinstance(self.name, str) or not self.name:
            raise ValueError("name must be a non-empty string")
            
        if not isinstance(self.magnitude, (int, float)) or self.magnitude < 0:
            raise ValueError("magnitude must be a non-negative number")
            
        if not isinstance(self.memory_limit, int) or self.memory_limit <= 0:
            raise ValueError("memory_limit must be a positive integer")
