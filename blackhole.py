# blackhole.py

import numpy as np
from typing import Tuple, Optional, Dict
from constants import G, c, h_bar, k_B, SOLAR_MASS

class BlackHole:
    def __init__(self, mass: float, position: np.ndarray, 
                 angular_momentum: Optional[np.ndarray] = None,
                 charge: float = 0.0):
        self.mass = mass
        self.position = np.array(position, dtype=np.float64)
        self.angular_momentum = np.array(angular_momentum, dtype=np.float64) if angular_momentum is not None else np.zeros(3)
        self.charge = charge
        self.accretion_disk_mass = 0.0
        self.time = 0.0
        self.history = []
        
    @property
    def schwarzschild_radius(self) -> float:
        """Calculate Schwarzschild radius"""
        return 2 * G * self.mass / c**2
    
    @property
    def hawking_temperature(self) -> float:
        """Calculate Hawking radiation temperature"""
        return h_bar * c**3 / (8 * np.pi * G * self.mass * k_B)
    
    @property
    def spin_parameter(self) -> float:
        """Calculate dimensionless spin parameter a = J/(Mc)"""
        return np.linalg.norm(self.angular_momentum) / (self.mass * c)
    
    @property
    def event_horizon_radius(self) -> float:
        """Calculate event horizon radius for Kerr black hole"""
        a = self.spin_parameter
        return self.mass * (1 + np.sqrt(1 - a**2))
    
    @property
    def ergosphere_radius(self) -> float:
        """Calculate ergosphere radius at equator"""
        return 2 * self.schwarzschild_radius
    
    def kerr_metric(self, r: float, theta: float) -> Tuple[float, float]:
        """Calculate Kerr metric components for rotating black hole"""
        a = self.spin_parameter
        if a > self.mass:  # Check for naked singularity
            a = self.mass  # Limit to maximum allowed rotation
            
        # Metric components
        rho2 = r**2 + (a * np.cos(theta))**2
        g_tt = -(1 - 2*G*self.mass*r/(c**2 * rho2))
        g_pp = (r**2 + a**2 + 2*G*self.mass*a**2*r*np.sin(theta)**2/(c**2 * rho2)) * np.sin(theta)**2
        
        return g_tt, g_pp
    
    def calculate_hawking_radiation(self, dt: float) -> float:
        """Calculate mass loss due to Hawking radiation"""
        # Stefan-Boltzmann constant
        sigma = 5.67e-8
        
        # Surface area of event horizon
        A = 4 * np.pi * self.schwarzschild_radius**2
        
        # Power radiated
        P = sigma * A * self.hawking_temperature**4
        
        # Mass loss
        dM = P * dt / c**2
        return dM
    
    def add_accretion_mass(self, dm: float):
        """Add mass through accretion with relativistic corrections"""
        self.accretion_disk_mass += dm
        
        # Calculate accretion efficiency (depends on spin)
        a = self.spin_parameter
        efficiency = 0.057 + 0.15 * a  # Efficiency increases with spin
        
        # Mass actually added to black hole
        dm_actual = dm * efficiency
        self.mass += dm_actual
        
    def update(self, dt: float):
        """Update black hole state"""
        # Calculate Hawking radiation
        P_hawking = self.hawking_radiation_power()
        dm_hawking = P_hawking * dt / c**2
        self.mass -= dm_hawking
        
        # Process accretion
        if self.accretion_disk_mass > 0:
            dm_accretion = min(self.accretion_disk_mass, 
                             self.eddington_limit() * dt / c**2)
            # Update angular momentum from accretion
            r_isco = 6 * G * self.mass / c**2
            dJ = dm_accretion * np.sqrt(G * self.mass * r_isco)
            self.angular_momentum += dJ * np.array([0, 0, 1])
            
            self.mass += dm_accretion
            self.accretion_disk_mass -= dm_accretion
        
        self.time += dt
        self.store_state()
    
    def store_state(self):
        """Store current state in history"""
        state = {
            'time': self.time,
            'mass': self.mass,
            'position': self.position.copy(),
            'angular_momentum': self.angular_momentum.copy(),
            'charge': self.charge,
            'temperature': self.hawking_temperature,
            'radius': self.schwarzschild_radius,
            'spin': self.spin_parameter,
            'accretion_rate': self.accretion_disk_mass
        }
        self.history.append(state)
        
    def eddington_limit(self) -> float:
        """Calculate Eddington luminosity limit"""
        sigma_T = 6.65e-29  # Thomson cross-section
        m_p = 1.67e-27  # Proton mass
        
        return 4 * np.pi * G * self.mass * m_p * c / sigma_T
    
    def frame_dragging_frequency(self, r: float) -> float:
        """Calculate frame dragging frequency"""
        a = self.spin_parameter
        return 2 * G * self.mass * a / (c * r**3)
    
    def gravitational_wave_power(self, orbital_frequency: float) -> float:
        """Calculate gravitational wave emission power"""
        # Quadrupole formula
        return 32 * G**4 * self.mass**5 * orbital_frequency**6 / (5 * c**5)
    
    def hawking_radiation_power(self) -> float:
        """Calculate power of Hawking radiation"""
        T = self.hawking_temperature()
        A = 4 * np.pi * self.schwarzschild_radius**2  # Horizon area
        sigma = 5.67e-8  # Stefan-Boltzmann constant
        return sigma * A * T**4
    
    def is_within_ergosphere(self, r: float, theta: float) -> bool:
        """Check if point is within ergosphere"""
        a = self.spin_parameter
        r_ergo = (G * self.mass / c**2) * (1 + np.sqrt(1 - (a * np.cos(theta))**2))
        return r < r_ergo
    
    def kerr_newman_metric(self, r: float, theta: float) -> Dict[str, float]:
        """Calculate Kerr-Newman metric components"""
        a = self.spin_parameter
        Q = self.charge
        
        # Metric parameters
        rho2 = r**2 + (a * np.cos(theta))**2
        rho2 = r**2 + (a * np.cos(theta))**2
        
        # Metric components
        g_tt = -(1 - 2*G*self.mass*r/(c**2 * rho2))
        g_pp = (r**2 + a**2 + 2*G*self.mass*a**2*r*np.sin(theta)**2/(c**2 * rho2)) * np.sin(theta)**2
        g_tp = -2*G*self.mass*a*r*np.sin(theta)**2/(c * rho2)
        
        return {'g_tt': g_tt, 'g_pp': g_pp, 'g_tp': g_tp}
    
    def gravitational_lensing(self, light_position: np.ndarray, 
                            light_direction: np.ndarray) -> Tuple[np.ndarray, float]:
        """Calculate gravitational lensing effect"""
        r = np.linalg.norm(light_position)
        impact_parameter = np.linalg.norm(np.cross(light_position, light_direction))
        
        # Einstein radius
        theta_E = np.sqrt(4 * G * self.mass * r / (c**2 * impact_parameter**2))
        
        # Deflection angle
        deflection = 4 * G * self.mass / (c**2 * impact_parameter)
        
        # New direction after lensing
        new_direction = light_direction + deflection * np.cross(
            np.cross(light_direction, light_position), light_direction
        )
        
        return new_direction / np.linalg.norm(new_direction), theta_E
    
    def accretion_disk_temperature(self, r: float) -> float:
        """Calculate temperature in accretion disk at radius r"""
        # Innermost stable circular orbit
        r_isco = 6 * G * self.mass / c**2
        
        if r < r_isco:
            return 0.0
        
        # Accretion efficiency
        efficiency = 0.1 * (1 - np.sqrt(r_isco/(3*r)))
        
        # Mass accretion rate (assuming Eddington limit)
        M_dot = self.eddington_limit() / c**2
        
        # Temperature profile (Shakura-Sunyaev model)
        T = (3 * G * self.mass * M_dot * (1 - np.sqrt(r_isco/r)) / 
             (8 * np.pi * 5.67e-8 * r**3))**(1/4)
        
        return T
    
    def jet_power(self) -> float:
        """Calculate relativistic jet power (Blandford-Znajek process)"""
        a = np.linalg.norm(self.angular_momentum) / (self.mass * c)
        B_field = 1e4  # Typical magnetic field strength in Tesla
        
        # Blandford-Znajek power
        P_BZ = 1e47 * (self.mass/SOLAR_MASS/1e9)**2 * B_field**2 * a**2  # erg/s
        
        return P_BZ
    
    def event_horizon_area(self) -> float:
        """Calculate event horizon area"""
        a = np.linalg.norm(self.angular_momentum) / (self.mass * c)
        r_plus = G * self.mass / c**2 * (1 + np.sqrt(1 - a**2))
        
        return 4 * np.pi * (r_plus**2 + a**2)
    
    def gravitational_time_dilation(self, r: float) -> float:
        """Calculate gravitational time dilation at distance r"""
        return np.sqrt(1 - self.schwarzschild_radius/r)
        
    def event_horizon_radius(self) -> float:
        """Get event horizon radius"""
        return self.schwarzschild_radius