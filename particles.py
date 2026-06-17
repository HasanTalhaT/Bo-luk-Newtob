# particles.py

import numpy as np
from typing import Optional, List
from dataclasses import dataclass
from constants import G, c, h_bar, ELECTRON_MASS, PROTON_MASS
from scipy.constants import e  # elementary charge

# Define constants
HYDROGEN_IONIZATION_ENERGY = 13.6 * e

@dataclass
class ParticleState:
    """Store particle state at a given time"""
    time: float
    position: np.ndarray
    velocity: np.ndarray
    energy: float
    momentum: np.ndarray

class Particle:
    _next_id = 0  # Class variable for generating unique IDs
    
    def __init__(self, 
                 mass: float, 
                 position: np.ndarray, 
                 velocity: np.ndarray,
                 charge: Optional[float] = None,
                 spin: float = 0,
                 radius: float = 1e-15,
                 name: str = ""):
        
        self.id = Particle._next_id  # Add unique ID
        Particle._next_id += 1
        
        self.mass = mass
        self.rest_mass = mass
        self.position = np.array(position, dtype=np.float64)
        self.velocity = np.array(velocity, dtype=np.float64)
        self.charge = charge
        self.spin = spin
        self.radius = radius
        self.force = np.zeros(3, dtype=np.float64)
        self.name = name
        
        # Tracking
        self.time = 0.0
        self.history: List[ParticleState] = []
        self.collision_count = 0
        self.max_history_length = 1000
        self._last_energy = None
        
    @property
    def gamma(self) -> float:
        """Relativistic gamma factor with safety checks"""
        v = np.linalg.norm(self.velocity)
        if v >= 0.99 * c:  # Limit velocity to 99% of c
            return 100.0  # Cap gamma factor
        return 1.0 / np.sqrt(1.0 - (v/c)**2)
    
    @property
    def relativistic_mass(self) -> float:
        """Relativistic mass with safety check"""
        return self.rest_mass * self.gamma if self.gamma != float('inf') else float('inf')
    
    @property
    def momentum(self) -> np.ndarray:
        """Relativistic momentum"""
        return self.relativistic_mass * self.velocity
    
    @property
    def kinetic_energy(self) -> float:
        """Relativistic kinetic energy"""
        return (self.gamma - 1) * self.rest_mass * c**2
    
    @property
    def total_energy(self) -> float:
        """Total relativistic energy"""
        return self.gamma * self.rest_mass * c**2
    
    @property
    def de_broglie_wavelength(self) -> float:
        momentum_norm = np.linalg.norm(self.momentum)
        if momentum_norm == 0:
            return np.inf
        return h_bar / momentum_norm
    
    def update(self, dt: float):
        """Update particle position and velocity"""
        self.position += self.velocity * dt
        self.velocity += (self.force / self.relativistic_mass) * dt
        self.time += dt
        self.store_state()
        
    def store_state(self):
        """Store current state in history"""
        state = ParticleState(
            time=self.time,
            position=self.position.copy(),
            velocity=self.velocity.copy(),
            energy=self.total_energy,
            momentum=self.momentum.copy()
        )
        self.history.append(state)
        if len(self.history) > self.max_history_length:
            self.history.pop(0)

    def interact_with_galaxy(self, galaxy):
        # Import here if needed to avoid circular import
        # from galaxy import Galaxy
        pass
        
    # ... [Other methods remain unchanged]