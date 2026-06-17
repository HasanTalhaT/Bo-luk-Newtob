# simulation.py

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import time
from dataclasses import dataclass
from typing import List, Dict, Optional
import forces
from particles import Particle
from blackhole import BlackHole
from scipy import stats
import h5py
from scipy.spatial import ConvexHull, qhull
import scipy.spatial
import logging

# Constants
G = 6.67430e-11         # Gravitational constant (m^3 kg^-1 s^-2)
SOLAR_MASS = 1.989e30   # Solar mass in kilograms
k_B = 1.380649e-23      # Boltzmann constant (m^2 kg s^-2 K^-1)
c = 2.99792458e8        # Speed of light in vacuum (m/s)

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class Galaxy:
    """Represents a galaxy with a given mass and position."""
    def __init__(self, mass: float, position: np.ndarray):
        self.mass = mass
        self.position = np.array(position, dtype=np.float64)

@dataclass
class SimulationConfig:
    """Configuration for universe simulation."""
    time_step: float  # in seconds
    total_steps: int
    dark_matter_density: float  # in kg/m^3
    dark_energy_density: float  # in kg/m^3
    hubble_constant: float  # in 1/s
    enable_relativity: bool = True
    enable_quantum_effects: bool = False
    enable_dark_energy: bool = True
    visualization_config: Dict = None
    output_config: Dict = None

    def __post_init__(self):
        """Initialize default configurations if None."""
        if self.visualization_config is None:
            self.visualization_config = {
                "enable": True,
                "trail_length": 100,
                "update_interval": 10,
                "save_frames": False,
                "frame_dir": "simulation_frames",
                "colormap": "viridis",
                "particle_size": 100,
                "show_trails": True,
                "show_forces": False,
                "show_velocities": False
            }
        if self.output_config is None:
            self.output_config = {
                "save_interval": 100,
                "output_dir": "simulation_results",
                "save_metrics": True,
                "save_trajectories": True,
                "compression": True
            }

class DataCollector:
    def __init__(self):
        self.particle_histories = {}
        self.system_metrics = {
            'total_energy': [],
            'kinetic_energy': [],
            'potential_energy': [],
            'temperature': [],
            'density': [],
            'entropy': [],
            'particle_count': [],
            'time': [],
            'step': []
        }

    def collect_particle_data(self, particles: List[Particle], step: int, time: float):
        """Collect particle data at each timestep."""
        for particle in particles:
            if particle.name not in self.particle_histories:
                self.particle_histories[particle.name] = []

            self.particle_histories[particle.name].append({
                'step': step,
                'time': time,
                'position': particle.position.copy(),
                'velocity': particle.velocity.copy(),
                'mass': particle.mass,
                'energy': 0.5 * particle.mass * np.sum(particle.velocity**2)
            })

    def collect_system_metrics(self, particles: List[Particle], time: float, step: int):
        """Collect system-wide metrics."""
        # Calculate energies
        ke = sum(0.5 * p.mass * np.sum(p.velocity**2) for p in particles)
        pe = self.calculate_potential_energy(particles)

        # Update metrics
        self.system_metrics['total_energy'].append(ke + pe)
        self.system_metrics['kinetic_energy'].append(ke)
        self.system_metrics['potential_energy'].append(pe)
        self.system_metrics['particle_count'].append(len(particles))
        self.system_metrics['time'].append(time)
        self.system_metrics['step'].append(step)

        # Calculate other metrics
        self.system_metrics['temperature'].append(self.calculate_temperature(particles))
        self.system_metrics['density'].append(self.calculate_density(particles))
        self.system_metrics['entropy'].append(self.calculate_entropy(particles))

    def calculate_potential_energy(self, particles: List[Particle]) -> float:
        """Calculate total gravitational potential energy."""
        pe = 0.0
        for i, p1 in enumerate(particles):
            for p2 in particles[i+1:]:
                r_vec = p2.position - p1.position
                r = np.linalg.norm(r_vec)
                if r > 0:  # Avoid division by zero
                    pe -= G * p1.mass * p2.mass / r
        return pe

    def calculate_temperature(self, particles: List[Particle]) -> float:
        """Calculate system temperature from kinetic energies."""
        if not particles:
            return 0.0
        avg_ke = np.mean([0.5 * p.mass * np.sum(p.velocity**2) for p in particles])
        return 2 * avg_ke / (3 * k_B)

    def calculate_density(self, particles: List[Particle]) -> float:
        """Calculate the density of the system using ConvexHull."""
        positions = np.array([particle.position for particle in particles])
        
        # Remove duplicate positions
        positions = np.unique(positions, axis=0)
        
        # Ensure there are enough unique points and they are 3D
        if positions.shape[0] < 4:
            print("Not enough points to form a ConvexHull.")
            return 0.0
        if positions.shape[1] != 3:
            raise ValueError("Positions must be 3-dimensional.")
        
        # Check for coplanarity
        try:
            hull = ConvexHull(positions)
            return hull.volume
        except scipy.spatial.qhull.QhullError as e:
            print(f"ConvexHull error: {e}")
            # Handle coplanar points by estimating density differently
            return self.estimate_density_coplanar(positions)
    
    def estimate_density_coplanar(self, positions: np.ndarray) -> float:
        """Estimate density when points are coplanar."""
        # Example: Use 2D ConvexHull
        try:
            hull = ConvexHull(positions[:, :2])
            area = hull.area
            return area  # Adjust based on how density is defined
        except scipy.spatial.qhull.QhullError as e:
            print(f"2D ConvexHull error: {e}")
            return 0.0

    def calculate_entropy(self, particles: List[Particle]) -> float:
        """Calculate system entropy (placeholder implementation)."""
        # Placeholder entropy calculation
        return 0.0

    def calculate_statistics(self) -> dict:
        """Calculate statistical metrics for the simulation data."""
        stats_dict = {}
        for metric, values in self.system_metrics.items():
            if values:
                stats_dict[metric] = {
                    'mean': np.mean(values),
                    'std': np.std(values),
                    'min': np.min(values),
                    'max': np.max(values),
                    'median': np.median(values),
                    'skewness': stats.skew(values),
                    'kurtosis': stats.kurtosis(values)
                }
        return stats_dict

class UniverseSimulation:
    def __init__(self, config: SimulationConfig):
        self.config = config
        self.time = 0.0
        self.data_collector = DataCollector()
        self.last_update_time = time.time()
        self.frame_count = 0
        self.max_fps = 30

        plt.ion()
        self.fig = plt.figure(figsize=(12, 8))
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.setup_visualization()

    def setup_visualization(self):
        """Setup the plotting environment."""
        self.ax.set_xlabel('X (m)')
        self.ax.set_ylabel('Y (m)')
        self.ax.set_zlabel('Z (m)')
        limit = 2e12
        self.ax.set_xlim([-limit, limit])
        self.ax.set_ylim([-limit, limit])
        self.ax.set_zlim([-limit, limit])
        self.ax.grid(True)
        self.ax.set_facecolor('black')
        self.fig.patch.set_facecolor('black')
        plt.title('Universe Simulation')
        plt.tight_layout()

    def calculate_forces(self, particle, black_hole, galaxy):
        """Calculate forces acting on a particle."""
        r_gal = galaxy.position - particle.position
        # Calculate other forces...
        return forces

    def update(self, step: int, particles: List[Particle], galaxy: Galaxy, black_hole: Optional[BlackHole] = None):
        """Update the simulation state."""
        self.time = step * self.config.time_step

        # Update particles
        for particle in particles:
            forces = self.calculate_forces(particle, black_hole, galaxy)
            particle.force = forces
            particle.update(self.config.time_step)

        # Collect data
        self.data_collector.collect_particle_data(particles, step, self.time)
        self.data_collector.collect_system_metrics(particles, self.time, step)

        # Visualization update
        if self.config.visualization_config["enable"]:
            if step % self.config.visualization_config["update_interval"] == 0:
                current_time = time.time()
                if current_time - self.last_update_time >= 1.0 / self.max_fps:
                    self.visualize(particles, galaxy, black_hole)
                    self.last_update_time = current_time

    def visualize(self, particles: List[Particle], galaxy: Galaxy, black_hole: Optional[BlackHole] = None):
        """Visualize the simulation state."""
        if not plt.fignum_exists(self.fig.number):
            return  # Exit if window was closed

        self.ax.cla()
        self.plot_particles(particles)
        self.plot_black_hole(black_hole)
        self.plot_galaxy_center(galaxy)
        self.set_labels_and_title()

        plt.draw()
        plt.pause(0.001)

    def plot_particles(self, particles: List[Particle]):
        """Plot particles and their trails."""
        for particle in particles:
            self.ax.scatter(*particle.position, c='blue',
                            s=self.config.visualization_config.get("particle_size", 100))
            self.plot_particle_trail(particle)

    def plot_particle_trail(self, particle: Particle):
        """Plot the trail of a particle."""
        if hasattr(particle, 'history') and self.config.visualization_config.get("show_trails", True):
            history = self.get_particle_history(particle)
            if len(history) > 1:
                self.ax.plot(history[:, 0], history[:, 1], history[:, 2], 'b-', alpha=0.3)

    def get_particle_history(self, particle: Particle) -> np.ndarray:
        """Get the history of a particle as a numpy array."""
        trail_length = self.config.visualization_config.get("trail_length", 50)
        history = self.data_collector.particle_histories.get(particle.name, [])
        return np.array([h['position'] for h in history[-trail_length:]])

    def plot_black_hole(self, black_hole: Optional[BlackHole]):
        """Plot the black hole if it exists."""
        if black_hole is not None:
            self.ax.scatter(*black_hole.position, c='red', s=200, label='Black Hole')

    def plot_galaxy_center(self, galaxy: Galaxy):
        """Plot the center of the galaxy."""
        self.ax.scatter(*galaxy.position, c='yellow', s=150, label='Galaxy Center')

    def set_labels_and_title(self):
        """Set the labels and title of the plot."""
        self.ax.set_xlabel('X (m)')
        self.ax.set_ylabel('Y (m)')
        self.ax.set_zlabel('Z (m)')
        self.ax.set_title(f'Time: {self.time:.2e} s')
        self.ax.legend()

def setup_simulation():
    """Setup simulation configuration."""
    return SimulationConfig(
        time_step=1e4,
        total_steps=1000,
        dark_matter_density=0.23,
        dark_energy_density=0.73,
        hubble_constant=70.0,
        visualization_config={"enable": True, "update_interval": 10, "show_trails": True},
        output_config={"save_interval": 100}
    )

def initialize_solar_system():
    """Initialize particles representing a simplified solar system."""
    SUN_MASS = 1.989e30
    MERCURY_MASS = 3.3011e23
    VENUS_MASS = 4.8675e24
    EARTH_MASS = 5.972e24
    MARS_MASS = 6.4171e23
    JUPITER_MASS = 1.8982e27

    particles = [
        Particle(SUN_MASS, [0, 0, 0], [0, 0, 0], "Sun"),
        Particle(MERCURY_MASS, [5.79e10, 0, 0], [0, 4.79e4, 0], "Mercury"),
        Particle(VENUS_MASS, [1.08e11, 0, 0], [0, 3.5e4, 0], "Venus"),
        Particle(EARTH_MASS, [1.496e11, 0, 0], [0, 2.978e4, 0], "Earth"),
        Particle(MARS_MASS, [2.279e11, 0, 0], [0, 2.41e4, 0], "Mars"),
        Particle(JUPITER_MASS, [7.786e11, 0, 0], [0, 1.307e4, 0], "Jupiter")
    ]
    return particles

def run_universe_simulation():
    """Run enhanced universe simulation."""
    try:
        print("Initializing simulation...")
        config = setup_simulation()
        sim = UniverseSimulation(config)

        # Initialize celestial bodies
        particles = initialize_solar_system()
        galaxy = Galaxy(1e12 * SOLAR_MASS, [0, 0, 0])
        black_hole = BlackHole(4e6 * SOLAR_MASS, [0, 0, 0],
                               angular_momentum=np.array([0, 0, 1e40]))

        # Run simulation
        print("\nRunning simulation...")
        for step in range(config.total_steps):
            sim.update(step, particles, galaxy, black_hole)

        # Finalize
        print("\nSimulation completed successfully!")
        plt.ioff()
        plt.show()

    except Exception as e:
        print(f"Simulation failed: {str(e)}")
        raise

if __name__ == "__main__":
    run_universe_simulation()
