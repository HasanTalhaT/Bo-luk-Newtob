import json
from datetime import datetime
import numpy as np
import os
from dataclasses import dataclass
from typing import Optional, Dict

# Custom JSON encoder for numpy types
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, datetime):
            return obj.isoformat()
        return super().default(obj)

# Import simulation components
from simulation import UniverseSimulation, SimulationConfig
from universe_analyzer import UniverseAnalyzer
from particles import Particle
from blackhole import BlackHole
from galaxy import Galaxy, GalaxyProperties
from constants import HUBBLE_CONSTANT, DARK_MATTER_DENSITY, CRITICAL_DENSITY, DARK_ENERGY_DENSITY  # Added DARK_ENERGY_DENSITY

@dataclass
class BaseForces:
    base_field1: int
    base_field2: float  # Non-default field

    def __post_init__(self):
        if self.base_field2 is None:
            self.base_field2 = 1.0  # Set default value here if needed

@dataclass
class Forces(BaseForces):
    # Non-default fields
    memory_limit: int
    field1: float

    # Default fields
    field2: Optional[Dict] = None
    default_field1: int = 0

# Define simulation constants
DAYS_IN_YEAR = 365
SIMULATION_YEARS = 10
TOTAL_STEPS = DAYS_IN_YEAR * SIMULATION_YEARS

# Define physical constants
SOLAR_MASS = 1.989e30  # kg
MERCURY_MASS = 3.285e23  # kg
VENUS_MASS = 4.867e24  # kg
EARTH_MASS = 5.972e24  # kg
MOON_MASS = 7.342e22  # kg
MARS_MASS = 6.39e23  # kg
JUPITER_MASS = 1.898e27  # kg

# Visualization configuration
visualization_config = {
    "enable": True,
    "trail_length": 100,
    "update_interval": 10,
    "save_frames": True,
    "frame_dir": "simulation_frames",
    "colormap": "viridis",
    "particle_size": 100,
    "show_forces": True,
    "show_velocities": True
}

# Optionally define DARK_ENERGY_DENSITY if not in constants.py
# DARK_ENERGY_DENSITY = 0.7  # Fraction of critical density

def setup_simulation():
    """Initialize simulation configuration"""
    config = SimulationConfig(
        time_step=60 * 60 * 24,  # 1 day in seconds
        total_steps=TOTAL_STEPS,  # 10 years
        dark_matter_density=DARK_MATTER_DENSITY * CRITICAL_DENSITY,
        hubble_constant=HUBBLE_CONSTANT,
        dark_energy_density=DARK_ENERGY_DENSITY * CRITICAL_DENSITY  # Added argument
    )
    
    config.output_config = {
        "save_interval": 100,
        "output_dir": "simulation_results",
        "save_metrics": True,
        "save_trajectories": True,
        "compression": True
    }
    return config

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

def initialize_simulation():
    """Initialize the simulation and celestial bodies"""
    config = setup_simulation()
    sim = UniverseSimulation(config)
    particles = initialize_solar_system()
    
    # Initialize Galaxy with default properties
    galaxy_properties = GalaxyProperties.default_properties()
    milky_way_galaxy = Galaxy(properties=galaxy_properties, scale_radius=galaxy_properties.scale_radius, position=np.array([0, 0, 0]))
    
    # Initialize BlackHole (assuming BlackHole takes mass, position, angular_momentum)
    sagittarius_a_black_hole = BlackHole(
        mass=4e6 * SOLAR_MASS, 
        position=np.array([0, 0, 0]), 
        angular_momentum=np.array([0, 0, 1e40])
    )
    
    return sim, config, particles, milky_way_galaxy, sagittarius_a_black_hole

def run_simulation_steps(sim, config, particles, galaxy, black_hole):
    """Run the simulation steps"""
    from tqdm import tqdm

    for step in tqdm(range(config.total_steps), desc="Simulation Progress"):
        sim.update(step, particles, galaxy, black_hole)
        if step % config.output_config["save_interval"] == 0:
            sim.data_collector.save_metrics(step)
            sim.data_collector.save_trajectories(step)
    # If visualization or frame saving is enabled, handle it here
    if visualization_config["save_frames"]:
        if not os.path.exists(visualization_config["frame_dir"]):
            os.makedirs(visualization_config["frame_dir"])
        # Implement frame saving logic if needed

def analyze_simulation(sim):
    """Analyze the simulation results"""
    analyzer = UniverseAnalyzer(sim.data_collector)
    results = analyzer.analyze_structure_formation()
    return results

def run_universe_simulation():
    """Run enhanced universe simulation"""
    try:
        sim, config, particles, milky_way_galaxy, black_hole = initialize_simulation()
        
        print("\nRunning simulation...")
        run_simulation_steps(sim, config, particles, milky_way_galaxy, black_hole)
        
        print("\nAnalyzing simulation results...")
        results = analyze_simulation(sim)
        
        print("\nSimulation completed successfully!")
        return results
            
    except Exception as e:
        print(f"Simulation failed: {str(e)}")
        raise

if __name__ == "__main__":
    run_universe_simulation()
