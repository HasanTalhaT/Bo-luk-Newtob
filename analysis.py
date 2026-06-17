# analysis.py

import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, Tuple
from scipy.fft import fftn
from scipy.spatial import cKDTree
import h5py

class UniverseAnalyzer:
    def __init__(self, simulation_data):
        self.data = simulation_data
        self.metrics = {}
        self.grid_size = 64  # For power spectrum calculations
        
    def _get_particle_positions(self) -> np.ndarray:
        """Retrieve particle positions from the simulation data"""
        return self.data['positions']
    
    def _get_particle_velocities(self) -> np.ndarray:
        """Retrieve particle velocities from the simulation data"""
        return self.data['velocities']
    
    def _get_particle_masses(self) -> np.ndarray:
        """Retrieve particle masses from the simulation data"""
        return self.data['masses']
    
    def analyze_structure_formation(self) -> Dict:
        """
        Analyze cosmic structure formation.
        """
        positions = self._get_particle_positions()
        velocities = self._get_particle_velocities()
        masses = self._get_particle_masses()
        
        return {
            'correlation_function': self._compute_correlation(positions),
            'power_spectrum': self._compute_power_spectrum(positions),
            'velocity_dispersion': self._compute_velocity_dispersion(velocities),
            'mass_function': self._compute_mass_function(masses),
            'clustering': self._compute_clustering(positions)
        }

        def _compute_correlation(self, positions: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
            radius_bins = np.logspace(-1, 2, 50)
            tree = cKDTree(positions)
            counts = []
            volumes = []
            
            # Count pairs in spherical shells
            for r1, r2 in zip(radius_bins[:-1], radius_bins[1:]):
                pairs = tree.count_neighbors(tree, r2) - tree.count_neighbors(tree, r1)
                volume = 4/3 * np.pi * (r2**3 - r1**3)
                counts.append(pairs)
                volumes.append(volume)
                
            counts = np.array(counts)
            volumes = np.array(volumes)
            correlation = counts / volumes / (len(positions) * (len(positions) - 1) / 2)
            return radius_bins[:-1], correlation
    
    def _compute_power_spectrum(self, positions: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Compute matter power spectrum"""
        # Create density field
        grid = np.zeros((self.grid_size,) * 3)
        positions_normalized = (positions - positions.min(axis=0)) / (positions.max(axis=0) - positions.min(axis=0))
        indices = (positions_normalized * (self.grid_size - 1)).astype(int)
        np.add.at(grid, tuple(indices.T), 1)
        # Compute FFT and power
        dk = fftn(grid)
        power = np.abs(dk)**2

        k = np.fft.fftfreq(self.grid_size) * 2 * np.pi
        k_mag = np.sqrt(np.sum(np.meshgrid(k, k, k, indexing='ij')**2, axis=0)).flatten()
        return k_mag, power.flatten()
    
    def visualize_results(self):
        """Create visualization dashboard"""
        fig = plt.figure(figsize=(15, 10))
        gs = fig.add_gridspec(3, 2)
        ax1 = fig.add_subplot(gs[0, 0])
        # Particle distribution
        ax2 = fig.add_subplot(gs[0, 1], projection='3d')
        # Plotting code here

        return fig

    def save_results(self, filename: str):
        """Save analysis results to an HDF5 file"""
        with h5py.File(filename, 'w') as f:
            results_group = f.create_group('analysis')
            results = self.analyze_structure_formation()
            for name, data in results.items():
                results_group.create_dataset(name, data=np.array(data))