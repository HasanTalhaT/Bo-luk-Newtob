# universe_analyzer.py
import numpy as np
from typing import Dict, Tuple
from scipy.fft import fftn
import warnings

class UniverseAnalyzer:
    def __init__(self, data_collector):
        self.data = data_collector
        self.metrics = {}
        self.grid_size = 64

    def _get_particle_positions(self) -> np.ndarray:
        """Extract particle positions from simulation data"""
        positions = []
        for particle_data in self.data.particle_histories.values():
            if particle_data:
                latest_data = particle_data[-1]
                positions.append(latest_data['position'])
        return np.array(positions)

    def _get_particle_velocities(self) -> np.ndarray:
        """Extract particle velocities from simulation data"""
        velocities = []
        for particle_data in self.data.particle_histories.values():
            if particle_data:
                latest_data = particle_data[-1]
                velocities.append(latest_data['velocity'])
        return np.array(velocities)

    def _compute_velocity_dispersion(self, velocities: np.ndarray) -> float:
        """Calculate velocity dispersion of particles"""
        if len(velocities) == 0:
            return 0.0

        # Suppress RuntimeWarnings for zero division or invalid operations
        with np.errstate(divide='ignore', invalid='ignore'):
            v_mag = np.linalg.norm(velocities, axis=1)
            dispersion = np.std(v_mag)
            if np.isnan(dispersion):
                return 0.0
            return dispersion

    def _compute_power_spectrum(self, positions: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Compute matter power spectrum with error handling"""
        if len(positions) == 0:
            return np.array([]), np.array([])

        with np.errstate(divide='ignore', invalid='ignore'):
            # Add small epsilon to avoid division by zero
            eps = 1e-10
            pos_min = np.min(positions, axis=0)
            pos_max = np.max(positions, axis=0)
            pos_range = pos_max - pos_min + eps

            # Safe normalization
            positions_normalized = (positions - pos_min) / pos_range
            positions_normalized = np.clip(positions_normalized, 0, 1)

            # Safe integer casting
            indices = np.floor(positions_normalized * (self.grid_size - 1)).astype(int)
            indices = np.clip(indices, 0, self.grid_size - 1)

        # Create density field
        grid = np.zeros((self.grid_size,) * 3)
        for idx in indices:
            grid[tuple(idx)] += 1

        # Compute FFT
        dk = fftn(grid)
        power = np.abs(dk) ** 2

        # Get k-modes
        k = np.fft.fftfreq(self.grid_size) * 2 * np.pi
        kx, ky, kz = np.meshgrid(k, k, k, indexing='ij')
        k_mag = np.sqrt(kx**2 + ky**2 + kz**2)

        return k_mag.flatten(), power.flatten()

    def analyze_structure_formation(self) -> Dict:
        """Analyze cosmic structure formation"""
        positions = self._get_particle_positions()
        velocities = self._get_particle_velocities()

        return {
            'power_spectrum': self._compute_power_spectrum(positions),
            'velocity_dispersion': self._compute_velocity_dispersion(velocities),
        }