# test_integration.py
import unittest
from main import run_universe_simulation
from simulation import UniverseSimulation
from particles import Particle
from galaxy import Galaxy
from blackhole import BlackHole

class TestIntegration(unittest.TestCase):
    def test_full_simulation(self):
        try:
            run_universe_simulation()
            self.assertTrue(True)
        except Exception as e:
            self.fail(f"Simulation failed with error: {str(e)}")