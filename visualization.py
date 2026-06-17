# visualization.py
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider, Button, CheckButtons
import numpy as np

class UniverseVisualizer:
    def __init__(self, data_collector):
        self.data = data_collector
        self.paused = False
        self.show_trails = True
        self.trail_length = 50
        self.particle_size = 100
        self.view_angle = 0
        
    def setup_controls(self, _):
        """Setup interactive controls"""
        # Sliders
        ax_time = plt.axes([0.1, 0.02, 0.3, 0.02])
        self.time_slider = Slider(ax_time, 'Time', 0, 1, valinit=0.0)
        
        ax_trail = plt.axes([0.1, 0.06, 0.3, 0.02])
        self.trail_slider = Slider(ax_trail, 'Trail Length', 0, 100, valinit=50.0)
        
        # Buttons
        ax_pause = plt.axes([0.5, 0.02, 0.1, 0.04])
        self.pause_button = Button(ax_pause, 'Pause/Resume')
        self.pause_button.on_clicked(self.toggle_pause)
        
        # Checkboxes
        ax_check = plt.axes([0.7, 0.02, 0.2, 0.1])
        self.check_buttons = CheckButtons(ax_check, ['Trails', 'Grid', 'Labels'], 
                                       [True, True, True])
        self.check_buttons.on_clicked(self.update_display_options)
        
    def toggle_pause(self, _):
        self.paused = not self.paused
        
    def update_display_options(self, label):
        if label == 'Trails':
            self.show_trails = not self.show_trails
            
    def plot_evolution_3d(self):
        """Enhanced 3D evolution animation"""
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection='3d')
        self.setup_controls(fig)
        
        def update(frame):
            if self.paused:
                return
                
            ax.clear()
            ax.grid(self.check_buttons.get_status()[1])
            
            # Get particle data
            positions = [h[frame]['position'] for h in self.data.particle_histories.values()]
            energies = [h[frame]['energy'] for h in self.data.particle_histories.values()]
            
            # Normalize energies for color mapping
            norm_energies = np.array(energies)
            if len(energies) > 0:
                norm_energies = (norm_energies - min(energies)) / (max(energies) - min(energies))
            
            # Plot particles
            scatter = ax.scatter(*zip(*positions), 
                               c=norm_energies,
                               cmap='plasma',
                               s=self.particle_size,
                               alpha=0.8)
            
            # Plot trails
            if self.show_trails:
                for i, hist in enumerate(self.data.particle_histories.values()):
                    trail = [h['position'] for h in hist[max(0, frame-self.trail_length):frame]]
                    if len(trail) > 1:
                        trail = np.array(trail)
                        ax.plot3D(trail[:,0], trail[:,1], trail[:,2], 
                                alpha=0.3, c=plt.cm.plasma(norm_energies[i]))
            
            # Labels
            if self.check_buttons.get_status()[2]:
                ax.set_xlabel('X (m)')
                ax.set_ylabel('Y (m)')
                ax.set_zlabel('Z (m)')
                ax.set_title(f'Time: {frame * self.data.time_step:.2e} s')
            
            # Rotate view
            ax.view_init(elev=30, azim=self.view_angle)
            self.view_angle += 0.5
            
            # Add colorbar
            plt.colorbar(scatter, label='Kinetic Energy')
            
        anim = FuncAnimation(fig, update, 
                           frames=len(next(iter(self.data.particle_histories.values()))),
                           interval=50,
                           blit=False)
        
        plt.show()
        return anim
        
    def plot_metrics_dashboard(self):
        """Enhanced metrics dashboard"""
        fig = plt.figure(figsize=(15, 10))
        gs = fig.add_gridspec(3, 3)
        
        # Energy plot
        ax_energy = fig.add_subplot(gs[0, :])
        times = [h['time'] for h in next(iter(self.data.particle_histories.values()))]
        energies = [sum(h['energy'] for h in frame) for frame in zip(*self.data.particle_histories.values())]
        ax_energy.plot(times, energies, 'b-', label='Total Energy')
        ax_energy.set_ylabel('Energy (J)')
        ax_energy.set_title('System Energy Evolution')
        ax_energy.grid(True)
        
        # Trajectory plot
        ax_traj = fig.add_subplot(gs[1:, :2], projection='3d')
        for hist in self.data.particle_histories.values():
            positions = np.array([h['position'] for h in hist])
            ax_traj.plot3D(positions[:,0], positions[:,1], positions[:,2])
        ax_traj.set_title('Particle Trajectories')
        
        # Statistics
        ax_stats = fig.add_subplot(gs[1:, 2])
        stats = self.calculate_statistics()
        ax_stats.axis('off')
        stats_text = '\n'.join([f"{k}: {v:.2e}" for k, v in stats.items()])
        ax_stats.text(0.1, 0.9, stats_text, transform=ax_stats.transAxes)
        
        plt.tight_layout()
        plt.show()
        
    def calculate_statistics(self):
        """Calculate system statistics"""
        latest_frame = -1
        positions = [h[latest_frame]['position'] for h in self.data.particle_histories.values()]
        velocities = [h['velocity'] for h in self.data.particle_histories.values()]
        
        return {
            'Average Speed': np.mean([np.linalg.norm(v) for v in velocities]),
            'Max Distance': np.max([np.linalg.norm(p) for p in positions]),
            'Total Energy': sum(h[latest_frame]['energy'] for h in self.data.particle_histories.values()),
            'Particle Count': len(self.data.particle_histories)
        }