import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.animation import FuncAnimation, PillowWriter
import colorsys

def hue_bright(th, phi): #auxiliary function
    #hue based phi, brightness based on th, saturation is costant
    hue = (phi % (2 * np.pi)) / (2 * np.pi)
    lig = 1. - (th % np.pi) / np.pi
    s = 1.
    return np.array([colorsys.hls_to_rgb(h,l,s) for h, l in zip(hue.flatten(), lig.flatten())]).reshape(hue.shape + (3,))

def animate_spin_sphere(s, dt, Title, save_as=None):
    # s is expected to be a NumPy array of shape (N, 3)
    n = s.shape[0]
    t = np.linspace(0., dt * n, n)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Set up initial empty scatter
    scatter = ax.scatter([], [], [], c=[], cmap='copper_r', marker='o', s=5)
    arrow = ax.quiver(0, 0, 0, 0, 0, 0, color='black', arrow_length_ratio=0.1)


    norm = plt.Normalize(t.min(), t.max())

    # Sphere wireframe for context
    phi = np.linspace(0, 2 * np.pi, 30)
    th = np.linspace(0, np.pi, 30)
    X = np.outer(np.cos(phi), np.sin(th))
    Y = np.outer(np.sin(phi), np.sin(th))
    Z = np.outer(np.ones_like(phi), np.cos(th))
    ax.plot_surface(X, Y, Z, color='gray', alpha=0.2, edgecolor='none')

    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.set_zlim([-1, 1])
    ax.set_xlabel('sx Projection', fontsize=16)
    ax.set_ylabel('sy Projection', fontsize=16)
    ax.set_zlabel('sz Projection', fontsize=16)
    ax.set_title(Title, fontsize=18)
    ax.tick_params(axis='both', labelsize=14)

    cmap = plt.get_cmap('copper_r')
    # Update function
    jump = 10
    def update(frame):
        nonlocal arrow
        arrow.remove()
        ax.collections.clear()
        ax.plot_surface(X, Y, Z, color='gray', alpha=0.2, edgecolor='none')
        ax.scatter(s[:frame*jump+1, 0], s[:frame*jump+1, 1], s[:frame*jump+1, 2], 
                   c=t[:frame*jump+1], cmap='copper_r', norm=norm, s=5)
        vx, vy, vz = s[frame*jump]
        arrow = ax.quiver(0, 0, 0, vx, vy, vz, color=cmap(norm(t[frame*jump+1])), arrow_length_ratio=0.1)

    ani = animation.FuncAnimation(fig, update, frames=int(n/jump), interval=100, repeat=False)

    if save_as:
        writer = animation.PillowWriter(fps=500)  
        ani.save(save_as, writer=writer)
    else:
        plt.show()

def animate_3D_traj_spin_trail(x_arr, y_arr, z_arr, Sx_arr, Sy_arr, Sz_arr, dt, Title="", Lims=0, save_path=None, show_axis=1):
    N = len(x_arr)       # Number of trajectories
    n = len(x_arr[0])    # Number of time steps

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Containers for markers and colors
    scatters = []
    colors_list = []

    # Create scatter objects for the whole trajectory
    for i in range(N):
        Sx, Sy, Sz = np.array(Sx_arr[i]), np.array(Sy_arr[i]), np.array(Sz_arr[i])
        phi = np.arctan2(Sy, Sx)
        theta = np.arccos(Sz / np.sqrt(Sx**2 + Sy**2 + Sz**2 + 1e-10))
        colors = hue_bright(theta, phi)
        colors_list.append(colors)

        # Create scatter for the entire trajectory, updating color by each point
        scatter = ax.scatter([], [], [], s=2, c=colors[0], marker='o')  # Initial color
        scatters.append(scatter)
        
    if show_axis == 1:
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False
        ax.set_axis_off()

    def init():
        if show_axis == 1:
            ax.set_xlabel('x [$c/\omega_0$]')
            ax.set_ylabel('y [$c/\omega_0$]')
            ax.set_zlabel('z [$c/\omega_0$]')

        ax.set_title(Title)
        if isinstance(Lims, (list, tuple)) and len(Lims) == 3:
            ax.set_xlim(Lims[0])
            ax.set_ylim(Lims[1])
            ax.set_zlim(Lims[2])
        
        return scatters

    def update(frame):
        for i in range(N):
            xi = x_arr[i][:frame+1]
            yi = y_arr[i][:frame+1]
            zi = z_arr[i][:frame+1]
            colors = colors_list[i][:frame+1]

            # Update scatter plot: x, y, z, and color at the current frame
            scatters[i]._offsets3d = (xi, yi, zi)
            scatters[i].set_color(colors)  # Update the color of the trajectory up to the current point

        # Rescale the axes after each frame update
        ax.relim()
        ax.autoscale_view(True, True, True)

        return scatters

    # Create the animation
    anim = FuncAnimation(fig, update, frames=n, init_func=init, blit=False, interval=30)

    # Save or show the animation
    if save_path:
        anim.save(save_path, writer="pillow", dpi=200)
    else:
        plt.show()

    return anim