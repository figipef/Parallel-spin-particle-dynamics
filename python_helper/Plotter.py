import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os
import re
from matplotlib.widgets import Slider
from python_helper import Reader as R


def plot_2D_traj(x, y, dt, Title = "", Grid = True, Lims = 0): 
    
    # Plots a 2D trajectory with colormap given by dt

    n = len(x)
    t = np.linspace(dt, dt*n, n)

    plt.figure(figsize=(8,5))
    sc = plt.scatter(x, y, c=t, cmap='copper_r', marker='o', s=1)
    cbar = plt.colorbar(sc)
    cbar.set_label('Time [$1/\omega_0$]', fontsize=16)
    cbar.ax.tick_params(labelsize=14)    

    plt.xlabel('x [$c/\omega_0$]', fontsize=16)
    plt.ylabel('y [$c/\omega_0$]', fontsize=16)
    plt.title(Title, fontsize=18)

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    if isinstance(Lims, (list, tuple)):
        plt.xlim((Lims[0][0], Lims[0][1]))
        plt.ylim((Lims[1][0], Lims[1][1]))
    
    cmap = plt.get_cmap("copper_r")
    norm = plt.Normalize(t.min(), t.max())
    for i in range(1, 10):
        ind = i*int(n/10)
        arrow_color = cmap(norm(t[ind]))
        plt.annotate("", xy=(x[ind], y[ind]), xytext=(x[ind - 1], y[ind - 1]), arrowprops=dict(arrowstyle="->", color=arrow_color, linewidth=2,mutation_scale=25))

    plt.grid(Grid)
    plt.show()

def plot_lots_2D_traj(x_arr, y_arr, dt, Title = "", Grid = True, Lims = 0):

    # Plots N trajectories based on x_arr and y_arr in a 2D Plot

    N = len(x_arr)
    n = len(x_arr[0])
    t = np.linspace(dt, dt*n, n)

    plt.figure(figsize=(8,5))
    plt.xlabel('x [$c/\omega_0$]', fontsize=16)
    plt.ylabel('y [$c/\omega_0$]', fontsize=16)
    plt.title(Title, fontsize=18)

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    if isinstance(Lims, (list, tuple)):
        plt.xlim((Lims[0][0], Lims[0][1]))
        plt.ylim((Lims[1][0], Lims[1][1]))

    plt.grid(Grid)
    norm = plt.Normalize(t.min(), t.max())
    cmaps = ['copper_r','Blues', 'Oranges', 'Greens', 'Purples', 'Reds', 'Greys', 'YlOrBr', 'OrRd']
    for j in range(N):
        plt.scatter(x_arr[j], y_arr[j], c=t, cmap=cmaps[j % len(cmaps)], marker='o', s=1)
        cmap = plt.get_cmap(cmaps[j])
        for i in range(1, 10):
            ind = i*int(n/10)
            arrow_color = cmap(norm(t[ind]))
            plt.annotate("", xy=(x_arr[j][ind], y_arr[ind]), xytext=(x_arr[ind - 1], y_arr[ind - 1]), arrowprops=dict(arrowstyle="->", color=arrow_color, linewidth=2,mutation_scale=25))


    plt.show()

def plot_3D_traj(x, y, z, dt, Title = "", Lims = 0):

    #Plots a 3D scatter plot of 3 lists x,y,z with a color map given by dt

    n = len(x)
    t = np.linspace(dt, dt*n, n)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    sc = ax.scatter(x, y, z, c=t, cmap='copper_r', marker='o', s=1)
    cbar = plt.colorbar(sc)
    cbar.set_label('t [$1/\omega_0$]')

    ax.set_xlabel('x [$c/\omega_0$]', fontsize=16)
    ax.set_ylabel('y [$c/\omega_0$]', fontsize=16)
    ax.set_zlabel('z [$c/\omega_0$]', fontsize=16)
    ax.set_title(Title, fontsize=18)

    ax.tick_params(axis='both', labelsize=14)

    if isinstance(Lims, (list, tuple)) and len(Lims) == 3:
        ax.set_xlim(Lims[0])
        ax.set_ylim(Lims[1])
        ax.set_zlim(Lims[2])
    else:
        if Lims != 0:
            print("Ignored limit indications due to incorrect data type being provided\n")

    plt.show()

def plot_lots_3D_traj(x_arr, y_arr, z_arr, dt, Title = "", Lims = 0):

    # Plots N trajectories in a 3D scatter plot with x_arr, y_arr, z_arr with colormap given by dt

    N = len(x_arr)  # Number of trajectories
    n = len(x_arr[0])  # Number of data points per trajectory
    t = np.linspace(dt, dt * n, n)  # Time array

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    cmaps = ['copper_r', 'Blues', 'Oranges', 'Greens', 'Purples', 'Reds', 'Greys', 'YlOrBr', 'OrRd']
    
    for i in range(N):
        ax.scatter(x_arr[i], y_arr[i], z_arr[i], c=t, cmap=cmaps[i % len(cmaps)], marker='o', s=1)
      
    ax.set_xlabel('x [$c/\omega_0$]', fontsize=16)
    ax.set_ylabel('y [$c/\omega_0$]', fontsize=16)
    ax.set_zlabel('z [$c/\omega_0$]', fontsize=16)
    ax.set_title(Title, fontsize=18)

    ax.tick_params(axis='both', labelsize=14)

    if isinstance(Lims, (list, tuple)) and len(Lims) == 3:
        ax.set_xlim(Lims[0])
        ax.set_ylim(Lims[1])
        ax.set_zlim(Lims[2])
    else:
        if Lims != 0:
            print("Ignored limit indications due to incorrect data type being provided\n")

    plt.show()

def plot_spin_sphere(s, dt, Title):

    # Plots spin lists [(Sx1,Sy1,Sz1),....] on a sphere

    n = len(s.transpose()[0])
    t = np.linspace(0., dt*n, n)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    sc = ax.scatter(s[:,0], s[:,1], s[:,2], c=t, cmap='copper_r', marker='o', s=1)
    cbar = plt.colorbar(sc)
    cbar.set_label('Time [$1/\omega_0$]', fontsize=16)
    cbar.ax.tick_params(labelsize=14)    

    ax.set_xlabel('sx Projection', fontsize=16)
    ax.set_ylabel('sy Projection', fontsize=16)
    ax.set_zlabel('sz Projection', fontsize=16)
    ax.set_title(Title, fontsize=18)

    phi = np.linspace(0, 2 * np.pi, 30)
    th = np.linspace(0, np.pi, 30)
    X = np.outer(np.cos(phi), np.sin(th))
    Y = np.outer(np.sin(phi), np.sin(th))
    Z = np.outer(np.ones_like(phi), np.cos(th))
    ax.plot_surface(X, Y, Z, color='gray', alpha=0.2, edgecolor='none')

    ax.tick_params(axis='both', labelsize=14)
    plt.show()

def plot_lots_v_time(Q, dt, Title = "", ylabel = "Values", Grid = True, Lims = 0, labels = ["Quantity 1", "Quantity 2"]):
    
    # Plots N variables (can be positions, momenta, spin)
    # on a 2D plot of the type (quantity, time)

    N = len(Q)
    if N != len(labels):
        labels = ["Quantity " + str(i + 1) for i in range(N)]
    
    plt.figure(figsize=(8,5))
    
    for i in range(N):
        q = Q[i]
        n = len(q)
        t = np.linspace(0., dt*n, n)

        plt.plot(t, q, label = labels[i])

    plt.xlabel('Time ($\omega t$)', fontsize=16) #I will eventually uncover our dimessions
    plt.ylabel(ylabel, fontsize=16)
    plt.title(Title, fontsize=18)
    plt.legend(fontsize=15)

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    if isinstance(Lims, (list, tuple)):
        plt.xlim((Lims[0][0], Lims[0][1]))
        plt.ylim((Lims[1][0], Lims[1][1]))

    plt.grid(Grid)
    plt.show()

def plot_v_time(q, dt, Title = "", ylabel = "Values", Grid = True, Lims = 0):
    
    # Plots a given quantity with respect to time on a 2D plot

    n = len(q)
    t = np.linspace(0., dt*n, n)

    plt.figure(figsize=(8,5))
    plt.plot(t, q)

    plt.xlabel('Time ($\omega t$)', fontsize=16) #I will eventually uncover our dimessions
    plt.ylabel(ylabel, fontsize=16)
    plt.title(Title, fontsize=18)

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    if isinstance(Lims, (list, tuple)):
        plt.xlim((Lims[0][0], Lims[0][1]))
        plt.ylim((Lims[1][0], Lims[1][1]))

    plt.grid(Grid)
    plt.show()

def plot_lots_v_space(Q, dx, Title = "", ylabel = "Values", Grid = True, Lims = 0, labels = ["Quantity 1", "Quantity 2"], space_axis="Space [$r_L$]"):
    
    # Plots N qantities in a 2D plot with respect to each other

    N = len(Q)
    if N != len(labels):
        labels = ["Quantity " + str(i + 1) for i in range(N)]
    
    plt.figure(figsize=(8,5))
    
    for i in range(N):
        q = Q[i]
        n = len(q)
        x = np.linspace(0., dx*n, n)

        plt.plot(x, q, label = labels[i])

    plt.xlabel(space_axis, fontsize=16) #I will eventually uncover our dimessions
    plt.ylabel(ylabel, fontsize=16)
    plt.title(Title, fontsize=18)
    plt.legend(fontsize=15)

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    if isinstance(Lims, (list, tuple)):
        plt.xlim((Lims[0][0], Lims[0][1]))
        plt.ylim((Lims[1][0], Lims[1][1]))

    plt.grid(Grid)
    plt.show()

## histogram plotters ## 

def plot_hists_through_time(variables, time_step):

    n_vars = len(variables)
    if n_vars not in [1, 2]:
        raise ValueError("Variables list must have 1 or 2 elements")
    
    parent_path = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
    folder = os.path.join(parent_path, "output")

    # Build regex pattern
    if n_vars ==1:
        if variables[0] in ["p1","p2","p3"]:
            pattern = re.compile(r"^position\.txt$")
        elif variables[0] in ["m1","m2","m3"]:
            pattern = re.compile(r"^momentum\.txt$")
        elif variables[0] in ["s1","s2","s3"]:
            pattern = re.compile(r"^spin\.txt$")
        
        else:
            raise ValueError("For 1D histograms, variable must be one of p1, p2, p3, m1, m2, m3, s1, s2, or s3.")
    elif n_vars == 2 and variables[0] in ["p1","p2","m1","m2"] and variables[1] in ["p1","p2","m1","m2"]:
        pattern = re.compile(rf"histogram(\d+)\.txt$")
    else:
        raise ValueError("Variables list must have 1 or 2 elements and meet the criteria.")


    # list all files in the folder.
    try:
        files = os.listdir(folder)
    except FileNotFoundError:
        print(f"The folder '{folder}' does not exist.")
        exit()

    matching_files = []
    
    # files that match the pattern
    for filename in files:
        match = pattern.match(filename)
        if match:
            if n_vars == 1:
                # For 1D cases we don't extract time from the filename.
                time_val = 0
            else:
                time_val = int(match.group(1))*time_step
            matching_files.append((time_val, filename))
    
    if not matching_files:
        print("No matching files found.")
        return

    # sort files by time
    matching_files.sort(key=lambda x: x[0])
    time_values = [time for time, _ in matching_files]
    file_paths = [os.path.join(folder, filename) for _, filename in matching_files]

    # Create one figure and a slider
    fig, ax = plt.subplots(figsize=(8, 6))
    plt.subplots_adjust(bottom=0.2)  # make space for slider

    # Plot the first histogram as default
    initial_index = 0
    cbar = plot_hists_txt(ax, file_paths[initial_index], variables, time_values[initial_index])
    
    # Create a slider axis below the main plot.
    slider_ax = plt.axes([0.15, 0.1, 0.7, 0.03])
    slider = Slider(slider_ax, 'Index', 0, len(file_paths)-1, valinit=initial_index, valstep=1, 
                    valfmt='%0.0f')

    current_cbar = [cbar]
    
    def update(val):
        index = int(slider.val)

        if current_cbar[0] is not None:
            current_cbar[0].remove()
            current_cbar[0] = None
        new_cbar = plot_hists_txt(ax, file_paths[index], variables, time_values[index])
        current_cbar[0] = new_cbar
        fig.canvas.draw_idle()
    
    slider.on_changed(update)
    plt.show()


def load_histogram_data(file_path):
    """Read histogram data from file and return as a list of lists."""
    data = []
    with open(file_path, "r") as f:
        first_line = f.readline().strip()  # Read the first line
        
        # Check if the first line contains "time", and skip it if true
        if "time" in first_line.lower():
            print("Skipping header:", first_line)
        
        # Process the remaining lines
        for line in f:
            line = line.strip()
            if line:  # Ignore empty lines
                row = [float(val) for val in line.split()]
                data.append(row)
    return data

def plot_hists_txt(ax, file_name, variables, time_val):
    
    full_path = os.path.abspath(file_name)
    data = load_histogram_data(full_path)
    
    if not data:
        raise ValueError("No data found in the file.")

    ax.clear()
    cbar = None
    
    # 1D or 2D hist
    if len(variables) == 1:
        
        if len(data) == 1:

            histogram = np.array(data[0])

        elif all(len(row) == 3 for row in data):

            if any("1" in s for s in variables):

                histogram = np.array([row[0] for row in data])

            elif any("2" in s for s in variables):

                histogram = np.array([row[1] for row in data])
            
            elif any("3" in s for s in variables):

                histogram = np.array([row[2] for row in data])
        else:

            raise ValueError("Data in file does not represent a 1D histogram as expected.")
        
        ax.bar(np.arange(len(histogram)), histogram, align='center')
        ax.set_xlabel("Count")
        ax.set_ylabel(variables[0])
        ax.set_title("1D Histogram")
        ax.legend()
    
    elif len(variables) == 2:

        # Ensure all rows have the same number of columns.
        row_lengths = [len(row) for row in data]

        if len(set(row_lengths)) != 1:

            raise ValueError("Rows in the file have an inconsistent number of columns.")
        
        # If the data appears to be 1D but two variable names were provided, raise an error.
        if all(len(row) == 1 for row in data):

            raise ValueError("Data in file represents a 1D histogram, but two variables were provided.")
        
        histogram = np.array(data)
        im = ax.imshow(histogram, origin='lower', interpolation='nearest', aspect='auto')
        ax.set_xlabel(variables[0])
        ax.set_ylabel(variables[1])
        ax.set_title("2D Histogram Heatmap")
        # Add time annotation on the plot
        ax.text(0.95, 0.95, f"Time: {time_val}", transform=ax.transAxes,
                fontsize=12, verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.5))
        cbar = plt.colorbar(im, ax=ax, label="Count")
    
    else:
        raise ValueError("The 'variables' parameter must be a list with either 1 or 2 elements.")
    
    return cbar

def plot_field_hists(variables):
    """
    Plot 1D histograms for field data for variables in E1, E2, E3, B1, B2, B3.
    The data file is assumed to have rows corresponding to time steps and columns as positions.
    A slider is provided to scroll through time steps.
    
    Parameters:
        variable (str): One of "E1", "E2", "E3", "B1", "B2", "B3".
        folder (str): Folder path containing the field data text files.
        bins (int): Number of bins for the histogram (default is 50).
    """
    # Map the variable to its corresponding filename.
    if variables[0] in ["E1", "E2", "E3"]:
        filename = f"e_field{variables[0][-1]}.txt"
    elif variables[0] in ["B1", "B2", "B3"]:
        filename = f"b_field{variables[0][-1]}.txt"
    else:
        raise ValueError("Variable must be one of E1, E2, E3, B1, B2, B3.")

    parent_path = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
    folder = os.path.join(parent_path, "output")
    
    file_path = os.path.join(folder, filename)
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File {file_path} not found.")
    
    bins = 50
    dt = R.read_input("build/input.txt")["TIME_STEP"]
    
    # Load the data.
    # Data is assumed to be in a text file where each row is a time step and each column is a position value.
    data = np.loadtxt(file_path)
    n_times = data.shape[0]
    
    # Create the main figure and axis.
    fig, ax = plt.subplots(figsize=(8, 6))
    plt.subplots_adjust(bottom=0.25)  # Adjust bottom space for the slider.
    
    # Plot the initial histogram for the first time step.
    initial_index = 0
    row_data = data[initial_index, :]
    current_time = initial_index * dt
    ax.hist(row_data, bins=bins, color='blue', alpha=0.7)
    ax.set_title(f"{variables[0]} Histogram at Time {current_time:.2f}")
    ax.set_xlabel("Field Value")
    ax.set_ylabel("Frequency")
    
    # Create a slider to select the time step.
    slider_ax = plt.axes([0.15, 0.1, 0.7, 0.03])
    time_slider = Slider(slider_ax, 'Time', 0, n_times-1, valinit=initial_index,
                         valstep=dt, valfmt='%0.0f')
    
    def update(val):
        time_index = int(time_slider.val)
        current_time = time_index * dt
        ax.clear()  # Clear previous histogram.
        row_data = data[time_index, :]
        ax.hist(row_data, bins=bins, color='blue', alpha=0.7)
        ax.set_title(f"{variables[0]} Histogram at Time {current_time:.2f}")
        ax.set_xlabel("Field Value")
        ax.set_ylabel("Frequency")
        fig.canvas.draw_idle()
    
    time_slider.on_changed(update)
    plt.show()

