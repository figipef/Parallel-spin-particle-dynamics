import numpy as np
import matplotlib.pyplot as plt

def plot_v_time(q, dt, Title = "", ylabel = "Values", Grid = True, Lims = 0): #takes and array q of some quantity and the time step between values and plots the data
    n = len(q)
    t = np.linspace(0., dt*n, n)

    plt.figure(figsize=(8,5))
    plt.plot(t, q)

    plt.xlabel('Time (units are overrated)', fontsize=16) #I will eventually uncover our dimessions
    plt.ylabel(ylabel, fontsize=16)
    plt.title(Title, fontsize=18)

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    if isinstance(Lims, (list, tuple)):
        plt.xlim((Lims[0][0], Lims[0][1]))
        plt.ylim((Lims[1][0], Lims[1][1]))

    plt.grid(Grid)
    plt.show()

def plot_lots_v_time(Q, dt, Title = "", ylabel = "Values", Grid = True, Lims = 0, labels = ["Quantity 1", "Quantity 2"]): #takes and array of arrays q of some quantities and the time step between values and plots the data
    N = len(Q)
    if N != len(labels):
        labels = ["Quantity " + str(i + 1) for i in range(N)]
    
    plt.figure(figsize=(8,5))
    
    for i in range(N):
        q = Q[i]
        n = len(q)
        t = np.linspace(0., dt*n, n)

        plt.plot(t, q, label = labels[i])

    plt.xlabel('Time (units are overrated)', fontsize=16) #I will eventually uncover our dimessions
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

#num primeiro aproach precisamos de plots tipo quantidade vs tempo assumindo que recebe a lista de outputs e o dicionario de input do reader
    #done

#Depois tambem haveremos de precisar de histogramas para espaço de fases

def hist_Phase_Space(q, bins=10, Title="", xlabel="Values", ylabel="Frequency", grid=True, density=False):

    plt.figure(figsize=(8, 5))
    plt.hist(q, bins=bins, density=density, edgecolor='black')
    plt.xlabel(xlabel, fontsize=16)
    plt.ylabel(ylabel, fontsize=16)
    plt.title(Title, fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    if grid:
        plt.grid(True)
    plt.show()