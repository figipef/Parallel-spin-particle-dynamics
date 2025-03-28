import numpy as np
import matplotlib.pyplot as plt

def plot_parts(parts, title, name):
    fig, ax = plt.subplots()
    ax.scatter(parts[:, 0], parts[:, 1], color='blue', marker='o')
    circle = plt.Circle((0, 0), R, color='red', fill=False, linewidth=2)
    ax.add_patch(circle)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(title)
    plt.savefig(name, dpi=300)
    plt.show()


iter = 500
k = 0.01
N = 8
R = 1.
jump = R/3

r = R* np.sqrt(np.random.rand(N))
th = 2.*np.pi*np.random.rand(N)

parts = np.column_stack([r*np.cos(th), r*np.sin(th)])
E_old = np.inf

plot_parts(parts, "Initial Configuration", r'C:\Users\afons\Desktop\Tópicos FC\homeworks\Annealing\1.png')

for i in range(iter):
    if i == 100:
        print(i)
        plot_parts(parts, "Configuration at t = 100", r'C:\Users\afons\Desktop\Tópicos FC\homeworks\Annealing\2.png')
    if i == 200:
        print(i)
        plot_parts(parts, "Configuration at t = 500", r'C:\Users\afons\Desktop\Tópicos FC\homeworks\Annealing\3.png')
    E_new = 0.
    ind = np.random.randint(0, N-1)
    accept = False
    while not accept:
        r_jump = jump* np.sqrt(np.random.rand(1))[0]
        th_jump = 2.*np.pi*np.random.rand(1)[0]
        new_pos = parts[ind] + np.array([r_jump*np.cos(th_jump), r_jump*np.sin(th_jump)])
        accept = True if new_pos[0]**2 + new_pos[1]**2 <= R*R else False

    old_pos = parts[ind].copy()
    parts[ind] = new_pos
    for j in range(N):
        for k in range(j+1, N):
            dist = np.sqrt((parts[k,0] - parts[j,0])**2 + (parts[k,1] - parts[j,1])**2)
            E_new += 1./dist

    if E_new > E_old:
        chance = k*np.exp(-k*i)
        prob = np.random.rand(1)[0]
        if prob < chance:
            E_old = E_new
        else:
            parts[ind] = old_pos
    else:
        E_old = E_new

plot_parts(parts, "Configuration at t = 1500", r'C:\Users\afons\Desktop\Tópicos FC\homeworks\Annealing\4.png')
