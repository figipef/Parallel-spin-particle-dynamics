import numpy as np
import matplotlib.pyplot as plt


iter = 50000
factor = .999
k = 1.
N = 10
R = 1.
jump = 0.05

r = R* np.sqrt(np.random.rand(N))
th = 2.*np.pi*np.random.rand(N)

parts = np.column_stack([r*np.cos(th), r*np.sin(th)])
E_old = 0.
for j in range(N):
    for l in range(j+1, N):
        dist = np.linalg.norm(parts[l] - parts[j])
        if dist > 0:
            E_old += 1./dist

fig, ax = plt.subplots()

ax.scatter(parts[:,0], parts[:,1], color='blue', marker='o')
circle = plt.Circle((0,0), R, color='red', fill=False, linewidth=2)
ax.add_patch(circle)

# Labels and title
plt.xlabel("X")
plt.ylabel("Y")
plt.title("Initial Particle positions")

# Show plot
plt.show()

for i in range(iter):
    ind = np.random.randint(0, N)
    accept = False
    while not accept:
        r_jump = jump #* np.sqrt(np.random.rand(1))[0]
        th_jump = 2.*np.pi*np.random.rand(1)[0]
        new_pos = parts[ind] + np.array([r_jump*np.cos(th_jump), r_jump*np.sin(th_jump)])
        accept = True if new_pos[0]**2 + new_pos[1]**2 <= R*R else False
    #print(parts[ind])
    #print(new_pos)
    #print('................')
    old_pos = parts[ind].copy()
    parts[ind] = new_pos
    #print(old_pos)
    #print(new_pos)
    #print("--------------")
    E_new = 0.
    for j in range(N):
        for l in range(j+1, N):
            dist = np.linalg.norm(parts[l] - parts[j])
            if dist > 0:
                E_new += 1./dist

    if E_new > E_old:
        chance = np.exp(-(E_new - E_old)/k)
        prob = np.random.rand(1)[0]
        if prob < chance:
            E_old = E_new
        else:
            parts[ind] = old_pos
    else:
        E_old = E_new
    k *= factor 

fig, ax = plt.subplots()

ax.scatter(parts[:,0], parts[:,1], color='blue', marker='o')
circle = plt.Circle((0,0), R, color='red', fill=False, linewidth=2)
ax.add_patch(circle)

# Labels and title
plt.xlabel("X")
plt.ylabel("Y")
plt.title("Final particle positions")

# Show plot
plt.show()

