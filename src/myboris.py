import numpy as np
import matplotlib.pyplot as plt

E = np.array([0., 1., 0.])
B = np.array([0., 0., 1.])

p = [np.array([0., 0., 0.])]

x = [np.array([0., 0., 0.])]
dt = 0.01
n = 10000
q = 1.
m = 1.

for i in range(n):
    pi = p[i]
    pm = pi + q*dt*0.5*E
    
    gam = np.sqrt(1. + np.inner(pi,pi)/m/m)
    t = B*q*dt*0.5/gam/m
    s = 2*t/(1. + np.inner(t,t))

    ppr = pm + np.cross(pm, t)
    ppl = pm + np.cross(ppr, s)
    pf = ppl + q*dt*0.5*E
    
    p.append(pf)  
    xf = x[i] + dt*p[i + 1] / m /np.sqrt(1. + np.inner(pf,pf)/m/m)
    x.append(xf)  

x_array = np.array(x)

plt.figure(figsize=(6, 4.5))
plt.scatter(x_array[:, 0], x_array[:, 1], c=np.arange(len(x)), cmap='viridis', marker='o', s=1)

plt.xlabel("X Position ")
plt.ylabel("Y Position ")
plt.title("Particle Trajectory (X vs Y)")
plt.colorbar(label="Time Step")
plt.grid()

plt.show()