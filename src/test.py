import numpy as np
import matplotlib.pyplot as plt

N = 10000
R = np.random.uniform(0., 1., 2*N)

phi = 2.*np.pi*R[:N]
costh = 2.*R[N:] - 1.
sinth = np.sqrt(1. - costh**2)

x = sinth*np.cos(phi)
y = sinth*np.sin(phi)
z = costh

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, z, s=2, color='blue')  # Small points to visualize distribution

# Set axis properties
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.set_title("Random Points on a Sphere")

plt.show()