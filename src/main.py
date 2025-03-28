import numpy as np
import Reader as R
import Plotter as P
import matplotlib.pyplot as plt

#v = R.read_output("../output/e_field.txt")
#w = R.read_output("../output/position.txt")
#z = R.read_output("../output/spin.txt")
p = R.read_output("../output/momentum.txt")
#s = R.read_output("../output/lots_spin.txt")


#print(v)
#print('-----------------------------')
#print(v[:,1])
print('-----------------------------')
d = R.read_input("../build/input.txt")
#print(d)

#P.plot_v_time(v[:,0], d['TIME_STEP'], Title="TESTING", ylabel="I honestly forgor", Grid=False, Lims=((0.,15.),(-2,3)))

#P.plot_lots_v_time([w[:,0]-w[:,0][0], w[:,1], w[:,2]], d['TIME_STEP'], Title="Position vs time", ylabel="Change in Position ($\Delta x$)", Grid=True, Lims=((0.,d['TOTAL_TIME']),(np.min([w[:, 0], w[:, 1], w[:, 2]])*1.1,np.max([w[:, 0]-w[:,0][0], w[:, 1], w[:, 2]])*1.1)), labels=["x", "y", "z"])

P.plot_lots_v_time([p[:,0], p[:,1], p[:,2]], d['TIME_STEP'], Title="Momentum vs time", ylabel="Momentum", Grid=True, Lims=((0.,d['TOTAL_TIME']),(np.min([p[:, 0], p[:, 1], p[:, 2]])*1.1,np.max([p[:, 0], p[:, 1], p[:, 2]])*1.1)), labels=["px", "py","pz"])

#P.plot_lots_v_time([z[:,0], z[:,1], z[:,2]], d['TIME_STEP'], Title="Spin vs time", ylabel="Spin", Grid=False, Lims=((0.,d['TOTAL_TIME']),(-1.5,1.5)), labels=["sx", "sy", "sz"])

#P.plot_lots_v_time([v[:,0], v[:,1], v[:,2]], d['TIME_STEP'], Title="TESTING", ylabel="I honestly forgor", Grid=False, Lims=((0.,15.),(-2,3)), labels=[d["PAR1"], "spin y", d["PAR3"]])

#P.plot_lots_v_space(v, 0.2, Title="TESTING", ylabel="Field Amplitude", Grid=False, Lims=((0.,4.),(-2,3)), space_axis="x [$r_L$]")

#P.plot_hists([v[:,0], v[:,1]], bins=10, Title="Teste", xlabel="Values", ylabel="Frequency", grid=True, density=False, labels=None, alpha=0.5)

#P.plot_hists_seaborn([v[:,0], v[:,1]], bins=10, title="Teste Seaborn", xlabel="Values", ylabel="Frequency", grid=True, density=False, labels=None, alpha=0.5)

#PLOT ANDRE PARA momentum em funcao de kx # bruh eu criei o lots_v_space no plotter for this reason

#plt.plot(w[:,0], w[:,1], color="Blue", label="$p_y$")
#plt.plot(w[:,0], p[:,0], color="Green", ls="--", label="$p_x$")

#plt.xlabel('Distance ($k x$)', fontsize=16) #I will eventually uncover our dimessions
#plt.ylabel("Momentum", fontsize=16)
#plt.legend(fontsize=15)
#plt.grid(True)
#plt.xlim(0,5)

#plt.show()

#sx = s[:,0]
#sy = s[:,1]
#sz = s[:,2]
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')
#ax.scatter(sx,sy,sz, s=1, c='blue', marker='o')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Spin')

u = np.linspace(0, 2 * np.pi, 30)  # Azimuth angles
v = np.linspace(0, np.pi, 30)  # Elevation angles
X = np.outer(np.cos(u), np.sin(v))  # X-coordinates
Y = np.outer(np.sin(u), np.sin(v))  # Y-coordinates
Z = np.outer(np.ones_like(u), np.cos(v))  # Z-coordinates

# Plot transparent sphere
ax.plot_surface(X, Y, Z, color='gray', alpha=0.2, edgecolor='none')

# Show plot
#plt.show()

#d = R.read_input("build/input.txt")
#print(d)

def histogram_caller(input_file_path): #input file path

    VAR1 = np.array([d["PAR1"],d["PAR2"]])
    time_step_1 = d["TIME_STEP"]

    P.plot_hists_through_time(VAR1, time_step_1)

    VAR2 = np.array(['p1'])
    #P.plot_hists_through_time(VAR2,time_step_1)

    VAR3 = np.array(["E2"])
    #P.plot_field_hists(VAR3)

histogram_caller(d)