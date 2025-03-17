import numpy as np
import Reader as R
import Plotter as P

v = R.read_output("output/e_field.txt")
w = R.read_output("output/position.txt")
z = R.read_output("output/spin.txt")
p = R.read_output("output/momentum.txt")


print(v)
print('-----------------------------')
print(v[:,1])
print('-----------------------------')
d = R.read_input("build/input.txt")
print(d)

#P.plot_v_time(v[:,0], d['TIME_STEP'], Title="TESTING", ylabel="I honestly forgor", Grid=False, Lims=((0.,15.),(-2,3)))

P.plot_lots_v_time([w[:,0]-w[:,0][0], w[:,1], w[:,2]], d['TIME_STEP'], Title="Position vs time", ylabel="pos", Grid=False, Lims=((0.,d['TOTAL_TIME']),(np.min([w[:, 0], w[:, 1], w[:, 2]])*1.1,np.max([w[:, 0]-w[:,0][0], w[:, 1], w[:, 2]])*1.1)), labels=["x", "y", "z"])

P.plot_lots_v_time([p[:,0], p[:,1], p[:,2]], d['TIME_STEP'], Title="Momentum vs time", ylabel="Mom", Grid=True, Lims=((0.,d['TOTAL_TIME']),(np.min([p[:, 0], p[:, 1], p[:, 2]])*1.1,np.max([p[:, 0], p[:, 1], p[:, 2]])*1.1)), labels=["px", "py","pz"])

#P.plot_lots_v_time([z[:,0], z[:,1], z[:,2]], d['TIME_STEP'], Title="Spin vs time", ylabel="Spin", Grid=False, Lims=((0.,d['TOTAL_TIME']),(-1.5,1.5)), labels=["sx", "sy", "sz"])

#P.plot_lots_v_time([v[:,0], v[:,1], v[:,2]], d['TIME_STEP'], Title="TESTING", ylabel="I honestly forgor", Grid=False, Lims=((0.,15.),(-2,3)), labels=[d["PAR1"], "spin y", d["PAR3"]])

P.plot_lots_v_space(v, 0.2, Title="TESTING", ylabel="Field Amplitude", Grid=False, Lims=((0.,4.),(-2,3)), space_axis="x [$r_L$]")

#P.plot_hists([v[:,0], v[:,1]], bins=10, Title="Teste", xlabel="Values", ylabel="Frequency", grid=True, density=False, labels=None, alpha=0.5)

#P.plot_hists_seaborn([v[:,0], v[:,1]], bins=10, title="Teste Seaborn", xlabel="Values", ylabel="Frequency", grid=True, density=False, labels=None, alpha=0.5)

