import numpy as np
import Reader as R
import Plotter as P

v = R.read_output("src/output_example.txt")
print(v)
print('-----------------------------')
print(v[:,1])
print('-----------------------------')
d = R.read_input("src/input_example.txt")
print(d)

#P.plot_v_time(v[:,0], d['TIME_STEP'], Title="TESTING", ylabel="I honestly forgor", Grid=False, Lims=((0.,15.),(-2,3)))

#P.plot_lots_v_time([v[:,0], v[:,1], v[:,2]], d['TIME_STEP'], Title="TESTING", ylabel="I honestly forgor", Grid=False, Lims=((0.,15.),(-2,3)), labels=[d["PAR1"], "spin y", d["PAR3"]])

P.hist_Phase_Space(v[:,0], bins=50, Title="Teste", xlabel="Values", ylabel="Frequency", grid=True, density=False)