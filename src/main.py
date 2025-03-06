import numpy as np
import Reader as R
import Plotter as P

v = R.read_output("src/output_example.txt")
print(v)
print('-----------------------------')
print(v[:,1])
print('-----------------------------')
v = R.read_input("src/input_example.txt")
print(v)

