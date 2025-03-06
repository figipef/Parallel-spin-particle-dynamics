import numpy as np

def read_output(file_name): #takes in a string with output file name and outputs numpy array of shape (#lines, 3)
    vect = []
    with open(file_name, "r") as file:
        for line in file:
            content = np.array(line.strip().split())
            vect += [float(val) for val in content],
    return np.array(vect)

def read_input(file_name): #takes in string with input file name and outputs dictionary with input parameters (will need to be updated consistently to include more parameters)
    dic = {}
    i = 1
    with open(file_name, "r") as file:
        for line in file:
            if ("INPUT" not in line) and ("DIAGNOSTICS" not in line):
                if "NUMBER" in line:
                    content = line.strip().split()
                    dic[content[0][:-1]] = int(content[1])
                else:
                    if "PAR" in line:
                        content = line.strip().split()
                        dic[content[0][:-1]] = content[1]
                    else:
                        content = line.strip().split()
                        dic[content[0][:-1]] = float(content[1])
            i += 1
    return dic