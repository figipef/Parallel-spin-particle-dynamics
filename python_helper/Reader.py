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
    with open(file_name, "r") as file:
        for line in file:
            content = line.strip().split()
            if "=" in line and len(content) > 1: 
                if "NUMBER" in line:
                    dic[content[0][:-1]] = int(content[1])
                else:
                    if "PAR" in line:
                        dic[content[0][:-1]] = content[1]
                    else:
                        if "E1=" in line or "B1=" in line or ',' in line:
                            dic[content[0][:-1]] = [float(val) for val in content[1].split(",")]
                        else:
                            dic[content[0][:-1]] = float(content[1])
    return dic