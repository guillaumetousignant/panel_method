# MCG 5136 Assignment 2
# Guillaume Tousignant, 0300151859
# February 3rd, 2020

import matplotlib.pyplot as plt
import numpy as np
import os
import re

x_arrays = []
y_arrays = []
cp_arrays = []
ue_arrays = []

N_finder = re.compile(r"I= \d*")

# Input from all the cp-X.dat files
filenames = [f for f in os.listdir(os.curdir) if os.path.isfile(f) and "cp-" in f and f.endswith(".dat")]
for filename in filenames:
    with open(filename, 'r') as file:
        lines = file.readlines()
        N_match = N_finder.search(lines[2])
        N = int(N_match.group(0)[3:])
        x_arrays.append(np.zeros(N))
        y_arrays.append(np.zeros(N))
        cp_arrays.append(np.zeros(N))
        ue_arrays.append(np.zeros(N))

        for i in range(0, N):
            numbers = lines[i+3].split()
            x_arrays[-1][i] = float(numbers[0])
            y_arrays[-1][i] = float(numbers[1])
            cp_arrays[-1][i] = float(numbers[2])
            ue_arrays[-1][i] = float(numbers[3])

        



#y = np.zeros(n_steps)

# Plotting
#plt.plot(t, y)
#plt.ylabel('y [m]')
#plt.xlabel('t [d]')
#plt.title(f"Height of fluid with time, delta_t = {delta_t}")
#plt.show()