# MCG 5136 Assignment 2
# Guillaume Tousignant, 0300151859
# February 3rd, 2020

import matplotlib.pyplot as plt
import numpy as np
import os
import re

alphas = []
x_arrays = []
y_arrays = []
cp_arrays = []
ue_arrays = []

alpha_finder = re.compile(r"alpha= \d*")
N_finder = re.compile(r"I= \d*")

# Input from all the cp-X.dat files
filenames = [f for f in os.listdir(os.curdir) if os.path.isfile(f) and "cp-" in f and f.endswith(".dat")]
for filename in filenames:
    with open(filename, 'r') as file:
        lines = file.readlines()
        alpha_match = alpha_finder.search(lines[0])
        alphas.append(float(alpha_match.group(0)[7:]))
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
            cp_arrays[-1][i] = -float(numbers[2])
            ue_arrays[-1][i] = float(numbers[3])

# Input from cl vs alpha file
with open("clalpha.dat", 'r') as file:
    lines = file.readlines()
    N_runs_match = N_finder.search(lines[2])
    N_runs = int(N_runs_match.group(0)[3:])
    alpha_array = np.zeros(N_runs)
    cl_array = np.zeros(N_runs)
    cm_array = np.zeros(N_runs)

    for i in range(0, N_runs):
        numbers = lines[i+3].split()
        alpha_array[i] = float(numbers[0])
        cl_array[i] = float(numbers[1])
        cm_array[i] = float(numbers[2])

# Plotting cp
legend_list = []
cp_fig, cp_ax = plt.subplots(1, 1)
for i in range(0, len(filenames)):
    cp_ax.plot(x_arrays[i], cp_arrays[i])
    legend_list.append(f"$\\alpha$ = {alphas[i]}°")

cp_ax.grid()
cp_ax.set_xlim(0, 1)
cp_ax.set_ylabel('-Cp')
cp_ax.set_xlabel('x/c')
cp_ax.set_title("Cp along chord")
cp_ax.legend(legend_list, loc='upper right')

# Plotting cl vs alpha
cl_legend_list = []
cl_fig, cl_ax = plt.subplots(1, 1)
cl_ax.plot(alpha_array, cl_array)

cl_ax.grid()
cl_ax.set_ylabel('CL')
cl_ax.set_xlabel('$\\alpha$ [°]')
cl_ax.set_title("CL vs $\\alpha$")

plt.show()