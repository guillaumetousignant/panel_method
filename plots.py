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
alphas_psi = []
x_arrays_psi = []
y_arrays_psi = []
psi_arrays = []
alphas_uxuy = []
x_arrays_uxuy = []
y_arrays_uxuy = []
ux_arrays = []
uy_arrays = []

alpha_finder = re.compile(r"alpha= \d*")
I_finder = re.compile(r"I= \d*")
J_finder = re.compile(r"J= \d*")

# Input from all the cp-X.dat files
filenames = [f for f in os.listdir(os.curdir) if os.path.isfile(f) and "cp-" in f and f.endswith(".dat")]
for filename in filenames:
    with open(filename, 'r') as file:
        lines = file.readlines()
        alpha_match = alpha_finder.search(lines[0])
        alphas.append(float(alpha_match.group(0)[7:]))
        N_match = I_finder.search(lines[2])
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

# Input from all the psi-X.dat files
filenames_psi = [f for f in os.listdir(os.curdir) if os.path.isfile(f) and "psi-" in f and f.endswith(".dat")]
for filename in filenames_psi:
    with open(filename, 'r') as file:
        lines = file.readlines()
        alpha_match = alpha_finder.search(lines[0])
        alphas_psi.append(float(alpha_match.group(0)[7:]))
        I_match = I_finder.search(lines[2])
        I = int(I_match.group(0)[3:])
        J_match = J_finder.search(lines[2])
        J = int(J_match.group(0)[3:])
        x_arrays_psi.append(np.zeros(I))
        y_arrays_psi.append(np.zeros(J))
        psi_arrays.append(np.zeros((J, I)))

        for i in range(0, I):
            for j in range(0, J):
                numbers = lines[i*J + j + 3].split()
                if j == 0:
                    x_arrays_psi[-1][i] = float(numbers[0])
                if i == 0:
                    y_arrays_psi[-1][j] = float(numbers[1])
                psi_arrays[-1][j, i] = float(numbers[2])

# Input from all the uxuy-X.dat files
filenames_uxuy = [f for f in os.listdir(os.curdir) if os.path.isfile(f) and "uxuy-" in f and f.endswith(".dat")]
for filename in filenames_uxuy:
    with open(filename, 'r') as file:
        lines = file.readlines()
        alpha_match = alpha_finder.search(lines[0])
        alphas_uxuy.append(float(alpha_match.group(0)[7:]))
        N_match = I_finder.search(lines[2])
        N = int(N_match.group(0)[3:])
        x_arrays_uxuy.append(np.zeros(N))
        y_arrays_uxuy.append(np.zeros(N))
        ux_arrays.append(np.zeros(N))
        uy_arrays.append(np.zeros(N))

        for i in range(0, N):
            numbers = lines[i+3].split()
            x_arrays_uxuy[-1][i] = float(numbers[0])
            y_arrays_uxuy[-1][i] = float(numbers[1])
            ux_arrays[-1][i] = float(numbers[2])
            uy_arrays[-1][i] = float(numbers[3])

# Input from cl vs alpha file
with open("clalpha.dat", 'r') as file:
    lines = file.readlines()
    N_runs_match = I_finder.search(lines[2])
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
cp_ax.set_ylabel('-$C_p$')
cp_ax.set_xlabel('x/c')
cp_ax.set_title("$C_p$ along chord")
cp_ax.legend(legend_list, loc='upper right')

# Plotting cl vs alpha
cl_fig, cl_ax = plt.subplots(1, 1)
cl_ax.plot(alpha_array, cl_array)

cl_ax.grid()
cl_ax.set_ylabel('$C_L$')
cl_ax.set_xlabel('$\\alpha$ [°]')
cl_ax.set_title("$C_L$ vs $\\alpha$")

# Plotting streamlines
for i in range(0, len(filenames_psi)):
    psi_fig, psi_ax = plt.subplots(1, 1)
    cp = psi_ax.contourf(x_arrays_psi[i], y_arrays_psi[i], psi_arrays[i])
    psi_ax.plot(x_arrays[i], y_arrays[i], color='red') # Assuming cp and psi files are the same airfoil

    psi_fig.colorbar(cp) # Add a colorbar to a plot
    psi_ax.set_xlabel('x/c')
    psi_ax.set_ylabel('y/c')
    psi_ax.set_title(f"Streamlines at $\\alpha$ = {alphas[i]}°")
    psi_ax.axis('scaled')

# Plotting u_x
legend_list_ux = []
ux_fig, ux_ax = plt.subplots(1, 1)
for i in range(0, len(filenames_uxuy)):
    ux_ax.plot(y_arrays_uxuy[i], ux_arrays[i])
    legend_list.append(f"$\\alpha$ = {alphas[i]}°")

ux_ax.grid()
ux_ax.set_ylabel('$U_x/V_\\infty$')
ux_ax.set_xlabel('$y/c$')
ux_ax.set_title("$U_x$ along c/4")
ux_ax.legend(legend_list, loc='upper right')

plt.show()