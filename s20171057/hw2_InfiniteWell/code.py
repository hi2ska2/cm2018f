import numpy as np
import matplotlib.pyplot as plt

###### input values ######
width = 5       # unit of nm
effective_mass = 0.19       # unit of electron mass
##########################

mass_electron = 9.10938356E-31               # unit of kg
mass_particle = effective_mass*mass_electron
h_bar = 1.054571800E-34         # unit of J*sec
Joule_to_eV = 6.242E+18
const_for_E = (h_bar**2 / (2*mass_particle)) *Joule_to_eV
E_history = []

N_grid = [5, 10, 20, 30, 50, 100, 200, 300, 500, 1000]
for N in N_grid:
    del_x = width*1E-9 / (N-1) 
    # define operator for (k * del_x)^2
    Operator = -2*np.eye(N-2, dtype=int) + np.eye(N-2, k=1, dtype=int) + np.eye(N-2, k=-1, dtype=int)
    eigValue , eigVector = np.linalg.eig(Operator)
    min_k_squere = -np.amax(eigValue)/(del_x**2)
    _E = const_for_E * min_k_squere
    E_history = np.append(E_history,_E)

plt.figure(tight_layout=True, figsize=(6,5))
plt.semilogx(N_grid,E_history, 'ro', markersize=8)
plt.xlabel("Number of grid", fontsize = 15)
plt.ylabel("Ground state energy (eV)",fontsize = 15)
plt.title("Infinite Well Convergence",fontsize=20)
plt.savefig("./infinitewell_result.png")
plt.show()
