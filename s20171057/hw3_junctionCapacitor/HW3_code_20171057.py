####
#
#	This code is implemented to calculate capacity(F/cm^2) of source-free hetero-juction.
#	User should input the width and dielectric constant of both materials.
#	The initial values of input are given for 10nm Bi2Te3(epsilon = 80) and 20nm PbTe(epsilon = 400).
#
####
###### input values ######
width_1 = 10       # unit of nm
epsilon_1 = 80
width_2 = 20       # unit of nm
epsilon_2 = 400
N_grid_per_nm = 1
##########################

import numpy as np
import matplotlib.pyplot as plt

total_width = width_1 + width_2
N_grid_1 = width_1*N_grid_per_nm
N_grid_2 = width_2*N_grid_per_nm
total_N_grid = N_grid_1 + N_grid_2 +1
print("Total width of the capacitor = ",total_width,"nm")
print("Total grid in the width = ",total_N_grid)
if(N_grid_1 != int(N_grid_1)):
	print("number of grid must be integer")
	exit()
if(N_grid_2 != int(N_grid_2)):
	print("number of grid must be integer")
	exit()

_matrix_1 = np.block([
        [np.array([1]),np.zeros(N_grid_1)],
        [np.transpose([np.block([np.array([epsilon_1]),np.zeros(N_grid_1-2)])]),-2*epsilon_1*np.eye(N_grid_1-1, dtype=int) + epsilon_1*np.eye(N_grid_1-1, k=1, dtype=int) + epsilon_1*np.eye(N_grid_1-1, k=-1, dtype=int),np.transpose([np.block([np.zeros(N_grid_1-2),np.array([epsilon_1])])])]
		])
_matrix_2 = np.block([
        [np.transpose([np.block([np.array([epsilon_2]),np.zeros(N_grid_2-2)])]),-2*epsilon_2*np.eye(N_grid_2-1, dtype=int) + epsilon_2*np.eye(N_grid_2-1, k=1, dtype=int) + epsilon_2*np.eye(N_grid_2-1, k=-1, dtype=int),np.transpose([np.block([np.zeros(N_grid_2-2),np.array([epsilon_2])])])], 
		[np.zeros(N_grid_2),np.array([1])]
		])
_matrix_juction = np.block([np.zeros(N_grid_1-1),np.array([epsilon_1,-epsilon_1-epsilon_2,epsilon_2]),np.zeros(N_grid_2-1)])
Matrix = np.block([
	    [_matrix_1, np.zeros((N_grid_1,N_grid_2))],
	    [_matrix_juction],
		[np.zeros((N_grid_2,N_grid_1)), _matrix_2]
		])

vector = np.block([np.zeros(total_N_grid-1),np.array([1])])

phi_numeric = np.matmul(np.linalg.inv(Matrix),vector)
print(phi_numeric[N_grid_1])
print("==============Analytic==============")
epsilon_zero = 8.85E-12
capacity_1 = epsilon_1*epsilon_zero/(width_1*1E-9) *1E-4
capacity_2 = epsilon_2*epsilon_zero/(width_2*1E-9) *1E-4
print("capacity of first material = ",capacity_1)
print("capacity of second material = ",capacity_2)
Capacity = 1/(1/capacity_1 + 1/capacity_2)
print("total capacity = ",Capacity)

print("potential at junction from analytic calculation")
print("= 1 - C*V / c_2 ")
phi_analytic = 1 - Capacity/capacity_2
print("=",phi_analytic) 
slope_1 = phi_analytic/width_1
slope_2 = (1-phi_analytic)/width_2
intercept_2 = phi_analytic - slope_2*width_1

x_gid = np.linspace(0, total_width, total_N_grid)
x_gid_1 = np.linspace(0, width_1, 10)
x_gid_2 = np.linspace(width_1, total_width, 10)

plt.figure(tight_layout=True, figsize=(6,5))
plt.plot(x_gid,phi_numeric, 'ro', markersize=8,label="Potential-numeric")
plt.plot(x_gid_1,x_gid_1*slope_1,color = 'blue')
plt.plot(x_gid_2,x_gid_2*slope_2+intercept_2,color = 'blue',label="Potential-analytic")
plt.xlabel("Width (nm)", fontsize = 15)
plt.ylabel("$\phi$ (V)",fontsize = 15)
plt.legend()
plt.title("source-free heterojuction",fontsize=20)
plt.savefig("./sourcefree_result.png")
plt.show()