import numpy as np
from numpy.linalg import inv

## In this code, we consider the two layers of the same thickness,
## a (nm) and b (nm), respectively. The thickness and the dielectric
## constant could be given by one's preference.
## In this homework, we consider GaAs / GaP heterostructure.

a = 10.0    # thickness of material 1 (nm)
b = 20.0    # thickness of material 2 (nm)
eps1 = 11.1 # dielectric constant of material 1
eps2 = 12.9 # dielectric constant of material 2
N_unit = 2.0     # Number of points per unit length (nm)
mid = int(a*N_unit)
N = int((a+b)*N_unit+1)

######## Analytic calculation of the capacaity ########

eps0 = 8.85E-12 # Permitivity at vacuum

cap1 = eps1*eps0/(a)*1e+9*1e-4  # 1e+9 : nm to meter, 1e-4 : per unit area!
cap2 = eps2*eps0/(b)*1e+9*1e-4

cap_tot = 1./(1./cap1 + 1./cap2)

print("cap1 = {}".format(cap1))
print("cap2 = {}".format(cap2))
print("cap_tot = {}".format(cap_tot))

A = np.zeros(shape=(N,N),dtype=np.float64)
A[0][0] = 1.0
for i in range(1,N-1):
    if i < mid:
        A[i][i-1] = eps1
        A[i][i] = -2.0*eps1
        A[i][i+1] = eps1
    elif i is mid:
        A[i][i-1] = eps1
        A[i][i] = -eps1-eps2
        A[i][i+1] = eps2
    else:
        A[i][i-1] = eps2
        A[i][i] = -eps2*2.0
        A[i][i+1] = eps2

A[N-1][N-1] = 1.0

b = np.zeros(shape=(N))
b[N-1] = 1.0

A_inv = inv(A)

x = np.dot(A_inv,b)

print("V at mat 1 = {}".format(x[mid]))
print("V at mat 2 = {}".format(1-x[mid]))

print("\nV1/V2 = {}.".format(x[mid]/(1-x[mid])))
print("C2/C1 = {}.".format(cap2/cap1))

f = open("vec.dat",'w')
for i in range(1,N):
    f.write("{}\t{}\n".format(float(i/N_unit),x[i]))
f.close()
