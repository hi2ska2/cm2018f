import numpy as np
from numpy.linalg import inv
import sys

hbar = 6.626176e-34 / (2*np.pi)
q = 1.602192e-19
m0 = 9.109534e-31
kB = 1.380662e-23
eps0 = 8.854187817e-12
T = 300.0
dx = 0.1e-9
N = 301    # 300 grids
tau = 1e-6  # Relaxation time in sec

x = np.linspace(0,30e-9,N)
mid1 = 100  # Interface 1
mid2 = 200  # Interface 2

f0_stack = np.zeros(shape=(91,N),dtype=np.float64)
f1_stack = np.zeros(shape=(91,N),dtype=np.float64)

H_grid = np.linspace(0.1,1.0,91)
#H = 0.1 # eV value
VD = 0.001  # V value
V = np.zeros(shape=(N),dtype=np.float64)
V[mid1:mid2] = np.linspace(0,1,100)
V[mid2:N] = VD

index = 0
for H in H_grid:
    fs = np.sqrt(2.0*np.pi)/(1+np.exp(q*H/(kB*T)))    #BCs at x=0
    fd = np.sqrt(2.0*np.pi)/(1+np.exp(q*(H+VD)/(kB*T))) # BCs at x=30nm

    V = np.zeros(shape=(N),dtype=np.float64)
    V[mid1:mid2] = np.linspace(0,1,100)
    V[mid2:N] = VD

    A = np.zeros(shape=(N,N),dtype=np.float64)

    A[0][0] = 1.0

    for i in range(1,N-1):
        c1 = H + 0.5*(V[i-1]+V[i])
        c2 = H + 0.5*(V[i]+V[i+1])
        A[i][i-1] = c1
        A[i][i] = -c1-c2
        A[i][i+1] = c2

    A[N-1][N-1] = 1.0

    b = np.zeros(shape=(N),dtype=np.float64)

    b[0] = fs
    b[N-1] = fd

    f0 = np.dot(inv(A),b)
    f1 = np.zeros_like(f0)


    f1[0] = (-3.0*f0[0]+4.0*f0[1]-f0[2])*0.5/dx
    f1[N-1] = (f0[N-3]-4.0*f0[N-2]+3.0*f0[N-1])*0.5/dx
    for i in range(1,N-1):
        df0 = (f0[i+1]-f0[i-1])*0.5/dx
        f1[i] = -tau*np.sqrt(2.0*q*(H+V[i])/m0)*df0

    f0_stack[index]=f0
    f1_stack[index]=f1
    index = index+1

f1 = open('f0.dat','w')
f2 = open('f1.dat','w')

for i in range(91):
    for j in range(N):
        f1.write("{}\t{}\t{}\n".format(j*0.1,H_grid[i],f0_stack[i][j]))
        f2.write("{}\t{}\t{}\n".format(j*0.1,H_grid[i],f1_stack[i][j]))

f1.close()
f2.close()
