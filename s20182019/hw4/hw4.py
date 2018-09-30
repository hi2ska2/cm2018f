import numpy as np
from numpy.linalg import inv

T = 300.0   # Temperature (K)
kB = 1.38066E-23    # Boltzmann constant (J/K)
N_acc = 1.0E+24  # Dopant density (# / cm^3)
Ni = 1.075E+16  # Intrinsic density
q = 1.602192E-19 # Elementary charge (C)

eps0 = 8.85E-12 # Vacuum Permitivity
eps1 = 11.7 # Si permitivity (relative)
eps2 = 3.9 # SiO2 permitivity (relative)

dx = 0.1e-9 # 0.1 nm Spacing
N = 61  # 6 nm thickness

mid1 = 5;   # interface at x=0.5 nm
mid2 = 55;  # interface at x=5.5 nm

A = np.zeros(shape=(N,N),dtype=np.float64)
A[0][0] = 1.0

for i in range(1,N-1):
    if i < mid1:
        A[i][i-1] = eps2
        A[i][i] = -2.0*eps2
        A[i][i+1] = eps2
    elif i is mid1:
        A[i][i-1] = eps2
        A[i][i] = -(eps1+eps2)
        A[i][i+1] = eps1
    elif i < mid2:
        A[i][i-1] = eps1
        A[i][i] = -2.0*eps1
        A[i][i+1] = eps1
    elif i is mid2:
        A[i][i-1] = eps1
        A[i][i] = -(eps1+eps2)
        A[i][i+1] = eps2
    else:
        A[i][i-1] = eps2
        A[i][i] = -2.0*eps2
        A[i][i+1] = eps2
A[N-1][N-1] = 1.0

A_inv = inv(A)
A_inv1 = inv(A)

b = np.zeros(shape=(N),dtype=np.float64)
bf = np.zeros(shape=(N),dtype=np.float64)
for i in range(11):
    b[0] = 0.33374-0.1*i
    b[N-1] = b[0]

    for j in range(mid1,mid2+1):
        if j is mid1:
            b[j] = dx*dx*q*N_acc*0.5/eps0
        elif j is mid2:
            b[j] = dx*dx*q*N_acc*0.5/eps0
        else:
            b[j] = dx*dx*q*N_acc/eps0

    bf[0] = 0.33374-0.1*i
    bf[N-1] = bf[0]



    x = np.zeros(shape=(N),dtype=np.float64)

    x = np.dot(A_inv,b) # Initial potential


    E_d = np.zeros(shape=(N),dtype=np.float64) # Electron density

    for j in range(mid1,mid2+1):
        E_d[j] = Ni*np.exp(q*x[j]/(kB*T))
        if j is mid1:
            bf[j] = dx*dx*q*(N_acc+E_d[j])*0.5/eps0
        elif i is mid2:
            bf[j] = dx*dx*q*(N_acc+E_d[j])*0.5/eps0
        else:
            bf[j] = dx*dx*q*(N_acc+E_d[j])/eps0

    xf = np.zeros(shape=(N),dtype=np.float64)
    xf = np.dot(A_inv,bf)

    f = open('data'+str(i)+'.dat','w')
    for j in range(N):
        f.write("{}\t{}\t{}\t{}\n".format(float(j*0.1),x[j],xf[j],E_d[j]))
    f.close()


