import numpy as np
from numpy.linalg import inv

N = 61

phi = np.zeros(shape=(101,N),dtype=np.float64)
phi[:][:] = 0.33374

T = 300.0   # Temperature (K)
kB = 1.38066E-23    # Boltzmann constant (J/K)
N_acc = 1.0E+24  # Dopant density (# / cm^3)
Ni = 1.075E+16  # Intrinsic density
q = 1.602192E-19 # Elementary charge (C)

Vt = kB*T/q    # Thermal voltage, V

eps0 = 8.85E-12 # Vacuum Permitivity
eps1 = 11.7 # Si permitivity (relative)
eps2 = 3.9 # SiO2 permitivity (relative)

dx = 0.1e-9 # 0.1 nm Spacing

mid1 = 5;   # interface at x=0.5 nm
mid2 = 55;  # interface at x=5.5 nm
coef = dx*dx*q/eps0
phi_now = np.zeros(shape=(N),dtype=np.float64)

e_d = np.zeros(shape=(101,N),dtype=np.float64)
e_d_int = np.zeros(shape=(101),dtype=np.float64)

for j in range(101):
    if j > 0:
        phi_now = phi[j-1]
        phi_now[0] = phi_now[0] + 0.01*j
        phi_now[N-1] = phi_now[0]
    else:
        phi_now = phi[j]
    for i in range(20):
        J = np.zeros(shape=(N,N),dtype=np.float64)  # Jacobian
        res = np.zeros(shape=(N),dtype=np.float64)
        res[0] = phi_now[0] - (0.33374+0.01*j)
        J[0][0] = 1.0
        res[N-1] = phi_now[N-1] - (0.33374+0.01*j)
        J[N-1][N-1] = 1.0
        for ii in range(1,N-1): # Laplacian part
            if ii < mid1 or ii > mid2:
                res[ii] = eps2*(phi_now[ii+1]-2.0*phi_now[ii]+phi_now[ii-1])
                J[ii][ii-1] = eps2
                J[ii][ii] = -2.0*eps2
                J[ii][ii+1] = eps2
            elif ii is mid1:
                res[ii] = eps1*(phi_now[ii+1]-phi_now[ii])+eps2*(phi_now[ii-1]-phi_now[ii])
                J[ii][ii-1] = eps2
                J[ii][ii] = -eps1-eps2
                J[ii][ii+1] = eps1
            elif ii is mid2:
                res[ii] = eps1*(phi_now[ii-1]-phi_now[ii])+eps2*(phi_now[ii+1]-phi_now[ii])
                J[ii][ii-1]=eps1
                J[ii][ii] = -eps1-eps2
                J[ii][ii+1] = eps2
            else:
                res[ii] = eps1*(phi_now[ii+1]-2*phi_now[ii]+phi_now[ii-1])
                J[ii][ii-1] = eps1
                J[ii][ii] = -2*eps1
                J[ii][ii+1] = eps1
        for ii in range(mid1,mid2+1):
            if ii is mid1:
                res[ii] = res[ii] - coef*(N_acc+Ni*np.exp(phi_now[ii]/Vt))*0.5
                J[ii][ii] = J[ii][ii] - coef*Ni*np.exp(phi_now[ii]/Vt)/Vt*0.5
            elif ii is mid2:
                res[ii] = res[ii] - coef*(N_acc+Ni*np.exp(phi_now[ii]/Vt))*0.5
                J[ii][ii] = J[ii][ii] - coef*Ni*np.exp(phi_now[ii]/Vt)/Vt*0.5
            else:
                res[ii] = res[ii] - coef*(N_acc+Ni*np.exp(phi_now[ii]/Vt))
                J[ii][ii] = J[ii][ii] - coef*Ni*np.exp(phi_now[ii]/Vt)/Vt

        
        update = np.dot(inv(J),-res)
        phi_now = phi_now + update
    
    phi[j] = phi_now

    for i in range(mid1,mid2+1):
        e_d[j][i] = Ni * np.exp(phi_now[i]/Vt)
        if i is mid1 or i is mid2:
            e_d_int[j] += dx*0.5*e_d[j][i]
        else:
            e_d_int[j] += dx*e_d[j][i]

#    f = open("d"+str(j)+".dat",'w')
#    for k in range(N):
#        f.write("{}\t{}\n".format(k,phi_now[k]))

f = open("integrated_density.dat",'w')
for j in range(101):
    f.write("{}\t{}\n".format(0.01*j,e_d_int[j]))
f.close()



