import numpy as np
from numpy.linalg import inv
import copy

nx = 121 # 120 nm : spacing with 1 nm
nz = 61 # 6 nm : spacing with 0.1 nm

T = 300.0
kB = 1.38066E-23
N_acc_1 = 5.0E+25
N_acc_2 = 2.0E+23

Ni = 1.075E+16
q = 1.602192E-19

Vt = kB*T/q

eps_ox = 3.9
eps_si = 11.7
eps0 = 8.85E-12

midx_1 = 40
midx_2 = 80

midz_1 = 5
midz_2 = 55

Nd = np.ones(shape=(nx))* N_acc_1
Nd[midx_1:midx_2] = N_acc_2

dx = 120e-9/(nx-1)
dz = 6e-9/(nz-1)

rt = int(round(dx/dz))

V_array = np.linspace(0,1,11)

phi_now = np.zeros(shape=(nx*nz),dtype=np.float64)
phi_now += 0.33374
ed = np.zeros(shape=(nx*nz),dtype=np.float64)
coef = (dx*dx*dz)*q/eps0

for V in V_array:
    phi = copy.deepcopy(phi_now)
    phi_bound = 0.33374 + V
    phi[midx_1:midx_2+1] = phi_bound
    phi[midx_1+(nz-1)*nx:(nz-1)*nx+midx_2+1] = phi_bound

    for iNewton in range(10):
        res = np.zeros(shape=(nx*nz),dtype=np.float64)
        J = np.zeros(shape=(nx*nz,nx*nz),dtype=np.float64)
        res[midx_1:midx_2+1] = 0
        res[midx_1+(nz-1)*nx:(nz-1)*nx+midx_2+1] = 0
        for k in range(midx_1,midx_2+1):
            J[k][k] = 1.0
            J[(nz-1)*nx+k][(nz-1)*nx+k] = 1.0
        for k in range(midz_1,midz_2+1):
            J[k*nx][k*nx] = 1.0
            J[k*nx+nx-1][k*nx+nx-1] = 1.0

        for j in range(1,nz-1):
            for k in range(1,nx-1):
                ind1 = j*nx + k
                if j < midz_1 or j > midz_2:
                    res[ind1] = eps_ox*(rt*phi[ind1-nx]+phi[ind1-1]-2.0*(1+rt)*phi[ind1]+phi[ind1+1]+rt*phi[ind1+nx])
                    J[ind1][ind1-nx] = rt*eps_ox
                    J[ind1][ind1-1] = eps_ox
                    J[ind1][ind1] = -2.0*(1+rt)*eps_ox
                    J[ind1][ind1+1] = eps_ox
                    J[ind1][ind1+nx] = rt*eps_ox
                elif j is midz_1:
                    res[ind1] = eps_ox*(rt*phi[ind1-nx]+0.5*(phi[ind1-1]+phi[ind1+1])-(1+rt)*phi[ind1]) + eps_si*(0.5*(phi[ind1-1]+phi[ind1+1])-(1+rt)*phi[ind1]+rt*phi[ind1+nx])
                    J[ind1][ind1-nx] = rt*eps_ox
                    J[ind1][ind1-1] = 0.5*(eps_ox+eps_si)
                    J[ind1][ind1] = -(1.0+rt)*(eps_ox+eps_si)
                    J[ind1][ind1+1] = 0.5*(eps_ox+eps_si)
                    J[ind1][ind1+nx] = rt*eps_si
                elif j is midz_2:
                    res[ind1] = eps_si*(rt*phi[ind1-nx]+0.5*(phi[ind1-1]+phi[ind1+1])-(1+rt)*phi[ind1]) + eps_ox*(0.5*(phi[ind1-1]+phi[ind1+1])-(1+rt)*phi[ind1]+rt*phi[ind1+nx])
                    J[ind1][ind1-nx] = rt*eps_si
                    J[ind1][ind1-1] = 0.5*(eps_ox+eps_si)
                    J[ind1][ind1] = -(1+rt)*(eps_ox+eps_si)
                    J[ind1][ind1+1] = 0.5*(eps_ox+eps_si)
                    J[ind1][ind1+nx] = rt*eps_ox
                else:
                    res[ind1] = eps_si*(rt*phi[ind1-nx]+phi[ind1-1]-2.0*(1+rt)*phi[ind1]+phi[ind1+1]+rt*phi[ind1+nx])
                    J[ind1][ind1-nx] = rt*eps_si
                    J[ind1][ind1-1] = eps_si
                    J[ind1][ind1] = -2.0*(1+rt)*eps_si
                    J[ind1][ind1+1] = eps_si
                    J[ind1][ind1+nx] = rt*eps_si

        for j in range(1,nz-1):
            ind1 = j*nx
            if j < midz_1 or j > midz_2:
                res[ind1] = eps_ox*(0.5*rt*phi[ind1-nx]-(1+rt)*phi[ind1]+phi[ind1+1]+rt*0.5*phi[ind1+nx])
                J[ind1][ind1-nx] = rt*0.5*eps_ox
                J[ind1][ind1] = -(1+rt)*eps_ox
                J[ind1][ind1+1] = eps_ox
                J[ind1][ind1+nx] = rt*0.5*eps_ox

        res[0] = eps_ox*(0.5*phi[1] + rt*0.5*phi[nx] - 0.5*(1+rt)*phi[0])
        J[0][0] = -0.5*(1+rt)*eps_ox
        J[0][1] = 0.5*eps_ox
        J[0][nx] = rt*0.5*eps_ox

        res[(nz-1)*nx] = eps_ox*(0.5*phi[(nz-1)*nx+1]+rt*0.5*phi[(nz-2)*nx]-0.5*(1+rt)*phi[(nz-1)*nx])
        J[(nz-1)*nx][(nz-1)*nx] = -0.5*(1+rt)*eps_ox
        J[(nz-1)*nx][(nz-2)*nx] = rt*0.5*eps_ox
        J[(nz-1)*nx][(nz-1)*nx+1] = 0.5*eps_ox

        for j in range(1,midx_1):
            ind1 = (nz-1)*nx+j
            res[j] = eps_ox*(0.5*phi[j-1]-(1+rt)*phi[j]+0.5*phi[j+1]+rt*phi[j+nx])
            J[j][j-1] = 0.5*eps_ox
            J[j][j] = -(1.0+rt)*eps_ox
            J[j][j+1] = 0.5*eps_ox
            J[j][j+nx] = rt*eps_ox
                    
            res[ind1] = eps_ox*(rt*phi[ind1-nx]+0.5*phi[ind1-1]-(1.0+rt)*phi[ind1]+0.5*phi[ind1+1])
            J[ind1][ind1-nx] = rt*eps_ox
            J[ind1][ind1-1] = 0.5*eps_ox
            J[ind1][ind1] = -(1.0+rt)*eps_ox
            J[ind1][ind1+1] = 0.5*eps_ox

        for j in range(midx_2+1,nx-1):
            ind1 = (nz-1)*nx+j
            res[j] = eps_ox*(0.5*phi[j-1]-(1.0+rt)*phi[j]+0.5*phi[j+1]+rt*phi[j+nx])
            J[j][j-1] = 0.5*eps_ox
            J[j][j] = -(1.0+rt)*eps_ox
            J[j][j+1] = 0.5*eps_ox
            J[j][j+nx] = rt*eps_ox
            res[ind1] = eps_ox*(rt*phi[ind1-nx]+0.5*phi[ind1-1]-(1.0+rt)*phi[ind1]+0.5*phi[ind1+1])
            J[ind1][ind1-nx] = rt*eps_ox
            J[ind1][ind1-1] = 0.5*eps_ox
            J[ind1][ind1] = -(1.0+rt)*eps_ox
            J[ind1][ind1+1] = 0.5*eps_ox

        res[nx-1] = eps_ox*(0.5*phi[nx-2]+rt*0.5*phi[nx-1+nx]-0.5*(1+rt)*phi[nx-1])
        J[nx-1][nx-1] = -0.5*(1+rt)*eps_ox
        J[nx-1][nx-2] = 0.5*eps_ox
        J[nx-1][nx-1+nx] = rt*0.5*eps_ox

        res[nx*nz-1] = eps_ox*(rt*0.5*phi[nx*nz-1-nx]+0.5*phi[nx*nz-2]-0.5*(1+rt)*phi[nx*nz-1])
        J[nx*nz-1][nx*nz-1] = -0.5*(1+rt)*eps_ox
        J[nx*nz-1][nx*nz-2] = 0.5*eps_ox
        J[nx*nz-1][nx*nz-1-nx] = rt*0.5*eps_ox

        for j in range(1,nz-1):
            ind1 = nx*j+nx-1
            if j < midz_1 or j > midz_2:
                res[ind1] = eps_ox*(rt*0.5*phi[ind1-nx]+phi[ind1-1]-(1.0+rt)*phi[ind1]+rt*0.5*phi[ind1+nx])
                J[ind1][ind1-nx] = rt*0.5*eps_ox
                J[ind1][ind1-1] = eps_ox
                J[ind1][ind1] = -(1.0+rt)*eps_ox
                J[ind1][ind1+nx] = rt*0.5*eps_ox

        for j in range(midz_1,midz_2+1):
            for k in range(1,nx-1):
                ind1 = j*nx+k
                if j is midz_1 or j is midz_2:
                    res[ind1] = res[ind1] - coef*(Nd[k]+Ni*np.exp(phi[ind1]/Vt))*0.5
                    J[ind1][ind1] = J[ind1][ind1] - coef*Ni*np.exp(phi[ind1]/Vt)/Vt*0.5
                else:
                    res[ind1] = res[ind1] - coef*(Nd[k]+Ni*np.exp(phi[ind1]/Vt))
                    J[ind1][ind1] = J[ind1][ind1] - coef*Ni*np.exp(phi[ind1]/Vt)/Vt

        update = np.dot(inv(J),-res)
        phi = phi + update
        print(np.linalg.norm(update,np.inf))

    for j in range(midz_1,midz_2+1):
        for k in range(nx):
            ind1 = j*nx+k
            ed[ind1] = Ni*np.exp(phi[ind1]/Vt)
                    
    if V < 20:
        f = open("phi"+str(int(V*10))+".dat",'w')
        g = open("ed"+str(int(V*10))+".dat",'w')
        for j in range(nz):
            for k in range(nx):
                f.write("{}\t{}\t{}\n".format(j*0.1,k,phi[j*nx+k]))
                g.write("{}\t{}\t{}\n".format(j*0.1,k,ed[j*nx+k]))
        f.close()
        g.close()

    phi_now = copy.deepcopy(phi)




