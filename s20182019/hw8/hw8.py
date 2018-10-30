import numpy as np
from numpy.linalg import inv
import sys

hbar = 6.626176e-34 / (2*np.pi)
q = 1.602192e-19
m0 = 9.109534e-31
kB = 1.380662e-23
eps0 = 8.854187817e-12
T = 300.0
Vt = kB*T/q
mxx = 0.19
myy = 0.19
mzz = 0.91
dz = 0.1e-9
Nz = 61


z = np.linspace(0,6e-9,Nz)
mid1 = 5
mid2 = 55
eps_si = 11.7
eps_ox = 3.9

N_acc = 1.0e24
Ni = 1.075e16

Lx = 100e-9
Ly = 100e-9
Lz = 5e-9

coef = dz*dz*q/eps0
Ec_Ei = 0.561004 # E_c-E_i, eV
ed_1 = []
ed_2 = []
delE = 0
Eb = 0

f_cl = open("integ_density_cl.dat",'w')
f_qm = open("integ_density_qm.dat",'w')


phi = np.zeros(shape=(Nz),dtype=np.float64)

Vg_grid = np.linspace(0,1.0,21)

for Vg in Vg_grid:
    phi_bound = 0.33374 + Vg
    phi[0] = phi_bound
    phi[Nz-1] = phi_bound
    ed1 = np.zeros(shape=(Nz),dtype=np.float64)

    for i in range(40):
        res = np.zeros(shape=(Nz),dtype=np.float64)
        J = np.zeros(shape=(Nz,Nz),dtype=np.float64)
        res[0] = 0
        res[Nz-1] = 0
        J[0][0] = 1.0
        J[Nz-1][Nz-1] = 1.0
        for j in range(1,Nz-1):
            if j < mid1 or j>mid2:
                res[j] = eps_ox*(phi[j+1] - 2.0*phi[j] + phi[j-1])
                J[j][j-1] = eps_ox
                J[j][j] = -2.0*eps_ox
                J[j][j+1] = eps_ox
            elif j is mid1:
                res[j] = eps_si*(phi[j+1]-phi[j]) + eps_ox*(phi[j-1]-phi[j])
                J[j][j-1] = eps_ox
                J[j][j] = -eps_ox-eps_si
                J[j][j+1] = eps_si
            elif j is mid2:
                res[j] = eps_ox*(phi[j+1]-phi[j]) + eps_si*(phi[j-1]-phi[j])
                J[j][j-1] = eps_si
                J[j][j] = -eps_ox-eps_si
                J[j][j+1] = eps_ox
            else:
                res[j] = eps_si*(phi[j+1] - 2.0*phi[j] + phi[j-1])
                J[j][j-1] = eps_si
                J[j][j] = -2.0*eps_si
                J[j][j+1] = eps_si

        for j in range(mid1,mid2+1):
            if j is mid1 or j is mid2:
                res[j] = res[j] - coef*(N_acc+Ni*np.exp(phi[j]/Vt))*0.5
                J[j][j] = J[j][j] - coef*Ni*np.exp(phi[j]/Vt)/Vt*0.5
            else:
                res[j] = res[j] - coef*(N_acc+Ni*np.exp(phi[j]/Vt))
                J[j][j] = J[j][j] - coef*Ni*np.exp(phi[j]/Vt)/Vt

        update = np.dot(inv(J),-res)
        phi = phi + update
    ed1[mid1:mid2+1] = Ni*np.exp(phi[mid1:mid2+1]/Vt)
    
    f_cl.write("{}\t{}\n".format(Vg,np.sum(ed1*dz)))

    f = open("ed_cl_"+str(int(Vg*20))+".dat",'w')
    for k in range(61):
        f.write("{}\t{}\n".format(0.1*k,ed1[k]))
    f.close()

    Eb = 100
    delE = 100
    iNewton = 0

#    for iNewton in range(int(sys.argv[3])):
    while(delE > 1e-5):
        ed = np.zeros(shape=(Nz),dtype=np.float64)
        for iValley in range(3):
            mass = np.ones(shape=(3),dtype=np.float64)*0.19
            mass[iValley] = 0.91
            coef_Sch = 2*Lx*Ly/(2*np.pi)*np.sqrt(mass[0]*mass[1])*m0/(hbar**2)*(kB*T)
        # Schrodinger solver starts here
            V = q*Ec_Ei - q*phi
            Nbulk = mid2-mid1-1
            H = np.zeros(shape=(Nbulk,Nbulk),dtype=np.float64)
            H[0][0] = -2.0
            H[0][1] = 1.0
            for j in range(1,Nbulk-1):
                H[j][j-1] = 1.0
                H[j][j] = -2.0
                H[j][j+1] = 1.0
            H[Nbulk-1][Nbulk-1] = -2.0
            H[Nbulk-1][Nbulk-2] = 1.0
            for j in range(Nbulk):
                H[j][j] = H[j][j] - 2*mass[2]*m0*(dz/hbar)**2*V[j+mid1+1]
            D, Vk = np.linalg.eig(H)
            Ez = D/(-2.0*mass[2]*m0*(dz/hbar)**2)
            Ezs = np.sort(Ez)
            Vs = Vk[:,Ez.argsort()]
            Vs = Vs.transpose()
            Vs = Vs**2

            for j in range(Nbulk):
                norm = np.sum(Vs[j])*dz
                Vs = Vs/norm
        
            value = coef_Sch*np.log(1+np.exp(-Ezs/(kB*T)))
            N_total = np.sum(value)
            E_total = np.sum(np.sum(value*Ezs))*6.242e+21
            ed_tmp = np.matmul(value,Vs)/Lx/Ly*2*2
            ed[mid1+1:mid2] += ed_tmp
#        for j in range(mid1+1,mid2):
#            ed[j] = ed[j] + ed_tmp[j-mid1-1]

        ed1 = ed
        for i in range(1):
            res = np.zeros(shape=(Nz),dtype=np.float64)
            J = np.zeros(shape=(Nz,Nz),dtype=np.float64)
            res[0] = 0
            res[Nz-1] = 0
            J[0][0] = 1.0
            J[Nz-1][Nz-1] = 1.0
            for j in range(1,Nz-1):
                if j < mid1 or j>mid2:
                    res[j] = eps_ox*(phi[j+1]-2.0*phi[j]+phi[j-1])
                    J[j][j-1] = eps_ox
                    J[j][j] = -2.0*eps_ox
                    J[j][j+1] = eps_ox
                elif j is mid1:
                    res[j] = eps_si*(phi[j+1]-phi[j]) + eps_ox*(phi[j-1]-phi[j])
                    J[j][j-1] = eps_ox
                    J[j][j] = -eps_ox-eps_si
                    J[j][j+1] = eps_si
                elif j is mid2:
                    res[j] = eps_ox*(phi[j+1]-phi[j]) + eps_si*(phi[j-1]-phi[j])
                    J[j][j-1] = eps_si
                    J[j][j] = -eps_ox-eps_si
                    J[j][j+1] = eps_ox
                else:
                    res[j] = eps_si*(phi[j+1] - 2.0*phi[j] + phi[j-1])
                    J[j][j-1] = eps_si
                    J[j][j] = -2.0*eps_si
                    J[j][j+1] = eps_si

            for j in range(mid1,mid2+1):
                if j is mid1 or j is mid2:
                    res[j] = res[j] - coef*(N_acc+ed[j])
                    J[j][j] = J[j][j] - coef*ed[j]/Vt
                else:
                    res[j] = res[j] - coef*(N_acc+ed[j])
                    J[j][j] = J[j][j] - coef*ed[j]/Vt

            update = np.dot(inv(J),-res)
            phi = phi + update

        delE = np.abs(E_total - Eb)
        Eb = E_total

    f = open("ed_cor_"+str(int(Vg*20))+".dat",'w')
    for k in range(61):
        f.write("{}\t{}\n".format(0.1*k,ed[k]))
    f.close()
    f_qm.write("{}\t{}\n".format(Vg,np.sum(ed1)*dz))


