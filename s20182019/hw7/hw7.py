import numpy as np
from numpy.linalg import inv

Lx = 100e-9
m_xx = 0.19
Ly = 100e-9
m_yy = 0.19
Lz = 5e-9
m_zz = 0.91

T = 300
kB = 1.38066e-23
hbar = 1.0545718e-34
q = 1.602192e-19
m0 = 9.109534e-31
Vt = kB*T/q
coef = 2*Lx*Ly/(2*np.pi)*np.sqrt(m_xx*m_yy)*m0/(hbar*hbar)*(kB*T)

dz = 0.1e-9

z = np.linspace(0,Lz,51)


n_values = np.arange(1,21)
E_F = q*np.linspace(-0.1,0.1,21)  # Fermi energy from -0.1 eV to 0.1 eV.
E_z = (hbar*hbar)/(2*m_zz*m0)*(np.pi*n_values/Lz)*(np.pi*n_values/Lz)

e_d = []
e_d_int = []
index = 0

for E_f in E_F:
    value = coef*np.log(1+np.exp(-(E_z-E_f)/(kB*T))) # 1*n_values matrix
    N_total = np.sum(value)
    subband = 2/(Lx*Ly*Lz)*(np.sin(np.outer(n_values,z)*np.pi/Lz))**2   # n_values * 51 matrix
    e_d = np.matmul(value,subband)  #1*51 matrix, electron density
    print(value.shape,subband.shape)
    e_d_int = np.append(e_d_int, np.sum(e_d)*dz)
    
    f = open('density_'+str(index)+'.dat','w')
    for i in range(51):
        f.write("{}\t{}\n".format(z[i]*1e9,e_d[i]*1e-6))
    f.close()
    index+=1

f = open('integrated_density.dat','w')
for i in range(21):
    f.write("{}\t{}\n".format((i-10)*0.01,e_d_int[i]*1e-4))
f.close()

