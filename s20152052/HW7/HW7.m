clear all
clc

h = 6.626176e-34; % Planck constant, J s 
hbar = h / (2*pi); % Reduced Planck constant, J s 
q = 1.602192e-19; % Elementary charge, C 
m0 = 9.109534e-31; % Electron rest mass, kg 
k_B = 1.380662e-23; % Boltzmann constant, J/K 
T = 300.0; % Temperature, K 
Lx = 100e-9; Ly = 100e-9; Lz = 5e-9; % Lenghs, m 
mxx = 0.19; myy = 0.19; mzz = 0.91; % Masses, m0  
nmax = 10; 
coef = 2*Lx*Ly/(2*pi)*sqrt(mxx*myy)*m0/(hbar^2)*(k_B*T);
E_F=transpose([-0.1:0.01:0.1]); % feri energy

for i=1:length(E_F)

totalNumber=0;
Nz = 51; 
z = transpose([0:Nz-1])*Lz/(Nz-1); 
elec = zeros(Nz,1); % Electron density, /m^3 
for n=1:nmax   
    Ez = (hbar^2)/(2*mzz*m0)*(pi*n/Lz)^2;   
    subbandNumber = coef*log(1+exp(-(Ez-E_F(i,1)*q)/(k_B*T)));   
    totalNumber = totalNumber + subbandNumber;   
    elec = elec + 2/(Lx*Ly*Lz)*(sin(n*pi*z/Lz).^2)*subbandNumber; 
end
elec1(:,i)=elec;
elec_int=sum(elec*Lz/(Nz-1),1); 
elec_ef(i,1)=elec_int; % fermi level¿¡ µû¸¥ electon density
end

f1=figure
plot(z/1e-9,elec1*1e-6)
xlabel('z (nm)')
ylabel('Electron density (cm^3)')

f2=figure
subplot(1,2,1)
semilogy(E_F,elec_ef*1e-4)
xlabel('Fermi energy(eV)')
ylabel('Electron density (cm^2)')
subplot(1,2,2)
plot(E_F,elec_ef*1e-4)
xlabel('Fermi energy(eV)')
ylabel('Electron density (cm^2)')
