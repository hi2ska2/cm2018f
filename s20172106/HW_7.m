clear;
clc;

h = 6.62617e-34; % Planck constant, J s
hbar = h/(2*pi); % Reduced planck constant, J s
q = 1.602192e-19; % Elementary charge, C
m0 = 9.109534e-31; % Electron rest mass, kg
k_B = 1.380662e-23; % Boltzmann constant, J/K
T = 300.0; % Temperature, K
Ef = linspace(-0.1,0.1,10)'; % Fermi energy, eV
Lx = 100e-9; Ly = 100e-9; Lz = 5e-9; % Length, m
mxx = 0.19; myy = 0.19; mzz = 0.91; % Effective masses, m_eff
Nmax = 10;
Coef = 2*Lx*Ly/(2*pi)*sqrt(mxx*myy)*m0/(hbar^2)*(k_B*T);

TotalNumber = 0 ;
Nz = 51;
z = transpose([0:Nz-1])*Lz/(Nz-1);
Elec = zeros(Nz,10); % Electron density, 1/m^3;

for e = 1:10
    Ef_current = Ef(e); 
    for n = 1:Nmax
        Ez = (hbar^2)/(2*mzz*m0)*(pi*n/Lz)^2;
        SubbandNumber = Coef*log(1 + exp(-(Ez-(Ef_current*q))/(k_B*T)));
        TotalNumber = TotalNumber + SubbandNumber;
        Elec(:,e) = Elec(:,e) + 2/(Lx*Ly*Lz)*(sin(n*pi*z/Lz).^2)*SubbandNumber;
    end
    TotalNumber = 0 ;
    
end

figure;
plot(z/1e-9,Elec(:,1)/1e6,z/1e-9,Elec(:,2)/1e6,z/1e-9,Elec(:,3)/1e6,z/1e-9,Elec(:,4)/1e6,z/1e-9,Elec(:,5)/1e6,z/1e-9,Elec(:,6)/1e6,z/1e-9,Elec(:,7)/1e6,z/1e-9,Elec(:,8)/1e6,z/1e-9,Elec(:,9)/1e6,z/1e-9,Elec(:,10)/1e6)
xlabel('z [nm]')
ylabel('Electron density [1/cm^3]')

Inte_Elec = sum(Elec./(1e14),1);
figure;
plot(Ef,Inte_Elec);
xlabel('Fermi Energy [eV]')
ylabel('Integrated electron density [1/cm^2]')
