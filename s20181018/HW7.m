close all;
clear all;
h = 6.626176e-34; % Planck constant, J s
hbar = h / (2*pi); % Reduced Planck constant, J s
q = 1.602192e-19; % Elementary charge, C
m0 = 9.109534e-31; % Electron rest mass, kg
k_B = 1.380662e-23; % Boltzmann constant, J/K
T = 300; % Temperature, K
Lx = 100e-9; Ly = 100e-9; Lz = 5e-9; % Lengths, m
mxx = 0.19; myy = 0.19; mzz = 0.91; % Masses, m0
nmax = 10;
coef = 2*Lx*Ly/(2*pi)*sqrt(mxx*myy)*m0/(hbar^2)*(k_B*T);
totalNumber = 0;
Nz = 51;
z = transpose([0:Nz-1])*Lz/(Nz-1);
elec = zeros(Nz,41); % Electron density, /m^3

for j = 1:41
    Fermi(1,j) = 0.005*(21-j);
    for n = 1:nmax
        Ez = (hbar^2)/(2*mzz*m0)*(pi*n/Lz)^2;
        subbandNumber = coef*log(1+exp((-Ez+(Fermi(1,j)*q))/(k_B*T)));
        totalNumber = totalNumber+subbandNumber;
        elec(:,j) = elec(:,j)+2/(Lx*Ly*Lz)*(sin(n*pi*z/Lz).^2)*subbandNumber;
    end
end

for j=1:41
    figure(1);
    plot(z/1e-9,elec(:,j)/1e6);
    hold all;
end
xlabel('Distance of z-direction (nm)');
ylabel('Electron density (cm^-^3)');
legend('0.100eV','0.095eV','0.090eV','0.085eV','0.080eV','0.075eV','0.070eV','0.065eV','0.060eV','0.055eV','0.050eV','0.045eV','0.040eV','0.035eV','0.030eV','0.025eV','0.020eV','0.015eV','0.010eV','0.005eV','0eV','-0.005eV','-0.010eV','-0.015eV','-0.020eV','-0.025eV','-0.030eV','-0.035eV','-0.040eV','-0.045eV','-0.050eV','-0.055eV','-0.060eV','-0.065eV','-0.070eV','-0.075eV','-0.080eV','-0.085eV','-0.090eV','-0.095eV','-0.100eV');
dz = Lz/(Nz-1);
elec_int = sum(elec)*dz;
figure(2);
plot(Fermi(1,:),elec_int(1,:)/1e4);
xlabel('Fermi enery (eV)');
ylabel('Integrated electron density (cm^-^2)');
