clear all;

%%
% Permittivity
% 
% Aluminum: 1.46
% mos2: 4
% silicon dioxide: 3.9
% silicon: 11.68
%%
h = 6.626176e-34; % Planck constant, J s
hbar = h / (2*pi); % Reduced Planck constant, J s
q = 1.602192e-19; % Elementary charge, C
e0 = 8.854187817e-12; % Vacuum permittivity, F/m
m0 = 9.109534e-31; % Electron rest mass, kg
k_B = 1.380662e-23; % Boltzmann constant, J/K
T = 300.0; % Temperature, K

thermal = k_B*T/q; % Thermal voltage, V


%%

N_point = 5; % N points in minimum thickness



mat_th = [100 100 5]*1e-9; % [Lx Ly Lz] material thickness (x0.1 nm)
mass = [0.19 0.19 0.91] %masses

% Num = mat_th.*N_point./d; % number of points in each material
% N = sum(Num)+1;
% 
% x = dx*transpose([0:sum(Num)]); % real space thickness

nmax = 1000;
coef = 2*mat_th(1)*mat_th(2)/(2*pi)*sqrt(mass(1)*mass(2))*m0/(hbar^2)*(k_B*T);
totalNumber = 0;
Nz = 51;
for Efr = 1:101
    Ef = -0.1 + (Efr-1)/100*0.2;
z = transpose([0:Nz-1])*mat_th(3)/(Nz-1);

elec = zeros(Nz,1); % Electron density, /m^3
for n=1:nmax
 Ez = (hbar^2)/(2*mass(3)*m0)*(pi*n/mat_th(3))^2;
 subbandNumber = coef*log(1+exp(-(Ez-q*Ef)/(k_B*T)));
totalNumber = totalNumber + subbandNumber;
elec = elec + 2/(mat_th(1)*mat_th(2)*mat_th(3))*(sin(n*pi*z/mat_th(3)).^2)*subbandNumber;
end
int_elec=sum(elec*1e-10);
elec_ef(Efr,:) = elec;
int_elec_ef(Efr) = int_elec;
end
figure(1)
plot(z/1e-9,elec_ef(1,:)/1e6,z/1e-9,elec_ef(11,:)/1e6,z/1e-9,elec_ef(21,:)/1e6,z/1e-9,elec_ef(31,:)/1e6,z/1e-9,elec_ef(41,:)/1e6,z/1e-9,elec_ef(51,:)/1e6,z/1e-9,elec_ef(61,:)/1e6,z/1e-9,elec_ef(71,:)/1e6,z/1e-9,elec_ef(81,:)/1e6,z/1e-9,elec_ef(91,:)/1e6,z/1e-9,elec_ef(101,:)/1e6)
legend('-0.1 eV', '-0.08 eV', '-0.06 eV', '-0.04 eV', '-0.02 eV', '-0.00 eV', '0.02 eV', '0.04 eV', '0.06 eV', '0.08 eV', '0.10 eV');
xlabel('position (m)')
ylabel('electron density cm^-^3')
figure(2)
plot(-0.1+[0:100]/100*0.2,int_elec_ef/1e4);
xlabel('fermi energy (eV)')
ylabel('integrated electron density cm^-^2')