close all;
clear all;



h = 6.626176e-34; %Planck constant (Js)
h_bar = h / (2*pi); %(Js)
Charge_q = 1.602192e-19; %charge Q
m0 = 9.109534e-31; %Electron reset mass (Kg)
Boltz = 1.380662e-23; %Boltzmann constant (J/K)
Temp = 300.0; %Temperature(K)
Lx = 100e-9; %Length of x-direction (m)
Ly = 100e-9;  %Length of y-direction (m)
Lz = 5e-9;   %Length of z-direction (m)
mxx = 0.19;  %Mass
myy = 0.19;  %Mass
mzz = 0.91;  %Mass

iter_max = 10;

coef = 2*Lx*Ly/(2*pi)*sqrt(mxx*myy)*m0/(h_bar^2)*(Boltz*Temp);
Fermi = [-0.1:0.01:0.1] * Charge_q ;%%-0.1 ~ 0.1ev(21Point Fermi_Level)
%Fermi = [-0.1:0.01:0.1] * 1 ;%%-0.1 ~ 0.1ev(21Point Fermi_Level)

totalNumber = 0;
Nz = 51;
z =  transpose([0:Nz-1])*Lz/(Nz-1);
del_z = Lz/(Nz-1);
elec = zeros(Nz,1); %Electron density (m^3)
Elec_density = cell(size(Fermi,2),1);

for iteriter=1:1:size(Fermi,2)
    for iter=1:1:iter_max
        Ez = (h_bar^2) / (2*mzz*m0) * (pi*iter/Lz)^2;
        subbandNumber = coef*log(1+exp(-(Ez-Fermi(1,iteriter))/(Boltz*Temp)));

        totalNumber = totalNumber + subbandNumber;
        elec = elec + 2/(Lx*Ly*Lz)*(sin(iter*pi*z/Lz).^2)*subbandNumber;
    end
    Elec_density{iteriter,1} = elec;
    elec = zeros(Nz,1); %Electron density (1/m^3)
end

for iter=1:1:size(Elec_density,1)
    Integrated_elec_density(iter,1) = sum(Elec_density{iter,1}) * (del_z); %(cm^2);
end

figure(1)
for iter=1:1:size(Fermi,2)
    plot(z/1e-9,Elec_density{iter,1}/1e6,'LineWidth',2); hold all;
end
xlabel('Distance of Z-direction[nm]');
ylabel('Electron density [cm^-3]');
grid on;
legend('-0.1eV','-0.09eV','-0.08eV','-0.07eV','-0.06eV','-0.05eV','-0.04eV','-0.03eV','-0.02eV','-0.01eV','0eV','0.01eV','0.02eV','0.03eV','0.04eV','0.05eV','0.06eV','0.07eV','0.08eV','0.09eV','0.1eV','Location','best');


figure(2)
plot(Fermi/Charge_q,Integrated_elec_density*1e-4)
xlabel('Fermi Level [eV]');
ylabel('Integrated Electron density [cm^-2]');
grid on;

figure(3)
semilogy(Fermi/Charge_q,Integrated_elec_density*1e-4)
xlabel('Fermi Level [eV]');
ylabel('Integrated Electron density [cm^-2]');
grid on;
