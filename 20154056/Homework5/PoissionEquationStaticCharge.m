%% Define the geometry
Deltax = 0.1e-9; % 0,1 nm spacing
N = 61; % 6 nm thick
interface1 = 6; % At x = 0.5 nm
interface2 = 56; % At x = 5.5 nm
epsSi = 11.7; epsOx = 3.9;
Nacc = 1e24;
ni = 1.075e16;
x = Deltax * transpose([0:N-1]);
GateVoltage = 1;
Voltage = 0.33374 + GateVoltage;

%% Constants
q = 1.6022e-19;
eps0 = 8.854e-12;
kB = 1.38e-23;
T = 300;

%% Define poission matrix
M = zeros(N, N);
M(1, 1) = 1;
for ii = 2 : N-1
    if ii < interface1 
        M(ii, ii-1) = epsOx; M(ii, ii) = -2 * epsOx; M(ii, ii+1) = epsOx;
    elseif ii == interface1
        M(ii, ii-1) = epsOx; M(ii, ii) = -epsOx-epsSi; M(ii, ii+1) = epsSi;
    elseif ii < interface2
        M(ii, ii-1) = epsSi; M(ii, ii) = -2 * epsSi; M(ii, ii+1) = epsSi;
    elseif ii == interface2
        M(ii, ii-1) = epsSi; M(ii, ii) = -epsSi-epsOx; M(ii, ii+1) = epsOx;
    elseif ii > interface2
        M(ii, ii-1) = epsOx; M(ii, ii) = -2 * epsOx; M(ii, ii+1) = epsOx;
    end
end
M(N, N) = 1;

%% Define boundary conditions
b = zeros(N, 1);
b(1, 1) = Voltage;
for ii = interface1 : interface2
    if ii == interface1
        b(ii, 1) = Deltax*Deltax*q*Nacc/eps0*0.5;
    elseif ii == interface2
        b(ii, 1) = Deltax*Deltax*q*Nacc/eps0*0.5;
    else
        b(ii, 1) = Deltax*Deltax*q*Nacc/eps0;
    end
end
b(N, 1) = Voltage;

%% Get the potential phi and electron density
phi = M \ b;
eleD = zeros(N, 1);
for ii = interface1 : interface2
    eleD(ii, 1) = ni * exp(q*phi(ii, 1)/(kB*T));
end
figure(1)
subplot(2, 2, 1)
plot(x/1e-9, eleD*1e-6, 'or-');
title(strcat('Results for : ', sprintf('%2.2f', Voltage), 'Volt'))
ylabel('Electron Density [#/cm^3]')

subplot(2, 2, 3)
plot(x/1e-9, phi, 'ob-')
xlabel('Position [nm]')
ylabel('Potential energy [eV]')

%% Recalculating of the potential
% boundary
b2 = zeros(N, 1);
b2(1, 1) = Voltage;
for ii = interface1 : interface2
    if ii == interface1
        b2(ii, 1) = Deltax*Deltax*q*eleD(ii, 1)/eps0*0.5;
    elseif ii == interface2
        b2(ii, 1) = Deltax*Deltax*q*eleD(ii, 1)/eps0*0.5;
    else
        b2(ii, 1) = Deltax*Deltax*q*eleD(ii, 1)/eps0;
    end
end
b2(N, 1) = Voltage;

% potential and electron density
phi2 = M \ b2;
eleD2 = zeros(N, 1);
for ii = interface1 : interface2
    eleD2(ii, 1) = ni * exp(q*phi2(ii, 1)/(kB*T));
end
figure(1)
subplot(2, 2, 2)
plot(x/1e-9, eleD2*1e-6, 'or-');
title('Reuslts with adapted potential')

subplot(2, 2, 4)
plot(x/1e-9, phi2, 'ob-')
xlabel('Position [nm]')