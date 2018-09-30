close all;
clear all;

q = 1.602192e-19; % Elementary charge (C)
eps0 = 8.854187817e-12; % Vacuum permittivity, (F/m)
k_B = 1.380662e-23; % Boltzmann constant (J/K) 
T = 300; % Temperature (K)
Delta_x = 0.1e-9; % 0.1 nm spacing
N = 61; % 6 nm thick
x = Delta_x*transpose([0:N-1]); % real space (m)
Inter1 = 6; % At x = 0.5 nm
Inter2 = 56; % At x = 5.5 nm
eps_si = 11.7; % Silicon relative permittivity
eps_ox = 3.9; % Silicon dioxide relative permittivity
Nacc = 1e24; % 1e18/cm^3
ni = 1.075e16; % 1.075e10/cm^3
Voltage_sweep = 0:0.1:1; % Votage step (V)

for i = 1:11;
    gate_potential = 0.33374-Voltage_sweep(i);
    A = zeros(N,N);
    A(1,1) = 1; A(N,N) = 1;
    for j = 2:N-1
        if     (j < Inter1)   A(j,j-1) = eps_ox;    A(j,j) = -2*eps_ox;         A(j,j+1) = eps_ox;
        elseif (j == Inter1)  A(j,j-1) = eps_ox;    A(j,j) = -eps_ox-eps_si;    A(j,j+1) = eps_si;
        elseif (j < Inter2)   A(j,j-1) = eps_si;    A(j,j) = -2*eps_si;         A(j,j+1) = eps_si;
        elseif (j == Inter2)  A(j,j-1) = eps_si;    A(j,j) = -eps_si-eps_ox;    A(j,j+1) = eps_ox;
        elseif (j > Inter2)   A(j,j-1) = eps_ox;    A(j,j) = -2*eps_ox;         A(j,j+1) = eps_ox;
        end
    end
    
    % Boundary value of the electrostatic potential %
    B = zeros(N,1);
    B(1,1) = gate_potential; B(N,1) = gate_potential;
    for j = Inter1:Inter2
        if     (j == Inter1)  B(j,1) = Delta_x*Delta_x*q*Nacc/eps0*0.5;
        elseif (j == Inter2)  B(j,1) = Delta_x*Delta_x*q*Nacc/eps0*0.5;
        else                  B(j,1) = Delta_x*Delta_x*q*Nacc/eps0;
        end
    end
    phi = inv(A)*B;

    % Electron density %
    n = zeros(N,1);
    for j = Inter1:Inter2
        n(j,1) = ni*exp(q*phi(j,1)/(k_B*T));
    end
    
    % Re-calculation of potential by using calculated carrier density %
    B_re = zeros(N,1);
    B_re(1,1) = gate_potential; B_re(N,1) = gate_potential;
    for j = Inter1:Inter2
        if     (j == Inter1)  B_re(j,1) = Delta_x*Delta_x*q*(Nacc+n(j,1))/eps0*0.5;
        elseif (j == Inter2)  B_re(j,1) = Delta_x*Delta_x*q*(Nacc+n(j,1))/eps0*0.5;
        else                  B_re(j,1) = Delta_x*Delta_x*q*(Nacc+n(j,1))/eps0;
        end
    end
    
    phi_re = inv(A)*B_re;
    diff(:,i) = phi-phi_re;
    
    plot(x/1e-9,diff);
    xlabel('Distance (nm)');
    ylabel('Potential (V)');
    legend('0 V','0.1 V','0.2 V','0.3 V','0.4 V','0.5 V','0.6 V','0.7 V','0.8 V','0.9 V','1.0 V');
end

% Initial Condition (i = 1) %
% figure(1);
% plot(x/1e-9,phi(:,i));
% xlabel('Distance (nm)');
% ylabel('Potential (V)');
% figure(2);
% plot(x/1e-9,n(:,i)*1e-6);
% xlabel('Position (nm)');
% ylabel('Electron density (cm^-^3)'); 
% figure(3);
% plot(x/1e-9,phi_re(:,i));
% xlabel('Distance (nm)');
% ylabel('Potential (V)');
% figure(4);
% plot(x/1e-9,diff(:,i));
% xlabel('Distance (nm)');
% ylabel('Potential (V)');
