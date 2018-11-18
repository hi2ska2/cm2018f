
clear all;
clc;

% HW 11

q = 1.602192e-19; % Elementary charge, C
k_B = 1.380662e-23; % Boltzmann constant, J/K
T = 300.0; % Temperature, K
tau = 1e-12; % Relaxation time
m_e = 9.109534e-31; % Electron mass, Kg
dx = 0.1e-9;

N_x = 301; % The number of X position grid
N_H = 300; % The number of H level grid

% H = 0.1; % (eV) % Test
H = (linspace(1,100,N_H)*0.01); % 0.01 ~ 0.1 eV H transformation
VD = 0.01; % (V) Drain voltage

X_posi = (linspace(0,300,N_x)*dx); % X position [nm]

interface1 = 101;
interface2 = 201;

f0 = zeros(N_x,N_H);
f1 = zeros(N_x,N_H);

for jj = 1:N_H
    
    fs = sqrt(2*pi)/(1 + exp(q*H(jj)/(k_B*T))); % At source
    fd = sqrt(2*pi)/(1 + exp(q*(H(jj)+VD)/(k_B*T))); % At drain

    V = zeros(N_x,1);
    V(interface1:interface2,1) = [0:1/(interface2-interface1):1]*VD;
    V(interface2:N_x,1) = VD;

    A = zeros(N_x,N_x);
    A(1,1) = 1.0;
    
    for ii=2:N_x-1
        c1 = H(jj) + 0.5*(V(ii,1)+V(ii-1,1));
        c2 = H(jj) + 0.5*(V(ii+1,1)+V(ii,1));
        A(ii,ii-1) = c1; A(ii,ii) = -c1-c2; A(ii,ii+1) = c2;
    end
    A(N_x,N_x) = 1.0;

    b = zeros(N_x,1);
    b(1,1) = fs;
    b(N_x,1) = fd;

    f0(:,jj) = A \ b;
    
    for i = 1:N_x
        df0 = gradient(f0(:,jj))./dx;
        f1(i,jj) = -tau*sqrt(1/pi)*sqrt(2.0*q*(H(jj) + V(i))./m_e)*df0(i);
    end
    
end

figure
mesh(H,X_posi,f0)
xlabel('H [eV]');
ylabel('Position [nm]');
zlabel('f0')
title('f0 Result')
colorbar

figure
mesh(H,X_posi,f1)
xlabel('H [eV]');
ylabel('Position [nm]');
zlabel('f1')
title('f1 Result')
colorbar

