close all;
clear all;


q = 1.602192e-19; % Elementary charge, C
eps0 = 8.854187817e-12; % Vacuum permittivity, F/m
k_B = 1.380662e-23; % Boltzmann constant, J/K
T = 300.0; % Temperature, K
thermal = k_B*T/q; % Thermal voltage, V
sel_size = 1; % 0 => long & 1 => short
sel_spacing = 2; % 0 => 0.5 nm or 0.2  & 1 => 1 nm & 2 => 10 nm or 5 nm
N_large = 600e-9; % [m] 600 nm
N_short = 120e-9; % [m] 120 nm
% Selection the specfic conditions

if (sel_size == 0) % long sturcture
    if (sel_spacing == 0) % 0.5 nm
        Deltax = 0.5e-9; %  spacing
        N = N_large/Deltax +1;
        N = round(N);
        x_12 = 201; % At x=200 nm
        x_23 = 1001; % At x=1000 nm
    elseif (sel_spacing == 1)
        Deltax = 1e-9; %  spacing
        N = N_large/Deltax +1;
        N = round(N);
        x_12 = 101; % At x=100 nm
        x_23 = 501; % At x=500 nm
    elseif (sel_spacing == 2)
        Deltax = 10e-9; %  spacing
        N = N_large/Deltax +1;
        N = round(N);
        x_12 = 11; % At x=10 nm
        x_23 = 51; % At x=50 nm
    end
        Ndon = 2e21*ones(N,1); % 2e15 /cm^3
        Ndon(1:x_12,1) = 5e23; % 5e17 /cm^3
        Ndon(x_23:N,1) = 5e23; % 5e17 /cm^3
elseif (sel_size == 1) % short sturcture
    if (sel_spacing == 0) % 0.5 nm
        Deltax = 0.2e-9; %  spacing
        N = N_short/Deltax +1;
        N = round(N);
        x_12 = 201; % At x=40 nm
        x_23 = 401; % At x=80 nm
    elseif (sel_spacing == 1)
        Deltax = 1e-9; %  spacing
        N = N_short/Deltax +1;
        N = round(N);
        x_12 = 41; % At x=40 nm
        x_23 = 81; % At x=80 nm
    elseif (sel_spacing == 2)
        Deltax = 5e-9; %  spacing
        N = N_short/Deltax +1;
        N = round(N);
        x_12 = 9; % At x=40 nm
        x_23 = 17; % At x=80 nm
    end
        Ndon = 2e23*ones(N,1); % 2e15 /cm^3
        Ndon(1:x_12,1) = 5e25; % 5e17 /cm^3
        Ndon(x_23:N,1) = 5e25; % 5e17 /cm^3
end
        

x = Deltax*transpose([0:N-1]); % real space, m

eps_si = 11.7; eps_ox = 3.9; % Relative permittivity
% Ndon = 2e21*ones(N,1); % 2e15 /cm^3
% Ndon(1:x_12,1) = 5e23; % 5e17 /cm^3
% Ndon(x_23:N,1) = 5e23; % 5e17 /cm^3



ni = 1.075e16; % 1.075e10 /cm^3
coef = Deltax*Deltax*q/eps0;

phi = zeros(N,1);
phi(:,1) = thermal*log(Ndon(:,1)/ni);


for newton = 1:10
    %%% Jaco and res are constructed here
    res = zeros(N,1);
    Jaco = sparse(N,N);
    res(1,1) = phi(1,1) - thermal*log(Ndon(1,1)/ni);
    Jaco(1,1) = 1.0;
    
    for ii = 2:N-1
        res(ii,1) = eps_si*(phi(ii+1,1) -2*phi(ii,1) + phi(ii-1,1));
        Jaco(ii,ii-1) = eps_si;
        Jaco(ii,ii) = -2*eps_si;
        Jaco(ii,ii+1) = eps_si;
    end

    res(N,1) = phi(N,1) - thermal*log(Ndon(N,1)/ni);
    Jaco(N,N) = 1.0;

    for ii = 2:N-1
        res(ii,1) = res(ii,1) + coef*(-ni*exp(phi(ii,1)/thermal) + Ndon(ii,1));
        Jaco(ii,ii) = Jaco(ii,ii) - coef*ni*exp(phi(ii,1)/thermal)/thermal;
    end

    update = Jaco\(-res);
    phi = phi + update;
%     norm(update,inf);
end


% figure;
% plot(x,phi);


% figure;
elec = zeros(N,1);
elec = ni*exp(phi/thermal);
% plot(x,elec,'r');
% hold on

phi_poi = phi;
elec_poi = elec;

%% coupled equation with continuity equation and Nonlinear poisson's equation

phi = zeros(N,1);
phi(:,1) = thermal*log(Ndon(:,1)/ni);

for newton = 1:10
    res = zeros(2*N,1);
    Jaco = sparse(2*N,2*N);
    res(1,1) = phi(1,1) - thermal*log(Ndon(1,1)/ni);
    Jaco(1,1) = 1.0;
    
    % Nonlinear Poisson's eqn
    
    for ii = 2:N-1
        res(2*ii-1,1) = eps_si*(phi(ii+1,1) - 2*phi(ii,1) + phi(ii-1,1)) + coef*(Ndon(ii,1) - elec(ii,1));
        Jaco(2*ii -1,2*ii +1) = eps_si;
        Jaco(2*ii -1,2*ii -1) = -2*eps_si;
        Jaco(2*ii -1,2*ii -3) = eps_si;
        Jaco(2*ii -1,2*ii ) = -coef;
    end
    
    res(2*N -1,1) = phi(N,1) - thermal*log(Ndon(N,1)/ni);
    Jaco(2*N -1, 2*N -1) = 1.0;
    
    % The continuity equation
    
    for ii = 1:N-1 % edge-wise construction
        n_av = 0.5*(elec(ii+1,1) + elec(ii,1));
        dphidx = (phi(ii+1,1) - phi(ii,1))/Deltax;
        delecdx = (elec(ii+1,1) - elec(ii,1))/Deltax;
        Jn = n_av*dphidx - thermal*delecdx;
        res(2*ii,1) = res(2*ii,1) + Jn;
        Jaco(2*ii,2*ii+2) = Jaco(2*ii,2*ii +2) + 0.5*dphidx - thermal/Deltax;
        Jaco(2*ii,2*ii) = Jaco(2*ii,2*ii ) + 0.5*dphidx + thermal/Deltax;
        Jaco(2*ii,2*ii+1) = Jaco(2*ii,2*ii +1) + n_av/Deltax;
        Jaco(2*ii,2*ii-1) = Jaco(2*ii,2*ii -1) - n_av/Deltax;
        res(2*ii+2,1) = res(2*ii+2,1) - Jn;
        Jaco(2*ii +2,2*ii+2) = Jaco(2*ii+2,2*ii +2) - 0.5*dphidx + thermal/Deltax;
        Jaco(2*ii +2,2*ii) = Jaco(2*ii+2,2*ii ) - 0.5*dphidx - thermal/Deltax;
        Jaco(2*ii +2,2*ii+1) = Jaco(2*ii+2,2*ii +1) - n_av/Deltax;
        Jaco(2*ii +2,2*ii-1) = Jaco(2*ii+2,2*ii -1) + n_av/Deltax;
    end
    
    % BCs
    
    res(2,1) = elec(1,1) - Ndon(1,1);
    Jaco(2,:) = 0.0;
    Jaco(2,2) = 1.0;
    res(2*N,1) = elec(N,1) - Ndon(N,1);
    Jaco(2*N,:) = 0.0;
    Jaco(2*N,2*N) = 1.0;
    
    % Scaling : Set R and C
    
    Cvector = zeros(2*N,1);
    Cvector(1:2:2*N-1,1) = thermal;
    Cvector(2:2:2*N,1) = max(abs(Ndon));
    Cmatrix = spdiags(Cvector,0,2*N,2*N);
    Jaco_scaled = Jaco * Cmatrix;
    Rvector = 1./sum(abs(Jaco_scaled),2);
    Rmatrix = spdiags(Rvector,0,2*N,2*N);
    Jaco_scaled = Rmatrix * Jaco_scaled;
    res_scaled = Rmatrix*res;
    update_scaled = Jaco_scaled \ (-res_scaled);
    update = Cmatrix * update_scaled;
    
    phi = phi + update(1:2:2*N-1,1);
    elec = elec + update(2:2:2*N,1);
end

% plot(x,elec,'o');

figure;
semilogy(x,elec_poi)
hold
semilogy(x,elec,'o')

legend('Nonlinear Poisson','Self-consistent')
title('Short structure')
xlabel('X-position [m]')
ylabel('Electronc density  [cm^-2]')
