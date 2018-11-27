close all;
clear all;

q = 1.602192e-19;
eps0 = 8.854187817e-12;
k_B = 1.380662e-23;
T = 300.0;
thermal = k_B*T/q;
eps_si = 11.7; eps_ox = 3.9;
ni = 1.075e16;

case_num = 6;

if case_num == 1 % For the long structure, spacing value is 0.5 nm
    N = 1201; Deltax = 0.5e-9; x_12 = 201; x_23 = 1001; 
elseif case_num == 2 % For the long structure, spacing value is 1 nm
    N = 601; Deltax = 1e-9; x_12 = 101; x_23 = 501;
elseif case_num == 3 % For the long structure, spacing value is 10 nm
    N = 61; Deltax = 1e-8; x_12 = 11; x_23 = 51;
elseif case_num == 4 % For the shor structure, spacing value is 0.2 nm
    N = 601; Deltax = 0.2e-9; x_12 = 201; x_23 = 401;
elseif case_num == 5 % For the shor structure, spacing value is 1nm
    N = 121; Deltax = 1e-9; x_12 = 41; x_23 = 81;
elseif case_num == 6 % For the shor structure, spacing value is 5 nm
    N = 25; Deltax = 5e-9; x_12 = 9; x_23 = 17;
end

x = Deltax*([0:N-1])';

if case_num <= 3
    Ndon = 2e21*ones(N,1);
    Ndon(1:x_12,1) = 5e23;
    Ndon(x_23:N,1) = 5e23;
else
    Ndon = 2e23*ones(N,1);
    Ndon(1:x_12,1) = 5e25;
    Ndon(x_23:N,1) = 5e25;
end

coef = Deltax*Deltax*q/eps0;

phi = zeros(N,1);
phi(:,1) = thermal*log(Ndon(:,1)/ni);

for newton = 1:10
    res = zeros(N,1);
    Jaco = sparse(N,N);
    res(1,1) = phi(1,1)-thermal*log(Ndon(N,1)/ni);
    Jaco(1,1) = 1;
    for  i = 2:N-1
        res(i,1) = eps_si*(phi(i+1,1)-2*phi(i,1)+phi(i-1,1));
        Jaco(i,i-1) = eps_si;
        Jaco(i,i) = -2*eps_si;
        Jaco(i,i+1) = eps_si;
    end
    res(N,1) - phi(N,1) - thermal*log(Ndon(N,1)/ni);
    Jaco(N,N) = 1;
    for i = 2:N-1
        res(i,1) = res(i,1) - coef*(-Ndon(i,1)+ni*exp(phi(i,1)/thermal));
        Jaco(i,i) = Jaco(i,i) - coef*ni*exp(phi(i,1)/thermal)/thermal;
    end
        
    update = Jaco\(-res);
    phi = phi+update;
end

n = zeros(N,1);
n = ni*exp(phi/thermal);
n_nonlinear = n;
subplot(2,1,1);
plot(x,n,'r');
hold on;

phi = zeros(N,1);
phi(:,1) = thermal*log(Ndon(:,1)/ni);
for newton=1:10
    res = zeros(2*N,1);
    Jaco = sparse(2*N,2*N);
    res(1,1) = phi(1,1) - thermal*log(Ndon(1,1)/ni);
    Jaco(1,1) = 1;
    
    for i=2:N-1
        res(2*i-1,1) = eps_si*(phi(i+1,1) - 2*phi(i,1) + phi(i-1,1)) + coef*(Ndon(i,1)-n(i,1));
        Jaco(2*i-1,2*i+1) = eps_si;
        Jaco(2*i-1,2*i-1) = -2*eps_si;
        Jaco(2*i-1,2*i-3) = eps_si;
        Jaco(2*i-1,2*i) = -coef;
    end
    
    res(2*N-1,1) = phi(N,1) - thermal*log(Ndon(N,1)/ni);
    Jaco(2*N-1,2*N-1) = 1.0;
    
    for i=1:N-1
        n_av = 0.5*(n(i+1,1)+n(i,1));
        dphidx = (phi(i+1,1) - phi(i,1))/Deltax;
        dndx = (n(i+1,1) - n(i,1))/Deltax;
        Jn = n_av * dphidx - thermal * dndx;
        res(2*i,1) = res(2*i,1) + Jn;
        Jaco(2*i,2*i+2) = Jaco(2*i,2*i+2) + 0.5*dphidx - thermal / Deltax;
        Jaco(2*i,2*i) = Jaco(2*i,2*i) + 0.5*dphidx + thermal / Deltax;
        Jaco(2*i,2*i+1) = Jaco(2*i,2*i+1) + n_av / Deltax;
        Jaco(2*i,2*i-1) = Jaco(2*i,2*i-1) - n_av / Deltax;
        res(2*i+2,1) = res(2*i+2,1) - Jn;
        Jaco(2*i+2,2*i+2) = Jaco(2*i+2,2*i+2) - 0.5*dphidx + thermal /Deltax;
        Jaco(2*i+2,2*i) = Jaco(2*i+2,2*i) - 0.5*dphidx - thermal / Deltax;
        Jaco(2*i+2,2*i+1) = Jaco(2*i+2,2*i+1) - n_av / Deltax;
        Jaco(2*i+2,2*i-1) = Jaco(2*i+2,2*i-1) + n_av / Deltax;
    end
    
    res(2,1) = n(1,1) - Ndon(1,1);
    Jaco(2,:) = 0;
    Jaco(2,2) = 1.0;
    res(2*N,1) = n(N,1) - Ndon(N,1);
    Jaco(2*N,:) = 0;
    Jaco(2*N,2*N) = 1.0;
    
    Cvector = zeros(2*N,1);
    Cvector(1:2:2*N-1,1) = thermal;
    Cvector(2:2:2*N,1) = max(abs(Ndon));
    Cmatrix = spdiags(Cvector,0,2*N,2*N);
    Jaco_scaled = Jaco * Cmatrix;
    Rvector = 1./sum(abs(Jaco_scaled),2);
    Rmatrix = spdiags(Rvector,0,2*N,2*N);
    Jaco_scaled = Rmatrix * Jaco_scaled;
    res_scaled = Rmatrix * res;
    update_scaled = Jaco_scaled \ (-res_scaled);
    update = Cmatrix * update_scaled;
        
    update = Jaco \ (-res);
    phi = phi + update(1:2:2*N-1,1);
    n = n + update(2:2:2*N,1);
    norm(update(1:2:2*N-1,1),inf)
end

n_self = n;
plot(x,n,'o');
legend('Nonlinear Poisson','Self-consistent');
xlabel('Distance [m]');
ylabel('Electron Density');

subplot(2,1,2);
error = (n_self-n_nonlinear)./n_self*100;
plot(x,error);
xlabel('Distance [m]');
ylabel('Error [%]')
title('Error of Nonlinear Poisson/Self-consistent Electron');
