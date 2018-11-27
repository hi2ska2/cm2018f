close all;
clear all;

Structure_sel = 0; % 0 is long, 1 is short structure.
Space_sel = 2; %1 is 0.5/0.2nm, 2 is 1/1nm, 3is 10nm/5nm;
q = 1.602192e-19;
k_B = 1.380662e-23;
eps0 = 8.854187817e-12; % Vaccum permittivity (F/m)
T = 300.0;
thermal = k_B*T/q;

%Deltax = 1e-9; % 1nm spacing


if(Structure_sel==0)%% Long structure
    if(Space_sel==1)
        N = 1201;
        Deltax = 0.5e-9;
        x_12 = 201; %At x=100nm
        x_23 = 1001; %At x=500nm
    elseif(Space_sel==2)
        N = 601;
        Deltax = 1e-9;
        x_12 = 101; %At x=100nm
        x_23 = 501; %At x=500nm
    elseif(Space_sel==3)
        N = 61;
        Deltax = 10e-9;
        x_12 = 11; %At x=100nm
        x_23 = 51; %At x=500nm
    end
   
else%%Short structure
    if(Space_sel==1)
        N = 601;
        Deltax = 0.2e-9;
        xs_12 = 201; %At x=40nm
        xs_23 = 401; %At x=80nm
    elseif(Space_sel==2)
        N = 121;
        Deltax = 1e-9;
        xs_12 = 41; %At x=40nm
        xs_23 = 81; %At x=80nm
    elseif(Space_sel==3)
        N = 13;
        Deltax = 5e-9;
        xs_12 = 5; %At x=40nm
        xs_23 = 9; %At x=80nm
    end
end
x = Deltax*transpose([0:N-1]); %real space (m)
eps_si = 11.7;
eps_ox = 3.9;

if(Structure_sel==0)
    Ndon = 2e21*ones(N,1); %2e15 /cm^3;
    Ndon(1:x_12,1) = 5e23; %5e17 /cm^3;
    Ndon(x_23:N,1) = 5e23; %5e17 /cm^3;
else
    Ndon = 2e23*ones(N,1); %2e17 /cm^3;
    Ndon(1:xs_12,1) = 5e25; %5e19 /cm^3;
    Ndon(xs_23:N,1) = 5e25; %5e19 /cm^3;
end
ni = 1.075e16; %1.075e10 /cm^3
coef = Deltax*Deltax*q/eps0;


phi = zeros(N,1);
phi(:,1) = thermal*log(Ndon(:,1)/ni);
for newton=1:10
    res = zeros(N,1);
    Jaco = sparse(N,N);
    res(1,1) = phi(1,1) - thermal*log(Ndon(1,1)/ni);
    Jaco(1,1) = 1.0;
    for iter=2:N-1
        res(iter,1) = eps_si*(phi(iter+1,1) - 2*phi(iter,1) + phi(iter-1,1));
        Jaco(iter,iter-1) = eps_si;
        Jaco(iter,iter) = -2*eps_si;
        Jaco(iter,iter+1) = eps_si;
    end
    res(N,1) = phi(N,1) - thermal*log(Ndon(N,1)/ni);
    Jaco(N,N) = 1.0;
    for iter=2:N-1
        res(iter,1) = res(iter,1) - coef*(-Ndon(iter,1) + ni*exp(phi(iter,1)/thermal));
        Jaco(iter,iter) = Jaco(iter,iter) - coef*ni*exp(phi(iter,1)/thermal)/thermal;
    end
    update = Jaco \ (-res);
    phi = phi + update;
    %norm(update,inf)
end
figure(1)
plot(x,phi)

elec = zeros(N,1);
elec = ni*exp(phi/thermal);
figure(2)
%plot(x,elec,'r');
semilogy(x,elec,'r');
hold on;



phi = zeros(N,1);
phi(:,1) = thermal*log(Ndon(:,1)/ni);
for newton=1:10
    res = zeros(2*N,1);
    Jaco = sparse(2*N,2*N);
    res(1,1) = phi(1,1) - thermal*log(Ndon(1,1)/ni);
    Jaco(1,1) = 1.0;
    for iter=2:N-1
        res(2*iter-1,1) = eps_si*(phi(iter+1,1) - 2*phi(iter,1) + phi(iter-1,1)) + coef*(Ndon(iter,1)-elec(iter,1));
        Jaco(2*iter-1,2*iter+1) = eps_si;
        Jaco(2*iter-1,2*iter-1) = -2*eps_si;
        Jaco(2*iter-1,2*iter-3) = eps_si;
        Jaco(2*iter-1,2*iter) = -coef;
    end
    res(2*N-1,1) = phi(N,1) - thermal*log(Ndon(N,1)/ni);
    Jaco(2*N-1,2*N-1) = 1.0;
    
    for iter=1:N-1
        n_av = 0.5*(elec(iter+1,1)+elec(iter,1));
        dphidx = (phi(iter+1,1) - phi(iter,1))/Deltax;
        delecdx = (elec(iter+1,1) - elec(iter,1))/Deltax;
        Jn = n_av * dphidx - thermal * delecdx;
        res(2*iter,1) = res(2*iter,1) + Jn;
        Jaco(2*iter,2*iter+2) = Jaco(2*iter,2*iter+2) + 0.5*dphidx - thermal / Deltax;
        Jaco(2*iter,2*iter) = Jaco(2*iter,2*iter) + 0.5*dphidx + thermal / Deltax;
        Jaco(2*iter,2*iter+1) = Jaco(2*iter,2*iter+1) + n_av / Deltax;
        Jaco(2*iter,2*iter-1) = Jaco(2*iter,2*iter-1) - n_av / Deltax;
        res(2*iter+2,1) = res(2*iter+2,1) - Jn;
        Jaco(2*iter+2,2*iter+2) = Jaco(2*iter+2,2*iter+2) - 0.5*dphidx + thermal /Deltax;
        Jaco(2*iter+2,2*iter) = Jaco(2*iter+2,2*iter) - 0.5*dphidx - thermal / Deltax;
        Jaco(2*iter+2,2*iter+1) = Jaco(2*iter+2,2*iter+1) - n_av / Deltax;
        Jaco(2*iter+2,2*iter-1) = Jaco(2*iter+2,2*iter-1) + n_av / Deltax;
    end
    
    res(2,1) = elec(1,1) - Ndon(1,1);
    Jaco(2,:) = 0;
    Jaco(2,2) = 1.0;
    res(2*N,1) = elec(N,1) - Ndon(N,1);
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
    update =Cmatrix * update_scaled;
    
    
    
    %update = Jaco \ (-res);
    phi = phi + update(1:2:2*N-1,1);
    elec = elec + update(2:2:2*N,1);
    %norm(update(1:2:2*N-1,1),inf)
end





%plot(x,elec,'o');
semilogy(x,elec,'o');

legend('Nonlinear Poisson','Self-consistent')
title('Electron Density')
xlabel('Distance [m]')
%figure(3)
%semilogy(x,elec)
