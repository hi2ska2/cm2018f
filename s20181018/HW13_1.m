close all;
clear all;

q = 1.602192e-19;
eps0 = 8.854187817e-12;
eps_si = 11.7; eps_ox = 3.9;
k_B = 1.380662e-23;
T = 300;
thermal = k_B*T/q;
mu = 0.14;
D_n = thermal*mu;
ni = 1.075e16;

case_num = 0;

if(case_num == 0) %% Long structure
    Deltax = 1e-9;
    N = 601;
    x_12 = 101;
    x_23 = 501;
else %% Short structure
    Deltax = 0.1e-9;
    N = 1201;
    xs_12 = 401;
    xs_23 = 801;
end

x = Deltax*([0:N-1])';
coef = Deltax*Deltax*q/eps0;

if(case_num == 0)
    Ndon = 2e21*ones(N,1);
    Ndon(1:x_12,1) = 5e23;
    Ndon(x_23:N,1) = 5e23;
else
    Ndon = 2e23*ones(N,1);
    Ndon(1:xs_12,1) = 5e25;
    Ndon(xs_23:N,1) = 5e25;
end
 
N_bias = 20;
phi_bias = zeros(N,N_bias);
n_bias = zeros(N,N_bias);

for bias = 0:N_bias;
    
    phi=zeros(N,1);
    phi(:,1)=thermal*log(Ndon(:,1)/ni);
    n=zeros(N,1);
    n=ni*exp(phi/thermal);
    V_applied = 0.02*bias;
    
    for newton = 1:10
        res = zeros(2*N,1);
        Jaco = sparse(2*N,2*N);
        res(1,1) = phi(1,1) - thermal*log(Ndon(1,1)/ni);
        Jaco(1,1) = 1.0;
        
        for i = 2:N-1
            res(2*i-1,1) = eps_si*(phi(i+1,1)-2*phi(i,1)+phi(i-1,1))+coef*(Ndon(i,1)-n(i,1));
            Jaco(2*i-1,2*i +1) = eps_si;
            Jaco(2*i-1,2*i -1) = -2*eps_si;
            Jaco(2*i-1,2*i -3) = eps_si;
            Jaco(2*i-1,2*i ) = -coef;
        end
        
        res(2*N-1,1) = phi(N,1)-thermal*log(Ndon(N,1)/ni)-V_applied;
        Jaco(2*N-1, 2*N-1) = 1.0;
        
        for i = 1:N-1
            n_av = 0.5*(n(i+1,1)+n(i,1));
            dphidx = (phi(i+1,1)-phi(i,1))/Deltax;
            delecdx = (n(i+1,1)-n(i,1))/Deltax;
            Jn = n_av*dphidx - thermal*delecdx;
            res(2*i,1) = res(2*i,1)+Jn;
            Jaco(2*i,2*i+2) = Jaco(2*i,2*i+2)+0.5*dphidx-thermal/Deltax;
            Jaco(2*i,2*i) = Jaco(2*i,2*i)+0.5*dphidx+thermal/Deltax;
            Jaco(2*i,2*i+1) = Jaco(2*i,2*i+1)+n_av/Deltax;
            Jaco(2*i,2*i-1) = Jaco(2*i,2*i-1)-n_av/Deltax;
            res(2*i+2,1) = res(2*i+2,1) - Jn;
            Jaco(2*i +2,2*i+2) = Jaco(2*i+2,2*i+2)-0.5*dphidx+thermal/Deltax;
            Jaco(2*i +2,2*i) = Jaco(2*i+2,2*i)-0.5*dphidx-thermal/Deltax;
            Jaco(2*i +2,2*i+1) = Jaco(2*i+2,2*i+1)-n_av/Deltax;
            Jaco(2*i +2,2*i-1) = Jaco(2*i+2,2*i-1)+n_av/Deltax;
        end
        
        res(2,1) = n(1,1) - Ndon(1,1);
        Jaco(2,:) = 0.0;
        Jaco(2,2) = 1.0;
        res(2*N,1) = n(N,1) - Ndon(N,1);
        Jaco(2*N,:) = 0.0;
        Jaco(2*N,2*N) = 1.0;
        
        Cvector = zeros(2*N,1);
        Cvector(1:2:2*N-1,1) = thermal;
        Cvector(2:2:2*N,1) = max(abs(Ndon));
        Cmatrix = spdiags(Cvector,0,2*N,2*N);
        Jaco_scaled = Jaco*Cmatrix;
        Rvector = 1./sum(abs(Jaco_scaled),2);
        Rmatrix = spdiags(Rvector,0,2*N,2*N);
        Jaco_scaled = Rmatrix*Jaco_scaled;
        res_scaled = Rmatrix*res;
        update_scaled = Jaco_scaled\(-res_scaled);
        update = Cmatrix*update_scaled;
        phi = phi+update(1:2:2*N-1,1);
        n = n+update(2:2:2*N,1);
    end
    phi_bias(:,bias+1) = phi;
    n_bias(:,bias+1) = n;
end

J = zeros(N-1,N_bias);

for j = 0:N_bias
    phi_J = phi_bias(:,j+1);
    n_J = n_bias(:,j+1);
    
    for i = 2:N
        x_1 = (phi_J(i) - phi_J(i-1) )./ thermal ;
        J(i,j+1) = -(q*D_n/Deltax)*( n_J(i,1)*(x_1/(exp(x_1)-1)) - n_J(i-1,1)* (-x_1/(exp(-x_1) -1)) );
    end
end

bias_tot = 0:20;
bias_tot = 0.1*bias_tot;

figure(1);
plot(bias_tot,J(N-1,:));
xlabel('Drain Voltage [V]'); ylabel('Current Density [A/m^2]');
