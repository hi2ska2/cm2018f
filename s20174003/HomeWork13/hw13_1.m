close all;
clear all;

Structure_sel = 1; % 0 is long, 1 is short structure.
q = 1.602192e-19;
k_B = 1.380662e-23;
eps0 = 8.854187817e-12; % Vaccum permittivity (F/m)
T = 300.0;
thermal = k_B*T/q;
Dn = 0.1;

%Deltax = 1e-9; % 1nm spacing


if(Structure_sel==0)%% Long structure
    N = 601;
    Deltax = 1e-9;
    x_12 = 101; %At x=100nm
    x_23 = 501; %At x=500nm
else%%Short structure
    N = 121;
    Deltax =1e-9;
    xs_12 = 41; %At x=40nm
    xs_23 = 81; %At x=80nm
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



elec = zeros(N,1);
elec = ni*exp(phi/thermal);
J = zeros(N-1,1);
sub_iter = 0;
phi = zeros(N,1);
phi(:,1) = thermal*log(Ndon(:,1)/ni);
for bias=0:10
    V_applied = 0.05*bias;
    sub_iter = sub_iter +1;
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
        res(2*N-1,1) = phi(N,1) - thermal*log(Ndon(N,1)/ni) - V_applied;
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


    for iteriter=1:N-1
        Bernoulli_x = (phi(iteriter+1,1)-phi(iteriter,1))/thermal;
        J(iteriter,1) = ((-q*Dn)/(Deltax))*(elec(iteriter+1,1)*(Bernoulli_x/(exp(Bernoulli_x)-1))-elec(iteriter,1)*(-Bernoulli_x/(exp(-Bernoulli_x)-1))); 
    end
    Current(sub_iter,1) = J(end,1);
end

V_sweep = 0:0.05:0.5;


figure(1)

semilogy(V_sweep,Current);
xlabel('Drain Voltage');
ylabel('Current Density at Drain[A/m^2]')
grid on;
