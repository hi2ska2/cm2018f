%% clear paramenter %%
clear all
clc
close all

%% parameter %%
N = 25;
newton_max = 10;
Layer_1 = [40e-9 11.7 5e25];  %[length (m) relative permittivity]
Layer_2 = [40e-9 11.7 2e23];
Layer_3 = [40e-9 11.7 5e25];

Width = Layer_1(1, 1) + Layer_2(1, 1) + Layer_3(1, 1);
Deltax = Width/(N-1);
x = Deltax*transpose([0:N-1]);
eps_0 = 8.8541878176e-12;
q = 1.602192e-19;
ni = 1.075e16;
coef = Deltax*Deltax*q/eps_0;

K_B = 1.380662e-23;
T = 300;
thermal = K_B*T/q;


Ndon = zeros(N, 1);

for num = 1 : N
    if Deltax*(num-1) < Layer_1(1, 1)
        Ndon(num, 1) = Layer_1(1, 3);
    elseif Deltax*(num-1) > Layer_1(1, 1) + Layer_2(1, 1)
        Ndon(num, 1) = Layer_3(1, 3);
    else 
        Ndon(num, 1) = Layer_2(1, 3);
    end
end

elec_CP = zeros(N, 1);
phi_CP = zeros(N, 1);
error = zeros(newton_max, 1); 

%% continuity-Poisson solver %%
for newton = 1:newton_max
    

    res = zeros(2*N, 1);
    Jaco = sparse(2*N, 2*N);

    for ii = 2:N-1
        res(2*ii-1, 1) = Layer_1(1, 2)*(phi_CP(ii+1, 1) - 2*phi_CP(ii, 1) + phi_CP(ii-1, 1)) + coef*(Ndon(ii, 1) - elec_CP(ii, 1));

        Jaco(2*ii-1, 2*ii+1) = Layer_1(1, 2);
        Jaco(2*ii-1, 2*ii-1) = -2* Layer_1(1, 2);
        Jaco(2*ii-1, 2*ii-3) = Layer_1(1, 2);
        Jaco(2*ii-1, 2*ii  ) = -coef;
    end

    res(1, 1) = phi_CP(1,1) - thermal * log(Ndon(1, 1)/ni);
    Jaco(1, 1) = 1;
    res(2*N-1, 1) = phi_CP(N, 1) - thermal*log(Ndon(N, 1)/ni);
    Jaco(2*N-1, 2*N-1) = 1;

    for jj = 1:N-1
        n_av = 0.5*(elec_CP(jj+1, 1) + elec_CP(jj, 1));
        dphidx = (phi_CP(jj+1, 1) - phi_CP(jj, 1))/Deltax;
        delecdx = (elec_CP(jj+1, 1) - elec_CP(jj, 1))/Deltax;
        Jn = n_av * dphidx - thermal * delecdx;

        res(2*jj, 1) = res(2*jj, 1) + Jn;

        Jaco(2*jj, 2*jj+2) = Jaco(2*jj, 2*jj+2) + 0.5*dphidx - thermal/Deltax;
        Jaco(2*jj, 2*jj  ) = Jaco(2*jj, 2*jj  ) + 0.5*dphidx + thermal/Deltax;
        Jaco(2*jj, 2*jj+1) = Jaco(2*jj, 2*jj+1) + n_av/Deltax;
        Jaco(2*jj, 2*jj-1) = Jaco(2*jj, 2*jj-1) - n_av/Deltax;

        res(2*jj+2, 1) = res(2*jj+2, 1) - Jn;

        Jaco(2*jj+2, 2*jj+2) = Jaco(2*jj+2, 2*jj+2) - 0.5*dphidx + thermal/Deltax;
        Jaco(2*jj+2, 2*jj  ) = Jaco(2*jj+2, 2*jj  ) - 0.5*dphidx - thermal/Deltax;
        Jaco(2*jj+2, 2*jj+1) = Jaco(2*jj+2, 2*jj+1) - n_av/Deltax;
        Jaco(2*jj+2, 2*jj-1) = Jaco(2*jj+2, 2*jj-1) + n_av/Deltax;
    end

    res(2, 1) = elec_CP(1, 1) - Ndon(1, 1);
    Jaco(2, :) = 0;
    Jaco(2, 2) = 1;
    res(2*N, 1) = elec_CP(N, 1) - Ndon(N, 1);
    Jaco(2*N, :) = 0;
    Jaco(2*N, 2*N) = 1;
   
   Cvector = zeros(2*N, 1);
   Cvector(1:2:2*N-1, 1) = thermal;
   Cvector(2:2:2*N  , 1) = max(abs(Ndon));
   Cmatrix = spdiags(Cvector, 0, 2*N, 2*N);
   Jaco_scaled = Jaco*Cmatrix;
   Rvector = 1./sum(abs(Jaco_scaled), 2);
   Rmatrix = spdiags(Rvector, 0, 2*N, 2*N);
   Jaco_scaled = Rmatrix * Jaco_scaled;
   res_scaled = Rmatrix * res;
   update_scaled = Jaco_scaled\(-res_scaled);
   update = Cmatrix * update_scaled;
   
   phi_CP = phi_CP + update(1:2:2*N-1, 1);
   elec_CP = elec_CP + update(2:2:2*N, 1);
   norm(update(1:2:2*N-1, 1), inf);
   
   
   error(newton, 1) = abs(sum(res(1:2:2*N-1, 1)))/N;
end

%% Non-linear Poisson solver %%

phi_NP = zeros(N, 1);
phi_NP(:, 1) = thermal*log(Ndon(:,1)/ni);

for newton = 1:newton_max
    

    res = zeros(N, 1);
    Jaco = sparse(N, N);

    for ii = 2:N-1
        res(ii, 1) = Layer_1(1, 2)*(phi_NP(ii+1, 1) - 2*phi_NP(ii, 1) + phi_NP(ii-1, 1));

        Jaco(ii, ii-1) = Layer_1(1, 2);
        Jaco(ii, ii  ) = -2* Layer_1(1, 2);
        Jaco(ii, ii+1) = Layer_1(1, 2);
    
    end

    res(1, 1) = phi_NP(1,1) - thermal * log(Ndon(1, 1)/ni);
    Jaco(1, 1) = 1;
    res(N, 1) = phi_NP(N, 1) - thermal*log(Ndon(N, 1)/ni);
    Jaco(N, N) = 1;
    
    for ii = 2:N-1
        res(ii, 1) = res(ii, 1) - coef*(-Ndon(ii, 1) + ni*exp(phi_NP(ii, 1)/thermal));
        Jaco(ii, ii) = Jaco(ii, ii) - coef*ni*exp(phi_NP(ii, 1)/thermal)/thermal;
    end
    
   update = Jaco\(-res);
   phi_NP = phi_NP + update;
   norm(update, inf);
end

elec_NP = ni*exp(phi_NP/thermal);

plot(x,elec_CP,'o');
hold on;
plot(x,elec_NP,'r');