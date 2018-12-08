clear all;


q = 1.602192e-19; % Elementary charge, C
eps0 = 8.854187817e-12; % Vacuum permittivity, F/m
k_B = 1.380662e-23; % Boltzmann constant, J/K
T = 300.0; % Temperature, K
thermal = k_B*T/q; % Thermal voltage, V
eps_si = 11.7; eps_ox = 3.9; % Relative permittivity
Dn = 0.01; %Diffusion coefficient
ni = 1.075e16; % 1.075e10 /cm^3

Jsum = zeros(26,2);
for jj=1:2
    
    if jj==1 % 600-nm-long structure, Spacing 0.5 nm
        N = 601; Deltax = 1.0e-9; x_12 = 101; x_23 = 501;
    elseif jj==2 % 120-nm-long structure, Spacing 0.2 nm
        N = 601; Deltax = 0.2e-9; x_12 = 201; x_23 = 401;
    end
    x = Deltax*transpose([0:N-1]); % real space, m
    if jj==1
        Ndon = 2e21*ones(N,1); % 2e15 /cm^3
        Ndon(1:x_12,1) = 5e23; % 5e17 /cm^3
        Ndon(x_23:N,1) = 5e23; % 5e17 /cm^3
    elseif jj==2
        Ndon = 2e23*ones(N,1); % 2e17 /cm^3
        Ndon(1:x_12,1) = 5e25; % 5e19 /cm^3
        Ndon(x_23:N,1) = 5e25; % 5e19 /cm^3
    end
    coef = Deltax*Deltax*q/eps0;
    phi=zeros(N,1);
    
    
    phi(:,1) = thermal*log(Ndon(:,1)/ni);
    elec = ni*exp(phi/thermal);
    res = zeros(2*N,1);
    Jaco = sparse(2*N,2*N);
    err = 1;
    
    for bias=0:25
        J = zeros(N,1);
        V_applied = 0.02*bias;  % Applied bias voltage
        V(bias+1,1) = V_applied;
        for newton=1:10
            res = zeros(2*N,1);
            Jaco = sparse(2*N,2*N);
            res(1,1) = phi(1,1) - thermal*log(Ndon(1,1)/ni) - V_applied; % Applied bias voltage
            Jaco(1,1) = 1.0;
            for ii=2:N-1
                res(2*ii-1,1) = eps_si*(phi(ii+1,1)-2*phi(ii,1)+phi(ii-1,1)) + coef*(Ndon(ii,1)-elec(ii,1));
                Jaco(2*ii-1,2*ii+1) = eps_si;
                Jaco(2*ii-1,2*ii-1) = -2*eps_si;
                Jaco(2*ii-1,2*ii-3) = eps_si;
                Jaco(2*ii-1,2*ii ) = -coef;
            end
            res(2*N-1,1) = phi(N,1) - thermal*log(Ndon(N,1)/ni);
            Jaco(2*N-1,2*N-1) = 1.0;
            
            
            res(2,1) = elec(1,1) - Ndon(1,1);
            Jaco(2,:) = 0.0;
            Jaco(2,2) = 1.0;
            res(2*N,1) = elec(N,1) - Ndon(N,1);
            Jaco(2*N,:) = 0.0;
            Jaco(2*N,2*N) = 1.0;
            
            for ii=1:N-1 % edge-wise construction
                n_av = 0.5*(elec(ii+1,1)+elec(ii,1));
                dphidx = (phi(ii+1,1)-phi(ii,1))/Deltax;
                delecdx = (elec(ii+1,1)-elec(ii,1))/Deltax;
                Jn = n_av * dphidx - thermal * delecdx;
                res(2*ii,1) = res(2*ii,1) + Jn;
                Jaco(2*ii,2*ii+2) = Jaco(2*ii,2*ii+2) + 0.5*dphidx - thermal / Deltax;
                Jaco(2*ii,2*ii ) = Jaco(2*ii,2*ii ) + 0.5*dphidx + thermal / Deltax;
                Jaco(2*ii,2*ii+1) = Jaco(2*ii,2*ii+1) + n_av / Deltax;
                Jaco(2*ii,2*ii-1) = Jaco(2*ii,2*ii-1) - n_av / Deltax;
                res(2*ii+2,1) = res(2*ii+2,1) - Jn;
                Jaco(2*ii+2,2*ii+2) = Jaco(2*ii+2,2*ii+2) - 0.5*dphidx + thermal / Deltax;
                Jaco(2*ii+2,2*ii ) = Jaco(2*ii+2,2*ii ) - 0.5*dphidx - thermal / Deltax;
                Jaco(2*ii+2,2*ii+1) = Jaco(2*ii+2,2*ii+1) - n_av / Deltax;
                Jaco(2*ii+2,2*ii-1) = Jaco(2*ii+2,2*ii-1) + n_av / Deltax;
            end
            
            Cvector = zeros(2*N,1);
            Cvector(1:2:2*N-1,1) = thermal;
            Cvector(2:2:2*N  ,1) = max(abs(Ndon));
            Cmatrix = spdiags(Cvector, 0, 2*N,2*N);
            Jaco_scaled = Jaco * Cmatrix;
            Rvector = 1./sum(abs(Jaco_scaled),2);
            Rmatrix = spdiags(Rvector,0,2*N,2*N);
            Jaco_scaled = Rmatrix * Jaco_scaled;
            res_scaled = Rmatrix * res;
            update_scaled = Jaco_scaled \ (-res_scaled);
            update = Cmatrix * update_scaled;
            
            phi = phi + update(1:2:2*N-1,1);
            elec = elec + update(2:2:2*N,1);
            err = norm(update(1:2:2*N-1,1),inf);
        end
        for kk=1:N-1
            t = (phi(kk+1,1)-phi(kk,1))/thermal;
            J(kk,1) = (q*Dn/Deltax)*(elec(kk+1,1)*t/(exp(t)-1) + elec(kk,1)*t/(exp(-t)-1));
            Jsum(bias+1,jj) = Jsum(bias+1,jj)+J(kk,1);
        end
    end
end
figure (1)
plot(V(:,1),-Jsum(:,1),'o');
xlabel('V_{\itD} (V)'); ylabel('Current density (A/m^2)');
figure (2)
plot(V(:,1),-Jsum(:,2),'o');
xlabel('V_{\itD} (V)'); ylabel('Current density (A/m^2)');