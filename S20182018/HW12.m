clear all;
q = 1.602192e-19; % Elementary charge, C
eps0 = 8.854187817e-12; % Vacuum permittivity, F/m
k_B = 1.380662e-23; % Boltzmann constant, J/K
T = 300.0; % Temperature, K
thermal = k_B*T/q; % Thermal voltage, V
eps_si = 11.7; eps_ox = 3.9; % Relative permittivity
ni = 1.075e16; % 1.075e10 /cm^3
count = zeros(6,1);
for jj=1:6
    
    if jj==1 % 600-nm-long structure, Spacing 0.5 nm
        N = 1201; Deltax = 0.5e-9; x_12 = 201; x_23 = 1001;
    elseif jj==2 % 600-nm-long structure, Spacing 1 nm
        N = 601; Deltax = 1e-9; x_12 = 101; x_23 = 501;
    elseif jj==3 % 600-nm-long structure, Spacing 10 nm
        N = 61; Deltax = 1e-8; x_12 = 11; x_23 = 51;
    elseif jj==4 % 120-nm-long structure, Spacing 0.2 nm
        N = 601; Deltax = 0.2e-9; x_12 = 201; x_23 = 401;
    elseif jj==5 % 120-nm-long structure, Spacing 1 nm
        N = 121; Deltax = 1e-9; x_12 = 41; x_23 = 81;
    elseif jj==6 % 120-nm-long structure, Spacing 10 nm
        N = 61; Deltax = 2e-9; x_12 = 21; x_23 = 41;
    end
    x = Deltax*transpose([0:N-1]); % real space, m
    Ndon = 2e21*ones(N,1); % 2e15 /cm^3
    Ndon(1:x_12,1) = 5e23; % 5e17 /cm^3
    Ndon(x_23:N,1) = 5e23; % 5e17 /cm^3
    coef = Deltax*Deltax*q/eps0;
    phi=zeros(N,1);
    
    phi(:,1) = thermal*log(Ndon(:,1)/ni);
    elec = ni*exp(phi/thermal);
    res = zeros(2*N,1);
    Jaco = sparse(2*N,2*N);
    err = 1;
    
    
    while err > 5e-17
        count(jj,1) = count(jj,1)+1;
        res = zeros(2*N,1);
        Jaco = sparse(2*N,2*N);
        res(1,1) = phi(1,1) - thermal*log(Ndon(1,1)/ni);
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
    elec_poisson = ni*exp(phi/thermal);
    for kk=1:N
        differ(kk,1) = (elec(kk,1)-elec_poisson(kk,1));
    end
    
    figure(jj)
    plot(x*1e+9,elec,'o'); hold on;
    plot(x*1e+9,elec_poisson,'r');
    xlabel('Position (nm)');
    ylabel('Electron density (cm^{-3})');
    hold off;
    figure(7)
    if jj<4
        for i=1:3
            plot(x*1e+9,differ(1:N,1)); hold on;
        end
    end
    xlabel('Position (nm)');
    ylabel('n_{Self-consistent} - n_{Nonlinear Poisson} (cm^{-3})');
    legend('0.5 nm spacing','1 nm spacing', '10 nm spacing');
    if jj>3
        figure(8)
        for i=4:6
            plot(x*1e+9,differ(1:N,1)); hold on;
        end
        xlabel('Position (nm)');
        ylabel('n_{Self-consistent} - n_{Nonlinear Poisson} (cm^{-3})');
           legend('0.2 nm spacing','1 nm spacing', '5 nm spacing');
    end
    
end


