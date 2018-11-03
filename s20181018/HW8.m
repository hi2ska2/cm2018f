close all;
clear all;
 
h = 6.626176e-34; % Planck constant (J s)
hbar = h/(2*pi);
q = 1.602192e-19; % Elementary charge (C)
m0 = 9.109534e-31;
eps0 = 8.854187817e-12; % Vacuum permittivity, (F/m)
k_B = 1.380662e-23; % Boltzmann constant (J/K) 
T = 300; % Temperature (K)
thermal = k_B*T/q;
Lx = 100e-9; Ly = 100e-9; Lz = 5e-9; % Lengths, m
mxx = 0.19; myy = 0.19; mzz = 0.91; % Masses,m0
Delta_z = 0.1e-9; % 0.1 nm spacing
N = 61; % 6 nm thick
z = Delta_z*transpose([0:N-1]); % real space (m)
Inter1 = 6; % At x = 0.5 nm
Inter2 = 56; % At x = 5.5 nm
eps_si = 11.7; % Silicon relative permittivity
eps_ox = 3.9; % Silicon dioxide relative permittivity
Nacc = 1e24; % 1e18/cm^3
ni = 1.075e16; % 1.075e10/cm^-4
coef_Poi = Delta_z*Delta_z*q/eps0;
coef_Sch = 2*Lx*Ly/(2*pi)*sqrt(mxx*myy)*m0/(hbar^2)*(k_B*T);
Vg_sweep = transpose(linspace(0,1,21));
Ec_Ei = 0.561004; % E_c-E_i (eV)
 
Vg = Vg_sweep+0.33374;
n_sweep_Poi = zeros(N,length(Vg));
n_sweep_Sch = zeros(N,length(Vg));
 
for i = 1:length(Vg)
    
    phi = zeros(N,1);
    phi(:,1) = Vg(i,1);
    res = zeros(N,1);
    Jaco = sparse(N,N);
    Jaco(1,1) = 1;
    Jaco(N,N) = 1;
    res(1,1) = phi(1,1)-Vg(i,1);
    res(N,1) = phi(N,1)-Vg(i,1);
    
    for newton = 1:100
        % Laplacian %
        for ii = 2:N-1
            if      (ii<Inter1 || ii>Inter2)
                res(ii,1) = eps_ox*phi(ii+1,1)-2*eps_ox*phi(ii,1)+eps_ox*phi(ii-1,1);
                Jaco(ii,ii-1) = eps_ox; Jaco(ii,ii) = -2*eps_ox; Jaco (ii,ii+1) = eps_ox;
            elseif  (ii == Inter1)
                res(ii,1) = eps_si*phi(ii+1,1)-(eps_si+eps_ox)*phi(ii,1)+eps_ox*phi(ii-1,1);
                Jaco(ii,ii-1) = eps_ox; Jaco(ii,ii) = -(eps_si+eps_ox); Jaco(ii,ii+1) = eps_si;
            elseif  (ii == Inter2)
                res(ii,1) = eps_ox*phi(ii+1,1)-(eps_ox+eps_si)*phi(ii,1)+eps_si*phi(ii-1,1);
                Jaco(ii,ii-1) = eps_si; Jaco(ii,ii) = -(eps_ox+eps_si); Jaco(ii,ii+1) = eps_ox;
            else
                res(ii,1) = eps_si*phi(ii+1,1)-2*eps_si*phi(ii,1)+eps_si*phi(ii-1,1);
                Jaco(ii,ii-1) = eps_si; Jaco(ii,ii) = -2*eps_si; Jaco(ii,ii+1) = eps_si;
            end
        end
        % Charge part %
        for ii = Inter1:Inter2
            if      (ii == Inter1)
                res(ii,1) = res(ii,1)-coef_Poi*(Nacc+ni*exp(phi(ii,1)/thermal))*0.5;
                Jaco(ii,ii) = Jaco(ii,ii)-coef_Poi*ni*exp(phi(ii,1)/thermal)/thermal*0.5;
            elseif  (ii==Inter2)
                res(ii,1) = res(ii,1)-coef_Poi*(Nacc+ni*exp(phi(ii,1)/thermal))*0.5;
                Jaco(ii,ii) = Jaco(ii,ii)-coef_Poi*ni*exp(phi(ii,1)/thermal)/thermal*0.5;
            else
                res(ii,1) = res(ii,1)-coef_Poi*(Nacc+ni*exp(phi(ii,1)/thermal));
                Jaco(ii,ii) = Jaco(ii,ii)-coef_Poi*ni*exp(phi(ii,1)/thermal)/thermal;
            end
        end
        update = inv(Jaco)*(-res);
        phi = phi+update;
        phi_Poi(:,i) = phi;
        
        for ii=Inter1:Inter2
            n_sweep_Poi(ii,i) = ni*exp(phi(ii,1)/thermal);
        end
    end
    
    
    % Schrodinger solver %
    for iNewton = 1:100
        totalNumber = 0;
        for iValley = 1:3
            mass = ones(3)*0.19;
            mass(iValley) = 0.91;
            coef_Sch = 2*Lx*Ly/(2*pi)*sqrt(mass(1)*mass(2))*m0/(hbar^2)*(k_B*T);
            
        V = q*Ec_Ei -q*phi; % Potential energy, J
        Nbulk = Inter2-Inter1-1;
        Hamil = zeros(Nbulk,Nbulk);
        Hamil(1,1) = -2; Hamil(1,2) = 1;
        Hamil(Nbulk,Nbulk) = -2; Hamil(Nbulk,Nbulk-1) = 1;
        for ii = 2:Nbulk-1
            Hamil(ii,ii+1) =  1;
            Hamil(ii,ii  ) = -2;
            Hamil(ii,ii-1) =  1;
        end
        for ii = 1:Nbulk
            Hamil(ii,ii) = Hamil(ii,ii)-2*mass(3)*m0*(Delta_z/hbar)^2*V(ii+Inter1,1);
        end
        [eigenvectors,eigenvalues] = eig(Hamil);
        scaledEz = diag(eigenvalues)/(-2*mass(3)*m0*(Delta_z/hbar)^2); % Eigenenergy, J
        [sortedEz,sortedIndex] = sort(scaledEz);
        nSubband = 10;
        n_Sch = zeros(N,1);
        totalNumber = 0;
        for n = 1:nSubband
            Ez = sortedEz(n,1);
            wavefunction2 = eigenvectors(:,sortedIndex(n)).^2;
            normalization = sum(wavefunction2)*Delta_z;
            wavefunction2 = wavefunction2 / normalization;
            subbandNumber = coef_Sch*log(1+exp(-Ez/(k_B*T)));
            totalNumber = totalNumber+subbandNumber;
            n_Sch(Inter1+1:Inter2-1,1) = n_Sch(Inter1+1:Inter2-1,1)+1/(Lx*Ly)*wavefunction2*subbandNumber;
        end
        end
        % Laplacian %
        for ii = 2:N-1
            if      (ii<Inter1 || ii>Inter2)
                res(ii,1) = eps_ox*phi(ii+1,1)-2*eps_ox*phi(ii,1)+eps_ox*phi(ii-1,1);
                Jaco(ii,ii-1) = eps_ox; Jaco(ii,ii) = -2*eps_ox; Jaco (ii,ii+1) = eps_ox;
            elseif  (ii == Inter1)
                res(ii,1) = eps_si*phi(ii+1,1)-(eps_si+eps_ox)*phi(ii,1)+eps_ox*phi(ii-1,1);
                Jaco(ii,ii-1) = eps_ox; Jaco(ii,ii) = -(eps_si+eps_ox); Jaco(ii,ii+1) = eps_si;
            elseif  (ii == Inter2)
                res(ii,1) = eps_ox*phi(ii+1,1)-(eps_ox+eps_si)*phi(ii,1)+eps_si*phi(ii-1,1);
                Jaco(ii,ii-1) = eps_si; Jaco(ii,ii) = -(eps_ox+eps_si); Jaco(ii,ii+1) = eps_ox;
            else
                res(ii,1) = eps_si*phi(ii+1,1)-2*eps_si*phi(ii,1)+eps_si*phi(ii-1,1);
                Jaco(ii,ii-1) = eps_si; Jaco(ii,ii) = -2*eps_si; Jaco(ii,ii+1) = eps_si;
            end
        end
        % Charge part %
        for ii = Inter1:Inter2
            if      (ii == Inter1)
                res(ii,1) = res(ii,1)-coef_Poi*(Nacc+n_Sch(ii,1))*0.5;
                Jaco(ii,ii) = Jaco(ii,ii)-coef_Poi*n_Sch(ii,1)/thermal*0.5;
            elseif  (ii==Inter2)
                res(ii,1) = res(ii,1)-coef_Poi*(Nacc+n_Sch(ii,1))*0.5;
                Jaco(ii,ii) = Jaco(ii,ii)-coef_Poi*n_Sch(ii,1)/thermal*0.5;
            else
                res(ii,1) = res(ii,1)-coef_Poi*(Nacc+n_Sch(ii,1));
                Jaco(ii,ii) = Jaco(ii,ii)-coef_Poi*n_Sch(ii,1)/thermal;
            end
        end
        update = inv(Jaco)*(-res);
        phi = phi+update;
        phi_Poi(:,i) = phi;
        n_sweep_Sch(:,i) = n_Sch;
    end
end
  
n_integrated_Poi = sum(n_sweep_Poi*0.1e-9,1);
n_integrated_Sch = sum(n_sweep_Sch*0.1e-9,1);
 
figure(1);
semilogy(z/1e-9,n_sweep_Poi*1e-6); % position vs electron density, *1e-6 : m^-3->cm^-3
xlabel('Position (nm)');
ylabel('Electron density (cm^-^3)');
legend('0V','0.05V','0.10V','0.15V','0.20V','0.25V','0.30V','0.35V','0.40V','0.45V','0.50V','0.55V','0.60V','0.65V','0.70V','0.75V','0.80V','0.85V','0.90V','0.95V','1.00V','Location','Best');

figure(2);
semilogy(z/1e-9,n_sweep_Sch*1e-6); % position vs electron density, *1e-6 : m^-3->cm^-3
xlabel('Position (nm)');
ylabel('Electron density (cm^-^3)');
legend('0V','0.05V','0.10V','0.15V','0.20V','0.25V','0.30V','0.35V','0.40V','0.45V','0.50V','0.55V','0.60V','0.65V','0.70V','0.75V','0.80V','0.85V','0.90V','0.95V','1.00V','Location','Best');

figure(3);
semilogy(Vg_sweep,n_integrated_Poi*1e-4,Vg_sweep,n_integrated_Sch*1e-4); % Gate voltage vs integrated electron density
xlabel('Gate Voltage (V)');
ylabel('Electron density (cm^-^2)');
legend('Nonlinear poisson equation','Schrodinger-Poisson equation','Location','Best');
