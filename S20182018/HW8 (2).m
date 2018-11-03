warning('off','all')
clear all;
close all;

h = 6.626176e-34; % Planck constant, J s
hbar = h / (2*pi); % Reduced Planck constant, J s
q = 1.602192e-19; % Elementary charge, C
m0 = 9.109534e-31; % Electron rest mass, kg
k_B = 1.380662e-23; % Boltzmann constant, J/K
eps0 = 8.854187817e-12; % Vacuum permittivity, F/m
T = 300.0; % Temperature, K
thermal = k_B*T/q; % Thermal voltage, V
mxx = 0.19; myy = 0.19; mzz = 0.91; % Masses, m0
Deltaz = 0.1e-9; % 0.1 nm spacing
Nz = 61; % 6 nm thick
z = Deltaz*transpose([0:Nz-1]); % real space, m
interface1 = 6; % At z=0.5 nm
interface2 = 56; % At z=5.5 nm
eps_si = 11.7; eps_ox = 3.9; % Relative permittivity
Nacc = 1e24; % 1e18 /cm^3
ni = 1.075e16; % 1.075e10 /cm^3
coef = Deltaz*Deltaz*q/eps0;
Ec_Ei = 0.561004; % E_c ? E_i, eV
Lx = 100e-9; Ly = 100e-9; % Lenghs, m
coef_Sch = 2*Lx*Ly/(2*pi)*sqrt(mxx*myy)*m0/(hbar^2)*(k_B*T);
M=101;

integ_elec(M,1)=zeros;
integ_elec_poi(M,1)=zeros;
gatev(M,1)=zeros;
elec = zeros(Nz,M); % Electron density, /m^3
phi = zeros(Nz,M);
phi(:,1) = 0.33374;
elec_poi=zeros(Nz,M);

for gate=1:M;

    gatev(gate,1) = (gate-1.0)/(101.0);
    phi(1,gate)=0.33374+((gate-1.0)/101.0);
    phi(Nz,gate)=phi(1,gate);
    Jaco(1,1) = 1.0;
    Jaco(Nz,Nz) = 1.0;
    res(1,1) = phi(1,gate) - 0.33374-((gate-1.0)/101.0);
    res(Nz,1) = res(1,1);
    
    for newton=1:50
        for ii=2:Nz-1
            if (ii< interface1 || ii>interface2)
                res(ii,1) = eps_ox*phi(ii+1,gate) - 2*eps_ox*phi(ii,gate) + eps_ox*phi(ii-1,gate);
                Jaco(ii,ii-1) = eps_ox;
                Jaco(ii,ii) = -2*eps_ox;
                Jaco(ii,ii+1) = eps_ox;
            elseif (ii==interface1)
                res(ii,1) = eps_si*phi(ii+1,gate) - (eps_si+eps_ox)*phi(ii,gate) + eps_ox*phi(ii-1,gate);
                Jaco(ii,ii-1) = eps_ox;
                Jaco(ii,ii) = -(eps_si+eps_ox);
                Jaco(ii,ii+1) = eps_si;
            elseif (ii==interface2)
                res(ii,1) = eps_ox*phi(ii+1,gate) - (eps_ox+eps_si)*phi(ii,gate) + eps_si*phi(ii-1,gate);
                Jaco(ii,ii-1) = eps_si;
                Jaco(ii,ii) = -(eps_ox+eps_si);
                Jaco(ii,ii+1) = eps_ox;
            else
                res(ii,1) = eps_si*phi(ii+1,gate) - 2*eps_si*phi(ii,gate) + eps_si*phi(ii-1,gate);
                Jaco(ii,ii-1) = eps_si;
                Jaco(ii,ii) = -2*eps_si;
                Jaco(ii,ii+1) = eps_si;
            end
        end   
        for ii=interface1:interface2
            if (ii==interface1)
                res(ii,1) = res(ii,1) - coef*(Nacc+ni*exp(phi(ii,gate)/thermal))*0.5;
                Jaco(ii,ii) = Jaco(ii,ii) - coef*ni*exp(phi(ii,gate)/thermal)/thermal*0.5;
            elseif(ii==interface2)
                res(ii,1) = res(ii,1) - coef*(Nacc+ni*exp(phi(ii,gate)/thermal))*0.5;
                Jaco(ii,ii) = Jaco(ii,ii) - coef*ni*exp(phi(ii,gate)/thermal)/thermal*0.5;
            else
                res(ii,1) = res(ii,1) - coef*(Nacc+ni*exp(phi(ii,gate)/thermal));
                Jaco(ii,ii) = Jaco(ii,ii) - coef*ni*exp(phi(ii,gate)/thermal)/thermal;
            end
        end   
        update= Jaco \ (-res);
        phi(:,gate) = phi(:,gate) + update;
    end
    
    for i=1:Nz
        elec_poi(i,gate) = ni*exp(phi(i,gate)/thermal);
    end
    
    for newton2=1:10
        V = q*Ec_Ei - q*phi(:,gate); % Potential energy, J
        Nbulk = interface2-interface1-1; % Number of bulk silicon nodes
        Hamil = zeros(Nbulk,Nbulk);
        Hamil(1,1) = -2; Hamil(1,2) = 1;
        for ii=2:Nbulk-1
            Hamil(ii,ii+1) = 1;
            Hamil(ii,ii ) = -2;
            Hamil(ii,ii-1) = 1;
        end
        Hamil(Nbulk,Nbulk) = -2; Hamil(Nbulk,Nbulk-1) = 1;
        for ii=1:Nbulk
            Hamil(ii,ii) = Hamil(ii,ii) -2*mzz*m0*(Deltaz/hbar)^2*V(ii+interface1,1);
        end
        
        [eigenvectors,eigenvalues] = eig(Hamil);
        scaledEz = diag(eigenvalues)/(-2*mzz*m0*(Deltaz/hbar)^2); % Eigenenergy, J
        [sortedEz,sortedIndex] = sort(scaledEz);
        nSubband = 10;
        totalNumber = 0;
        for n=1:nSubband
            Ez = sortedEz(n,1);
            wavefunction2 = eigenvectors(:,sortedIndex(n)).^2;
            normalization = sum(wavefunction2)*Deltaz;
            wavefunction2 = wavefunction2 / normalization;
            subbandNumber = coef_Sch*log(1+exp(-Ez/(k_B*T)));
            totalNumber = totalNumber + subbandNumber;
            elec(interface1+1:interface2-1,gate) = elec(interface1+1:interface2-1,gate) + 1/(Lx*Ly)*wavefunction2*subbandNumber;
        end
        for new3=1:10

        for ii=2:Nz-1
            if (ii< interface1 || ii>interface2)
                res(ii,1) = eps_ox*phi(ii+1,gate) - 2*eps_ox*phi(ii,gate) + eps_ox*phi(ii-1,gate);
                Jaco(ii,ii-1) = eps_ox;
                Jaco(ii,ii) = -2*eps_ox;
                Jaco(ii,ii+1) = eps_ox;
            elseif (ii==interface1)
                res(ii,1) = eps_si*phi(ii+1,gate) - (eps_si+eps_ox)*phi(ii,gate) + eps_ox*phi(ii-1,gate);
                Jaco(ii,ii-1) = eps_ox;
                Jaco(ii,ii) = -(eps_si+eps_ox);
                Jaco(ii,ii+1) = eps_si;
            elseif (ii==interface2)
                res(ii,1) = eps_ox*phi(ii+1,gate) - (eps_ox+eps_si)*phi(ii,gate) + eps_si*phi(ii-1,gate);
                Jaco(ii,ii-1) = eps_si;
                Jaco(ii,ii) = -(eps_ox+eps_si);
                Jaco(ii,ii+1) = eps_ox;
            else
                res(ii,1) = eps_si*phi(ii+1,gate) - 2*eps_si*phi(ii,gate) + eps_si*phi(ii-1,gate);
                Jaco(ii,ii-1) = eps_si;
                Jaco(ii,ii) = -2*eps_si;
                Jaco(ii,ii+1) = eps_si;
            end
        end   
        for ii=interface1:interface2
            if (ii==interface1)
                res(ii,1) = res(ii,1) - coef*(Nacc+elec(ii,gate))*0.5;
                Jaco(ii,ii) = Jaco(ii,ii) - coef*(elec(ii+1,gate)-elec(ii,gate))/(phi(ii+1,gate)-phi(ii,gate))*0.5;
            elseif(ii==interface2)
                res(ii,1) = res(ii,1) - coef*(Nacc+elec(ii,gate))*0.5;
                Jaco(ii,ii) = Jaco(ii,ii) - coef*(elec(ii,gate)-elec(ii-1,gate))/(phi(ii,gate)-phi(ii-1,gate))*0.5;
            else
                res(ii,1) = res(ii,1) - coef*(Nacc+elec(ii,gate));
                Jaco(ii,ii) = Jaco(ii,ii) - coef*(elec(ii+1,gate)-elec(ii-1,gate))/(phi(ii+1,gate)-phi(ii-1,gate));
            end
        end   
        update= Jaco \ (-res);
        phi(:,gate) = phi(:,gate) + update;
        end
    end
        
    for ll=interface1:interface2
         if (ll== interface1 || ll==interface2)
             integ_elec(gate,1) = integ_elec(gate,1)+ Deltaz/2 * elec(ll,gate);
             integ_elec_poi(gate,1) = integ_elec_poi(gate,1)+ Deltaz/2 * elec_poi(ll,gate);
         else
             integ_elec(gate,1) = integ_elec(gate,1)+ Deltaz * elec(ll,gate);
             integ_elec_poi(gate,1) = integ_elec_poi(gate,1)+ Deltaz * elec_poi(ll,gate);
         end       
    end
end



figure(4)
semilogy(gatev(:,1),integ_elec_poi(:,1)*1e-4); hold on;
semilogy(gatev(:,1),integ_elec(:,1)*1e-4); hold on;
xlabel ('Gate Voltage (V)');
ylabel ('electron density (cm^-^2)');
    
    