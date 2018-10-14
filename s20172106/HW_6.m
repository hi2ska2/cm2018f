clear all;
clc

%% Defining variables

q = 1.602192e-19; % Elementary charge, C
eps0 = 8.854187817e-12; % Vacuum permittivity, F/m
k_B = 1.380662e-23; % Boltzmann constant, J/K
T = 300.0; % Temperature, K
thermal = k_B*T/q; % Thermal voltage, V
Deltax = 0.1e-9; % 0.1 nm spacing
N = 61; % 6 nm thick
x = Deltax*transpose([0:N-1]); % real space, m
interface1 = 6; % At x=0.5 nm
interface2 = 56; % At x=5.5 nm
eps_si = 11.7; eps_ox = 3.9; % Relative permittivity
Nacc = 1e24; % 1e18 /cm^3
ni = 1.075e16; % 1.075e10 /cm^3
coef = Deltax*Deltax*q/eps0;
Work_Ftn=-4.30; % vacuum level-fermi level [eV]
Int_Fermi=-4.63374; % vacuum level-intrinsic fermi level [eV]
Volt_Gate = [-1:0.1:1]'; % set: fermi level= 0 V at equilibrium 
kk = size(Volt_Gate); % The number of the different Volt_Gate's

phi = zeros(N,kk(1));
%% Potential by using the Newton method
Contact_Volt = zeros(kk(1),1);

for k = 1:kk(1) %% or kk(1) for Gate voltage from 0 V to 1 V 

Contact_Volt(k) = ((Work_Ftn - Volt_Gate(k,1))- Int_Fermi);

for newton = 1:1000
  

    res = zeros(N,1);
    Jaco = sparse(N,N);
    res(1,1) = phi(1,k) - Contact_Volt(k);
    Jaco(1,1) = 1.0;
    res(N,1) = phi(N,k) - Contact_Volt(k);
    Jaco(N,N) = 1.0;
    %%%%%% Laplacian part %%%%%%%%%%%%%%%%%%%%%%%%
    
    for ii = 2:N-1
        if (ii< interface1 || ii> interface2)
            res(ii,1) = eps_ox*phi(ii+1,k) - 2*eps_ox*phi(ii,k) + eps_ox*phi(ii-1,k);
            Jaco(ii,ii-1) = eps_ox; Jaco(ii,ii) = -2*eps_ox; Jaco(ii,ii+1) = eps_ox;
        elseif (ii == interface1)
            res(ii,1) = eps_si*phi(ii+1,k) - (eps_ox + eps_si)*phi(ii,k) + eps_ox*phi(ii-1,k);
            Jaco(ii,ii-1) = eps_ox; Jaco(ii,ii) = - (eps_ox + eps_si); Jaco(ii,ii+1) = eps_si;
        elseif (ii == interface2)
             res(ii,1) = eps_ox*phi(ii+1,k) - (eps_ox + eps_si)*phi(ii,k) + eps_si*phi(ii-1,k);
            Jaco(ii,ii-1) = eps_si; Jaco(ii,ii) = - (eps_ox + eps_si); Jaco(ii,ii+1) = eps_ox;
        else
             res(ii,1) = eps_si*phi(ii-1,k) - 2*eps_si*phi(ii,k) + eps_si*phi(ii+1,k);
            Jaco(ii,ii-1) = eps_si; Jaco(ii,ii) = -2*eps_si; Jaco(ii,ii+1) = eps_si;
        end
    end
    
    %%%%%%% If there are charges (==> Poisson's eqn) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for ii = interface1:interface2
         if (ii == interface1)
             res(ii,1) = res(ii,1) - coef*(Nacc + ni*exp(phi(ii,k)/thermal))*0.5;
             Jaco(ii,ii) = Jaco(ii,ii) - coef*ni*exp(phi(ii,k)/thermal)/thermal*0.5;
         elseif (ii == interface2)
             res(ii,1) = res(ii,1) - coef*(Nacc + ni*exp(phi(ii,k)/thermal))*0.5;
             Jaco(ii,ii) = Jaco(ii,ii) - coef*ni*exp(phi(ii,k)/thermal)/thermal*0.5;
         else
             res(ii,1) = res(ii,1) - coef*(Nacc + ni*exp(phi(ii,k)/thermal));
             Jaco(ii,ii) = Jaco(ii,ii) - coef*ni*exp(phi(ii,k)/thermal)/thermal;
         end
    end
    
    update = Jaco \ (-res);
    phi(:,k) = phi(:,k) + update;
end

end       
            
%% The integrated electron density
elec = zeros(N,kk(1));
Inte_elec = zeros(kk(1),1);
for s = 1:kk(1)

for ii=interface1:interface2
    elec(ii,s)=ni*exp(phi(ii,s)/thermal);
end
Inte_elec(s) = trapz(elec(:,s))/1e8; % Scaling for the integrated electron density 1/cm^2

end
figure
for tt = 1:kk(1)
semilogy(x,elec(:,tt)); hold on
end
xlabel('Position [nm]');
ylabel('n_e [1/cm^3]]');
figure
semilogy(Volt_Gate,Inte_elec); 
xlabel('Gate Voltage [V]');
ylabel('Integrated n_e [1/cm^2]]');


