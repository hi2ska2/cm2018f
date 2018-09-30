clear all
clc

% Input variables 

q=1.602192e-19; % elementary charge, [C]
eps0=8.854187817e-12; % vacuum permittivity, [F/m]
k_B=1.380662e-23; % Boltzmann constant, [J/K]
T=300.0; % temperature, [K]
Deltax=0.1e-9; % 0.1 nm [m]
N=61; % 6 nm thick
x=Deltax*[0:N-1]'; % real space, nm
interface1=6; % at x=0.5 nm
interface2=56; % at x=5.5 nm
eps_si=11.7; eps_ox=3.9; % relatvie permittivity
Nacc=1e24; %1e18/cm^3
ni=1.075e16; %1.075e10/cm^3 
Work_Ftn=-4.30; % vacuum level-fermi level [eV]
Int_Fermi=-4.63374; % vacuum level-intrinsic fermi level [eV]
Volt_Gate = [0:0.1:1]'; % set: fermi level= 0 V at equilibrium 
kk = size(Volt_Gate); % The number of the different Volt_Gate's
b_recal = zeros(N,kk(1));
phi_recal = zeros(N,kk(1));
diff_phi = zeros(N,kk(1));

Gate_volt_vs_diff_phi = zeros(kk(1),2); 

% Constructing of Poisson's eqn

for k = 1:kk(1) %% or kk(1) for Gate voltage from 0 V to 1 V 
    %%%%%%%%%%%%%%%%%%%%% Construction for the euation of Poisson's eqn%%%%%%%%%%%%%%%%%%%%%
    A = zeros(N,N);
    b_char = ((Work_Ftn - Volt_Gate(k,1))- Int_Fermi); % Charge density under Gate_voltage
    
 A(1,1) = 1.0;
 
   for ii=2:N-1
    if     (ii<interface1)  A(ii,ii-1)=eps_ox; A(ii,ii)=-2*eps_ox;      A(ii,ii+1)=eps_ox;
    elseif (ii==interface1) A(ii,ii-1)=eps_ox; A(ii,ii)=-eps_ox-eps_si; A(ii,ii+1)=eps_si;
    elseif (ii<interface2)  A(ii,ii-1)=eps_si; A(ii,ii)=-2*eps_si;      A(ii,ii+1)=eps_si;
    elseif (ii==interface2) A(ii,ii-1)=eps_si; A(ii,ii)=-eps_si-eps_ox; A(ii,ii+1)=eps_ox;
    elseif (ii>interface2)  A(ii,ii-1)=eps_ox; A(ii,ii)=-2*eps_ox;      A(ii,ii+1)=eps_ox;
    end
    end
A(N,N)=1.0;

%%%%%%%%%%%%%%%%%%%%%%%%%% Considering only Nacc as the charge density %%%%%%%%%%%%%% 

b = zeros(N,1);
b(1,1) = b_char;
for ii=interface1:interface2
    if     (ii==interface1) b(ii,1)=Deltax*Deltax*q*Nacc/eps0*0.5;
    elseif (ii==interface2) b(ii,1)=Deltax*Deltax*q*Nacc/eps0*0.5;
    else                    b(ii,1)=Deltax*Deltax*q*Nacc/eps0;
    end
end
b(N,1) = b_char;

phi = A\b;

% Plot %
% figure;
% plot(x*1e9,phi*1e3);
% xlabel('X distance [nm]');
% ylabel('Potential [mV]');
% title('Potential of Nacc vs x distance','fontsize',15);



%%%%%%%%%%%%%%%%%%%%%% Calculating for the electron density %%%%%%%%%%%%%%%%%%%%%%%

elec = zeros(N,1);
for ii=interface1:interface2
    elec(ii,1)=ni*exp(q*phi(ii,1)/(k_B*T));
end

% Plot %
% figure;
% plot(x*1e9,elec*1e-6);
% xlabel('X distance [nm]');
% ylabel('Electron density [cm^-^3]');
% title('Electron density vs X distance','fontsize',15);

%%%%%%%%%%%%%%%%%%%%%% Recalculating for the electrostatic potential under reconsidering for the updated electron density  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


b_recal(1,k) = b_char;
b_recal(N,k) = b_char;
for ii=interface1:interface2
    if     (ii==interface1) b_recal(ii,k)=Deltax*Deltax*q*(Nacc + elec(ii,1))/eps0*0.5;
    elseif (ii==interface2) b_recal(ii,k)=Deltax*Deltax*q*(Nacc + elec(ii,1))/eps0*0.5;
    else                    b_recal(ii,k)=Deltax*Deltax*q*(Nacc + elec(ii,1))/eps0;
    end
end

phi_recal(:,k) = A\b_recal(:,k);

% Plot %
% figure;
% plot(x*1e9,phi*1e3,'r',x*1e9,phi_recal(:,1)*1e3,'b');
% xlabel('X distance [nm]');
% ylabel('Potential [mV]');
% legend('Only Nacc','N acc + Elec den')
% title('Potentials vs x distance','fontsize',15);

%%%%%%%%%%%%%%%%%%%%%%% Check their difference for several gate voltage from 0 V to 1V %%%%%%%%%%%%%%%%%%%%%%%

diff_phi(:,k) = phi - phi_recal(:,k);

% Plot %
% figure;
% plot(x*1e9,diff_phi(:,1));
% xlabel('X distance [nm]');
% ylabel('Differnce btw two potentials [mV]');
% title('Diff btw two potential vs x distance','fontsize',15);

end

Gate_volt_vs_diff_phi(:,1) = Volt_Gate;
Gate_volt_vs_diff_phi(:,2) = diff_phi(30,:);
figure;
plot(Gate_volt_vs_diff_phi(:,1),Gate_volt_vs_diff_phi(:,2)*1e3)
xlabel('Gate Potential [V]');
ylabel('Difference of potential [mV]');
title('Gate Voltage vs Diff Potential','fontsize',15);
