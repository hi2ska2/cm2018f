clear all
clc

q=1.602192e-19; % elementary charge, [C]
eps0=8.854187817e-12; % vacuum permittivity, [F/m]
k_B=1.380662e-23; %Boltzmann constant, [J/K]
T=300.0; % temperature, [K]
Deltax=0.1e-9; % 0.1 nm [m]
N=61; % 6 nm thick
x=Deltax*transpose([0:N-1]); % real space, nm
interface1=6; % at x=0.5 nm
interface2=56; % at x=5.5 nm
eps_si=11.7; eps_ox=3.9; % relatvie permittivity
Nacc=1e24; %1e18/cm^3
ni=1.075e16; %1.075e10/cm^3 
WF=-4.30; % vacuum level-fermi level [eV]
IF=-4.63374; % vacuum level-intrinsic fermi level [eV]
V_g=[0:0.05:1]'; % set: fermi level= 0 V at equilibrium 
% b_cont=0-(IF-WF-V_g) % using elecmentary charge unit [V]

for j=1:21;
    b_cont(j,1)=-(IF-WF-V_g(j,1)); % using elecmentary charge unit
A=zeros(N,N);
A(1,1)=1.0;
for ii=2:N-1
    if     (ii<interface1)  A(ii,ii-1)=eps_ox; A(ii,ii)=-2*eps_ox;      A(ii,ii+1)=eps_ox;
    elseif (ii==interface1) A(ii,ii-1)=eps_ox; A(ii,ii)=-eps_ox-eps_si; A(ii,ii+1)=eps_si;
    elseif (ii<interface2)  A(ii,ii-1)=eps_si; A(ii,ii)=-2*eps_si;      A(ii,ii+1)=eps_si;
    elseif (ii==interface2) A(ii,ii-1)=eps_si; A(ii,ii)=-eps_si-eps_ox; A(ii,ii+1)=eps_ox;
    elseif (ii>interface2)  A(ii,ii-1)=eps_ox; A(ii,ii)=-2*eps_ox;      A(ii,ii+1)=eps_ox;
    end
end
A(N,N)=1.0;

%%%%%%%%%%%%%%%% first cacluation of Poisson's equation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b=zeros(N,1);
b(1,1)=b_cont(j,1);
for ii=interface1:interface2
    if     (ii==interface1) b(ii,1)=Deltax*Deltax*q*Nacc/eps0*0.5;
    elseif (ii==interface2) b(ii,1)=Deltax*Deltax*q*Nacc/eps0*0.5;
    else                    b(ii,1)=Deltax*Deltax*q*Nacc/eps0;
    end
end
b(N,1)=b_cont(j,1);

phi=A\b;


%%%%%%%%% carrier density calculated from Poisson's equation and Boltzmann distribution%%%%%%%%%%%%%%%%%%%%%%%%
elec=zeros(N,1);
for ii=interface1:interface2
    elec(ii,1)=ni*exp(q*phi(ii,1)/(k_B*T));
end
elec1(:,j)=elec; % V_g vs. elec

plot(x/1e-9,elec1*1e-6);
xlabel('Position (nm)')
ylabel('Electron density (cm^-^3)')


%%%%%%%%%% 2nd calculation of Poisson's equation by using calculated carrier density %%%%%%%%%%%%%%%%%%%%%%%%%%%
b_sc=zeros(N,1); % source term Nacc-> (Nacc+elec)
b_sc(1,1)=b_cont(j,1); 
for ii=interface1:interface2
    if     (ii==interface1) b_sc(ii,1)=Deltax*Deltax*q*(Nacc+elec(ii,1))/eps0*0.5;
    elseif (ii==interface2) b_sc(ii,1)=Deltax*Deltax*q*(Nacc+elec(ii,1))/eps0*0.5;
    else                    b_sc(ii,1)=Deltax*Deltax*q*(Nacc+elec(ii,1))/eps0;
    end
end
b_sc(N,1)=b_cont(j,1);

phi_sc=A\b_sc;

diff_potential(:,j)=phi-phi_sc; % difference between 1st and 2nd calculation of poisson's equation

% plot(x*1e9,phi,'r',x*1e9,phi_sc,'b')
%  plot(x*1e9,phi_sc,'b')
% plot(x*1e9,phi,'r')
% plot(x*1e9,diff_potential)
% xlabel('position [nm]')
% ylabel('potential [V]')
% legend('1st calculation','2nd calculation')


end

subplot(1,2,1)
semilogy(V_g,diff_potential(31,:))
xlabel('Applied voltage [V]')
ylabel('Difference in two potential [V]')
subplot(1,2,2)
plot(V_g,diff_potential(31,:))
xlabel('Applied voltage [V]')
ylabel('Difference in two potential [V]')