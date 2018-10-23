clear all
clc

q=1.602192e-19; % elementary charge, C
eps0=8.854187817e-12; % vacuum permittivity, F/m
k_B=1.380662e-23; % Boltzmann constatant, J/K
T=300.0; % temperature, K
thermal=k_B*T/q; % Thermal voltage, V
Deltax=0.1e-9; % 0.1 nm spacing
N=61; % 6 nm thick
x=Deltax*transpose([0:N-1]); % real space, m
interface1=6; % At x=0.5 nm
interface2=56; % At x=5.5 nm
eps_si=11.7; eps_ox=3.9; % Relative permittivity
Nacc=1e24; % 1e18 /cm^3
ni=1.075e16; %1.075e10/cm^3
coef=Deltax*Deltax*q/eps0;
IF=-4.63374; % vacuum level-intrinsic fermi level [eV]
WF=-4.30; % vacuum level-fermi level [eV]
V_appl=transpose([0:0.05:1]); % applied voltage set: fermi level= 0 V at equilibrium 
x=Deltax*transpose([0:N-1]); % real space, nm
elec=zeros(N,1);

for j=1:length(V_appl)
    
    V_g(j,1)=0-((IF-WF-V_appl(j,1))); % applied gate voltage  [V]
%% boundary condition
phi=zeros(N,1);
phi(:,1)=V_g(j,1); % initial value

    res=zeros(N,1);
    Jaco=sparse(N,N);
    res(1,1)=phi(1,1)-V_g(j,1);
    Jaco(1,1)=1.0;
    res(N,1)=phi(N,1)-V_g(j,1);
    Jaco(N,N)=1.0; 


count(j,1)=0;

for newton=1:50000


% Laplacian
for ii=2:N-1
    if (ii<interface1 || ii>interface2)
        res(ii,1)*eps_ox*phi(ii+1,1)-2*eps_ox*phi(ii,1)+eps_ox*phi(ii-1,1);
        Jaco(ii,ii-1)=eps_ox; Jaco(ii,ii)=-2*eps_ox; Jaco(ii,ii+1)=eps_ox;
    elseif (ii==interface1)
        res(ii,1)=eps_si*phi(ii+1,1)-(eps_si+eps_ox)*phi(ii,1)+eps_ox*phi(ii-1,1);
        Jaco(ii,ii-1)=eps_ox; Jaco(ii,ii)=-(eps_si+eps_ox); Jaco(ii,ii+1)=eps_si;
    elseif (ii==interface2);
        res(ii,1)=eps_ox*phi(ii+1,1)-(eps_ox+eps_si)*phi(ii,1)+eps_si*phi(ii-1,1);
        Jaco(ii,ii-1)=eps_si; Jaco(ii,ii)=-(eps_ox+eps_si); Jaco(ii,ii+1)=eps_ox;
    else
        res(ii,1)=eps_si*phi(ii+1,1)-2*eps_si*phi(ii,1)+eps_si*phi(ii-1,1);
        Jaco(ii,ii-1)=eps_si; Jaco(ii,ii)=-2*eps_si; Jaco(ii,ii+1)=eps_si;
    end
end


%% charge part
for ii=interface1:interface2
    if (ii==interface1)
        res(ii,1)=res(ii,1)-coef*(Nacc+ni*exp(phi(ii,1)/thermal))*0.5;
        Jaco(ii,ii)=Jaco(ii,ii)-coef*ni*exp(phi(ii,1)/thermal)/thermal*0.5;
    elseif (ii==interface2)
        res(ii,1)=res(ii,1)-coef*(Nacc+ni*exp(phi(ii,1)/thermal))*0.5;
        Jaco(ii,ii)=Jaco(ii,ii)-coef*ni*exp(phi(ii,1)/thermal)/thermal*0.5;
    else
        res(ii,1)=res(ii,1)-coef*(Nacc+ni*exp(phi(ii,1)/thermal));
        Jaco(ii,ii)=Jaco(ii,ii)-coef*ni*exp(phi(ii,1)/thermal)/thermal;
    end
end
update=Jaco\(-res);
phi=phi+update;
count(j,1)=count(j,1)+1;

if abs(update(interface1+1:interface2-1,1))<= 5e-16;
   break
end
phi1(:,j)=phi;

end
count(j,1)=count(j,1);


for ii=interface1:interface2
%     elec(ii,j)=ni*exp(q*phi(ii,1)/(k_B*T))*1e-6;
%     
    elec(ii,j)=ni*exp(q*phi(ii,1)/(k_B*T));
end
end

% elec_int=sum(elec,1)*1e-6*5e-7; % integrated carrier density, cm^-2
elec_int=sum(elec*0.1e-9,1); % integrated carrier density, m^-2

f1=figure
plot(x/1e-9,phi1) % position vs potential
xlabel('Position (nm)')
ylabel('Potential (V)')

f2=figure
plot(x/1e-9,elec*1e-6); % position vs electron density, *1e-6 : m^-3->cm^-3
xlabel('Position (nm)')
ylabel('Electron density (cm^-^3)')

f3=figure
semilogy(V_appl,elec_int*1e-4); % Gate voltage vs integrated electron density
xlabel('Gate Voltage')
ylabel('Electron density (cm^-^2)')





% %% Compare #4 potential % electron density
% for j=1:length(V_appl);
% A(1,1)=1.0;
% for ii=2:N-1
%     if     (ii<interface1)  A(ii,ii-1)=eps_ox; A(ii,ii)=-2*eps_ox;      A(ii,ii+1)=eps_ox;
%     elseif (ii==interface1) A(ii,ii-1)=eps_ox; A(ii,ii)=-eps_ox-eps_si; A(ii,ii+1)=eps_si;
%     elseif (ii<interface2)  A(ii,ii-1)=eps_si; A(ii,ii)=-2*eps_si;      A(ii,ii+1)=eps_si;
%     elseif (ii==interface2) A(ii,ii-1)=eps_si; A(ii,ii)=-eps_si-eps_ox; A(ii,ii+1)=eps_ox;
%     elseif (ii>interface2)  A(ii,ii-1)=eps_ox; A(ii,ii)=-2*eps_ox;      A(ii,ii+1)=eps_ox;
%     end
% end
% A(N,N)=1.0;
% 
% %%%%%%%%%%%%%%%% first cacluation of Poisson's equation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% b=zeros(N,1);
% b(1,j)=V_g(j,1);
% for ii=interface1:interface2
%     if     (ii==interface1) b(ii,j)=Deltax*Deltax*q*Nacc/eps0*0.5;
%     elseif (ii==interface2) b(ii,j)=Deltax*Deltax*q*Nacc/eps0*0.5;
%     else                    b(ii,j)=Deltax*Deltax*q*Nacc/eps0;
%     end
% end
% b(N,j)=V_g(j,1);
% 
% phi_pre(:,j)=A\b(:,j);
% 
% 
% %%%%%%%%% carrier density calculated from Poisson's equation and Boltzmann distribution%%%%%%%%%%%%%%%%%%%%%%%%
% elec_pre=zeros(N,1);
% for ii=interface1:interface2
%     elec_pre(ii,1)=ni*exp(q*phi_pre(ii,j)/(k_B*T));
% end
% elec_pre1(:,j)=elec_pre(:,1);
% 
% end
% 
% f4=figure 
% plot(x/1e-9,elec_pre1*1e-6);
% xlabel('Position (nm)')
% ylabel('Electron density (cm^-^3)')
% 
% 
% f5=figure
% plot(x/1e-9,(elec-elec_pre1)*1e-6)
% xlabel('Position (nm)')
% ylabel('Difference of electron density (cm^-^3)')
% 
% 
% f6=figure
% plot(x/1e-9,phi1-phi_pre)
% xlabel('Position (nm)')
% ylabel('Difference of potential (V)')
