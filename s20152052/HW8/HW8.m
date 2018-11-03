clear all
clc


h = 6.626176e-34; % Planck constant, J s 
hbar = h / (2*pi); % Reduced Planck constant, J s 
q = 1.602192e-19; % Elementary charge, C 
m0 = 9.109534e-31; % Electron rest mass, kg 
k_B = 1.380662e-23; % Boltzmann constant, J/K 
eps0 = 8.854187817e-12; % Vacuum permittivity, F/m 
T = 300.0; % Temperature, K 
thermal = k_B*T/q; % Thermal voltage, V 

Lx=100e-9; Ly=100e-9; % length, m
% mxx = 0.19; myy = 0.19; mzz = 0.91; % Masses, m0  
Deltaz = 0.1e-9; % 0.1 nm spacing 
Nz = 61; % 6 nm thick 
z = Deltaz*transpose([0:Nz-1]); % real space, m 
interface1 = 6; % At z=0.5 nm 
interface2 = 56; % At z=5.5 nm 
eps_si = 11.7; eps_ox = 3.9; % Relative permittivity 
Nacc = 1e24; % 1e18 /cm^3 
ni = 1.075e16; % 1.075e10 /cm^3 
coef_Poi=Deltaz*Deltaz*q/eps0;
% coef_Sch= 2*Lx*Ly/(2*pi)*sqrt(mxx*myy)*m0/(hbar^2)*(k_B*T); 
Ec_Ei = 0.561004; % E_c ? E_i, eV 
elec=zeros(Nz,1);

IF=-4.63374; % vacuum level-intrinsic fermi level [eV]
WF=-4.30; % vacuum level-fermi level [eV]
V_appl=transpose([0:0.1:1]) 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi = zeros(Nz,1);  

for j=1:length(V_appl)
    V_G(j,1)=-(IF-WF-V_appl(j,1)); % applied gate voltage  [V]
    
    %% semiclassical nonlinear poisson solver
for newton=1:100   
    res = zeros(Nz,1);     
    Jaco = sparse(Nz,Nz);    
    res(1,1) = phi(1,1) - V_G(j,1);  
    Jaco(1,1) = 1.0;     
    res(Nz,1) = phi(Nz,1) - V_G(j,1);   
    Jaco(Nz,Nz) = 1.0;     
    
    % laplacian part
    for ii=2:Nz-1       
        if     (ii< interface1 || ii> interface2)           
            res(ii,1) = eps_ox*phi(ii+1,1) - 2*eps_ox*phi(ii,1) + eps_ox*phi(ii-1,1);          
            Jaco(ii,ii-1) = eps_ox; 
            Jaco(ii,ii) = -2*eps_ox;        
            Jaco(ii,ii+1) = eps_ox;        
        elseif (ii==interface1)           
            res(ii,1) = eps_si*phi(ii+1,1) - (eps_si+eps_ox)*phi(ii,1) + eps_ox*phi(ii-1,1);          
            Jaco(ii,ii-1) = eps_ox; 
            Jaco(ii,ii) = -(eps_si+eps_ox); 
            Jaco(ii,ii+1) = eps_si;        
        elseif (ii==interface2)           
            res(ii,1) = eps_ox*phi(ii+1,1) - (eps_ox+eps_si)*phi(ii,1) + eps_si*phi(ii-1,1);          
            Jaco(ii,ii-1) = eps_si; 
            Jaco(ii,ii) = -(eps_ox+eps_si); 
            Jaco(ii,ii+1) = eps_ox;         
        else
            res(ii,1) = eps_si*phi(ii+1,1) - 2*eps_si*phi(ii,1) + eps_si*phi(ii-1,1);
            Jaco(ii,ii-1) = eps_si; 
            Jaco(ii,ii) = -2*eps_si;        
            Jaco(ii,ii+1) = eps_si;       
        end
    end
    
    % charge part
    for ii=interface1:interface2        
        if     (ii==interface1)           
            res(ii,1) = res(ii,1) - coef_Poi*(Nacc+ni*exp(phi(ii,1)/thermal))*0.5;           
            Jaco(ii,ii) = Jaco(ii,ii) - coef_Poi*ni*exp(phi(ii,1)/thermal)/thermal*0.5;        
        elseif (ii==interface2)           
            res(ii,1) = res(ii,1) - coef_Poi*(Nacc+ni*exp(phi(ii,1)/thermal))*0.5;           
            Jaco(ii,ii) = Jaco(ii,ii) - coef_Poi*ni*exp(phi(ii,1)/thermal)/thermal*0.5;       
        else
            res(ii,1) = res(ii,1) - coef_Poi*(Nacc+ni*exp(phi(ii,1)/thermal));           
            Jaco(ii,ii) = Jaco(ii,ii) - coef_Poi*ni*exp(phi(ii,1)/thermal)/thermal;           
        end
    end
    update = Jaco \ (-res);    
    phi = phi + update; 
end
phi_NP(:,j)=phi


for ii=interface1:interface2
    elec(ii,1)=ni*exp(q*phi(ii,1)/(k_B*T));
end    
    elec1(:,j)=elec(:,1);
integ_elec(j,1)=sum(elec*0.1e-9,1)
% f1=figure
% plot(z,phi)
% 
% f2=figure
% plot(z,elec)


count(j,1)=0;

for iNewton=1:100
    totalNumber=0;
    elec_Sch=zeros(Nz,1);
    for iValley=1:3 
        mass=ones(3)*0.19;
        mass(iValley)=0.91;
        coef_Sch=2*Lx*Ly/(2*pi)*sqrt(mass(1)*mass(2))*m0/(hbar^2)*(k_B*T);
       
        %% Schrodinger solver
        V=q*Ec_Ei-q*phi; % potential energy, J
        Nbulk=interface2-interface1-1; % number of bulk silicon nodes
        Hamil=zeros(Nbulk,Nbulk);
        Hamil(1,1)=-2;
        Hamil(1,2)=1;
        
        for ii=2:Nbulk-1
            Hamil(ii,ii+1)=1;
            Hamil(ii,ii)=-2;
            Hamil(ii,ii-1)=1;
        end
        Hamil(Nbulk,Nbulk)=-2;
        Hamil(Nbulk,Nbulk-1)=1;
        
        for ii=1:Nbulk
            Hamil(ii,ii)=Hamil(ii,ii)-2*mass(3)*m0*(Deltaz/hbar)^2*V(ii+interface1,1);
        end
        [eigenvectors,eigenvalues]=eig(Hamil);
        Ez=diag(eigenvalues)/(-2*mass(3)*m0*(Deltaz/hbar)^2); 
        scaledEz=diag(eigenvalues)/(-2*mass(3)*m0*(Deltaz/hbar)^2); % Eigenenergy, J
        [sortedEz,sortedIndex]=sort(scaledEz)
        
        % f3=figure
        % plot(1:Nbulk,Ez)
        % f4=figure
        % plot(1:Nbulk,scaledEz)
        
       %% electron density
       nSubband=10;
%        totalNumber=0;
%        elec=zeros(Nz,1);
       for n=1:nSubband
           Ez=sortedEz(n,1);
           wavefunction2=eigenvectors(:,sortedIndex(n)).^2;
           normalization=sum(wavefunction2)*Deltaz;
           wavefunction2=wavefunction2/normalization;
           subbandNumber=coef_Sch*log(1+exp(-Ez/(k_B*T)));
           totalNumber=totalNumber+subbandNumber;
           elec_Sch(interface1+1:interface2-1,1)=elec_Sch(interface1+1:interface2-1,1)+...
               1/(Lx*Ly)*wavefunction2*subbandNumber;
       end
    end
    
%     plot(z,elec)
    

%% Laplacian %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=2:Nz-1
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
        res(ii,1)=res(ii,1)-coef_Poi*(Nacc+elec_Sch(ii,1))*0.5;
        Jaco(ii,ii)=Jaco(ii,ii)-coef_Poi*elec_Sch(ii,1)/thermal*0.5;
    elseif (ii==interface2)
        res(ii,1)=res(ii,1)-coef_Poi*(Nacc+elec_Sch(ii,1))*0.5;
        Jaco(ii,ii)=Jaco(ii,ii)-coef_Poi*elec_Sch(ii,1)/thermal*0.5;
    else
        res(ii,1)=res(ii,1)-coef_Poi*(Nacc+elec_Sch(ii,1));
        Jaco(ii,ii)=Jaco(ii,ii)-coef_Poi*elec_Sch(ii,1)/thermal;
    end
end

update=Jaco\(-res);
phi=phi+update;
count(j,1)=count(j,1)+1;

%    elec(ii,1)=ni*exp(q*phi(ii,1)/(k_B*T));
if abs(update(interface1+1:interface2-1,1))<= 5e-16;
   break
end
phi1(:,j)=phi;
elec_Sch1(:,j)=elec_Sch;
end

  integ_elec_Sch(j,1)=sum(elec_Sch*0.1e-9,1)
end
% 


f1=figure
subplot(1,2,1)
semilogy(V_appl,integ_elec*1e-4,'b',V_appl,integ_elec_Sch*1e-4,'r')
xlabel('Gate Voltage (V)')
ylabel('Electron density (cm^2)')
legend('Nonlinear Poisson','Schrodinger- Poisson') 
    
subplot(1,2,2)
plot(V_appl,integ_elec*1e-4,'b',V_appl,integ_elec_Sch*1e-4,'r')
xlabel('Gate Voltage (V)')
ylabel('Electron density (cm^2)')   
legend('Nonlinear Poisson','Schrodinger- Poisson')
    
f2=figure
plot(z,phi1)
xlabel('z (nm)')
ylabel('Potential (V)')   


f3=figure
plot(z,elec1*1e-6)
xlabel('z (nm)')
ylabel('Electron density (cm^-^3)')   

f4=figure
plot(z,elec_Sch1*1e-6)
xlabel('z (nm)')
ylabel('Electron density  (cm^-^3)')   