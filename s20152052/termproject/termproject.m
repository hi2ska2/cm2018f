%% term project

% device mesh

clear all
clc


q = 1.602192e-19; % Elementary charge, C 
eps0 = 8.854187817e-12; % Vacuum permittivity, F/m 
k_B = 1.380662e-23; % Boltzmann constant, J/K 
T = 300.0; % Temperature, K 
thermal = k_B*T/q; % Thermal voltage, V 

ny=25;
nz=25;
Nnode=ny*nz;

Deltay = 5e-9; % 1 nm spacing 
y = Deltay*transpose([0:ny-1]); % real space, m 
y_12 = 9; % At y= 40 nm
y_23 = 17; % At y= 80 nm (25x25)
eps_si = 11.7; eps_ox = 3.9; % Relative permittivity 
Ndon13 = 2e23*ones(ny,1); % 2e17/cm^3 
Ndon2 = 5e25; % 5e19 /cm^3 
ni = 1.075e16; % 1.075e10 /cm^3 
coef = Deltay*Deltay*q/eps0; 
% mob_si=0.14; % Si mobility m^2/Vs
% Dn_Si=thermal*mob_si;
IF=-4.63374; % vacuum level-intrinsic fermi level [eV]
WF=-4.30; % vacuum level-fermi level [eV]
% V_appl=transpose([0:0.1:1]); % applied voltage set: fermi level= 0 V at equilibrium 
V_appl=0;
% V_g(j,1)=-(IF-WF-V_appl(j,1)); % applied gate voltage  [V]
V_g=-(IF-WF-V_appl); % applied gate voltage  [V]
% phi(:,1)=V_g; % initial value
phi=zeros(Nnode,1)+V_g;
A=zeros(Nnode,Nnode);

ydivz=20;
zdivy=0.05;

% ind=ii+ny*(jj-1); %index
for newton=1:10
    
%% Laplacian
for jj=1 % z-axis 
    for ii=1
        res(ii+ny*(jj-1),1)=0.5*eps_ox*(zdivy)*phi(ii+1+ny*(jj-1),1)... % ii+1,jj
                            -0.5*eps_ox*(ydivz+zdivy)*phi(ii+ny*(jj-1),1)... % ii,jj
                            +0.5*eps_ox*(ydivz)*phi(ii+ny*(jj),1); % ii,jj+1
        Jaco(ii+ny*(jj-1),ii+1+ny*(jj-1))=0.5*eps_ox*(zdivy);               % ii+1,jj 
        Jaco(ii+ny*(jj-1),ii+ny*(jj-1))=-0.5*eps_ox*(ydivz+zdivy);      % ii,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj))=0.5*eps_ox*(ydivz);                 % ii,jj+1  
    end
    
    for ii=2:ny-1
        res(ii+ny*(jj-1),1)=0.5*eps_ox*(zdivy)*phi(ii-1+ny*(jj-1))... % ii-1,jj
                            +0.5*eps_ox*(zdivy)*phi(ii+1+ny*(jj-1))... % ii+1,jj
                            -0.5*eps_ox*(ydivz+zdivy)*phi(ii+ny*(jj-1))... % ii,jj
                            +eps_ox*(ydivz)*phi(ii+ny*(jj)); % ii,jj+1
        Jaco(ii+ny*(jj-1),ii-1+ny*(jj-1))=0.5*eps_ox*(zdivy); % ii-1,jj
        Jaco(ii+ny*(jj-1),ii+1+ny*(jj-1))=0.5*eps_ox*(zdivy); % ii+1,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj-1))=-0.5*eps_ox*(ydivz+zdivy); % ii,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj))=eps_ox*(ydivz); % ii,jj+1
    end
        
    for ii=ny
        res(ii+ny*(jj-1),1)=0.5*eps_ox*(zdivy)*phi(ii-1+ny*(jj-1),1)... % ii-1,jj
                            -0.5*eps_ox*(ydivz+zdivy)*phi(ii+ny*(jj-1),1)... % ii,jj
                            +eps_ox*(ydivz)*phi(ii+ny*(jj),1); % ii,jj+1
        Jaco(ii+ny*(jj-1),ii-1+ny*(jj-1))=0.5*eps_ox*(zdivy); 
        Jaco(ii+ny*(jj-1),ii+ny*(jj-1))=-0.5*eps_ox*(ydivz+zdivy);  
        Jaco(ii+ny*(jj-1),ii+ny*(jj))=eps_ox*(ydivz);
    end
            
end
% %%%%%%%%%%%%%%%%%%%%%%%%
for jj=2
    for ii=1
        res(ii+ny*(jj-1),1)=eps_ox*(zdivy)*phi(ii+1+ny*(jj-1),1)... % ii+1,jj
                          -eps_ox*(ydivz+zdivy)*phi(ii+ny*(jj-1),1)... % ii,jj
                          +0.5*eps_ox*(ydivz)*phi(ii+ny*(jj-2),1)... % ii,jj-1
                          +0.5*eps_ox*(ydivz)*phi(ii+ny*(jj),1);  % ii,jj+1
        Jaco(ii+ny*(jj-1),ii+1+ny*(jj-1))=eps_ox*(zdivy); % ii+1,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj-1))=-eps_ox*(ydivz+zdivy);      % ii,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj-2))=0.5*eps_ox*(ydivz);  % ii,jj-1
        Jaco(ii+ny*(jj-1),ii+ny*(jj))=0.5*eps_ox*(ydivz);    % ii,jj+1
    end
    
    for ii=2:ny-1
        res(ii+ny*(jj-1),1)=eps_ox*(zdivy)*phi(ii-1+ny*(jj-1),1)... % ii-1,jj...
                            +eps_ox*(zdivy)*phi(ii+1+ny*(jj-1),1)... % ii+1,jj
                            -2*eps_ox*(ydivz+zdivy)*phi(ii+ny*(jj-1),1)... % ii,jj
                            +eps_ox*(ydivz)*phi(ii+ny*(jj-2),1)... % ii,jj-1
                            +eps_ox*(ydivz)*phi(ii+ny*(jj),1);  % ii,jj+1
        Jaco(ii+ny*(jj-1),ii-1+ny*(jj-1))=eps_ox*(zdivy); % ii-1,jj
        Jaco(ii+ny*(jj-1),ii+1+ny*(jj-1))=eps_ox*(zdivy); % ii+1,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj-1))=-2*eps_ox*(ydivz+zdivy);      % ii,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj-2))=eps_ox*(ydivz);  % ii,jj-1
        Jaco(ii+ny*(jj-1),ii+ny*(jj))=eps_ox*(ydivz);    % ii,jj+1
    end
      
    for ii=ny
        res(ii+ny*(jj-1),1)=eps_ox*zdivy*phi(ii-1+ny*(jj-1),1)... % ii-1,jj
                            -eps_ox*(ydivz+zdivy)*phi(ii+ny*(jj-1),1)... % ii,jj
                            +0.5*eps_ox*(ydivz)*phi(ii+ny*(jj-2),1)... % ii,jj-1
                            +0.5*eps_ox*(ydivz)*phi(ii+ny*(jj),1);  % ii,jj+1
        Jaco(ii+ny*(jj-1),ii-1+ny*(jj-1))=eps_ox*zdivy; % ii-1,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj-1))=-eps_ox*(ydivz+zdivy);      % ii,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj-2))=0.5*eps_ox*(ydivz);  % ii,jj-1
        Jaco(ii+ny*(jj-1),ii+ny*(jj))=0.5*eps_ox*(ydivz);    % ii,jj+1
    end  
end

for jj=3  %% interface1 SiO -> Si
    for ii=1
        res(ii+ny*(jj-1),1)=0.5*(eps_ox+eps_si)*(zdivy)*phi(ii+1+ny*(jj-1),1)... % ii+1,jj
                          -0.5*(eps_ox+eps_si)*(ydivz+zdivy)*phi(ii+ny*(jj-1),1)... % ii,jj
                          +0.5*eps_ox*(ydivz)*phi(ii+ny*(jj-2),1)... % ii,jj-1
                          +0.5*eps_si*(ydivz)*phi(ii+ny*(jj),1);  % ii,jj+1
        Jaco(ii+ny*(jj-1),ii+1+ny*(jj-1))=0.5*(eps_ox+eps_si)*(zdivy); % ii+1,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj-1))=-0.5*(eps_ox+eps_si)*(ydivz+zdivy);      % ii,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj-2))=0.5*eps_ox*(ydivz);  % ii,jj-1
        Jaco(ii+ny*(jj-1),ii+ny*(jj))=0.5*eps_si*(ydivz);    % ii,jj+1
    end
    
    for ii=2:ny-1
        res(ii+ny*(jj-1),1)=0.5*(eps_ox+eps_si)*(zdivy)*phi(ii-1+ny*(jj-1),1)... % ii-1,jj...
                            +0.5*(eps_ox+eps_si)*(zdivy)*phi(ii+1+ny*(jj-1),1)... % ii+1,jj
                            -(eps_ox+eps_si)*(ydivz+zdivy)*phi(ii+ny*(jj-1),1)... % ii,jj
                            +eps_ox*(ydivz)*phi(ii+ny*(jj-2),1)... % ii,jj-1
                            +eps_si*(ydivz)*phi(ii+ny*(jj),1);  % ii,jj+1
        Jaco(ii+ny*(jj-1),ii-1+ny*(jj-1))=0.5*(eps_ox+eps_si)*(zdivy); % ii-1,jj
        Jaco(ii+ny*(jj-1),ii+1+ny*(jj-1))=0.5*(eps_ox+eps_si)*(zdivy); % ii+1,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj-1))=-(eps_ox+eps_si)*(ydivz+zdivy);      % ii,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj-2))=eps_ox*(ydivz);  % ii,jj-1
        Jaco(ii+ny*(jj-1),ii+ny*(jj))=eps_si*(ydivz);    % ii,jj+1
    end
            
    for ii=ny
        res(ii+ny*(jj-1),1)=0.5*(eps_ox+eps_si)*zdivy*phi(ii-1+ny*(jj-1),1)... % ii-1,jj
                            -0.5*(eps_ox+eps_si)*(ydivz+zdivy)*phi(ii+ny*(jj-1),1)... % ii,jj
                            +0.5*eps_ox*(ydivz)*phi(ii+ny*(jj-2),1)... % ii,jj-1
                            +0.5*eps_si*(ydivz)*phi(ii+ny*(jj),1);  % ii,jj+1
        Jaco(ii+ny*(jj-1),ii-1+ny*(jj-1))=0.5*(eps_ox+eps_si)*zdivy; % ii-1,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj-1))=-0.5*(eps_ox+eps_si)*(ydivz+zdivy);      % ii,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj-2))=+0.5*eps_ox*(ydivz)*phi(ii+ny*(jj-2),1);  % ii,jj-1
        Jaco(ii+ny*(jj-1),ii+ny*(jj))=+0.5*eps_si*(ydivz)*phi(ii+ny*(jj),1);    % ii,jj+1
    end  
end
    

for jj=4:22
    for ii=1
        res(ii+ny*(jj-1),1)=eps_si*(zdivy)*phi(ii+1+ny*(jj-1),1)... % ii+1,jj
                          -eps_si*(ydivz+zdivy)*phi(ii+ny*(jj-1),1)... % ii,jj
                          +0.5*eps_si*(ydivz)*phi(ii+ny*(jj-2),1)... % ii,jj-1
                          +0.5*eps_si*(ydivz)*phi(ii+ny*(jj),1);  % ii,jj+1
        Jaco(ii+ny*(jj-1),ii+1+ny*(jj-1))=eps_si*(zdivy); % ii+1,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj-1))=-eps_si*(ydivz+zdivy);      % ii,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj-2))=0.5*eps_si*(ydivz);  % ii,jj-1
        Jaco(ii+ny*(jj-1),ii+ny*(jj))=0.5*eps_si*(ydivz);    % ii,jj+1
    end
    
    for ii=2:ny-1
        res(ii+ny*(jj-1),1)=eps_si*(zdivy)*phi(ii-1+ny*(jj-1),1)... % ii-1,jj...
                            +eps_si*(zdivy)*phi(ii+1+ny*(jj-1),1)... % ii+1,jj
                            -2*eps_si*(ydivz+zdivy)*phi(ii+ny*(jj-1),1)... % ii,jj
                            +eps_si*(ydivz)*phi(ii+ny*(jj-2),1)... % ii,jj-1
                            +eps_si*(ydivz)*phi(ii+ny*(jj),1);  % ii,jj+1
        Jaco(ii+ny*(jj-1),ii-1+ny*(jj-1))=eps_si*(zdivy); % ii-1,jj
        Jaco(ii+ny*(jj-1),ii+1+ny*(jj-1))=eps_si*(zdivy); % ii+1,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj-1))=-2*eps_si*(ydivz+zdivy);      % ii,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj-2))=eps_si*(ydivz);  % ii,jj-1
        Jaco(ii+ny*(jj-1),ii+ny*(jj))=eps_si*(ydivz);    % ii,jj+1
    end
      
    for ii=ny
        res(ii+ny*(jj-1),1)=eps_si*zdivy*phi(ii-1+ny*(jj-1),1)... % ii-1,jj
                            -eps_si*(ydivz+zdivy)*phi(ii+ny*(jj-1),1)... % ii,jj
                            +0.5*eps_si*(ydivz)*phi(ii+ny*(jj-2),1)... % ii,jj-1
                            +0.5*eps_si*(ydivz)*phi(ii+ny*(jj),1);  % ii,jj+1
        Jaco(ii+ny*(jj-1),ii-1+ny*(jj-1))=eps_si*zdivy; % ii-1,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj-1))=-eps_si*(ydivz+zdivy);      % ii,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj-2))=0.5*eps_si*(ydivz);  % ii,jj-1
        Jaco(ii+ny*(jj-1),ii+ny*(jj))=0.5*eps_si*(ydivz);    % ii,jj+1
    end  
end

for jj=23  %% interface1 SiO -> Si
    for ii=1
        res(ii+ny*(jj-1),1)=0.5*(eps_ox+eps_si)*(zdivy)*phi(ii+1+ny*(jj-1),1)... % ii+1,jj
                          -0.5*(eps_ox+eps_si)*(ydivz+zdivy)*phi(ii+ny*(jj-1),1)... % ii,jj
                          +0.5*eps_si*(ydivz)*phi(ii+ny*(jj-2),1)... % ii,jj-1
                          +0.5*eps_ox*(ydivz)*phi(ii+ny*(jj),1);  % ii,jj+1
        Jaco(ii+ny*(jj-1),ii+1+ny*(jj-1))=0.5*(eps_ox+eps_si)*(zdivy); % ii+1,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj-1))=-0.5*(eps_ox+eps_si)*(ydivz+zdivy);      % ii,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj-2))=0.5*eps_si*(ydivz);  % ii,jj-1
        Jaco(ii+ny*(jj-1),ii+ny*(jj))=0.5*eps_ox*(ydivz);    % ii,jj+1
    end
    
    for ii=2:ny-1
        res(ii+ny*(jj-1),1)=0.5*(eps_ox+eps_si)*(zdivy)*phi(ii-1+ny*(jj-1),1)... % ii-1,jj...
                            +0.5*(eps_ox+eps_si)*(zdivy)*phi(ii+1+ny*(jj-1),1)... % ii+1,jj
                            -(eps_ox+eps_si)*(ydivz+zdivy)*phi(ii+ny*(jj-1),1)... % ii,jj
                            +eps_si*(ydivz)*phi(ii+ny*(jj-2),1)... % ii,jj-1
                            +eps_ox*(ydivz)*phi(ii+ny*(jj),1);  % ii,jj+1
        Jaco(ii+ny*(jj-1),ii-1+ny*(jj-1))=0.5*(eps_ox+eps_si)*(zdivy); % ii-1,jj
        Jaco(ii+ny*(jj-1),ii+1+ny*(jj-1))=0.5*(eps_ox+eps_si)*(zdivy); % ii+1,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj-1))=-(eps_ox+eps_si)*(ydivz+zdivy);      % ii,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj-2))=eps_si*(ydivz);  % ii,jj-1
        Jaco(ii+ny*(jj-1),ii+ny*(jj))=eps_ox*(ydivz);    % ii,jj+1
    end
            
    for ii=ny
        res(ii+ny*(jj-1),1)=0.5*(eps_ox+eps_si)*zdivy*phi(ii-1+ny*(jj-1),1)... % ii-1,jj
                            -0.5*(eps_ox+eps_si)*(ydivz+zdivy)*phi(ii+ny*(jj-1),1)... % ii,jj
                            +0.5*eps_si*(ydivz)*phi(ii+ny*(jj-2),1)... % ii,jj-1
                            +0.5*eps_ox*(ydivz)*phi(ii+ny*(jj),1);  % ii,jj+1
        Jaco(ii+ny*(jj-1),ii-1+ny*(jj-1))=0.5*(eps_ox+eps_si)*zdivy; % ii-1,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj-1))=-0.5*(eps_ox+eps_si)*(ydivz+zdivy);      % ii,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj-2))=+0.5*eps_si*(ydivz)*phi(ii+ny*(jj-2),1);  % ii,jj-1
        Jaco(ii+ny*(jj-1),ii+ny*(jj))=+0.5*eps_ox*(ydivz)*phi(ii+ny*(jj),1);    % ii,jj+1
    end  
end

for jj=24
    for ii=1
        res(ii+ny*(jj-1),1)=eps_ox*(zdivy)*phi(ii+1+ny*(jj-1),1)... % ii+1,jj
                          -eps_ox*(ydivz+zdivy)*phi(ii+ny*(jj-1),1)... % ii,jj
                          +0.5*eps_ox*(ydivz)*phi(ii+ny*(jj-2),1)... % ii,jj-1
                          +0.5*eps_ox*(ydivz)*phi(ii+ny*(jj),1);  % ii,jj+1
        Jaco(ii+ny*(jj-1),ii+1+ny*(jj-1))=eps_ox*(zdivy); % ii+1,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj-1))=-eps_ox*(ydivz+zdivy);      % ii,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj-2))=0.5*eps_ox*(ydivz);  % ii,jj-1
        Jaco(ii+ny*(jj-1),ii+ny*(jj))=0.5*eps_ox*(ydivz);    % ii,jj+1
    end
    
    for ii=2:ny-1
        res(ii+ny*(jj-1),1)=eps_ox*(zdivy)*phi(ii-1+ny*(jj-1),1)... % ii-1,jj...
                            +eps_ox*(zdivy)*phi(ii+1+ny*(jj-1),1)... % ii+1,jj
                            -2*eps_ox*(ydivz+zdivy)*phi(ii+ny*(jj-1),1)... % ii,jj
                            +eps_ox*(ydivz)*phi(ii+ny*(jj-2),1)... % ii,jj-1
                            +eps_ox*(ydivz)*phi(ii+ny*(jj),1);  % ii,jj+1
        Jaco(ii+ny*(jj-1),ii-1+ny*(jj-1))=eps_ox*(zdivy); % ii-1,jj
        Jaco(ii+ny*(jj-1),ii+1+ny*(jj-1))=eps_ox*(zdivy); % ii+1,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj-1))=-2*eps_ox*(ydivz+zdivy);      % ii,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj-2))=eps_ox*(ydivz);  % ii,jj-1
        Jaco(ii+ny*(jj-1),ii+ny*(jj))=eps_ox*(ydivz);    % ii,jj+1
    end
      
    for ii=ny
        res(ii+ny*(jj-1),1)=eps_ox*zdivy*phi(ii-1+ny*(jj-1),1)... % ii-1,jj
                            -eps_ox*(ydivz+zdivy)*phi(ii+ny*(jj-1),1)... % ii,jj
                            +0.5*eps_ox*(ydivz)*phi(ii+ny*(jj-2),1)... % ii,jj-1
                            +0.5*eps_ox*(ydivz)*phi(ii+ny*(jj),1);  % ii,jj+1
        Jaco(ii+ny*(jj-1),ii-1+ny*(jj-1))=eps_ox*zdivy; % ii-1,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj-1))=-eps_ox*(ydivz+zdivy);      % ii,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj-2))=0.5*eps_ox*(ydivz);  % ii,jj-1
        Jaco(ii+ny*(jj-1),ii+ny*(jj))=0.5*eps_ox*(ydivz);    % ii,jj+1
    end  
end

for jj=25 % z-axis 
    for ii=1
        res(ii+ny*(jj-1),1)=0.5*eps_ox*(zdivy)*phi(ii+1+ny*(jj-1),1)... % ii+1,jj
                            -0.5*eps_ox*(ydivz+zdivy)*phi(ii+ny*(jj-1),1)... % ii,jj
                            +0.5*eps_ox*(ydivz)*phi(ii+ny*(jj-2),1); % ii,jj-1
        Jaco(ii+ny*(jj-1),ii+1+ny*(jj-1))=0.5*eps_ox*(zdivy);               % ii+1,jj 
        Jaco(ii+ny*(jj-1),ii+ny*(jj-1))=-0.5*eps_ox*(ydivz+zdivy);      % ii,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj-2))=0.5*eps_ox*(ydivz);                 % ii,jj-1  
    end
    
    for ii=2:ny-1
        res(ii+ny*(jj-1),1)=0.5*eps_ox*(zdivy)*phi(ii-1+ny*(jj-1))... % ii-1,jj
                            +0.5*eps_ox*(zdivy)*phi(ii+1+ny*(jj-1))... % ii+1,jj
                            -0.5*eps_ox*(ydivz+zdivy)*phi(ii+ny*(jj-1))... % ii,jj
                            +eps_ox*(ydivz)*phi(ii+ny*(jj-2)); % ii,jj-1
        Jaco(ii+ny*(jj-1),ii-1+ny*(jj-1))=0.5*eps_ox*(zdivy); % ii-1,jj
        Jaco(ii+ny*(jj-1),ii+1+ny*(jj-1))=0.5*eps_ox*(zdivy); % ii+1,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj-1))=-0.5*eps_ox*(ydivz+zdivy); % ii,jj
        Jaco(ii+ny*(jj-1),ii+ny*(jj-2))=eps_ox*(ydivz); % ii,jj-1
    end
        
    for ii=ny
        res(ii+ny*(jj-1),1)=0.5*eps_ox*(zdivy)*phi(ii-1+ny*(jj-1),1)... % ii-1,jj
                            -0.5*eps_ox*(ydivz+zdivy)*phi(ii+ny*(jj-1),1)... % ii,jj
                            +eps_ox*(ydivz)*phi(ii+ny*(jj-2),1); % ii,jj-1
        Jaco(ii+ny*(jj-1),ii-1+ny*(jj-1))=0.5*eps_ox*(zdivy); 
        Jaco(ii+ny*(jj-1),ii+ny*(jj-1))=-0.5*eps_ox*(ydivz+zdivy);  
        Jaco(ii+ny*(jj-1),ii+ny*(jj-2))=eps_ox*(ydivz);
    end
            
end

%%
for jj=1
    for ii=y_12:y_23
        res(ii+ny*(jj-1),1)=phi(ii+ny*(jj-1),1)-V_g
        Jaco(ii+ny*(jj-1),ii+ny*(jj-1))=1;
    end
end

for jj=nz
    for ii=y_12:y_23
        res(ii+ny*(jj-1),1)=phi(ii+ny*(jj-1),1)-V_g
        Jaco(ii+ny*(jj-1),ii+ny*(jj-1))=1;
    end
end

    
for jj=3:23
    if (ii==1)
        res(ii+ny*(jj-1),1)=res(ii+ny*(jj-1),1)-coef*(-Ndon13+ni*exp(phi(ii+ny*(jj-1),1)/thermal))*0.5;
        Jaco(ii+ny*(jj-1),ii+ny*(jj-1))=Jaco(ii+ny*(jj-1),ii+ny*(jj-1))-coef*ni*exp(phi(ii+ny*(jj-1),1)/thermal)/termal*0.5;
    elseif (ii==y_12)
        res(ii+ny*(jj-1),1)=res(ii+ny*(jj-1),1)-coef*(-Ndon13+ni*exp(phi(ii+ny*(jj-1),1)/thermal))*0.5;
    for ii=1:y_12
        n=coef*(ni*exp(phi(ii,1)/thermal))