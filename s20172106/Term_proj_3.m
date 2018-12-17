
clc;

q=1.602192e-19; % elementary charge, [C]
eps0=8.854187817e-12; % vacuum permittivity, [F/m]
k_B=1.380662e-23; % Boltzmann constant, [J/K]
T=300.0; % temperature, [K]
eps_si=11.7; eps_ox=3.9; % relatvie permittivity
thermal = k_B*T/q; % Thermal voltage, V
mu = 1430; % [cm/V*s] mobility
D_n = thermal*mu; % [] Diffusion constant by Einstein relation

Nacc_h = 5e25; % 1/m^3 5e19 [1/cm^3] Doping density [high]
Nacc_l = 2e23; % 1/m^3 2e17 [1/cm^3] Doping density [low]

ni=1.075e16; %1.075e10/cm^3 
Work_Ftn=-4.30; % vacuum level-fermi level [eV]
Int_Fermi=-4.63374; % vacuum level-intrinsic fermi level for oxide [eV]
% Volt_Gate = [0:0.1:1]'; % set: fermi level= 0 V at equilibrium

Volt_gate = 0:0.1:1; % [V]
Volt_source = 0.0; % [V]
Volt_drain = 0.0; % [V]
x_fullsize = 120 ; % [nm]
y_fullsize = 6 ; %[nm]
nx = 60; 
ny = 60;
nx = nx + 1;
ny = ny + 1;

interface_x_1 = ((nx-1)/3) +1; 
interface_x_2 = (2*(nx-1)/3) +1;

interface_y_1 = 6; 
interface_y_2 = 56;

N = nx*ny;

delta_x = x_fullsize/(nx-1)*1e-9; % 2 nm
delta_y = 0.1e-9; % 0.1 nm
ratio = delta_x/delta_y;
ratio_sq = ratio^2;

len_gate = length(Volt_gate);
b_gate = zeros(len_gate,1);
b_gate(:) = ((Work_Ftn - Volt_gate(:))- Int_Fermi); % Charge density under Gate_voltage
b_source = ((Work_Ftn - Volt_source)- Int_Fermi); % Charge density under Source_voltage
b_drain = ((Work_Ftn - Volt_drain)- Int_Fermi); % Charge density under Drain_voltage


coef = delta_x*delta_x*q/eps0;
phi_gate = zeros(nx,ny,len_gate);
elec_gate = zeros(nx,ny,len_gate);

for ss =  1: len_gate
    b_gate_state = b_gate(ss,1);

    phi = zeros(N,1);
% phi(:,1) = 

for newton = 1:10;
   

            res = zeros(N,1);
            Jaco = zeros(N,N);
        
    %% Bulk part
        
    % Laplacian part
        
        for ii = 1:nx
            for jj = 1:ny
                index = ii + nx*(jj-1);
                
                % Coner boundary
                    
                if ( (ii == 1 && jj == 1) || (ii == nx && jj == 1) || (ii == 1 && jj == ny) || (ii == nx && jj == ny) )
                    if (ii == 1 && jj == 1) % Left_bottom in Oxide
                    
                        res(index,1) = -(eps_ox + ratio_sq*eps_ox)*phi(index,1) + ratio_sq*eps_ox*phi(index + nx,1) + eps_ox*phi(index +1,1);
                        Jaco(index,index) = -(eps_ox + ratio_sq*eps_ox);
                        Jaco(index,index + 1) = eps_ox;
                        Jaco(index,index + nx) = ratio_sq*eps_ox;
                    
                    elseif (ii == nx && jj == 1 ) % Right_bottom in Oxide
                        
                        res(index,1) = -(eps_ox + ratio_sq*eps_ox)*phi(index,1) + ratio_sq*eps_ox*phi(index + nx,1) + eps_ox*phi(index -1,1);
                        Jaco(index,index) = -(eps_ox + ratio_sq*eps_ox);
                        Jaco(index,index - 1) = eps_ox;
                        Jaco(index,index + nx) = ratio_sq*eps_ox;
                        
                    elseif (ii == 1 && jj == ny ) % Left_Top in Oxide
                        
                        res(index,1) = -(eps_ox + ratio_sq*eps_ox)*phi(index,1) + ratio_sq*eps_ox*phi(index - nx,1) + eps_ox*phi(index + 1,1);
                        Jaco(index,index) = -(eps_ox + ratio_sq*eps_ox);
                        Jaco(index,index + 1) = eps_ox;
                        Jaco(index,index - nx) = ratio_sq*eps_ox;
                    
                    elseif ( ii == nx && jj == ny ) % Right_Top in Oxide
                        
                        res(index,1) = -(eps_ox + ratio_sq*eps_ox)*phi(index,1) + ratio_sq*eps_ox*phi(index - nx,1) + eps_ox*phi(index -1,1);
                        Jaco(index,index) = -(eps_ox + ratio_sq*eps_ox);
                        Jaco(index,index - 1) = eps_ox;
                        Jaco(index,index - nx) = ratio_sq*eps_ox;
                  
                    end
                
                
                % Side boundaries
                
                elseif ( ii == 1 || jj == 1 || ii == nx || jj == ny) 
                    
                    if ( jj ==  1) % Bottom side
                        
                        if ( ( index < interface_x_1 || index > interface_x_2 ) && ( index ~= 1 || index ~= nx ) ) % no gate
                             
                            res(index, 1) = 0.5*eps_ox*phi(index - 1,1) -(1.0*eps_ox + 1.0*ratio_sq*eps_ox)*phi(index,1) + 0.5*eps_ox*phi(index + 1,1) + 1.0*ratio_sq*eps_ox*phi(index + nx,1);
                            Jaco(index,index - 1) = 0.5*eps_ox; % x = -1
                            Jaco(index,index) =  -(1.0*eps_ox + 1.0*ratio_sq*eps_ox); % x,y = 0
                            Jaco(index,index + 1) = 0.5*eps_ox; % x = +1
                            Jaco(index,index + nx) = 1.0*ratio_sq*eps_ox; % y = + 1
                            
                            
                        elseif ( index >= interface_x_1 && index <= interface_x_2 ) % Bottom gate
                            
                            res(index,1) = phi(index,1) + b_gate_state;
                            Jaco(index,index) = 1.0 ;
                        end
                        
                    elseif ( jj == ny ) %  Top side
                        
                        if ( ( index < nx*( ny -1 ) + interface_x_1 || index > nx*( ny -1 ) + interface_x_2 ) && ( index ~= nx*(ny -1) +1 || index ~= nx*ny ) )
                            
                            res(index,1) = 0.5*eps_ox*phi(index - 1,1) -(1.0*eps_ox + 1.0*ratio_sq*eps_ox)*phi(index,1) + 0.5*eps_ox*phi(index + 1,1) + 1.0*ratio_sq*eps_ox*phi(index -nx,1 );
                            Jaco(index,index-1) = 0.5*eps_ox;
                            Jaco(index,index) = -(1.0*eps_ox + 1.0*ratio_sq*eps_ox);
                            Jaco(index,index+1) = 0.5*eps_ox;
                            Jaco(index,index-nx) = 1.0*ratio_sq*eps_ox;
                            
                        elseif ( index >= nx*( ny -1 ) + interface_x_1  && index <= nx*( ny -1 ) + interface_x_2 ) % Top gate
                            
                            res(index,1) = phi(index,1) + b_gate_state;
                            Jaco(index,index) = 1.0 ;
                        end
                        
                    elseif ( ii == 1 ) % Left side
                        
                        if ( jj < interface_y_1 && jj ~= 1 ) % Ox
                            
                            res(index,1) = 1.0*eps_ox*phi(index + 1,1)  + 0.5*ratio_sq*eps_ox*phi(index - nx,1) -(1.0*eps_ox + 1.0*ratio_sq*eps_ox)*phi(index,1) +  0.5*ratio_sq*eps_ox*phi(index + nx,1);
                            Jaco(index,index+1) = 1.0*eps_ox;
                            Jaco(index,index - nx) = 0.5*ratio_sq*eps_ox;
                            Jaco(index,index) = -(1.0*eps_ox + 1.0*ratio_sq*eps_ox);
                            Jaco(index,index+nx) = 0.5*ratio_sq*eps_ox;
                            
%                         
%                         elseif ( jj == interface_y_1 )
%                             
%                             res(index,1) = 1.0*(eps_ox + eps_si)*0.5*phi(index+1,1) + 0.5*ratio_sq*eps_si*phi(index + nx,1) -(1.0*(eps_ox + eps_si)*0.5 + 1.0*ratio_sq*(eps_ox + eps_si)*0.5)*phi(index,1) + 0.5*ratio_sq*eps_ox*phi(index - nx,1);
%                             Jaco(index,index + 1) = 1.0*(eps_ox + eps_si)*0.5; % x = +1
%                             Jaco(index,index - nx) = 0.5*ratio_sq*eps_ox; % y = -1
%                             Jaco(index,index) = -(1.0*(eps_ox + eps_si)*0.5 + 1.0*ratio_sq*(eps_ox + eps_si)*0.5); % x,y = 0,0
%                             Jaco(index,index + nx) = 0.5*ratio_sq*eps_si; % y = +1
%                             
%                             
                        elseif ( jj >= interface_y_1 && jj <= interface_y_2 ) % Si  % Left source 
                            
                            res(index,1) = phi(index,1) + b_source;
                            Jaco(index,index) = 1.0 ;
                            
%                         elseif ( jj == interface_y_2 )
%                             
%                             res(index,1) = 1.0*(eps_ox + eps_si)*0.5*phi(index+1,1) + 0.5*ratio_sq*eps_ox*phi(index + nx,1) -(1.0*(eps_ox + eps_si)*0.5 + 1.0*ratio_sq*(eps_ox + eps_si)*0.5)*phi(index,1) + 0.5*ratio_sq*eps_si*phi(index - nx,1);
%                             Jaco(index,index + 1) = 1.0*(eps_ox + eps_si)*0.5; % x = +1
%                             Jaco(index,index - nx) = 0.5*ratio_sq*eps_si; % y = -1
%                             Jaco(index,index) = -(1.0*(eps_ox + eps_si)*0.5 + 1.0*ratio_sq*(eps_ox + eps_si)*0.5); % x,y = 0,0
%                             Jaco(index,index + nx) = 0.5*ratio_sq*eps_ox; % y = +1
                            
                        elseif ( jj > interface_y_2 && jj ~= ny) % Ox
                            
                             res(index,1) = 1.0*eps_ox*phi(index + 1,1)  + 0.5*ratio_sq*eps_ox*phi(index - nx,1) -(1.0*eps_ox + 1.0*ratio_sq*eps_ox)*phi(index,1) +  0.5*ratio_sq*eps_ox*phi(index + nx,1);
                            Jaco(index,index+1) = 1.0*eps_ox;
                            Jaco(index,index - nx) = 0.5*ratio_sq*eps_ox;
                            Jaco(index,index) = -(1.0*eps_ox + 1.0*ratio_sq*eps_ox);
                            Jaco(index,index + nx) = 0.5*ratio_sq*eps_ox;
                            
                            
                        end
                        
                    elseif ( ii == nx ) % right_side
                        
                        if ( jj < interface_y_1 && jj ~= 1 ) % Ox
                            
                            res(index,1) = 1.0*eps_ox*phi(index - 1,1)  + 0.5*ratio_sq*eps_ox*phi(index - nx,1) -(1.0*eps_ox + 1.0*ratio_sq*eps_ox)*phi(index,1) +  0.5*ratio_sq*eps_ox*phi(index + nx,1);
                            Jaco(index,index-1) = 1.0*eps_ox;
                            Jaco(index,index - nx) = 0.5*ratio_sq*eps_ox;
                            Jaco(index,index) = -(1.0*eps_ox + 1.0*ratio_sq*eps_ox);
                            Jaco(index,index+nx) = 0.5*ratio_sq*eps_ox;
                            
                         
%                          
%                         elseif ( jj == interface_y_1 )
%                             
%                             res(index,1) = 1.0*(eps_ox + eps_si)*0.5*phi(index-1,1) + 0.5*ratio_sq*eps_si*phi(index + nx,1) -(1.0*(eps_ox + eps_si)*0.5 + 1.0*ratio_sq*(eps_ox + eps_si)*0.5)*phi(index,1) + 0.5*ratio_sq*eps_ox*phi(index - nx,1);
%                             Jaco(index,index - 1) = 1.0*(eps_ox + eps_si)*0.5; % x = +1
%                             Jaco(index,index - nx) = 0.5*ratio_sq*eps_ox; % y = -1
%                             Jaco(index,index) = -(1.0*(eps_ox + eps_si)*0.5 + 1.0*ratio_sq*(eps_ox + eps_si)*0.5); % x,y = 0,0
%                             Jaco(index,index + nx) = 0.5*ratio_sq*eps_si; % y = +1
                            
                            
                        elseif ( jj >= interface_y_1 && jj <= interface_y_2 ) % Si % Right Drain
                            
                             res(index,1) = phi(index,1) + b_drain;
                            Jaco(index,index) = 1.0 ;
                            
%                         elseif ( jj == interface_y_2 )
%                             
%                             res(index,1) = 1.0*(eps_ox + eps_si)*0.5*phi(index - 1,1) + 0.5*ratio_sq*eps_ox*phi(index + nx,1) -(1.0*(eps_ox + eps_si)*0.5 + 1.0*ratio_sq*(eps_ox + eps_si)*0.5)*phi(index,1) + 0.5*ratio_sq*eps_si*phi(index - nx,1);
%                             Jaco(index,index - 1) = 1.0*(eps_ox + eps_si)*0.5; % x = +1
%                             Jaco(index,index - nx) = 0.5*ratio_sq*eps_si; % y = -1
%                             Jaco(index,index) = -(1.0*(eps_ox + eps_si)*0.5 + 1.0*ratio_sq*(eps_ox + eps_si)*0.5); % x,y = 0,0
%                             Jaco(index,index + nx) = 0.5*ratio_sq*eps_ox; % y = +1
                            
                        elseif ( jj > interface_y_2 && jj ~= ny) % Ox
                            
                             res(index,1) = 1.0*eps_ox*phi(index - 1,1)  + 0.5*ratio_sq*eps_ox*phi(index - nx,1) -(1.0*eps_ox + 1.0*ratio_sq*eps_ox)*phi(index,1) +  0.5*ratio_sq*eps_ox*phi(index + nx,1);
                            Jaco(index,index - 1) = 1.0*eps_ox;
                            Jaco(index,index - nx) = 0.5*ratio_sq*eps_ox;
                            Jaco(index,index) = -(1.0*eps_ox + 1.0*ratio_sq*eps_ox);
                            Jaco(index,index+nx) = 0.5*ratio_sq*eps_ox;
                            
                        end
                    end
                        
                elseif  (( ii >= 2 && ii <= nx -1 ) && ( jj >= 2 && jj <= ny-1 )) % bulk terms
                    
                    if ( jj < interface_y_1 ) % Oxide _ bottom

                        res(index,1) = ratio_sq*eps_ox*phi(index - nx,1)  + 1.0*eps_ox*phi(index - 1,1) -(2.0*eps_ox + ratio_sq*2*eps_ox)*phi(index, 1) + 1.0*eps_ox*phi(index +1,1) + ratio_sq*eps_ox*phi(index+nx,1);
                        Jaco(index,index - nx) =  ratio_sq*eps_ox;
                        Jaco(index,index - 1) = 1.0*eps_ox;
                        Jaco(index,index) =  -(2.0*eps_ox + ratio_sq*2*eps_ox);
                        Jaco(index,index + 1) = 1.0*eps_ox;
                        Jaco(index,index + nx) = ratio_sq*eps_ox;

                    elseif ( jj == interface_y_1 )  % interface_y_1

                        res(index,1) = ratio_sq*eps_ox*phi(index - nx,1)  + 1.0*(eps_si + eps_ox)/2*phi(index - 1,1) -(2.0*(eps_si + eps_ox)/2 + ratio_sq*2*(eps_si + eps_ox)/2)*phi(index, 1) + 1.0*(eps_si + eps_ox)/2*phi(index +1,1) + ratio_sq*eps_si*phi(index+nx,1);
                        Jaco(index ,index - nx) =  ratio_sq*eps_ox;
                        Jaco(index ,index - 1) = 1.0*(eps_si + eps_ox)/2;
                        Jaco(index ,index) =  -2.0*(eps_si + eps_ox)/2 - ratio_sq*2*(eps_si + eps_ox)/2;
                        Jaco(index ,index + 1) = 1.0*(eps_si + eps_ox)/2; % x = +1
                        Jaco(index ,index + nx) = ratio_sq*eps_si; % y = +1

                    elseif ( jj > interface_y_1 && jj < interface_y_2 ) % Si

                        res(index,1) = ratio_sq*eps_si*phi(index - nx,1)  + 1.0*eps_si*phi(index - 1,1) -(2.0*eps_si + ratio_sq*2*eps_si)*phi(index, 1) + 1.0*eps_si*phi(index +1,1) + ratio_sq*eps_si*phi(index+nx,1);
                        Jaco(index,index - nx) = ratio_sq*eps_si; % y = -1
                        Jaco(index,index - 1) = 1.0*eps_si; % x = -1
                        Jaco(index,index) = -2.0*eps_si -ratio_sq*2*eps_si; % x,y = 0,0
                        Jaco(index,index + 1) = 1.0*eps_si; % x = +1
                        Jaco(index,index + nx) = ratio_sq*eps_si; % y = +1

                    elseif ( jj == interface_y_2 ) % interface_y_2

                        res(index,1) = ratio_sq*eps_si*phi(index - nx,1)  + 1.0*(eps_si + eps_ox)/2*phi(index - 1,1) -(2.0*(eps_si + eps_ox)/2 + ratio_sq*2*(eps_si + eps_ox)/2)*phi(index, 1) + 1.0*(eps_si + eps_ox)/2*phi(index +1,1) + ratio_sq*eps_ox*phi(index+nx,1);
                        Jaco(index,index - nx) =  ratio_sq*eps_si;
                        Jaco(index,index - 1) = 1.0*(eps_si + eps_ox)/2;
                        Jaco(index,index) =  -2.0*(eps_si + eps_ox)/2 - ratio_sq*2*(eps_si + eps_ox)/2;
                        Jaco(index,index + 1) = 1.0*(eps_si + eps_ox)/2; % x = +1
                        Jaco(index,index + nx) = ratio_sq*eps_ox; % y = +1

                    elseif ( jj > interface_y_2 ) % Oxide _ Top

                        res(index,1) = ratio_sq*eps_ox*phi(index - nx,1)  + 1.0*eps_ox*phi(index - 1,1) -(2.0*eps_ox + ratio_sq*2*eps_ox)*phi(index, 1) + 1.0*eps_ox*phi(index +1,1) + ratio_sq*eps_ox*phi(index+nx,1);
                        Jaco(index,index - nx) =  ratio_sq*eps_ox;
                        Jaco(index,index - 1) = 1.0*eps_ox;
                        Jaco(index ,index) =  -(2.0*eps_ox + ratio_sq*2*eps_ox);
                        Jaco(index,index + 1) = 1.0*eps_ox;
                        Jaco(index,index + nx) = ratio_sq*eps_ox;
                        
                    end
                end
            end
        end
        
        %%%%%% If there are charges (==> Poisson's eqn) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
        for ii = 2 : nx-1
            for jj = interface_y_1:interface_y_2
                index = ii + nx*(jj - 1);
                
                if ( ii < interface_x_1 )
                    res(index,1) = res(index,1) - coef*(Nacc_h  + ni*exp(phi(index,1)/thermal));
                    Jaco(index,index) =  Jaco(index,index) - coef*ni*exp(phi(index,1)/thermal)/thermal;
                elseif ( ii == interface_x_1 )
                    res(index,1) = res(index,1) - coef*(Nacc_h*0.5 + ni*exp(phi(index,1)/thermal));
                    Jaco(index,index) = Jaco(index,index) - coef*ni*exp(phi(index,1)/thermal)/thermal;
                elseif ( ii > interface_x_1 && ii < interface_x_2 )
                    res(index,1) = res(index,1) - coef*(Nacc_l*0.5  + ni*exp(phi(index,1)/thermal));
                    Jaco(index,index) =  Jaco(index,index) - coef*ni*exp(phi(index,1)/thermal)/thermal;
                elseif ( ii == interface_x_2 )
                    res(index,1) = res(index,1) - coef*(Nacc_h*0.5 + ni*exp(phi(index,1)/thermal));
                    Jaco(index,index) = Jaco(index,index) - coef*ni*exp(phi(index,1)/thermal)/thermal;    
                elseif ( ii > interface_x_2 )
                    res(index,1) = res(index,1) - coef*(Nacc_h  + ni*exp(phi(index,1)/thermal));
                    Jaco(index,index) =  Jaco(index,index) - coef*ni*exp(phi(index,1)/thermal)/thermal;
                end
            end
        end
        
        
        update = Jaco \ (-res);
        phi_1(:,1) = phi(:,1) + update;
                
end
        

      
res_real = zeros(nx,ny);
    for ii = 1:nx
        for jj = 1:ny
            index = ii + nx*(jj-1);
            res_real(ii,jj) = res(index,1);
        end
    end
   
 %% Potential    

 phi_real = zeros(nx,ny);
    for ii = 1:nx
        for jj = 1:ny
            index = ii + nx*(jj-1);
            phi_real(ii,jj) = phi_1(index,1);
        end
    end
    
    phi_gate(:,:,ss) = phi_real;
    
%% The integrated electron density

elec = zeros(nx,ny);
% Inte_elec = zeros(kk(1),1);
% for s = 1:kk(1)

for ii=2:nx-1
    for jj = 2: ny -1
        
    elec(ii,jj)=ni*exp(phi_real(ii,jj)/thermal);
    end
end    
    
    elec_gate(:,:,ss) = elec;

end
    


%% plot

%%%%% Potential

x = linspace(0,60,61)*2;
y = linspace(0,6,61);

%%%% Position-dependent distribution at Gate_Volt = 0 V and 1 V



potential_gate_0V = phi_gate(:,:,1); % Gate_Volt = 0V

figure;
surf(y,x,potential_gate_0V);

legend('Gate = 0V')
title('Potential distribution')
% xlabel('Source_drain direction [nm]');
% ylabel('Gates direction  [nm]');

potential_gate_1V = phi_gate(:,:,11); % Gate_Volt = 0V

figure;
surf(y,x,potential_gate_1V);
legend('Gate = 1V')
title('Potential distribution')
% xlabel('Source_drain direction [nm]');
% ylabel('Gates direction  [nm]');

%%%% Gate dependence

%%% Source-drain direction

figure;
pot_gate_re = zeros(ny,11);
hold
for ii = 1: 11
    pot_gate_re(:,ii) =  phi_gate(:,30,ii);
    plot(x,pot_gate_re(:,ii))
end
% xlable('Source_drain direction [nm]');
% ylable('Potential  [V]');
% 








%%%%% Electron density

%%%% Position-dependent distribution at Gate_Volt = 0 V and 1 V

elec_gate_0V = elec_gate(:,:,1); % Gate_Volt = 0V

figure;
surf(y,x,elec_gate_0V);

legend('Gate = 0V')
title('Potential distribution')
% xlabel('Source_drain direction [nm]');
% ylabel('Gates direction  [nm]');

elec_gate_1V = elec_gate(:,:,11); % Gate_Volt = 0V

figure;
surf(y,x,elec_gate_1V);
legend('Gate = 1V')
title('Potential distribution')


%%%% Gate dependence

%%% Source-drain direction

figure;
elec_gate_re = zeros(nx,11);
hold
for ii = 1: 11
    elec_gate_re(:,ii) =  elec_gate(:,30,ii);
    plot(x,elec_gate_re(:,ii))
end

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    