clear all;

q = 1.602192e-19; % Elementary charge, C
eps0 = 8.854187817e-12; % Vacuum permittivity, F/m
k_B = 1.380662e-23; % Boltzmann constant, J/K
T = 300.0; % Temperature, K
thermal = k_B*T/q; % Thermal voltage, V
dy = 1e-10; % 1 nm spacing
dx = 1e-9;
dx2 = dx*dx;
dy2 = dy*dy;
Nx =121; % 120-nm-long structure
Ny =61; % 6 nm long width
x = dx*transpose([0:Nx-1]); % real space, m
x_12 = 41; % At x=40 nm
x_23 = 81; % At x=80 nm
y = dy*transpose([0:Ny-1]); % real space, m
y_12 = 6;
y_23 = 56;
eps_si = 11.7; eps_ox = 3.9; eps_avg = (eps_ox+eps_si)/2;% Relative permittivity
Ndon = 2e23*ones(Nx,1); % 2e17 /cm^3
Ndon(1:x_12,1) = 5e25; % 5e19 /cm^3
Ndon(x_23:Nx,1) = 5e25; % 5e19 /cm^3
ni = 1.075e16; % 1.075e10 /cm^3
u = 1430e-4; % Electron mobility, m2/Vs
Dn = thermal*u;
I = zeros(11,11);
 coef = dy*dx*q/eps0;
 
%% Step 1 
 for ii=1:Ny
    if (ii >y_12 && ii<y_23)
        Ndon1(:,ii) = Ndon(:,1);
        phi(:,ii)= thermal*log(Ndon1(:,ii)/ni);
    elseif(ii == y_12 || ii == y_23)
        Ndon1(:,ii) = Ndon(:,1)/2;
        phi(:,ii)= thermal*log(Ndon1(:,ii)/ni);
    elseif (ii<y_12 || ii> y_23)
        phi(:,ii) = 0.33374;
    end
end
elec = zeros(Nx*Ny,1);

phi=reshape(phi,Nx*Ny,1);
Ndon1=reshape(Ndon1,Nx*Ny,1);

for ii=1:Nx*Ny
    if (ii>Nx*(y_12-1) && ii<=Nx*y_23)
        elec(ii,1) = ni*exp(phi(ii,1)/thermal);
    else
        elec(ii,1) = 0;
    end
end
NewPhi = zeros(Nx*Ny,1);

count = zeros(11,1);
for gate=5:10
    err=1.0;
    while err > 5e-4
        count(gate+1,1) = count(gate+1,1) +1;
        count(gate+1,1)
        gate
        for ii=1:Nx*Ny;
            
            bound2 = ((ii)/Nx - fix(ii/Nx))*ii; % right boundary
            bound1 = ((ii-1)/Nx - fix((ii-1)/Nx))*ii; % left boundar
            NewPhi(1,1          ) = (dx2*phi(1+Nx,1) + dy2*phi(1+1,1))/(dx2 + dy2);
            NewPhi(Nx,1         ) = (dx2*phi(Nx+Nx,1) + dy2*phi(Nx-1,1))/(dx2 + dy2);
            NewPhi(Nx*(Ny-1)+1,1) = (dx2*phi(Nx*(Ny-1)+1-Nx,1) + dy2*phi(Nx*(Ny-1)+1+1,1))/(dx2 + dy2);
            NewPhi(Nx*Ny,1      ) = (dx2*phi(Nx*Ny-Nx,1) + dy2*phi(Nx*Ny-1,1))/(dx2+dy2);
            
            
            if((ii>1 && ii<=x_12) || (ii>=x_23 && ii<Nx)) %% Bottom surface of SiO2
                NewPhi(ii,1) = (2*dx2*phi(ii+Nx,1) + dy2*phi(ii-1,1) + dy2*phi(ii+1,1))/2/(dx2+dy2);
                
            elseif(ii>x_12 && ii<x_23) %% Bottom gate electrode
                NewPhi(ii,1) = 0.33374+gate*0.1;
                
            elseif(((ii>Nx && ii <=Nx*(y_12-1)) || (ii>Nx*y_23 && ii<=Nx*(Ny-1))) && bound1 ==0) %% left surface of SiO2
                NewPhi(ii,1) = (2*dy2*phi(ii+1,1) + dx2*phi(ii+Nx,1) + dx2*phi(ii-Nx,1))/2/(dx2+dy2);
                
            elseif(((ii>Nx && ii <=Nx*(y_12-1)) || (ii>Nx*y_23 && ii<=Nx*(Ny-1))) && bound2 ==0) %% righte surface of SiO2
                NewPhi(ii,1) = (2*dy2*phi(ii-1,1) + dx2*phi(ii+Nx,1) + dx2*phi(ii-Nx,1))/2/(dx2+dy2);
                
            elseif(((ii>Nx && ii <=Nx*(y_12-1)) || (ii>Nx*y_23 && ii<=Nx*(Ny-1))) && bound2 ~=0 && bound1 ~= 0) %% inside SiO2
                NewPhi(ii,1) = (dy2*phi(ii-1,1) + dy2*phi(ii+1,1) + dx2*phi(ii+Nx,1) + dx2*phi(ii-Nx,1))/2/(dx2+dy2);
                
            elseif((ii > Nx*(y_12-1) && ii <Nx*y_12) && (bound2 ~=0 && bound1 ~= 0)) %% Interface y_12
                NewPhi(ii,1) = (coef*(Ndon1(ii,1) - elec(ii,1))/2 + dy2*eps_avg*phi(ii-1,1) + dy2*eps_avg*phi(ii+1,1) + dx2*eps_si*phi(ii+Nx,1) + dx2*eps_ox*phi(ii-Nx,1))/2/(dx2+dy2)/eps_avg;
                
            elseif((ii > Nx*(y_12-1) && ii <=Nx*y_23) && (bound1 == 0)) %% Source and drain electrode
                NewPhi(ii,1) = thermal*log(Ndon1(ii,1)/ni);
                
            elseif((ii > Nx*(y_12-1) && ii <=Nx*y_23) && bound2 ==0 ) %% Source and drain electrode
                NewPhi(ii,1) =thermal*log(Ndon1(ii,1)/ni);
                
            elseif(ii>Nx*y_12 && ii<=Nx*(y_23-1)&&(bound2 ~=0 && bound1 ~= 0)) %%Inside the silicon
                NewPhi(ii,1) = (coef*(Ndon1(ii,1) - elec(ii,1)) + eps_si*dx2*phi(ii+Nx,1) + eps_si*dx2*phi(ii-Nx,1) + eps_si*dy2*phi(ii+1,1) + eps_si*dy2*phi(ii-1,1))/eps_si/2/(dx2+dy2);
                
            elseif((ii > Nx*(y_23-1) && ii <Nx*y_23) && (bound2 ~=0 && bound1 ~= 0)) %% Interface y_23
                NewPhi(ii,1) = (coef*(Ndon1(ii,1) - elec(ii,1))/2 + dy2*eps_avg*phi(ii-1,1) + dy2*eps_avg*phi(ii+1,1) + dx2*eps_ox*phi(ii+Nx,1) + dx2*eps_si*phi(ii-Nx,1))/2/(dx2+dy2)/eps_avg;
                
            elseif((ii>Nx*(Ny-1)+1&& ii<=Nx*(Ny-1)+x_12) || (ii>=Nx*(Ny-1)+x_23 && ii < Nx*Ny)) %% Top surface of SiO2
                NewPhi(ii,1) = (2*dx2*phi(ii-Nx,1) + dy2*phi(ii-1,1) + dy2*phi(ii+1,1))/2/(dx2+dy2);
                
            elseif(ii>Nx*(Ny-1)+x_12 && ii<Nx*(Ny-1)+x_23) % Top gate electrode
                NewPhi(ii,1) = 0.33374 + 0.1*gate;
            end
        end
        
        for kk=1:Nx*Ny
            if (kk>Nx*(y_12-1) && kk<=Nx*(y_23))
                elec(kk,1) = ni*exp(phi(kk,1)/thermal);
            elseif (kk<=Nx*(y_12-1) || kk>Nx*y_23)
                elec(kk,1) =0;
            end
        end
        err=sqrt(sum(sum((NewPhi-phi).^2)));
        phi = NewPhi;
    end
    figure(2*(gate+1))
    surf(reshape(phi,Nx,Ny)');
    xlabel('x (nm)'); ylabel('y (nm)'); zlabel('Potential (V)');
    figure(2*(gate+1)+1)
    surf(reshape(1e-6*elec,Nx,Ny)');
    set(gca,'ZScale','log');
    xlabel('x (nm)'); ylabel('y (nm)'); zlabel('Electron density(cm^{-3})');
end