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
 

        
%% Step 2
Ndon1 = zeros(Nx,Ny);
phi = zeros(Nx,Ny);
        for ii=1:Ny %Arrange of doping density and its potential
            if (ii >y_12 && ii<y_23)
                Ndon1(:,ii) = Ndon(:,1);
                phi(:,ii)= thermal*log(Ndon1(:,ii)/ni);
            elseif(ii == y_12 || ii == y_23)
                Ndon1(:,ii) = Ndon(:,1);
                phi(:,ii)= thermal*log(Ndon1(:,ii)/ni);
            elseif (ii<y_12 || ii> y_23)
                phi(:,ii) = 0.33374+0.1;
            end
        end
        
        phi=reshape(phi,Nx*Ny,1);
        Ndon1=reshape(Ndon1,Nx*Ny,1);
        elec = zeros(Nx*Ny,1);
        
        for ii=1:Nx*Ny % Electrond ensity inside the Si
            if (ii>Nx*(y_12-1) && ii<=Nx*y_23)
                elec(ii,1) = ni*exp(phi(ii,1)/thermal);
            else
                elec(ii,1) = 0;
            end
        end

A = zeros(1,4);

for bias=0:10
    for gate=0:10
        
        V_applied = 0.1*bias; % Drain voltage
        V_gate= 0.1*gate; % Gate voltage
        err=1;

        for Newton=1:10 % Nonlinear poisson solver using Newton's method
            res = zeros(2*Nx*Ny,1);
            Jaco = sparse(2*Nx*Ny,2*Nx*Ny);
            
            res(1,1) = 2*dx/dy*eps_ox*phi(1+Nx,1) + 2*dy/dx*eps_ox*phi(1+1,1) - 2*(dx/dy+dy/dx)*eps_ox*phi(1,1);
            Jaco(1,2*(1+Nx)-1) = 2*dx/dy*eps_ox;
            Jaco(1,2*(1+1)-1) = 2*dy/dx*eps_ox;
            Jaco(1,1) = - 2*(dx/dy+dy/dx)*eps_ox;
            
            res(2*Nx-1,1) = 2*dx/dy*eps_ox*phi(Nx+Nx,1) + 2*dy/dx*eps_ox*phi(Nx-1,1) - 2*(dx/dy+dy/dx)*eps_ox*phi(Nx,1);
            Jaco(2*Nx-1,2*(Nx+Nx)-1) = 2*dx/dy*eps_ox;
            Jaco(2*Nx-1,2*(Nx-1)-1) = 2*dy/dx*eps_ox;
            Jaco(2*Nx-1,2*Nx-1) = - 2*(dx/dy+dy/dx)*eps_ox;
            
            res(2*(Nx*(Ny-1)+1)-1,1) = 2*dx/dy*eps_ox*phi(Nx*(Ny-1)+1-Nx,1) + 2*dy/dx*eps_ox*phi(Nx*(Ny-1)+1+1,1) - 2*(dx/dy+dy/dx)*eps_ox*phi(Nx*(Ny-1)+1,1);
            Jaco(2*(Nx*(Ny-1)+1)-1,2*(Nx*(Ny-1)+1-Nx)-1) = 2*dx/dy*eps_ox;
            Jaco(2*(Nx*(Ny-1)+1)-1,2*(Nx*(Ny-1)+1+1)-1) = 2*dy/dx*eps_ox;
            Jaco(2*(Nx*(Ny-1)+1)-1,2*(Nx*(Ny-1)+1)-1) = -2*(dx/dy + dy/dx)*eps_ox;
            
            res(2*Nx*Ny-1,1)  =2*dx/dy*eps_ox*phi(Nx*Ny-Nx,1) + 2*dy/dx*eps_ox*phi(Nx*Ny-1,1) - 2*(dx/dy+dy/dx)*eps_ox*phi(Nx*Ny,1);
            Jaco(2*Nx*Ny-1,2*(Nx*Ny-Nx)-1) = 2*dx/dy*eps_ox;
            Jaco(2*Nx*Ny-1,2*(Nx*Ny-1)-1) = 2*dy/dx*eps_ox;
            Jaco(2*Nx*Ny-1,2*Nx*Ny-1) = - 2*(dx/dy+dy/dx)*eps_ox;
            for ii=1:Nx*Ny
                bound2 = ((ii)/Nx - fix(ii/Nx))*ii; % right boundary
                bound1 = ((ii-1)/Nx - fix((ii-1)/Nx))*ii; % left boundary
                if((ii>1 && ii<=x_12) || (ii>=x_23 && ii<Nx)) %% Bottom surface of SiO2
                    res(2*ii-1,1) = 2*(dx/dy)*eps_ox*phi(ii+Nx,1) + (dy/dx)*eps_ox*phi(ii+1,1) -  2*(dx/dy+dy/dx)*eps_ox*phi(ii,1) + (dy/dx)*eps_ox*phi(ii-1,1);
                    Jaco(2*ii-1,2*(ii+Nx)-1) = 2*dx/dy*eps_ox;
                    Jaco(2*ii-1,2*(ii+1)-1) = dy/dx*eps_ox;
                    Jaco(2*ii-1,2*ii-1) =-2*(dy/dx+dx/dy)*eps_ox;
                    Jaco(2*ii-1,2*(ii-1)-1) = dy/dx*eps_ox;
                    
                elseif(((ii>Nx && ii <=Nx*(y_12-1)) || (ii>Nx*y_23 && ii<=Nx*(Ny-1))) && bound1 ==0) %% left surface of SiO2
                    res(2*ii-1,1) = (dx/dy)*eps_ox*phi(ii+Nx,1) + (dy/dx)*2* eps_ox*phi(ii+1,1) - (dx/dy+dy/dx)*2*eps_ox*phi(ii,1) + (dx/dy)*eps_ox*phi(ii-Nx,1);
                    Jaco(2*ii-1,2*(ii+Nx)-1) = (dx/dy)*eps_ox;
                    Jaco(2*ii-1,2*(ii+1)-1) = (dy/dx)*2*eps_ox;
                    Jaco(2*ii-1,2*ii-1) = -(dx/dy+dy/dx)*2*eps_ox;
                    Jaco(2*ii-1,2*(ii-Nx)-1) = (dx/dy)*eps_ox;
                    
                elseif(((ii>Nx && ii <=Nx*(y_12-1)) || (ii>Nx*y_23 && ii<=Nx*(Ny-1))) && bound2 ==0) %% righte surface of SiO2
                    res(2*ii-1,1) = (dx/dy)*eps_ox*phi(ii+Nx,1) + (dy/dx)*2* eps_ox*phi(ii-1,1) - (dx/dy)*2*eps_ox*phi(ii,1)-(dy/dx)*2*eps_ox*phi(ii,1) + (dx/dy)*eps_ox*phi(ii-Nx,1);
                    Jaco(2*ii-1,2*(ii+Nx)-1) = (dx/dy)*eps_ox;
                    Jaco(2*ii-1,2*(ii-1)-1) = (dy/dx)*2*eps_ox;
                    Jaco(2*ii-1,2*ii-1) = -(dx/dy)*2*eps_ox - (dy/dx)*2*eps_ox;
                    Jaco(2*ii-1,2*(ii-Nx)-1) = (dx/dy)*eps_ox;
                    
                elseif(((ii>Nx && ii <=Nx*(y_12-1)) || (ii>Nx*y_23 && ii<=Nx*(Ny-1))) && bound2 ~=0 && bound1 ~= 0) %% inside SiO2
                    res(2*ii-1,1)      = (dx/dy)*eps_ox*phi(ii+Nx,1) + (dy/dx)*eps_ox*phi(ii+1,1) - 2*(dy/dx)*eps_ox*phi(ii,1)-(dx/dy)*2*eps_ox*phi(ii,1) + (dy/dx)*eps_ox*phi(ii-1,1) + (dx/dy)*eps_ox*phi(ii-Nx,1);
                    Jaco(2*ii-1,2*(ii+Nx)-1) = (dx/dy)*eps_ox;
                    Jaco(2*ii-1,2*(ii+1)-1 ) = (dy/dx)*eps_ox;
                    Jaco(2*ii-1,2*ii-1   ) = -(dy/dx)*2*eps_ox-(dx/dy)*2*eps_ox;
                    Jaco(2*ii-1,2*(ii-1)-1 ) = (dy/dx)*eps_ox;
                    Jaco(2*ii-1,2*(ii-Nx)-1) = (dx/dy)*eps_ox;
                    
                elseif((ii > Nx*(y_12-1) && ii <Nx*y_12) && (bound2 ~=0 && bound1 ~= 0)) %% Interface y_12
                    res(2*ii-1,1) = (dy/dx)*(eps_si+eps_ox)*phi(ii-1,1)/2 + (dx/dy)*eps_si*phi(ii+Nx,1) - (dy/dx)*(eps_si+eps_ox)*phi(ii,1) - (dx/dy)*(eps_si+eps_ox)*phi(ii,1) + (dx/dy)*eps_ox*phi(ii-Nx,1)+(dy/dx)*(eps_si+eps_ox)*phi(ii+1,1)/2+ coef*(Ndon1(ii,1)-elec(ii,1))/2;
                    Jaco(2*ii-1,2*(ii+Nx)-1) = (dx/dy)*eps_si;
                    Jaco(2*ii-1,2*(ii+1)-1) = (dy/dx)*(eps_si+eps_ox)/2;
                    Jaco(2*ii-1,2*ii-1) = -(dy/dx)*(eps_si + eps_ox)-(dx/dy)*(eps_si + eps_ox);%-coef*elec(ii,1)/thermal/2;
                    Jaco(2*ii-1,2*ii) = -coef/2;
                    Jaco(2*ii-1,2*(ii-1)-1) = (dy/dx)*(eps_si+eps_ox)/2;
                    Jaco(2*ii-1,2*(ii-Nx)-1) = (dx/dy)*eps_ox;
                    
                elseif(ii>Nx*y_12 && ii<=Nx*(y_23-1)&&(bound2 ~=0 && bound1 ~= 0)) %%Inside the silicon
                    res(2*ii-1,1) = (dx/dy)*eps_si*phi(ii+Nx,1)+ (dy/dx)*eps_si*phi(ii+1,1) - (dy/dx)*2*eps_si*phi(ii,1)- (dx/dy)*2*eps_si*phi(ii,1) + (dy/dx)*eps_si*phi(ii-1,1)+(dx/dy)*eps_si*phi(ii-Nx,1) + coef*(Ndon1(ii,1)-elec(ii,1));
                    Jaco(2*ii-1,2*(ii+Nx)-1) = (dx/dy)*eps_si;
                    Jaco(2*ii-1,2*(ii+1)-1) = (dy/dx)*eps_si;
                    Jaco(2*ii-1,2*ii-1) = -(dy/dx)*2*eps_si-(dx/dy)*2*eps_si;%-coef*elec(ii,1)/thermal;
                    Jaco(2*ii-1,2*ii) = -coef;
                    Jaco(2*ii-1,2*(ii-1)-1) = (dy/dx)*eps_si;
                    Jaco(2*ii-1,2*(ii-Nx)-1) = (dx/dy)*eps_si;
                    
                elseif((ii > Nx*(y_23-1) && ii <Nx*y_23) && (bound2 ~=0 && bound1 ~= 0)) %% Interface y_23
                    res(2*ii-1,1) = (dy/dx)*(eps_si+eps_ox)*phi(ii-1,1)/2 + (dx/dy)*eps_ox*phi(ii+Nx,1) - (dy/dx)*(eps_si+eps_ox)*phi(ii,1)-(dx/dy)* (eps_si+eps_ox)*phi(ii,1) + (dx/dy)*eps_si*phi(ii-Nx,1)+(dy/dx)*(eps_si+eps_ox)*phi(ii+1,1)/2+ coef*(Ndon1(ii,1)-elec(ii,1))/2;
                    Jaco(2*ii-1,2*(ii+Nx)-1) = (dx/dy)*eps_ox;
                    Jaco(2*ii-1,2*(ii+1)-1) = (dy/dx)*(eps_si+eps_ox)/2;
                    Jaco(2*ii-1,2*ii-1) = -(dy/dx)*(eps_si + eps_ox)-(dx/dy)*(eps_si + eps_ox);%coef*elec(ii,1)/thermal/2;
                    Jaco(2*ii-1,2*ii) = -coef/2;
                    Jaco(2*ii-1,2*(ii-1)-1) = (dy/dx)*(eps_si+eps_ox)/2;
                    Jaco(2*ii-1,2*(ii-Nx)-1) = (dx/dy)*eps_si;
                    
                elseif((ii>Nx*(Ny-1)+1&& ii<=Nx*(Ny-1)+x_12) || (ii>=Nx*(Ny-1)+x_23 && ii < Nx*Ny)) %% Top surface of SiO2
                    res(2*ii-1,1) = 2*(dx/dy)*eps_ox*phi(ii-Nx,1) + (dy/dx)*eps_ox*phi(ii+1,1) - 2*(dy/dx)*eps_ox*phi(ii,1) - 2*(dx/dy)*eps_ox*phi(ii,1) + (dy/dx)*eps_ox*phi(ii-1,1);
                    Jaco(2*ii-1,2*(ii-Nx)-1) = 2*dx/dy*eps_ox;
                    Jaco(2*ii-1,2*(ii+1)-1) = dy/dx*eps_ox;
                    Jaco(2*ii-1,2*ii-1) = -2*(dx/dy+dy/dx)* eps_ox;
                    Jaco(2*ii-1,2*(ii-1)-1) = dy/dx*eps_ox;
                    
                elseif(ii>=Nx*(Ny-1)+x_12 && ii<=Nx*(Ny-1)+x_23) % Top gate electrode
                    res(2*ii-1,1)= phi(ii,1) - 0.33374-V_gate;
                    Jaco(2*ii-1,2*ii-1) = 1.0;
                    
                elseif(ii>=x_12 && ii<=x_23) %% Bottom gate electrode
                    res(2*ii-1,1)= phi(ii,1) - 0.33374-V_gate;
                    Jaco(2*ii-1,2*ii-1) = 1.0;
                    
                elseif((2*ii-1 > Nx*(y_12-1) && ii <=Nx*y_23) && bound2 ==0 ) %% Source and drain electrode
                    res(2*ii-1,1) =phi(ii,1) - thermal*log(Ndon1(ii,1)/ni) - V_applied;
                    Jaco(2*ii-1,2*ii-1) =1.0;
                elseif((ii > Nx*(y_12-1) && ii <=Nx*y_23) && ( bound1 == 0)) %% drain silicon contact
                    res(2*ii-1,1) = phi(ii,1) - thermal*log(Ndon1(ii,1)/ni);
                    Jaco(2*ii-1,2*ii-1) =1.0;
                end
            end
            
            for jj=0:Ny-1
                for ii=1:Nx
                    if(jj<y_12-1 || jj>y_23-1)
                        res(2*(jj*Nx+ii),1) = 0;
                        Jaco(2*(jj*Nx+ii),:) = 0;
                        Jaco(2*(jj*Nx+ii),2*(jj*Nx+ii)) = 1.0;
                    end
                end
                
                for ii=1:Nx
                    if((jj>=y_12-1 && jj<=y_23-1) && (ii>=1 && ii<Nx))
                        n_avx = 0.5*(elec(jj*Nx+ii+1,1)+elec(jj*Nx+ii,1));
                        n_avy = 0.5*(elec(jj*Nx+ii+Nx,1)+elec(jj*Nx+ii,1));
                        dphidx = (phi(jj*Nx+ii+1,1)-phi(jj*Nx+ii,1))/dx;
                        dphidy = (phi(jj*Nx+ii+Nx,1)-phi(jj*Nx+ii,1))/dy;
                        delecdx = (elec(jj*Nx+ii+1,1)-elec(jj*Nx+ii,1))/dx;
                        delecdy = (elec(jj*Nx+ii+Nx,1)-elec(jj*Nx+ii,1))/dy;
                        Jnx = n_avx * dphidx - thermal * delecdx;
                        Jny = n_avy * dphidy - thermal * delecdy;
                        
                        res(2*(jj*Nx+ii),1) = res(2*(jj*Nx+ii),1) + Jnx;
                        Jaco(2*(jj*Nx+ii),2*(jj*Nx+ii+1)) = Jaco(2*(jj*Nx+ii),2*(jj*Nx+ii+1)) + 0.5*dphidx - thermal / dx;
                        Jaco(2*(jj*Nx+ii),2*(jj*Nx+ii )) = Jaco(2*(jj*Nx+ii),2*(jj*Nx+ii) ) + 0.5*dphidx + thermal / dx;
                        Jaco(2*(jj*Nx+ii),2*(jj*Nx+ii)+1) = Jaco(2*(jj*Nx+ii),2*(jj*Nx+ii)+1) + n_avx/dx;
                        Jaco(2*(jj*Nx+ii),2*(jj*Nx+ii)-1) = Jaco(2*(jj*Nx+ii),2*(jj*Nx+ii)-1) - n_avx/dx;
                        
                        res(2*(jj*Nx+ii+1),1) = res(2*(jj*Nx+ii+1),1) - Jnx;
                        Jaco(2*(jj*Nx+ii+1),2*(jj*Nx+ii+1)) = Jaco(2*(jj*Nx+ii+1),2*(jj*Nx+ii+1)) - 0.5*dphidx + thermal / dx;
                        Jaco(2*(jj*Nx+ii+1),2*(jj*Nx+ii) ) = Jaco(2*(jj*Nx+ii+1),2*(jj*Nx+ii )) - 0.5*dphidx - thermal / dx;
                        Jaco(2*(jj*Nx+ii+1),2*(jj*Nx+ii)+1) = Jaco(2*(jj*Nx+ii+1),2*(jj*Nx+ii)+1) - n_avx/dx;
                        Jaco(2*(jj*Nx+ii+1),2*(jj*Nx+ii)-1) = Jaco(2*(jj*Nx+ii+1),2*(jj*Nx+ii)-1) + n_avx/dx;
                        
                        %res(2*(jj*Nx+ii),1) = res(2*(jj*Nx+ii),1) + Jny;
                        %Jaco(2*(jj*Nx+ii),2*(jj*Nx+ii+Nx)) = Jaco(2*(jj*Nx+ii),2*(jj*Nx+ii+Nx)) + 0.5*dphidy - thermal / dy;
                        %Jaco(2*(jj*Nx+ii),2*(jj*Nx+ii )) = Jaco(2*(jj*Nx+ii),2*(jj*Nx+ii) ) + 0.5*dphidy + thermal / dy;
                        %Jaco(2*(jj*Nx+ii),2*(jj*Nx+ii+Nx)-1) = Jaco(2*(jj*Nx+ii),2*(jj*Nx+ii+Nx)-1) + n_avy/dy;
                        %Jaco(2*(jj*Nx+ii),2*(jj*Nx+ii)-1) = Jaco(2*(jj*Nx+ii),2*(jj*Nx+ii)-1) - n_avy/dy;
                        
                        %res(2*(jj*Nx+ii+Nx),1) = res(2*(jj*Nx+ii+Nx),1) + Jny;
                        %Jaco(2*(jj*Nx+ii+Nx),2*(jj*Nx+ii+Nx)) = Jaco(2*(jj*Nx+ii+Nx),2*(jj*Nx+ii+Nx)) + 0.5*dphidy - thermal / dy;
                        %Jaco(2*(jj*Nx+ii+Nx),2*(jj*Nx+ii )) = Jaco(2*(jj*Nx+ii+Nx),2*(jj*Nx+ii) ) + 0.5*dphidy + thermal / dy;
                        %Jaco(2*(jj*Nx+ii+Nx),2*(jj*Nx+ii+Nx)-1) = Jaco(2*(jj*Nx+ii+Nx),2*(jj*Nx+ii+Nx)-1) + n_avy/dy;
                        %Jaco(2*(jj*Nx+ii+Nx),2*(jj*Nx+ii)-1) = Jaco(2*(jj*Nx+ii+Nx),2*(jj*Nx+ii)-1) - n_avy/dy;
                    end
                end
                for ii=1:Nx
                    if((jj>=y_12-1 && jj<=y_23-1) && (ii==1))
                        res(2*(jj*Nx+ii),1) = elec(jj*Nx+ii,1) - Ndon1(jj*Nx+ii,1);
                        Jaco(2*(jj*Nx+ii),2*(jj*Nx+ii)) = 1.0;
                    elseif((jj>=y_12-1 && jj<=y_23-1) && (ii==Nx))
                        res(2*(jj*Nx+ii),1) = elec(jj*Nx+ii,1) - Ndon1(jj*Nx+ii,1);
                        Jaco(2*(jj*Nx+ii),2*(jj*Nx+ii)) = 1.0;
                    end
                end
            end
            
            Cvector = zeros(2*Nx*Ny,1);
            Cvector(1:2:2*Nx*Ny-1,1) = thermal;
            Cvector(2:2:2*Nx*Ny ,1) = max(abs(Ndon));
            Cmatrix = spdiags(Cvector,0,2*Nx*Ny,2*Nx*Ny);
            Jaco_scaled = Jaco * Cmatrix;
            Rvector = 1./sum(abs(Jaco_scaled),2);
            Rmatrix = spdiags(Rvector,0,2*Nx*Ny,2*Nx*Ny);
            Jaco_scaled = Rmatrix * Jaco_scaled;
            res_scaled = Rmatrix * res;
            update_scaled = Jaco_scaled \ (-res_scaled);
            update = Cmatrix * update_scaled;
            
            phi = phi + update(1:2:2*Nx*Ny-1,1);
            elec = elec + update(2:2:2*Nx*Ny,1);
            err = norm(update(1:2:2*Nx*Ny-1,1),inf);
            A(1,1) = bias;
            A(1,2) = gate;
            A(1,3) = Newton;
            A(1,4) = err;
            A
        end
        J = zeros(Nx-1,Ny);
        for jj=y_12-1:y_23
            for ii=1:Nx-1
                t(ii,jj+1) = (phi(Nx*jj+ii+1,1)-phi(Nx*jj+ii,1))/thermal;
                J(ii,jj+1) = (q*Dn/dx)*(elec(Nx*jj+ii+1,1)*t(ii,jj+1)/(exp(t(ii,jj+1))-1) + elec(Nx*jj+ii,1)*t(ii,jj+1)/(exp(-t(ii,jj+1))-1));
            end
            I(bias+1,gate+1) = I(bias+1,gate+1)+J(Nx-1,jj)*dy;
        end
    end
end
%% Step 3
h = 6.626176e-34; % Planck constant, J s
hbar = h / (2*pi); % Reduced Planck constant, J s
m0 = 9.109534e-31; % Electron rest mass, kg
mxx = 0.91; myy = 0.19; mzz = 0.19; % Masses, m0
%mxx = 0.19; myy = 0.91; mzz = 0.19; % Masses, m0
%mxx = 0.19; myy = 0.19; mzz = 0.91; % Masses, m0
coef_Sch = 2*Nx*dx/(2*pi)*sqrt(mxx*myy)*m0/(hbar^2)*(k_B*T);
Ec_Ei = 0.561004; % E_c ? E_i, eV
V = q*Ec_Ei - q*phi; % Potential energy, J
Nbulk = y_23-y_12-1; %Number of bulk silicon nodes
elec_Sch = zeros(Ny,Nx); % Electron density, /m^3
for jj=1:Nx
    for ii=1:Ny
        phi_y(ii,1) = phi((ii-1)*Nx+jj,1);
    end
    V = q*Ec_Ei - q*phi_y;
    Hamil = zeros(Nbulk,Nbulk);
    Hamil(1,1) = -2; Hamil(1,2) = 1;
    for ii=2:Nbulk-1
        Hamil(ii,ii+1) = 1;
        Hamil(ii,ii ) = -2;
        Hamil(ii,ii-1) = 1;
    end
    Hamil(Nbulk,Nbulk) = -2; Hamil(Nbulk,Nbulk-1) = 1;
    for ii=1:Nbulk
        Hamil(ii,ii) = Hamil(ii,ii) -2*mzz*m0*(dy/hbar)^2*V(ii+y_12,1);
    end
    
    [eigenvectors,eigenvalues] = eig(Hamil);
    scaledEz = diag(eigenvalues)/(-2*mzz*m0*(dy/hbar)^2); % Eigenenergy, J
    [sortedEz,sortedIndex] = sort(scaledEz);
    
    nSubband = 1; % The lowest subband number
    
    totalNumber = 0;
    for n=1:nSubband
        Ez = sortedEz(n,1);
        wavefunction2 = eigenvectors(:,sortedIndex(n)).^2;
        normalization = sum(wavefunction2)*dy;
        wavefunction2 = wavefunction2 / normalization;
        subbandNumber = coef_Sch*log(1+exp(-Ez/(k_B*T)));
        totalNumber = totalNumber + subbandNumber;
        elec_Sch(y_12+1:y_23-1,jj) = elec_Sch(y_12+1:y_23-1,jj) + 1/(Nx*dx)*wavefunction2*subbandNumber;
    end
end
figure(1)
surf(elec_Sch);
xlabel('x (nm)');
ylabel('y (nm)');
zlabel('n (cm^{-3})');
set(gca,'zscale','log');


%% Step 4

vf=1e+6;
t = 1e-12;
M=101;
Z = m0/hbar/hbar/2/2/pi/pi;
H = (linspace(0.1,1,M))'; % (eV)
f1 = zeros(Nx*Ny-1,M);
A= zeros(Nx*Ny,Nx*Ny);
b = zeros(Nx*Ny,1);
for kk=1:M
    for jj=0:Ny-1
        b(jj*Nx+1,1) = sqrt(2*pi)/(1 + exp(q*H(kk,1)/(k_B*T)));  %fs
        b((jj+1)*Nx,1) = sqrt(2*pi)/(1 + exp(q*(H(kk,1)+phi((jj+1)*Nx,1))/(k_B*T)));  %fd
        A(jj*Nx+1,jj*Nx+1) = 1.0;
        for ii=2:Nx-1
            c1 = H(kk,1) + 0.5*(phi(jj*Nx+ii,1)+phi(jj*Nx+ii-1,1));
            c2 = H(kk,1) + 0.5*(phi(jj*Nx+ii+1,1)+phi(jj*Nx+ii,1));
            A(jj*Nx+ii,jj*Nx+ii-1) = c1; A(jj*Nx+ii,jj*Nx+ii) = -c1-c2; A(jj*Nx+ii,jj*Nx+ii+1) = c2;
        end
        A((jj+1)*Nx,(jj+1)*Nx) = 1.0;
    end
    f0(:,kk) = A \ b;
    
    for  jj=0:Ny-1
        for ii=1:Nx-1
            f1(jj*Nx+ii,kk) = -t*vf*(1/sqrt(pi))*(f0(jj*Nx+ii+1,1) - f0(jj*Nx+ii,1))/dx;
        end
    end
end
f1_H = zeros(Nx*Ny-1,1);
f0_H = zeros(Nx*Ny-1,1);
for ii=1:Nx*Ny-1
    for kk=1:M
        f0_H(ii,1) = f0_H(ii,1)*(H(2,1)-H(1,1))+f0(ii,kk);
        f1_H(ii,1) = f1_H(ii,1)*(H(2,1)-H(1,1))+f1(ii,kk);
    end
end
f1_H(Nx*Ny,1)=0;
f0_H(Nx*Ny,1)=0;

surf(reshape(f0_H,Nx,Ny));
xlabel('x (nm)');
ylabel('y (nm)');
zlabel('f0');