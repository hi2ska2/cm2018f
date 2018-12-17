clear all;

q = 1.602192e-19; % Elementary charge, C
e0 = 8.854187817e-12; % Vacuum permittivity, F/m
k_B = 1.380662e-23; % Boltzmann constant, J/K
T = 300.0; % Temperature, K

thermal = k_B*T/q; % Thermal voltage, V


%%
ni = 1.075e16; % 1.075e10 /cm^3 intrinsic carrier density
Nacc = 1e24; % 1e18 /cm^3 doping density
dx = 0.1e-9;
dy = 0.1e-9;

coef = dx*dx*q/e0;

nx = 1201; % 120 nm
x_12 = 401;
x_23 = 801;
x = dx*transpose([0:nx-1]);

ny = 61; %6 nm
y_12 = 6;
y_23 = 56;

e_si = 11.7; e_ox = 3.9;
Ndon = zeros(nx*ny,1);
for iy = y_12:y_23
    Ndon((iy-1)*nx+(1:nx),1) = 5e23; %5e23 -> 5e17 /cm^3
    Ndon((iy-1)*nx+(1:x_12),1) = 5e25; % 5e25 -> 5e19 /cm^3
    Ndon((iy-1)*nx+(x_23:nx),1) = 5e25; % 5e25 -> 5e19 /cm^3
end
%% for checking
% for ix = 1:nx
%     for iy = 1: ny
% N_c(iy,ix) = phi((iy-1)*nx+ix,1);
%     end
% end
%%
for ii = 1:10
Vg = (ii-1)/10;
phi = zeros(nx*ny,1);
phi(:,1) = 0.33374+Vg;

for ix = 2:nx-1
for iy = y_12:y_23
phi((iy-1)*nx+ix,1) = phi((iy-1)*nx+ix,1) + thermal*log(Ndon((iy-1)*nx+ix,1)/ni);
end
end


for newton = 1:10
    
res = zeros(nx*ny,1);
Jaco = sparse(nx*ny,nx*ny);


for ix = 2:nx-1
for iy = 2:ny-1
    n = (iy-1)*nx+ix;
    
    % Jacobian & residue

     if (iy< y_12 || iy> y_23)
         res(n,1) = e_ox*phi(n+1,1) + e_ox*phi(n+nx,1) - 4*e_ox*phi(n,1)+ e_ox*phi(n-nx,1) + e_ox*phi(n-1,1);
         Jaco(n,n-1) = e_ox;
         Jaco(n,n-nx) = e_ox;
         Jaco(n,n) = -4*e_ox;
         Jaco(n,n+nx) = e_ox;
         Jaco(n,n+1) = e_ox;
     elseif (iy==y_12)
         res(n,1) = (e_si+e_ox)/2*phi(n+1,1) + e_si*phi(n+nx,1) - 2*(e_si+e_ox)*phi(n,1) +e_ox*phi(n-nx,1) + (e_si+e_ox)/2*phi(n-1,1);
         Jaco(n,n-1) = (e_si+e_ox)/2;
         Jaco(n,n-nx) = e_ox;
         Jaco(n,n) = -2*(e_si+e_ox);
         Jaco(n,n+nx) = e_si;
         Jaco(n,n+1) = (e_si+e_ox)/2;
     elseif (iy==y_23)
         res(n,1) = (e_ox+e_si)/2*phi(n+1,1) + e_ox*phi(n+1,1) - 2*(e_ox+e_si)*phi(n,1)+ e_si*phi(n-1,1) + (e_ox+e_si)/2*phi(n-1,1);
         Jaco(n,n-1) = (e_si+e_ox)/2;
         Jaco(n,n-nx) = e_si;
         Jaco(n,n) = -2*(e_si+e_ox);
         Jaco(n,n+nx) = e_ox;
         Jaco(n,n+1) = (e_si+e_ox)/2;
     else
         res(n,1) = e_si*phi(n+1,1) + e_si*phi(n+nx,1) - 2*e_si*phi(n,1)+ e_si*phi(n-nx,1) + e_si*phi(n-1,1);
         Jaco(n,n-1) = e_si;
         Jaco(n,n-nx) = e_si;
         Jaco(n,n) = -4*e_si;
         Jaco(n,n+nx) = e_si;
         Jaco(n,n+1) = e_si;
     end
     
     
end
end


for ix = 2:nx-1
for iy = 2:ny-1
    n = (iy-1)*nx+ix;
    for iy=y_12:y_23
     if (iy==y_12)
         res(n,1) = res(n,1) - coef*(-Ndon(n,1)+ni*exp(phi(n,1)/thermal))*0.5;
         Jaco(n,n) = Jaco(n,n) - coef*ni*exp(phi(n,1)/thermal)/thermal*0.5;
     elseif (iy==y_23)
         res(n,1) = res(n,1) - coef*(-Ndon(n,1)+ni*exp(phi(n,1)/thermal))*0.5;
         Jaco(n,n) = Jaco(n,n) - coef*ni*exp(phi(n,1)/thermal)/thermal*0.5;
     else
         res(n,1) = res(n,1) - coef*(-Ndon(n,1)+ni*exp(phi(n,1)/thermal));
         Jaco(n,n) = Jaco(n,n) - coef*ni*exp(phi(n,1)/thermal)/thermal;
     end
     end
end
end


%      %        %         %   boundary  %    %    %   %
%%%%%%%%%%%%     %bottom        %%%%%%%%%%%%
for ix = 2 : x_23 - 1
    n = ix;
    res(n,1) = phi(n,1) - (0.33374+Vg);
    Jaco(n,n + 1) = 0.5;
    Jaco(n,n + nx) = 1;
    Jaco(n,n) = -2;
    Jaco(n,n - 1) = 0.5;
end

for ix = x_12 + 1 : x_23 - 1
    n = ix;
    res(n,1) = phi(n,1) - (0.33374+Vg);
    Jaco(n,n) = 1;
end

for ix = x_23 : nx - 1
    n = ix;
    res(n,1) = phi(n,1) - (0.33374+Vg);
    Jaco(n,n + 1) = 0.5;
    Jaco(n,n + nx) = 1;
    Jaco(n,n) = -2;
    Jaco(n,n - 1) = 0.5;
end


    res(1,1) = phi(1,1) - (0.33374+Vg);
    Jaco(1,1 + 1) = 1;
    Jaco(1,1 + nx) = 1;
    Jaco(1,1) = -2;
    
    res(nx,1) = phi(nx,1) - (0.33374+Vg);
    Jaco(nx,nx + nx) = 1;
    Jaco(nx,nx) = -2;
    Jaco(nx,nx - 1) = 1;

%%%%%%%%        top        %%%%%%%%%
for ix = 2:x_12
    n = (ny-1)*nx+ix;
    res(n,1) = phi(n,1) - (0.33374+Vg);
    Jaco(n,n + 1) = 0.5;
    Jaco(n,n - nx) = 1;
    Jaco(n,n) = -2;
    Jaco(n,n - 1) = 0.5;
end

for ix = x_12 + 1 : x_23 - 1
    n = (ny-1)*nx+ix;
    res(n,1) = phi(n,1) - (0.33374+Vg);
    Jaco(n,n) = 1;
end

for ix = x_23 : nx-1
    n = (ny-1)*nx+ix;
    res(n,1) = phi(n,1) - (0.33374+Vg);
    Jaco(n,n + 1) = 0.5;
    Jaco(n,n - nx) = 1;
    Jaco(n,n) = -2;
    Jaco(n,n - 1) = 0.5;
end

    res((ny-1)*nx + 1, 1) = phi((ny-1)*nx+1,1) - (0.33374+Vg);
    Jaco((ny-1)*nx + 1, (ny-1)*nx + 1 + 1) = 1;
    Jaco((ny-1)*nx + 1, (ny-1)*nx + 1) = -2;
    Jaco((ny-1)*nx + 1, (ny-1-1)*nx + 1) = 1;
    
    res((ny-1)*nx+nx,1) = phi((ny-1)*nx+nx,1) - (0.33374+Vg);
    Jaco((ny-1)*nx + nx, (ny-1)*nx + nx - 1 ) = 1;
    Jaco((ny-1)*nx + nx, (ny-1-1)*nx + nx) = 1;
    Jaco((ny-1)*nx + nx, (ny-1)*nx + nx) = -2;

%%%%%%%%%%      left       %%%%%%%%%%%%

for iy = 2:y_12
    n = (iy-1)*nx+1;
    res(n,1) = phi(n,1) - (0.33374+Vg);
    Jaco(n,n + 1) = 1;
    Jaco(n,n + nx) = 0.5;
    Jaco(n,n) = -2;
    Jaco(n,n - nx) = 0.5;
end


for iy = y_12:y_23
    n = (iy-1)*nx+1;
    res(n,1) = phi(n,1) - (0.33374) - thermal*log(Ndon(n,1)/ni);
    Jaco(n,n) = 1;
end

for iy = y_23:ny-1
    n = (iy-1)*nx+1;
    res(n,1) = phi(n,1) - (0.33374+Vg);
    Jaco(n,n + 1) = 1;
    Jaco(n,n + nx) = 0.5;
    Jaco(n,n) = -2;
    Jaco(n,n - nx) = 0.5;
end

%%%%%%%%%%%%     right     %%%%%%%%%%%%%%%


for iy = 2:y_12
    n = (iy-1)*nx+nx;
    res(n,1) = phi(n,1) - (0.33374+Vg);
    
    Jaco(n,n + nx) = 0.5;
    Jaco(n,n) = -2;
    Jaco(n,n - nx) = 0.5;
    Jaco(n,n - 1) =1;
end


for iy = y_12:y_23
    n = (iy-1)*nx+nx;
    res(n,1) = phi(n,1) - (0.33374) - thermal*log(Ndon(n,1)/ni);
    Jaco(n,n) = 1;
end

for iy = y_23:ny-1
    n = (iy-1)*nx+nx;
    res(n,1) = phi(n,1) - (0.33374+Vg);
    Jaco(n,n + nx) = 0.5;
    Jaco(n,n) = -2;
    Jaco(n,n - nx) = 0.5;
    Jaco(n,n - 1) =1;
end




update = Jaco \ (-res);
phi = phi + update;

end

elec = ni*exp(phi/thermal);

for ix = 1:nx
     for iy = 1: ny
 phi_3D(iy,ix,ii) = phi((iy-1)*nx+ix,1);
     end
end

for ix = 1:nx
     for iy = 1: ny
 elec_3D(iy,ix,ii) = elec((iy-1)*nx+ix,1);
     end
end

phi_stack(:,1,ii) = phi(:,1);
elec_stack(:,1,ii) = phi(:,1);

end
figure(201)
surf(phi_3D(y_12:y_23,2:nx-1,1));
shading interp;
xlabel('x-axis(delta)');
ylabel('y-axis(delta)');
zlabel('phi (V)');
colorbar;


figure(202)
surf(elec_3D(y_12:y_23,2:nx-1,1));
shading interp;
xlabel('x-axis(delta)');
ylabel('y-axis(delta)');
zlabel('electron density (m^-^3)');
colorbar;