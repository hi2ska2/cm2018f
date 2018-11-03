clear all;

%%
% Permittivity
% 
% Aluminum: 1.46
% mos2: 4
% silicon dioxide: 3.9
% silicon: 11.68
%%
q = 1.602192e-19; % Elementary charge, C
e0 = 8.854187817e-12; % Vacuum permittivity, F/m
k_B = 1.380662e-23; % Boltzmann constant, J/K
T = 300.0; % Temperature, K
h = 6.626176e-34; % Planck constant, J s
hbar = h / (2*pi); % Reduced Planck constant, J s
m0 = 9.109534e-31; % Electron rest mass, kg

thermal = k_B*T/q; % Thermal voltage, V


%%
ni = 1.075e16; % 1.075e10 /cm^3 intrinsic carrier density
Nacc = 1e24; % 1e18 /cm^3 doping density
dx = 0.1e-9;
dy = 0.1e-9;
dz = 0.1e-9;
Lx = 100e-9; Ly = 100e-9; % Lenghs, m
mxx = 0.19; myy = 0.19; mzz = 0.91; % Masses, m0
coef = dz*dz*q/e0;
coef_Sch = 2*Lx*Ly/(2*pi)*sqrt(mxx*myy)*m0/(hbar^2)*(k_B*T);
Ec_Ei = 0.561004; % E_c - E_i, eV
N_point = 5; % N points in minimum thickness



mat_th = [5 50 5]; % material thickness (x0.1 nm)
e = [3.9 11.7 3.9]; % permittivity

d = mat_th(1);
for n = 2:numel(mat_th)
  d = gcd(d,mat_th(n)); % greatest common divisor
end

Num = mat_th.*N_point./d; % number of points in each material
N = sum(Num)+1;

for n = 1:numel(mat_th)-1
    interface(n) = sum(mat_th(1:n))+1;
    
end % Interface define
z = dz*transpose([0:sum(Num)]); % real space thickness

for Vg1 = 1 : 11
    Vg=(Vg1-1)/10;
phi = zeros(N,1);
phi(:,1) = 0.33374+Vg;

for n = 1:1000
    
Jaco=sparse(sum(Num)+1,sum(Num)+1);
Jaco(1,1)=1;
Jaco(sum(Num)+1,sum(Num)+1)=1;

res = zeros(N,1);
res(1,1) = phi(1,1) - (0.33374+Vg);
res(N,1) = phi(N,1) - (0.33374+Vg);
sum_num=1;
for i = 1:numel(mat_th)
    for j = sum_num+1:sum_num+Num(i)-1
        Jaco(j,j-1) = e(i);
        Jaco(j,j) = -2*e(i);
        Jaco(j,j+1) = e(i);
        
        res(j,1) = e(i)*phi(j+1,1) - 2*e(i)*phi(j,1) + e(i)*phi(j-1,1);
        
    end
    
    if i < numel(mat_th)
    sum_num=sum_num+Num(i);
    
    Jaco(sum_num,sum_num-1) = e(i);
    Jaco(sum_num,sum_num) = -e(i) -e(i+1);
    Jaco(sum_num,sum_num+1) = e(i+1);
    
    res(sum_num,1) = e(i+1)*phi(sum_num+1,1) -(e(i)+e(i+1))*phi(sum_num,1) + e(i)*phi(sum_num-1,1);
    
    end
end

% Charge
for i = Num(1)+1:Num(1)+Num(2)+1
    if i==Num(1)+1
        res(i,1) = res(i,1) - coef*(Nacc+ni*exp(phi(i,1)/thermal))*0.5;
        Jaco(i,i) = Jaco(i,i) - coef*ni*exp(phi(i,1)/thermal)/thermal*0.5;
    elseif i==Num(1)+Num(2)+1
        res(i,1)=res(i,1) - coef*(Nacc+ni*exp(phi(i,1)/thermal))*0.5;
        Jaco(i,i) = Jaco(i,i) - coef*ni*exp(phi(i,1)/thermal)/thermal*0.5;
    else
        res(i,1) = res(i,1) - coef*(Nacc+ni*exp(phi(i,1)/thermal));
        Jaco(i,i) = Jaco(i,i) - coef*ni*exp(phi(i,1)/thermal)/thermal;
    end
end

update = Jaco \ (-res);
phi = phi + update;
end

%% Electron density 

%%%%%%%%%%%Electron density simply calculated with semiclassical
%%%%%%%%%%%model%%%%%%%%%%%%%%%%%%
elec_poi = zeros(sum(Num)+1,1);
for i = Num(1)+1:Num(1)+Num(2)+1
elec_poi(i,1) = ni * exp(phi(i,1)/thermal);
end
elec2_poi(Vg1,:)=elec_poi;
elec_den_poi(Vg1) = sum(elec_poi*1E-14); % /cm^2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Electron density calculated with Schrodinger
%%%%%%%%%% equation%%%%%%%%%%%%
for loop = 1:100
V = q*Ec_Ei-q*phi; % Potential energy, J

Nbulk = interface(2)-interface(1)-1;
Hamil = zeros(Nbulk,Nbulk);

Hamil(1,1) = -2; Hamil(1,2) = 1;
for i = 2:Nbulk-1
    Hamil(i,i+1) = 1;
    Hamil(i,i) = -2;
    Hamil(i,i-1) = 1;
end
Hamil(Nbulk,Nbulk)=-2; Hamil(Nbulk,Nbulk-1) = 1;

for i=1:Nbulk
 Hamil(i,i) = Hamil(i,i) -2*mzz*m0*(dz/hbar)^2*V(interface(1)+i,1);
end

[eigenvectors,eigenvalues] = eig(Hamil);
scaledEz = diag(eigenvalues)/(-2*mzz*m0*(dz/hbar)^2); % Eigenenergy, J
[sortedEz,sortedIndex] = sort(scaledEz);

nSubband = 10;
elec = zeros(N,1); % Electron density, /m^3
totalNumber = 0;
for n=1:nSubband
    Ez = sortedEz(n,1);
    inn_wavefunc = eigenvectors(:,sortedIndex(n)).^2;
    norm = sum(inn_wavefunc)*dz;
    inn_wavefunc = inn_wavefunc / norm;
    subbandNumber = coef_Sch*log(1+exp(-(Ez)/(k_B*T)));
    totalNumber = totalNumber + subbandNumber;
    elec(interface(1)+1:interface(2)-1, 1) = elec(interface(1)+1:interface(2)-1, 1) + 1/(Lx*Ly)*inn_wavefunc.*subbandNumber;
end

    
Jaco=sparse(sum(Num)+1,sum(Num)+1);
Jaco(1,1)=1;
Jaco(sum(Num)+1,sum(Num)+1)=1;

res = zeros(N,1);
res(1,1) = phi(1,1) - (0.33374+Vg);
res(N,1) = phi(N,1) - (0.33374+Vg);
sum_num=1;
for i = 1:numel(mat_th)
    for j = sum_num+1:sum_num+Num(i)-1
        Jaco(j,j-1) = e(i);
        Jaco(j,j) = -2*e(i);
        Jaco(j,j+1) = e(i);
        
        res(j,1) = e(i)*phi(j+1,1) - 2*e(i)*phi(j,1) + e(i)*phi(j-1,1);
        
    end
    
    if i < numel(mat_th)
    sum_num=sum_num+Num(i);
    
    Jaco(sum_num,sum_num-1) = e(i);
    Jaco(sum_num,sum_num) = -e(i) -e(i+1);
    Jaco(sum_num,sum_num+1) = e(i+1);
    
    res(sum_num,1) = e(i+1)*phi(sum_num+1,1) -(e(i)+e(i+1))*phi(sum_num,1) + e(i)*phi(sum_num-1,1);
    
    end
end

% Charge
for i = interface(1):interface(2)
    if i==interface(1)
        res(i,1) = res(i,1) - coef*(Nacc+elec(i,1))*0.5;
        Jaco(i,i) = Jaco(i,i) - coef*elec(i,1)/thermal*0.5;
    elseif i==interface(2)
        res(i,1)=res(i,1) - coef*(Nacc+elec(i,1))*0.5;
        Jaco(i,i) = Jaco(i,i) - coef*elec(i,1)/thermal*0.5;
    else
        res(i,1) = res(i,1) - coef*(Nacc+elec(i,1));
        Jaco(i,i) = Jaco(i,i) - coef*elec(i,1)/thermal;
    end
end

update = Jaco \ (-res);
phi = phi + update;

end
elec2(Vg1,:)=elec;
elec_den(Vg1) = sum(elec*1E-14); % /cm^2
end

figure(201)
plot([0:10]/10, elec_den, [0:10]/10, elec_den_poi)
legend('Poisson+Schrodinger', 'Poisson')
set(gca, 'YScale', 'log'); xlabel('Gate voltage(V)'); ylabel('electron density (/cm^2)')

figure(202)
plot(z/1e-9,elec2(1,:)/1e6,z/1e-9,elec2(2,:)/1e6,z/1e-9,elec2(3,:)/1e6,z/1e-9,elec2(4,:)/1e6,z/1e-9,elec2(5,:)/1e6,z/1e-9,elec2(6,:)/1e6,z/1e-9,elec2(7,:)/1e6,z/1e-9,elec2(8,:)/1e6,z/1e-9,elec2(9,:)/1e6,z/1e-9,elec2(10,:)/1e6,z/1e-9,elec2(11,:)/1e6)
legend('0.0 V', '0.1 V', '0.2 V', '0.3 V', '0.4 V', '0.5 V', '0.6 V', '0.7 V', '0.8 V', '0.9 V', '1.0 V');
xlabel('position (nm)')
ylabel('electron density cm^-^3')

figure(203)
plot(z/1e-9,elec2_poi(1,:)/1e6,z/1e-9,elec2_poi(2,:)/1e6,z/1e-9,elec2_poi(3,:)/1e6,z/1e-9,elec2_poi(4,:)/1e6,z/1e-9,elec2_poi(5,:)/1e6,z/1e-9,elec2_poi(6,:)/1e6,z/1e-9,elec2_poi(7,:)/1e6,z/1e-9,elec2_poi(8,:)/1e6,z/1e-9,elec2_poi(9,:)/1e6,z/1e-9,elec2_poi(10,:)/1e6,z/1e-9,elec2_poi(11,:)/1e6)
legend('0.0 V', '0.1 V', '0.2 V', '0.3 V', '0.4 V', '0.5 V', '0.6 V', '0.7 V', '0.8 V', '0.9 V', '1.0 V');
xlabel('position (nm)')
ylabel('electron density cm^-^3')
