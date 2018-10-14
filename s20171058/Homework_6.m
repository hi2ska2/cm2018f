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

thermal = k_B*T/q; % Thermal voltage, V


%%
ni = 1.075e16; % 1.075e10 /cm^3 intrinsic carrier density
Nacc = 1e24; % 1e18 /cm^3 doping density
dx = 0.1e-9;
coef = dx*dx*q/e0;
N_point = 5; % N points in minimum thickness



mat_th = [5 50 5]; % material thickness (x0.1 nm)
e = [3.9 11.7 3.9]; % permittivity

d = mat_th(1);
for n = 2:numel(mat_th)
  d = gcd(d,mat_th(n)); % greatest common divisor
end

Num = mat_th.*N_point./d; % number of points in each material
N = sum(Num)+1;

x = dx*transpose([0:sum(Num)]); % real space thickness
for Vg1 = 1 : 101
    Vg=(Vg1-1)/100;
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

figure(1)
plot(x,phi)


elec = zeros(sum(Num)+1,1);
for i = Num(1)+1:Num(1)+Num(2)+1
elec(i,1) = ni * exp(phi(i,1)/thermal);
end
elec2(Vg1,:)=elec;
elec_den(Vg1) = sum(elec*1E-14); % /cm^2
end
figure(2)
plot([0:100]/100, elec_den)
set(gca, 'YScale', 'log')