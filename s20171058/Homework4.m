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


%%
ni = 1.075e16; % 1.075e10 /cm^3
Nacc = 1e24; % 1e18 /cm^3
dx = 0.1e-9;
N_point = 5; % N points in minimum thickness

mat_th = [5 50 5]; % material thickness (x0.1 nm)
e = [3.9 11.7 3.9]; % permittivity

d = mat_th(1);
for n = 2:numel(mat_th)
  d = gcd(d,mat_th(n)); % greatest common divisor
end

Num = mat_th.*N_point./d; % number of points in each material

x = dx*transpose([0:sum(Num)]); % real space thickness

A=zeros(sum(Num)+1,sum(Num)+1);
A(1,1)=1;
A(sum(Num)+1,sum(Num)+1)=1;

sum_num=1;
for i = 1:numel(mat_th)
    for j = sum_num+1:sum_num+Num(i)-1
        A(j,j-1)= e(i);
        A(j,j) = -2*e(i);
        A(j,j+1)= e(i);
    end
    
    if i<numel(mat_th)
    sum_num=sum_num+Num(i);
    A(sum_num,sum_num-1) = e(i);
    A(sum_num,sum_num) = -e(i) -e(i+1);
    A(sum_num,sum_num+1) = e(i+1);
    end
end

% Boundary condition
b = zeros(sum(Num)+1,1);
b(1,1) = 0.33374;
b(sum(Num)+1,1) = 1;
for i = Num(1)+1:Num(1)+Num(2)+1
    if i==Num(1)+1
        b(i,1) = dx*dx*q*Nacc/e0*0.5;
    elseif i==Num(1)+Num(2)+1
        b(i,1) = dx*dx*q*Nacc/e0*0.5;
    else b(i,1) = dx*dx*q*Nacc/e0;
    end
end
b(sum(Num)+1,1) = 0.33374;
phi = A \ b;

elec = zeros(sum(Num)+1,1);
for i = Num(1)+1:Num(1)+Num(2)+1
elec(i,1) = ni * exp(q*phi(i,1)/(k_B*T));
end

%%
b2 = zeros(sum(Num)+1,1);
b2(1,1) = 0.33374;
b2(sum(Num)+1,1) = 1;
for i = Num(1)+1:Num(1)+Num(2)+1
    if i==Num(1)+1
        b2(i,1) = dx*dx*q*(Nacc+elec(i,1))/e0*0.5;
    elseif i==Num(1)+Num(2)+1
        b2(i,1) = dx*dx*q*(Nacc+elec(i,1))/e0*0.5;
    else b2(i,1) = dx*dx*q*(Nacc+elec(i,1))/e0;
    end
end
b2(sum(Num)+1,1) = 0.33374;
phi2 = A \ b2;
pot_diff = phi2 - phi
%%
figure(1)
plot(x*1e9,phi,'r-',x*1e9,phi2,'b-')
set(gcf,'Color','w')
title('Potential','FontSize',15)
xlabel('position (nm)','FontSize',12)
ylabel('Potential (V)','FontSize',12)
% 
% figure(2)
% plot(x*1e9,elec*1e-6,'bo-')
% set(gcf,'Color','w')
% title('Electron density','FontSize',15)
% xlabel('position (nm)','FontSize',12)
% ylabel('Electron density (cm^-3)','FontSize',12)
% 
% figure(3)
% plot(x*1e9,phi2,'r-')
% set(gcf,'Color','w')
% title('Potential2','FontSize',15)
% xlabel('position (nm)','FontSize',12)
% ylabel('Potential (V)','FontSize',12)
% 
figure(4)
plot(x*1e9,pot_diff,'r-')
set(gcf,'Color','w')
title('Potential Difference','FontSize',15)
xlabel('position (nm)','FontSize',12)
ylabel('Potential (V)','FontSize',12)

%%

% Boundary condition
c = zeros(sum(Num)+1,11);
for j = 1:11
c(1,j) = 0.33374-0.1*(j-1);
for i = Num(1)+1:Num(1)+Num(2)+1
    if i==Num(1)+1
        c(i,j) = dx*dx*q*Nacc/e0*0.5;
    elseif i==Num(1)+Num(2)+1
        c(i,j) = dx*dx*q*Nacc/e0*0.5;
    else c(i,j) = dx*dx*q*Nacc/e0;
    end
end
c(sum(Num)+1,j) = 0.33374-0.1*(j-1);
end
phic = A \ c;

elecc = zeros(sum(Num)+1,1);
for i = Num(1)+1:Num(1)+Num(2)+1
    for j =1:11
elecc(i,j) = ni * exp(q*phic(i,j)/(k_B*T));
    end
end


c2 = zeros(sum(Num)+1,11);
for j = 1:11
c2(1,j) = 0.33374-0.1*(j-1);
for i = Num(1)+1:Num(1)+Num(2)+1
    if i==Num(1)+1
        c2(i,j) = dx*dx*q*(Nacc+elecc(i,j))/e0*0.5;
    elseif i==Num(1)+Num(2)+1
        c2(i,j) = dx*dx*q*(Nacc+elecc(i,j))/e0*0.5;
    else c2(i,j) = dx*dx*q*(Nacc+elecc(i,j))/e0;
    end
end
c2(sum(Num)+1,j) = 0.33374-0.1*(j-1);
end
phi2c = A \ c2;
pot_diffc = phi2c - phic
figure(5)
plot(x,pot_diffc(:,1),x,pot_diffc(:,2),x,pot_diffc(:,3),x,pot_diffc(:,4),x,pot_diffc(:,5),x,pot_diffc(:,6),x,pot_diffc(:,7),x,pot_diffc(:,8),x,pot_diffc(:,9),x,pot_diffc(:,10),x,pot_diffc(:,11))
legend('0 V','0.1 V','0.2 V','0.3 V','0.4 V','0.5 V','0.6 V','0.7 V','0.8 V','0.9 V','1.0 V')
set(gcf,'Color','w')
title('Potential Difference','FontSize',15)
xlabel('position (m)','FontSize',12)
ylabel('Potential difference (V)','FontSize',12)
