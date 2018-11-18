clear all
clc

q = 1.602192e-19; % Elementary charge, C 
k_B = 1.380662e-23; % Boltzmann constant, J/K 
T = 300.0; % Temperature, K
H = transpose(linspace(0.01,0.1,300));
VD = 0.001; % V
m0 = 9.109534e-31; %electron mass in vacuum, Kg
m_eff=0.5*m0; % effectvie mass
tau=1e-5; % decay time s

N = 301; 
interface1 = 101; 
interface2 = 201;
x=transpose([0:N-1]*1e-10);

for jj= 1:length(H)
% boundary condition
fs = sqrt(2*pi)/(1 + exp(q*H(jj,1)/(k_B*T))); 
fd = sqrt(2*pi)/(1 + exp(q*(H(jj,1)+VD)/(k_B*T))); 

% potential
V = zeros(N,1); 
V(interface1:interface2,1) = [0:1/(interface2-interface1):1]*VD; 
V(interface2:N,1) = VD; 

A = zeros(N,N); 
A(1,1) = 1.0;  
for ii=2:N-1       
    c1 = H(jj,1) + 0.5*(V(ii,1)+V(ii-1,1));     
    c2 = H(jj,1) + 0.5*(V(ii+1,1)+V(ii,1));     
    A(ii,ii-1) = c1; 
    A(ii,ii) = -c1-c2; 
    A(ii,ii+1) = c2; 
end
A(N,N) = 1.0; 

b = zeros(N,1); 
b(1,1) = fs; 
b(N,1) = fd;  

f0(:,jj) = A \ b;

df0(:,jj)=gradient(f0(:,jj));
f1(:,jj)=-sqrt(pi)*tau*sqrt(2.*(H(jj,1)+q*V)./m_eff).*df0(:,jj);
end

figure
mesh(H,x*1e9,f0)
xlabel('H (eV)')
ylabel('x (nm)')
zlabel('f0')

figure
mesh(H,x*1e9,f1)
xlabel('H (eV)')
ylabel('x (nm)')
zlabel('f1')