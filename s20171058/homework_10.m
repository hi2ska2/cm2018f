clear all;

q = 1.602192e-19; % Elementary charge, C
k_B = 1.380662e-23; % Boltzmann constant, J/K
T = 300.0; % Temperature, K
v = 1e6; % ~ Velocity of electron of graphite
tau = 1e-12; % ~ Relaxation time of graphite
VD = 0.01; % (V)
N = 301;
interface1 = 101;
interface2 = 201;
H = linspace(0.001, 0.1, 200);
x = linspace(0,300,301);
for j = 1:200
fs = sqrt(2*pi)./(1 + exp(q*H(j)/(k_B*T)));
fd = sqrt(2*pi)./(1 + exp(q*(H(j)+VD)/(k_B*T)));

V = zeros(N,1);
V(interface1:interface2,1) = [0:1/(interface2-interface1):1]*VD;
V(interface2:N,1) = VD; 
A = zeros(N,N);
A(1,1) = 1.0;
for ii=2:N-1
 c1 = H(j) + 0.5*(V(ii,1)+V(ii-1,1));
 c2 = H(j) + 0.5*(V(ii+1,1)+V(ii,1));
 A(ii,ii-1) = c1; A(ii,ii) = -c1-c2; A(ii,ii+1) = c2;
end
A(N,N) = 1.0;
b = zeros(N,1);
b(1,1) = fs;
b(N,1) = fd; 
f0 = A \ b;
f_0(j,:)=f0;

f1(j,:) = -tau*v*(1/sqrt(pi))*gradient(f0)*(1/sqrt(2*pi));

end

figure(201)
mesh(x/10,H,f_0);

xlabel('x (nm)');
ylabel('H (eV)');
zlabel('f_{0}');

figure(202)
mesh(x/10,H,f1);

xlabel('x (nm)');
ylabel('H (eV)');
zlabel('f_{1}');
