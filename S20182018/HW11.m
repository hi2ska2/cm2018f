clear all;

q = 1.602192e-19; % Elementary charge, C
k_B = 1.380662e-23; % Boltzmann constant, J/K
T = 300.0; % Temperature, K
v = 1e+6; % Velocity
t = 1e-12; % Relaxation time

N = 301;
M = 301;
H = (linspace(0.1,1,M))'; % (eV)
VD = 0.01; % (V)
interface1 = 101;
interface2 = 201;
f1 = zeros(N-2,M);
x = (linspace(0,N-1,N)*1e-1)'; % Unit nm
x2 = (linspace(1,N-2,N-2)*1e-1)'; % Unit nm


for jj=1:M

    fs(jj,1) = sqrt(2*pi)/(1 + exp(q*H(jj,1)/(k_B*T)));
    fd(jj,1) = sqrt(2*pi)/(1 + exp(q*(H(jj,1)+VD)/(k_B*T)));

    V = zeros(N,1);
    V(interface1:interface2,1) = [0:1/(interface2-interface1):1]*VD;
    V(interface2:N,1) = VD;

    A = zeros(N,N);
    A(1,1) = 1.0;
    
    for ii=2:N-1
        c1 = H(jj,1) + 0.5*(V(ii,1)+V(ii-1,1));
        c2 = H(jj,1) + 0.5*(V(ii+1,1)+V(ii,1));
        A(ii,ii-1) = c1; A(ii,ii) = -c1-c2; A(ii,ii+1) = c2;
    end
    A(N,N) = 1.0;

    b = zeros(N,1);
    b(1,1) = fs(jj,1);
    b(N,1) = fd(jj,1);

    f0(:,jj) = A \ b;

    for kk=1:N-2
        f1(kk,jj) = -t*v*(1/sqrt(pi))*(f0(kk+2,jj) - f0(kk,jj))/(2e-10);
    end
end
figure(1)
mesh(H*1000,x,f0)
xlabel('H (meV)'); ylabel('Position (nm)'); zlabel('f_{0}');
figure(2)
mesh(H*1000,x2,f1)
xlabel('H (meV)'); ylabel('Position (nm)'); zlabel('f_{1}');

