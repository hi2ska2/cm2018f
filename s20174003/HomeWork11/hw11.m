close all;
clear all;

q = 1.602192e-19;
k_B = 1.380662e-23;
T = 300.0;
tau = 1e-6; %Relaxation time
m0 = 9.109534e-31; %Electron reset mass (Kg)

%H = 0.1; %(eV)
H = [0.01:0.0025:0.0575]';%(0.01~0575eV)
%VD = 0.001; %(V) Default
VD = 0.001; %(V)
X = [0:0.1:30]; %Distance of X (nm)

df0 = zeros(301,1);

N = 301;
f0 = zeros(N,20);
f1 = zeros(N,20);

interface1 = 101;
interface2 = 201;
for iteriter=1:1:20
    fs = sqrt(2*pi)/(1 + exp(q*H(iteriter,1)/(k_B*T)));
    fd = sqrt(2*pi)/(1 + exp(q*(H(iteriter,1)+VD)/(k_B*T)));

    V = zeros(N,1);
    V(interface1:interface2,1) = [0:1/(interface2-interface1):1]*VD;
    V(interface2:N,1) = VD;

    A = zeros(N,N);
    A(1,1) = 1.0;
    for iter=2:N-1
        c1 = H(iteriter,1) + 0.5*(V(iter,1)+V(iter-1,1));
        c2 = H(iteriter,1) + 0.5*(V(iter+1,1)+V(iter,1));
        A(iter,iter-1)  = c1;
        A(iter,iter) = -c1-c2;
        A(iter,iter+1) = c2;
    end

    A(N,N) = 1.0;

    b = zeros(N,1);
    b(1,1) = fs;
    b(N,1) = fd; 

    f0(:,iteriter) = A \ b;
    for subiter=1:1:300
        df0(subiter+1,1) = (f0(subiter+1,iteriter) - f0(subiter,iteriter))/(0.1e-9);
        
    end
    for subiter=1:1:301
        f1(subiter,iteriter) = -tau * sqrt((2*(H(iteriter,1)+q*V(subiter,1)))/m0) * sqrt(pi) * df0(subiter,1);
    end
    
end
figure(1)
mesh(H,X,f0)
xlabel('H [eV]');
ylabel('Position [nm]');
zlabel('F0')
title('F0 Result')
colorbar

figure(2)
mesh(H,X,f1)
xlabel('H [eV]');
ylabel('Position [nm]');
zlabel('F1')
title('F1 Result')
colorbar

