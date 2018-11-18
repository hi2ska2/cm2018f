close all;
clear all;

q = 1.602192e-19;
k_B = 1.380662e-23;
T = 300;
tau = 0.1e-12;
m0 = 9.109534e-31;
m_eff = 0.8*m0;

VD = 0.01;
N = 301;
Inter1 = 101;
Inter2 = 201;
H = linspace(0.001, 0.1, 200);
x = linspace(0,300,301);

for j = 1:200
      fs = sqrt(2*pi)/(1+exp(q*H(j)/(k_B*T)));
      fd = sqrt(2*pi)/(1+exp(q*(H(j)+VD)/(k_B*T)));
      V = zeros(N,1);
      V(Inter1:Inter2,1) = [0:1/(Inter2-Inter1):1]*VD;
      V(Inter2:N,1) = VD;
      A = zeros(N,N); A(1,1) = 1;
      for i = 2:N-1
          c1 = H(j)+0.5*(V(i,1)+V(i-1,1));
          c2 = H(j)+0.5*(V(i+1,1)+V(i,1));
          A(i,i-1) = c1; A(i,i) = -c1-c2; A(i,i+1) = c2;
      end
      A(N,N) = 1;
      b = zeros(N,1);
      b(1,1) = fs;
      b(N,1) = fd;
      f0 = A\b;
      f0_tot(j,:) = f0;
      f1(j,:) = -sqrt(pi)*tau*sqrt(2.*(H(j)+q*V)/m_eff).*gradient(f0);
    end
      
figure(1);
mesh(x/10,log(H),f0_tot);
xlabel('Position (nm)'); ylabel('log_{10}(H) (eV)'); zlabel('f_{0}'); colorbar;
figure(2);
mesh(x/10,log(H),f1);
xlabel('Position (nm)'); ylabel('log_{10}(H) (eV)'); zlabel('f_{1}'); colorbar;      
