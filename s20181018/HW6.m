close all;
clear all;

q = 1.602192e-19; % Elementary charge (C)
eps0 = 8.854187817e-12; % Vacuum permittivity, (F/m)
k_B = 1.380662e-23; % Boltzmann constant (J/K) 
T = 300; % Temperature (K)
thermal = k_B*T/q;
Delta_x = 0.1e-9; % 0.1 nm spacing
N = 61; % 6 nm thick
x = Delta_x*transpose([0:N-1]); % real space (m)
Inter1 = 6; % At x = 0.5 nm
Inter2 = 56; % At x = 5.5 nm
eps_si = 11.7; % Silicon relative permittivity
eps_ox = 3.9; % Silicon dioxide relative permittivity
Nacc = 1e24; % 1e18/cm^3
ni = 1.075e16; % 1.075e10/cm^-4
coef = Delta_x*Delta_x*q/eps0;

res = zeros(N,11);
Jaco = sparse(N,N);
Jaco(1,1) = 1;
Jaco(N,N) = 1;
phi = zeros(N,11);
for i = 1:11
    phi(1,i) = 0.33374-(i-1)/10;
    phi(N,i) = 0.33374-(i-1)/10;
end


for j = 1:11
    for newton = 1:100
        for i = 2:N-1
            if      (i< Inter1 || i> Inter2)
                res(i,j) = eps_ox*phi(i+1,j)-2*eps_ox*phi(i,j)+eps_ox*phi(i-1,j);
                Jaco(i,i-1) = eps_ox; Jaco(i,i) = -2*eps_ox; Jaco (i,i+1) = eps_ox;
            elseif  (i == Inter1)
                res(i,j) = eps_si*phi(i+1,j)-(eps_si+eps_ox)*phi(i,j)+eps_ox*phi(i-1,j);
                Jaco(i,i-1) = eps_ox; Jaco(i,i) = -(eps_si+eps_ox); Jaco(i,i+1) = eps_si;
            elseif  (i == Inter2)
                res(i,j) = eps_ox*phi(i+1,j)-(eps_ox+eps_si)*phi(i,j)+eps_si*phi(i-1,j);
                Jaco(i,i-1) = eps_si; Jaco(i,i) = -(eps_ox+eps_si); Jaco(i,i+1) = eps_ox;
            else
                res(i,j) = eps_si*phi(i+1,j)-2*eps_si*phi(i,j)+eps_si*phi(i-1,j);
                Jaco(i,i-1) = eps_si; Jaco(i,i) = -2*eps_si; Jaco(i,i+1) = eps_si;
            end
        end
        for i = Inter1:Inter2
            if      (i == Inter1)
                res(i,j) = res(i,j)-coef*(Nacc+ni*exp(phi(i,j)/thermal))*0.5;
                Jaco(i,i) = Jaco(i,i)-coef*ni*exp(phi(i,j)/thermal)/thermal*0.5;
            elseif  (i==Inter2)
                res(i,j) = res(i,j)-coef*(Nacc+ni*exp(phi(i,j)/thermal))*0.5;
                Jaco(i,i) = Jaco(i,i)-coef*ni*exp(phi(i,j)/thermal)/thermal*0.5;
            else
                res(i,j) = res(i,j)-coef*(Nacc+ni*exp(phi(i,j)/thermal));
                Jaco(i,i) = Jaco(i,i)-coef*ni*exp(phi(i,j)/thermal)/thermal;
            end
        end
        update(:,j) = inv(Jaco)*(-res(:,j));
        phi(:,j) = phi(:,j)+update(:,j);        
    end
end

n = zeros(N,11);
for j = 1:11
    for i = Inter1:Inter2
        n(i,j) = ni*exp(phi(i,j)/thermal);
    end
end
n_integrated = sum(n);

figure(1);
for i = 1:11
    plot(x/1e-9,phi(:,i));
    hold all;
end
title ('Potential difference over gate voltages');
xlabel('Position (nm)');
ylabel('Potential (V)');
legend('0V','0.1V','0.2V','0.3V','0.4V','0.5V','0.6V','0.7V','0.8V','0.9V','1.0V','Location','Best');

figure(2);
Voltage_sweep = 0:0.1:1;
semilogy(Voltage_sweep,n_integrated);
title('Integrated electron concentration');
xlabel('Gate voltage (V)');
ylabel('Electron concentration (cm-2)');
