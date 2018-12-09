clear all;


R = 2e6; % Ohm
C = 5e-12; % F
freq = 10000000e0; % Hz
dt = 1/freq/100; % 0.01 of a period
w = 2*pi*freq
A = zeros(5,5);
A = [0 0 0 1 0;
     0 1 0 -C/dt C/dt;
     0 0 1 0 -1/R;
     1 1 0 0 0;
     0 -1 1 0 0;];

b = zeros(5,1);
sol = [0 0 0 1 0]';
N = 1000;
for i=1:N
 t = i*dt;
 sol_old = sol;
 b(1,1) = cos(2*pi*freq*t);
 b(2,1) = -C/dt*(sol_old(4,1)-sol_old(5,1));
 sol = A \ b;
 Curr(i) = sol(2);
 time(i) = t;
end

I_ana = ((w.^2.*R.*C.^2)./(1+(w.*R.*C).^2)).*cos(w.*time) - ((w.*C)./(1+(w.*R.*C).^2)).*sin(w.*time);

figure(345)
plot(time, Curr,'bo', time, I_ana,'r-')
xlabel('time (s)')
ylabel('Current (A)')
