close all;
clear all;

R = 2e6;
C = 5e-12;
freq = 1e6;
deltat = 1/freq/100;

A = zeros(5,5);
A(1,:) = [0 0 0 1 0];
A(2,:) = [0 1 0 -C/deltat C/deltat];
A(3,:) = [0 0 1 0 -1/R];
A(4,:) = [1 1 0 0 0];
A(5,:) = [0 -1 1 0 0];

b = zeros(5,1);
solution = [0 0 0 1 0]';

N = 1000;

for i = 1:N
    t = i*deltat;
    solution_old = solution;
    b(1,1) = cos(2*pi*freq*t);
    b(2,1) = -C/deltat*(solution_old(4,1)-solution_old(5,1));
    solution = A\b;
    
    time(i,1) = t;
    solution_num(:,i) = solution;
end

w=2*pi*freq;
t = i*deltat;
I=(w.^2*R*C.^2./(1+(w.*R.*C).^2)).*cos(w.*time)-(w.*C./(1+(w.*R.*C).^2)).*sin(w.*time);

figure
plot (time,solution_num(3,:),'r',time,I,'b');
xlabel('Time (sec)')
ylabel('Current (A)')
legend('Numerical','Analytical')
