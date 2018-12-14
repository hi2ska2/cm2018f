%% HW_13-2 RC filter
clear;
clc;


% define some constants

R = 2e6; % Ohm
C = 5e-12; % F
freq = 1e6; % Hz
delta_t = 1/freq/100; % 0.01 of a period

% The system matrix

% Numerical sol

A = zeros(5,5);
A(1,:) = [0 0 0 1 0];
A(2,:) = [0 1 0 -C/delta_t C/delta_t];
A(3,:) = [0 0 1 0 -1/R];
A(4,:) = [1 1 0 0 0];
A(5,:) = [0 -1 1 0 0];

b = zeros(5,1);
solution = [0 0 0 1 0]';

N = 1000;

for ii = 1:N
    t = ii*delta_t;
    solution_old = solution;
    b(1,1) = cos(2*pi*freq*t);
    b(2,1) = -C/delta_t*(solution_old(4,1) - solution_old(5,1));
    solution = A \ b;
    time(ii,1) = t;
    solution_num(:,ii) = solution;
end

% Analytic sol

w = 2*pi*freq;
I_anal = zeros(N,1);
for ii = 1:N
    I_anal(ii,1) = (w^2*R*C^2)*cos(w*time(ii,1))/(1 + (w*R*C)^2) - w*C*sin(w*time(ii,1))/(1 + (w*R*C)^2);
end

    



figure
plot(time,solution_num(3,:),'r',time,I_anal,'b')
xlabel('Time (s)')
ylabel('Current (A)')
legend('numerical','analytical')
title('Frequency = 1 MHz')




