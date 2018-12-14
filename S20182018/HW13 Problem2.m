
clear all;
R = 2e6; % Ohm
C = 5e-12; % F
for jj=0:5
freq = 10^(jj); % Hz
w=2*pi*freq;
deltat = 1/freq/100; % 0.01 of a period

A = zeros(5,5);
A(1,:) = [0 0 0 1 0];
A(2,:) = [0 1 0 -C/deltat C/deltat];
A(3,:) = [0 0 1 0 -1/R];
A(4,:) = [1 1 0 0 0];
A(5,:) = [0 -1 1 0 0];

b = zeros(5,1);
solution = [0 0 0 1 0]';
N = 1000;
for ii=1:N
    t = ii*deltat;
    x(ii,1)=t;
    solution_old = solution;
    b(1,1) = cos(w*t);
    b(2,1) = -C/deltat*(solution_old(4,1)-solution_old(5,1));
    
    solution = A \ b;
    I(ii,1) = solution(3,1);
    c(ii,1) = w*w*R*C*C*cos(w*t)/(1+(w*R*C)^2) - w*C*sin(w*t)/(1+(w*R*C)^2);
end
figure(jj+1)
plot(x(:,1),I(:,1),'o'); hold on;
plot(x(:,1),c(:,1));
xlabel('Time (sec)');
ylabel('Current (A)');
end
