%% Homework 13-2 RC filter

R = 2e6; % Ohm
C = 5e-12; % F
freq = 1e5; % Hz
deltat = 1/freq/100; % 0.01 of a period


% system matrix
A = zeros(5,5);
A(1,:) = [0 0 0 1 0];
A(2,:) = [0 1 0 -C/deltat C/deltat];
A(3,:) = [0 0 1 0 -1/R];
A(4,:) = [1 1 0 0 0];
A(5,:) = [0 -1 1 0 0];

b = zeros(5,1);
solution = [0 0 0 1 0]';
N = 1000;

% Numerical solution
for ii=1:N
 t = ii*deltat;
 solution_old = solution;
 b(1,1) = cos(2*pi*freq*t);
 b(2,1) = -C/deltat*(solution_old(4,1)-solution_old(5,1));
 solution = A \ b;
 
 time(ii,1)=t;
 solution3(:,ii)=solution;
end

% Analytical solution
w=2*pi*freq;
t = ii*deltat;
I=(w.^2*R*C.^2./(1+(w.*R.*C).^2)).*cos(w.*time)-(w.*C./(1+(w.*R.*C).^2)).*sin(w.*time);

% plot(time,I)

figure
plot (time,solution3(3,:),'r',time,I,'b');
xlabel('Time (s)')
ylabel('Current (A)')
legend('numerical','analytical')