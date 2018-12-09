clear all;
close all;

R = 2e6; %Ohm
C = 5e-12; %F
%freq = 1e0; %1Hz
freq = 15e3; %15kHz
deltat = 1/freq/100; %0.01 of a period

A = zeros(5,5);
A(1,:) = [ 0 0 0 1 0];
A(2,:) = [ 0 1 0 -C/deltat C/deltat];
A(3,:) = [ 0 0 1 0 -1/R];
A(4,:) = [ 1 1 0 0 0 ];
A(5,:) = [ 0 -1 1 0 0 ];

b = zeros(5,1);

solution = [ 0 0 0 1 0 ]';
N = 1000;
%N= 10000;
Vin_t = zeros(1,N);
Vout_t = zeros(1,N);
Ir_t = zeros(1,N);
Ic_t = zeros(1,N);
Iin_t = zeros(1,N);

for iter=1:N
    t = iter*deltat;
    solution_old = solution;
    b(1,1) = cos(2*pi*freq*t);
    b(2,1) = -C/deltat*(solution_old(4,1)-solution_old(5,1));
    solution =A\b;
    %plot(t,solution(4,1),'-o');hold on; % Vin
    Vin_t(1,iter) = solution(4,1);
    Vout_t(1,iter) = solution(5,1);
    Ir_t(1,iter) = solution(3,1);
    Ic_t(1,iter) = solution(2,1);
    Iin_t(1,iter) = solution(1,1);
    %VoutVout
end

time = [deltat:deltat:N*deltat];
time_2 = [deltat:deltat:N*deltat-deltat];
Vin = cos(2*pi*freq*time);



Vout(1,1) = 0;
for iter=1:1:N-1
    Vout(1,iter+1) = (R*C)*(Vin(1,iter+1)-Vin(1,iter))/deltat;
    
end


Ir = Vout/R;
Ic = Ir;
Iin = -Ic;

figure(1)
plot(time,Vin_t,'Color','b','Marker','o')
hold on;
plot(time,Vin,'r')
title('Vin f=1Hz')
legend('Discrete','Analytic')
xlabel('time [sec]')
ylabel('Voltage')


figure(2)
plot(time,Vout_t,'Color','b','Marker','o')
hold on;
plot(time,Vout,'r')
%title('Vout f=1Hz')
title('Vout f=15kHz')
legend('Discrete','Analytic')
xlabel('time [sec]')
ylabel('Voltage')

figure(3)
plot(time,Iin_t,'Color','b','Marker','o')
hold on;
plot(time,Iin,'r')
%title('Input Current Iin f=1Hz')
title('Input Current Iin f=15kHz')
legend('Discrete','Analytic')
xlabel('time [sec]')
ylabel('Voltage')




