function [E] = HW2(N)

m0=0.511*10^(6)/(2.99*10^8)^2; % Electron rest mass 
h=6.582*10^(-16); % h-bar [eV*s]
a=5*10^(-9); % Well width
dx=a/(N-1); % Delta x

A=[];
A(1,1)=-2;
A(1,2)=1;
for n=2:N-3
    A(n,n+1)=1;
    A(n,n)=-2;
    A(n,n-1)=1;
end
A(N-2,N-3)=1;
A(N-2,N-2)=-2;

e = eig(A);
ksquare=-e(N-2)/(dx)^2;
E=h*h*ksquare/2/m0/0.19; % Unit: eV
