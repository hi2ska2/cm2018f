function [E] = HW2(N) % Point values (5, 50, 500)

x = 5*10^-9/(N-1); % Length of the box (unit: m)
m0 = (0.511*10^6)/(3*10^8)^2; % Electron rest mass
h = 6.582*10^-16; % h-bar (unit: eV*s)

% Matrix form with the boundary conditions %
A = zeros(N-2,N-2);
for i = 2:N-3
    A(i,i-1) = 1;
    A(i,i) = -2;
    A(i,i+1) = 1;
end
A(1,1) = -2; A(1,2) = 1;
A(N-2,N-3)=1; A(N-2,N-2)= -2;
[V,D] = eig(A);

a = D(1,1);
for k=1:N-2
    if abs(a) > abs(D(k,k))
        a=abs(D(k,k));
    end
end

ksquare = a/x^2;
E = h^2*ksquare/(2*m0*0.19); % Ground state energy (unit: eV)
