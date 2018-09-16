t1 = 2*10^-9; % GaN thickness (unit: m)
t2 = 27*10^-9; % AlGaN thickness (unit: m)
e1 = 8.9*8.854*10^-12; % GaN permittivity (unit: F/m)
e2 = 8.8*8.854*10^-12; % AlGaN permittivity (unit: F/m)

T1 = t1/10^-9;
T2 = t2/10^-9;
N = T1+T2+1; % Point values

% Matrix form for approximation %
A = zeros(N,N);
A(1,1) = 1; A(N,N) = 1;
A(2,1) = e1; A(2,2) = -2*e1; A(2,3) = e1;
A(3,2) = e1; A(3,3) = -e2-e1; A(3,4) = e2;
for i = 4:N-1
    A(i,i-1) = e2;
    A(i,i) = -2*e2;
    A(i,i+1) = e2;
end

b = zeros(N,1); b(N,1) = 1;
x = inv(A)*b;
Vmatrix = x(3,1); % Result

a = [0:N-1];
stem(a,x);
xlabel('Distance (nm)');
ylabel('Electric potential (V)');

C1 = e1/t1; 
C2 = e2/t2;
Vformular = C2/(C1+C2); % Analytic solution

Error = (abs(Vmatrix - Vformular))/Vformular*100; % Comparision (error)
