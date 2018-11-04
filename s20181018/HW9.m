close all;
clear all;

Nz = 5; Ny = 9;

A1 = zeros(Nz,Ny); A2 = zeros(Nz,Ny); A3 = zeros(Nz,Ny); A4 = zeros(Nz,Ny);

A1(Nz,Ny) = 1; A1(Nz,Ny-1) = 1;
A2(Nz,1) = 1; A2(Nz,2) = 1;
A3(1,:) = 1;
A4(Nz,Ny) = 1; A4(Nz,Ny-1) = 1; A4(Nz,1) = 1; A4(Nz,2) = 1; A4(1,:) = 1;

A = A4; % (1)-A1 (2)-A2 (3)-A3 (4)-A4

B = A;

for ii = 1:1000
    for j = 1 : Nz
        for i = 1: Ny
            if     (1<i && i<Ny && 1<j && j<Nz) % Not boundary 
                B(j,i) = 0.25*(A(j+1,i)+A(j,i+1)+A(j-1,i)+A(j,i-1));
            elseif (i==1 && 1<j && j<Nz) % Left condition
                B(j,i) = 0.25*(2*A(j,i+1)+A(j+1,i)+A(j-1,i));
            elseif (i==Ny && 1<j && j<Nz) % Right condition
                B(j,i) = 0.25*(2*A(j,i-1)+A(j+1,i)+A(j-1,i));
            elseif (j==Nz && 2<i && i<Ny-1) % Top condition
                B(j,i) = 0.25*(2*A(j-1,i)+A(j,i+1)+A(j,i-1));
            end
        end
    end
    A = B;
end

surf(A);
xlabel('Ny');
ylabel('Nz');
colorbar;
