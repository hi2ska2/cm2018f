clear all;

ny = 9;
nz = 5;

A = zeros(ny*nz,ny*nz);

for iy = 2:ny-1
    for iz = 2:nz-1
        n = (iz-1)*ny+iy;
        A(n,n+1) = 1;
        A(n,n+9) = 1;
        A(n,n) = -4;
        A(n,n-1) = 1;
        A(n,n-9) = 1;
    end
end


% bottom boundary
for iy = 1:ny
    A(iy,iy) = 1;
end

% left boundary
for iz = 2:nz-1
    n = (iz-1)*ny+1;
    A(n,n + 1) = 1;
    A(n,n + ny) = 0.5;
    A(n,n) = -2;
    A(n,n - ny) = 0.5;
end

% right boundary
for iz = 2:nz-1
    n = (iz-1)*ny+ny;
    
    A(n,n + ny) = 0.5;
    A(n,n) = -2;
    A(n,n - ny) = 0.5;
    A(n,n - 1) = 1;
end

% top boundary
for iy = 3:ny-2
    n = (nz-1)*ny+iy;
    A(n,n + 1) = 0.5;
    A(n,n) = -2;
    A(n,n - ny) = 1;
    A(n,n - 1) = 0.5;
end

% top left
for iy = 1:2
    n = (nz-1)*ny+iy;
    A(n,n) = 1;
end

% top right
for iy = ny-1:ny
    n = (nz-1)*ny+iy;
    A(n,n) = 1;
end

b = zeros(ny*nz,1);


%% Boundary condition
% Bottom
for iy = 1:ny
    b(iy,1) = 1;
end

% top left
for iy = 1:2
    n = (nz-1)*ny+iy;
    b(n,1) = 1;
end

% top right
for iy = ny-1:ny
    n = (nz-1)*ny+iy;
    b(n,1) = 1;
end

x = A\b;
psi = zeros(ny,nz);
for i = 1:nz
psi(1:9,i)=x((i-1)*ny+1:(i-1)*ny+9);
end

figure(201)
surf(psi');
xlabel('y-axis(delta)');
ylabel('z-axis(delta)');
colorbar;