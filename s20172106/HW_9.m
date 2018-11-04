clear;
clc;

ny = 9;
nz = 5;
N = ny*nz;


psi_1 = zeros(ny,nz); %% Case 1
psi_2 = zeros(ny,nz); %% Case 2
psi_3 = zeros(ny,nz); %% Case 3
psi_4 = zeros(ny,nz); %% Case 4

b_init = zeros(N,1);

A = zeros(N,N);


% Bulk terms

for i = 2:ny-1 %% y index
    for j = 2:nz-1 %% z index
        index = i + (ny)*(j-1);
        A(index,index - ny) = 1.0;
        A(index,index - 1) = 1.0;
        A(index,index) = -4.0;
        A(index,index + 1) = 1.0;
        A(index,index + ny) = 1.0;
    end
end

% Bottom boundary
% 
% for i = 2:ny-1
%     A(i,i - 1) = 0.5;
%     A(i,i) = -2.0;
%     A(i,i + 1) = 0.5;
%     A(i,i + ny) = 1.0;
% end

% Left boundary
% 
for j = 2:nz-1
    index = (j-1)*ny + 1;
    A(index,index + 1) = 1.0;
    A(index,index - ny) = 0.5;
    A(index,index) = -2.0;
    A(index,index + ny) = 0.5;
end
% % 
% % % Right boundary
% % 
for j = 2:nz-1
    index = (j)*ny;
    A(index,index - ny) = 0.5;
    A(index,index) = -2.0;
    A(index,index + ny) = 0.5;
    A(index,index - 1) = 1.0;
end
% % 
% % % Top boundary
% % 
for i = 3:ny-2
    index = (nz-1)*ny + i;
    A(index,index - 1) = 0.5;
    A(index,index) = -2.0;
    A(index,index + 1) = 0.5;
    A(index,index - ny) = 1.0;
end
% % 
% % Boundary conditions
% 
A((nz-1)*ny + 1, (nz-1)*ny + 1) = 1.0; % 
A((nz-1)*ny + 2, (nz-1)*ny + 2) = 1.0; % 
A(nz*ny - 1, nz*ny - 1) = 1.0; % 
A(nz*ny, nz*ny) = 1.0; % 
% % 

% Bottom boundary condition for blue circle

%% Case 1

for i = 1: ny
    A(i,i) = 1.0;
end

% % Original position(top/right position)
% 
b = b_init;
b(ny*nz-1,1) = 1.0;
b(ny*nz,1) = 1.0;
A_re = A(:,:);
psi = inv(A_re)*b;


for jj = 1:nz
    for ii = 1:ny
        psi_1(ii,jj) = psi(ii + 9*(jj-1));
    end
end



%% Case 2

% The top/left position

b = b_init;
b((nz-1)*ny + 1,1) = 1.0;
b((nz-1)*ny + 2,1) = 1.0;
A_re = A(:,:);
psi = inv(A_re)*b;


for jj = 1:nz
    for ii = 1:ny
        psi_2(ii,jj) = psi(ii + 9*(jj-1));
    end
end

%% Case 3


% 
% % the bottom position
% 
b = b_init;
for ii = 1:ny
    b(ii,1) = 1.0;
end
A_re = A(:,:);
psi = inv(A_re)*b;


for jj = 1:nz
    for ii = 1:ny
        psi_3(ii,jj) = psi(ii + 9*(jj-1));
    end
end

%% Case 4


% 
% % All three positions above
% 
b = b_init;

b(ny*nz-1,1) = 1.0;
b(ny*nz,1) = 1.0;

b((nz-1)*ny + 1,1) = 1.0;
b((nz-1)*ny + 2,1) = 1.0;


for ii = 1:ny
    b(ii,1) = 1.0;
end

A_re = A(:,:);
psi = inv(A_re)*b;


for jj = 1:nz
    for ii = 1:ny
        psi_4(ii,jj) = psi(ii + 9*(jj-1));
    end
end




        

