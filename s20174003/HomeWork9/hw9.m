close all;
clear all;

Ny = 9;
Nz = 5;
Situation =4;

Phi = zeros(Nz,Ny);
%Phi = ones(Nz,Ny);
if(Situation ==1)
    Phi(1,Ny) = 1;
    Phi(1,Ny-1) = 1;
elseif(Situation ==2)
    Phi(1,2) = 1;
    Phi(1,1) = 1;
elseif(Situation == 3)
    Phi(Nz,:) = 1;
else
    Phi(1,Ny) = 1;
    Phi(1,Ny-1) = 1;
    Phi(1,2) = 1;
    Phi(1,1) = 1;
    Phi(Nz,:) = 1;
end

for iter=1:1:1
    for i=1:1:Nz
       for j=1:1:Ny 
           if( i == 1 && 2<j && j<Ny-1)
               Phi(i,j) = (2*Phi(i+1,j) + Phi(i,j+1) + Phi(i,j-1))/4;
           elseif( (1<i && i<Nz ) && (j==1))
               Phi(i,j) = (2*Phi(i,j+1) + Phi(i-1,j) + Phi(i+1,j))/4;
           elseif( (1<i && i<Nz ) && (j==Ny))
               Phi(i,j) = (2*Phi(i,j-1) + Phi(i-1,j) + Phi(i+1,j))/4;    
           elseif( (1<i && i<Nz ) && ( 1<j && j<Ny))
               Phi(i,j) = (Phi(i,j+1) + Phi(i,j-1) + Phi(i+1,j)+Phi(i-1,j))/4;
           end 
       end
    end
end

figure
surf(Phi)
xlabel('Y');
ylabel('Z');
colorbar
