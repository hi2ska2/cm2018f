Nz = 5;
Ny = 9;
count=zeros(4,1);

for k=1:4
    err=1;
    Psi = zeros(Nz,Ny);
    if k==1 %Case1
        Psi(1,Ny-1)=1; Psi(1,Ny)=1;
    elseif k==2 %case2
        Psi(1,1)=1; Psi(1,2)=1; 
    elseif k==3 %case3
        Psi(Nz,:)=1;
    elseif k==4%case4
        Psi(1,Ny-1)=1; Psi(1,Ny)=1;
        Psi(1,1)=1; Psi(1,2)=1;
        Psi(Nz,:)=1; 
    end
    NewPsi = Psi;
    while err>1e-8
        count(k) = count(k)+1;
        for j = 1 : Nz
            for i = 1: Ny
                if 1<i && i<Ny && 1<j && j<Nz
                    NewPsi(j,i) = 0.25*(Psi(j+1,i) + Psi(j,i+1) + Psi(j-1,i) + Psi(j,i-1));
                elseif i==1 && 1<j && j<Nz
                    NewPsi(j,i) = 0.25*(2*Psi(j,i+1)+Psi(j+1,i)+Psi(j-1,i));
                elseif i==Ny && 1<j && j<Nz
                    NewPsi(j,i) = 0.25*(2*Psi(j,i-1)+Psi(j+1,i)+Psi(j-1,i));
                elseif j==1 && 2<i && i<Ny-1
                    NewPsi(j,i) = 0.25*(2*Psi(j+1,i)+Psi(j,i+1)+Psi(j,i-1));

                end
            end
        end
        err=sqrt(sum(sum((NewPsi-Psi).^2)));
        Psi= NewPsi;
    end
    figure(k)
    surf(Psi)
    xlabel('Ny'); ylabel('Nz');
end