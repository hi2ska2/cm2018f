close all;
clear all;

set(0, 'DefaultFigureRenderer', 'painters');

q = 1.602192e-19;
eps0 = 8.854187817e-12;
k_B = 1.380662e-23;
T = 300;
thermal = k_B*T/q;
Inter1 = 2;
Inter2 = 12;
eps_si = 11.7;
eps_ox = 3.9;
Nacc = 2e23;
ni = 1.075e16;
Vg_sweep = transpose(linspace(0,1,11));
Nx = 13; Deltax = 0.5e-9;
Ny = 241; Deltay = 0.5e-9; y_12 = 81; y_23 = 161;
x = Deltax*([0:Nx-1])';
y = Deltay*([0:Ny-1])';
Ndon = 2e23*ones(1,Ny);
Ndon(1,1:y_12) = 5e25;
Ndon(1,y_23:Ny) = 5e25;
coef_Poi = Deltax*Deltax*q/eps0;

n_tot = zeros(Nx,Ny);
 
Vg = Vg_sweep+0.33374;
 
for i = 1:length(Vg)
    A = zeros(Nx,Ny);
    for kk = 81:161
        A(1,kk) = Vg(i,1);
        A(11,kk) = Vg(i,1);
    end
    B = A;
    for ii = 1:1000
        for kk = 1 : Ny
            for jj = 1: Nx
                if     (1<kk && kk<Ny && 1<jj && jj<Nx) % Not boundary 
                    B(jj,kk) = 0.25*(A(jj+1,kk)+A(jj,kk+1)+A(jj-1,kk)+A(jj,kk-1));
                elseif (kk==1 && 1<jj && jj<Nx) % Left condition
                    B(jj,kk) = 0.25*(2*A(jj,kk+1)+A(jj+1,kk)+A(jj-1,kk));
                elseif (kk==Ny && 1<jj && jj<Nx) % Right condition
                    B(jj,kk) = 0.25*(2*A(jj,kk-1)+A(jj+1,kk)+A(jj-1,kk));
                elseif (jj==Nx && 1<kk && kk<Ny) % Bottom condition
                    B(jj,kk) = 0.25*(2*A(jj-1,kk)+A(jj,kk+1)+A(jj,kk-1));
                elseif (jj==1 && 1<kk && kk<Ny) % Top condition
                    B(jj,kk) = 0.25*(2*A(jj+1,kk)+A(jj,kk+1)+A(jj,kk-1));
                end
            end
        end
        A = B;
    end
    phi = zeros(Nx,1);
    phi(:,1) = Vg(i,1);
    res = zeros(Nx,1);
    Jaco = sparse(Nx,Nx);
    Jaco(1,1) = 1;
    Jaco(Nx,Nx) = 1;
    res(1,1) = phi(1,1)-Vg(i,1);
    res(Nx,1) = phi(Nx,1)-Vg(i,1);
    for newton = 1:100
        for ii = 2:Nx-1
            if      (ii<Inter1 || ii>Inter2)
                res(ii,1) = eps_ox*phi(ii+1,1)-2*eps_ox*phi(ii,1)+eps_ox*phi(ii-1,1);
                Jaco(ii,ii-1) = eps_ox; Jaco(ii,ii) = -2*eps_ox; Jaco (ii,ii+1) = eps_ox;
            elseif  (ii == Inter1)
                res(ii,1) = eps_si*phi(ii+1,1)-(eps_si+eps_ox)*phi(ii,1)+eps_ox*phi(ii-1,1);
                Jaco(ii,ii-1) = eps_ox; Jaco(ii,ii) = -(eps_si+eps_ox); Jaco(ii,ii+1) = eps_si;
            elseif  (ii == Inter2)
                res(ii,1) = eps_ox*phi(ii+1,1)-(eps_ox+eps_si)*phi(ii,1)+eps_si*phi(ii-1,1);
                Jaco(ii,ii-1) = eps_si; Jaco(ii,ii) = -(eps_ox+eps_si); Jaco(ii,ii+1) = eps_ox;
            else
                res(ii,1) = eps_si*phi(ii+1,1)-2*eps_si*phi(ii,1)+eps_si*phi(ii-1,1);
                Jaco(ii,ii-1) = eps_si; Jaco(ii,ii) = -2*eps_si; Jaco(ii,ii+1) = eps_si;
            end
        end
        for ii = Inter1:Inter2
            if      (ii == Inter1)
                res(ii,1) = res(ii,1)-coef_Poi*(Nacc+ni*exp(phi(ii,1)/thermal))*0.5;
                Jaco(ii,ii) = Jaco(ii,ii)-coef_Poi*ni*exp(phi(ii,1)/thermal)/thermal*0.5;
            elseif  (ii == Inter2)
                res(ii,1) = res(ii,1)-coef_Poi*(Nacc+ni*exp(phi(ii,1)/thermal))*0.5;
                Jaco(ii,ii) = Jaco(ii,ii)-coef_Poi*ni*exp(phi(ii,1)/thermal)/thermal*0.5;
            else
                res(ii,1) = res(ii,1)-coef_Poi*(Nacc+ni*exp(phi(ii,1)/thermal));
                Jaco(ii,ii) = Jaco(ii,ii)-coef_Poi*ni*exp(phi(ii,1)/thermal)/thermal;
            end
        end
        update = inv(Jaco)*(-res);
        phi = phi+update;
        phi_Poi(:,i) = phi;
        for l = 1:241
            phi_tot(:,l) = phi;
        end
        phi_tot = phi_tot+A;
        for ii = Inter1:Inter2
            n_tot(ii,:) = ni*exp(phi_tot(ii,:)/thermal);    
        end
    end
    figure(1);
    C =jet;
    CC = C(ceil(i/length(Vg_sweep)*64),:);
    mesh(y,x,phi_tot,'FaceColor',CC,'EdgeColor','none');
    hold on;
    figure(i+1);
    CC = C(ceil(i/length(Vg_sweep)*64),:);
    mesh(y,x,n_tot*1e-6,'FaceColor',CC,'EdgeColor','none');
    xlabel('Length of Structure [nm]'); ylabel('Thickness of Structure [nm]'); zlabel('Electron Concentration [cm^-^3]');
    if (i == 1) 
        title('Gate Voltage : 0 V');
    elseif (i == 2)
        title('Gate Voltage : 0.1 V');
    elseif (i == 3)
        title('Gate Voltage : 0.2 V');
    elseif (i == 4)
        title('Gate Voltage : 0.3 V');
    elseif (i == 5)
        title('Gate Voltage : 0.4 V');
    elseif (i == 6)
        title('Gate Voltage : 0.5 V');
    elseif (i == 7)
        title('Gate Voltage : 0.6 V');
    elseif (i == 8)
        title('Gate Voltage : 0.7 V');
    elseif (i == 9)
        title('Gate Voltage : 0.8 V');
    elseif (i == 10)
        title('Gate Voltage : 0.9 V');
    elseif (i == 11)
        title('Gate Voltage : 1.0 V');
    end
end
figure(1);
xlabel('Length of Structure [nm]'); ylabel('Thickness of Structure [nm]'); zlabel('Potential [V]');
legend('0 V','0.1 V','0.2 V','0.3 V','0.4 V','0.5 V','0.6 V','0.7 V','0.8 V','0.9 V','1.0 V','Location','Best');
