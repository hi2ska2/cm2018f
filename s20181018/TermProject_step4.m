close all;
clear all;

set(0, 'defaultFigureRenderer', 'opengl');

h = 6.626176e-34;
hbar = h/(2*pi);
q = 1.602192e-19;
m0 = 9.109534e-31;
eps0 = 8.854187817e-12;
k_B = 1.380662e-23;
T = 300;
thermal = k_B*T/q;
Lx = 100e-9; Ly = 100e-9; Lz = 5e-9;
mxx = 0.19; myy = 0.19; mzz = 0.91;
Inter1 = 2;
Inter2 = 12;
eps_si = 11.7;
eps_ox = 3.9;
Nacc = 1e24;
ni = 1.075e16;
coef_Sch = 2*Lx*Ly/(2*pi)*sqrt(mxx*myy)*m0/(hbar^2)*(k_B*T);
Vg_sweep = [0.1, 0.5]';
Vd = [0.1, 0.5]';
Ec_Ei = 0.561004;

Nx = 13; Deltax = 0.5e-9;
Ny = 241; Deltay = 0.5e-9; y_12 = 81; y_23 = 161;
x = Deltax*([0:Nx-1])';
y = Deltay*([0:Ny-1])';
Ndon = 2e23*ones(1,Ny);
Ndon(1,1:y_12) = 5e25;
Ndon(1,y_23:Ny) = 5e25;
coef_Poi = Deltax*Deltax*q/eps0;

H = transpose(linspace(0.01,0.1,Ny));
m0 = 9.109534e-31;
m_eff=0.5*m0;
tau=0.1e-12;

n_tot = zeros(Nx,Ny);
mobility = 1430;
 
Vg = Vg_sweep+0.33374;

n_sweep_Poi = zeros(Nx,length(Vg));
n_sweep_Sch = zeros(Nx,length(Vg));
 
for i = 1:length(Vg)
    for s = 1:length(Vd)
    A = zeros(Nx,Ny);
    for kk = 81:161
        A(1,kk) = Vg(i,1);
        A(11,kk) = Vg(i,1);
    end
    for kkk = 1:13
        A(kkk,241) = Vd(s,1);
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
    for ll = 1:241
    phi = phi_tot(:,ll);
    for iNewton = 1:10
        totalNumber = 0;
        for iValley = 1:3
            mass = ones(3)*0.19;
            mass(iValley) = 0.91;
            coef_Sch = 2*Lx*Ly/(2*pi)*sqrt(mass(1)*mass(2))*m0/(hbar^2)*(k_B*T);
                V = q*Ec_Ei -q*phi;
                Nbulk = Inter2-Inter1-1;
                Hamil = zeros(Nbulk,Nbulk);
                Hamil(1,1) = -2; Hamil(1,2) = 1;
                Hamil(Nbulk,Nbulk) = -2; Hamil(Nbulk,Nbulk-1) = 1;
                for ii = 2:Nbulk-1
                    Hamil(ii,ii+1) =  1;
                    Hamil(ii,ii  ) = -2;
                    Hamil(ii,ii-1) =  1;
                end
                for ii = 1:Nbulk
                    Hamil(ii,ii) = Hamil(ii,ii)-2*mass(3)*m0*(Deltax/hbar)^2*V(ii+Inter1,1);
                end
                [eigenvectors,eigenvalues] = eig(Hamil);
                scaledEz = diag(eigenvalues)/(-2*mass(3)*m0*(Deltax/hbar)^2);
                [sortedEz,sortedIndex] = sort(scaledEz);
                nSubband = 9;
                n_Sch = zeros(Nx,1);
                totalNumber = 0;        
                for n = 1:nSubband
                    Ez = sortedEz(n,1);
                    wavefunction2 = eigenvectors(:,sortedIndex(n)).^2;
                    normalization = sum(wavefunction2)*Deltax;
                    wavefunction2 = wavefunction2 / normalization;
                    subbandNumber = coef_Sch*log(1+exp(-Ez/(k_B*T)));
                    totalNumber = totalNumber+subbandNumber;
                    n_Sch(Inter1+1:Inter2-1,1) = n_Sch(Inter1+1:Inter2-1,1)+1/(Lx*Ly)*wavefunction2*subbandNumber;
                end
        end
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
    end
    end
    for jj= 1:length(H)
        fs = sqrt(2*pi)/(1 + exp(q*H(jj,1)/(k_B*T))); 
        fd = sqrt(2*pi)/(1 + exp(q*(H(jj,1)+Vd(s,1))/(k_B*T))); 
        V = phi;
        V(Inter1:Inter2,1) = [0:1/(Inter2-Inter1):1]*Vd(s,1); 
        V(Inter2:Nx,1) = Vd(s,1); 
        C = zeros(Nx,Nx); 
        C(1,1) = 1.0;  
        for ii=2:Nx-1
            c1 = H(jj,1) + 0.5*(V(ii,1)+V(ii-1,1));     
            c2 = H(jj,1) + 0.5*(V(ii+1,1)+V(ii,1));     
            C(ii,ii-1) = c1;
            C(ii,ii) = -c1-c2; 
            C(ii,ii+1) = c2; 
        end
        C(Nx,Nx) = 1.0;
        b = zeros(Nx,1); 
        b(1,1) = fs; 
        b(Nx,1) = fd;  
        f0(:,jj) = C \ b;
    end
    figure
    mesh(y,x,f0);
    xlabel('Length of Structure [nm]'); ylabel('Thickness of Structure [nm]'); zlabel('f0');
    end
end
