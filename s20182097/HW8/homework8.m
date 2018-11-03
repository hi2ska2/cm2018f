%% clear paramenter %%
clear all
clc
close all

%% parameter %%
N = 500;
Layer_1 = [0.5e-9 3.9];  %[length (m) relative permittivity]
Layer_2 = [5e-9 11.7];
Layer_3 = [0.5e-9 3.9];



A = zeros(N);
Jaco = zeros(N);
res = zeros(N,1);
phi = zeros(N,1);
mass = zeros(3,1);



Width = Layer_1(1, 1) + Layer_2(1, 1) + Layer_3(1, 1);
Lx = 100e-9;
Ly = 100e-9;
eps_0 = 8.8541878176e-12;
hbar = 6.626176e-34 / (2*pi);
q = 1.602192e-19;
m0 = 9.109534e-31;
Nacc = 1e24;
ni = 1.075e16;
ni_before = 0;
ni_after = 0;
K_B = 1.380662e-23;
T = 300;
del_z = Width/(N-1);
coef_Poi = del_z*del_z*q/eps_0;

Ec_Ei = 0.561004;

P_1 = 1;
P_2 = N;

while del_z*(P_1 - 1) <= Layer_1(1, 1)
    P_1 = P_1 + 1;
end

while del_z*(P_2 - 1) >= Layer_1(1, 1) + Layer_2(1, 1)
    P_2 = P_2 - 1;
end



phi(:,1) = 0.33374;
Jaco(1,1) = 1;
Jaco(N, N) = 1;

%% Schrodinger-Poisson solver %%

% Voltage change %
for voltage = 1:30
   
    gate_voltage = 1/29*(voltage-1);
    phi(1,1) = 0.33374 + gate_voltage;
    phi(N,1) = 0.33374 + gate_voltage;
    res(1,1) = phi(1,1) - 0.33374 - gate_voltage;
    res(N,1) = phi(N,1) - 0.33374 - gate_voltage;
    
% Nonlinear Poisson solver %
    for Newton = 1:100

        for i = 2 : N-1

                Position_before = del_z*((i - 1) - 0.5);
                Position_after = del_z*((i - 1) + 0.5);

                if (Position_before < Layer_1(1, 1)) || (Position_before > Layer_1(1, 1)+ Layer_2(1, 1)) 
                    ni_before = 0;

                elseif (Position_after < Layer_1(1, 1)) || (Position_after > Layer_1(1, 1)+ Layer_2(1, 1)) 
                    ni_after = 0;

                else 
                    ni_before = 0.5*ni;
                    ni_after = 0.5*ni;
                end

                ni_eff = ni_before + ni_after;

            for j = 1 : N

                if i == j
                    Jaco(i, j) = -(r_permitivity(Position_before, Layer_1, Layer_2, Layer_3) + r_permitivity(Position_after, Layer_1, Layer_2, Layer_3)) - coef_Poi*ni_eff*q/(K_B*T)*exp(q*phi(i,1)/(K_B*T));

                elseif i == j + 1 
                    Jaco(i, j) = r_permitivity(Position_before, Layer_1, Layer_2, Layer_3);

                elseif i == j - 1  
                    Jaco(i, j) = r_permitivity(Position_after, Layer_1, Layer_2, Layer_3);

                end
            end

                res(i,1) = r_permitivity(Position_after, Layer_1, Layer_2, Layer_3)*(phi(i+1, 1) - phi(i, 1)) - r_permitivity(Position_before, Layer_1, Layer_2, Layer_3)*(phi(i,1) - phi(i-1,1)) - coef_Poi*(Nacc + ni_eff*exp(q*phi(i,1)/(K_B*T)));  
        end

    update = Jaco\(-res);
    phi = phi + update;
    
    end    




    
% Schrodinger - Poisson Newton method %
    for Newton2 = 1:10 
        elec = zeros(N,1);
        total_electron = 0;
    
    


% 6-vally %
        for massii = 1 : 3
            mass = [0.19; 0.19; 0.19];
            mass(massii,1) = 0.91;
            coef_Sch = 2*Lx*Ly/(2*pi)*sqrt(mass(1,1)*mass(2,1))*m0/(hbar)^2*K_B*T;

        % Schrodinger solver %

            V = q*Ec_Ei - q*phi;
            Hamil = zeros(P_2 - P_1 +1);

            for ii = 2 : P_2 - P_1

                Hamil(ii, ii) = -2;
                Hamil(ii, ii - 1) = 1;
                Hamil(ii, ii + 1) = 1;  
                
            end

            Hamil(1, 1) = -2;
            Hamil(1, 2) = 1;
            Hamil(P_2 - P_1 + 1, P_2 - P_1 + 1) = -2;
            Hamil(P_2 - P_1 + 1, P_2 - P_1 ) = 1;

            for jj = 1 : P_2 - P_1 + 1

                Hamil(jj, jj) = Hamil(jj, jj) - 2*mass(3,1)*m0*(del_z/hbar)^2*V(jj + P_1 - 1, 1);

            end

            [eigenvector, eigenvalue] = eig(Hamil);
            ScaledEz = diag(eigenvalue)/ (-2*mass(3,1)*m0*(del_z/hbar)^2);
            [SortedEz, Sortedindex] = sort(ScaledEz);

        % Electron density %

            nSubband = 10;


            for n = 1:nSubband

                Ez = SortedEz(n,1);
                Probability = eigenvector(:, Sortedindex(n)).^2;
                normalization = sum(Probability)*del_z;
                Probability = Probability/normalization;
                Subbandnumber = coef_Sch*log(1+exp(-Ez/(K_B*T)));
                elec(P_1:P_2, 1) = elec(P_1:P_2, 1) + 1/(Lx*Ly)*Probability*Subbandnumber;
                total_electron = total_electron + Subbandnumber;
            end
        end

     % Poisson solver %

        for m = 1 : N
            for l = 1 : N
                Position_before = Width/(N-1)*((m - 1) - 0.5);
                Position_after = Width/(N-1)*((m - 1) + 0.5);

                if (m == 1) || (m == N)
                    A(m, l) = 0;

                elseif m == l
                    A(m, l) = -(r_permitivity(Position_before, Layer_1, Layer_2, Layer_3) + r_permitivity(Position_after, Layer_1, Layer_2, Layer_3));

                elseif m == l + 1 
                    A(m, l) = r_permitivity(Position_before, Layer_1, Layer_2, Layer_3);

                elseif m == l - 1  
                    A(m, l) = r_permitivity(Position_after, Layer_1, Layer_2, Layer_3);

                end
            end
        end

        A(1,1) = 1;
        A(N, N) = 1;

        b = q*(Nacc + elec)*(del_z)^2/eps_0;
        b(1, 1) = 0.33374 + gate_voltage;
        b(N, 1) = 0.33374 + gate_voltage;

        phi = A\b;        

    end

Potential(:,voltage) = phi;
integrated_electron(:,voltage) = total_electron/(Lx*Ly);
electron_density(:,voltage) = elec;


end        



        
        
