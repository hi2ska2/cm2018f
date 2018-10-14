%% clear paramenter %%
clear all
clc
close all

%% parameter %%
N = 500;
Layer_1 = [0.5e-9 3.9];  %[length (m) relative permittivity]
Layer_2 = [5e-9 11.7];
Layer_3 = [0.5e-9 3.9];


Jaco = zeros(N);
res = zeros(N,1);
phi = zeros(N,1);




Width = Layer_1(1, 1) + Layer_2(1, 1) + Layer_3(1, 1);
eps_0 = 8.8541878176e-12;
q = 1.602192e-19;
Nacc = 1e24;
ni = 1.075e16;
ni_before = 0;
ni_after = 0;
K_B = 1.380662e-23;
T = 300;
del_x = Width/(N-1);
coef = del_x*del_x*q/eps_0;

phi(:,1) = 0.33374;
Jaco(1,1) = 1;
Jaco(N, N) = 1;


%% Newton method %%
for voltage = 1:50
    gate_voltage = 1/49*(voltage-1);
    phi(1,1) = 0.33374 + gate_voltage;
    phi(N,1) = 0.33374 + gate_voltage;
    res(1,1) = phi(1,1) - 0.33374 - gate_voltage;
    res(N,1) = phi(N,1) - 0.33374 - gate_voltage;
    
for Newton = 1:100

    for i = 2 : N-1

            Position_before = del_x*((i - 1) - 0.5);
            Position_after = del_x*((i - 1) + 0.5);

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
                Jaco(i, j) = -(r_permitivity(Position_before, Layer_1, Layer_2, Layer_3) + r_permitivity(Position_after, Layer_1, Layer_2, Layer_3)) - coef*ni_eff*q/(K_B*T)*exp(q*phi(i,1)/(K_B*T));

            elseif i == j + 1 
                Jaco(i, j) = r_permitivity(Position_before, Layer_1, Layer_2, Layer_3);

            elseif i == j - 1  
                Jaco(i, j) = r_permitivity(Position_after, Layer_1, Layer_2, Layer_3);

            end
        end

            res(i,1) = r_permitivity(Position_after, Layer_1, Layer_2, Layer_3)*(phi(i+1, 1) - phi(i, 1)) - r_permitivity(Position_before, Layer_1, Layer_2, Layer_3)*(phi(i,1) - phi(i-1,1)) - coef*(Nacc + ni_eff*exp(q*phi(i,1)/(K_B*T)));  
    end
    
update = Jaco\(-res);
phi = phi + update;
    
end    

result(:,voltage) = phi;

end
%% result %%

electron_density = zeros(voltage,1);

for ii = 1:voltage
    for jj = 1:N
        if (del_x*(jj-1) < Layer_1(1, 1)) || (del_x*(jj-1) > Layer_1(1, 1)+Layer_2(1, 1))
            electron_density(ii,1) = electron_density(ii) + 0;
        else
            electron_density(ii,1) = electron_density(ii) + del_x*ni*exp(q*result(jj,ii)/(K_B*T));
        end
    end
end

