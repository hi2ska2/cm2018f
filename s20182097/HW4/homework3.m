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
b = zeros(N,1);
Electron_density = zeros(N,1);

gate_voltage = 1;
Width = Layer_1(1, 1) + Layer_2(1, 1) + Layer_3(1, 1);
eps_0 = 8.8541878176e-12;
q = 1.602192e-19;
Nacc = 1e24;
ni = 1.075e16;
K_B = 1.380662e-23;
T = 300;


%% define matrix %%
for i = 1 : N
    for j = 1 : N
        Position_before = Width/(N-1)*((i - 1) - 0.5);
        Position_after = Width/(N-1)*((i - 1) + 0.5);
        
        if (i == 1) || (i == N)
            A(i, j) = 0;
            
        elseif i == j
            A(i, j) = -(r_permitivity(Position_before, Layer_1, Layer_2, Layer_3) + r_permitivity(Position_after, Layer_1, Layer_2, Layer_3));
            
        elseif i == j + 1 
            A(i, j) = r_permitivity(Position_before, Layer_1, Layer_2, Layer_3);
            
        elseif i == j - 1  
            A(i, j) = r_permitivity(Position_after, Layer_1, Layer_2, Layer_3);
            
        end
       
    end
end

A(1,1) = 1;
A(N, N) = 1;

for k = 1 : N
   
    Position = Width/(N-1)*(k - 1);
    
    if (Position == Layer_1(1, 1)) || (Position == Layer_1(1, 1) + Layer_2(1, 1))
        b(k, 1) = q*Nacc*power(Width/(N-1), 2)/eps_0*0.5;
    
    elseif Position < Layer_1(1, 1)
        b(k, 1) = 0;
        
    elseif Position < Layer_1(1, 1) + Layer_2(1, 1)
        b(k, 1) = q*Nacc*power(Width/(N-1), 2)/eps_0;
        
    end
end

b(1, 1) = 0.33374 - gate_voltage;
b(N, 1) = 0.33374 - gate_voltage;

%% solution %%
x = Width/(N-1)*transpose([0:N-1]);
Sample_potential = [x A\b];


for m = 1 : N
   
    Position = Width/(N-1)*(m - 1);
     
    if (Position >= Layer_1(1, 1)) && (Position <= Layer_1(1, 1) + Layer_2(1, 1))
        Electron_density(m, 1) = ni*exp(q*Sample_potential(m, 2)/(K_B*T));
        
    end
end

b_new = b + q*Electron_density*power(Width/(N-1), 2)/eps_0;
Sample_potential_new = [x A\b_new];
difference = Sample_potential - Sample_potential_new;

