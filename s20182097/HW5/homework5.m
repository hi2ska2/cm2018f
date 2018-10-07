%% clear paramenter %%
clear all
clc
close all

%% parameter %%
rep = 1000;  % The number of Newton-Raphson method repetition 
num_Ncc = 1000; % The number of Ncc step 

ini_phi = 10; % initial value
n_int = 1e16; % m^(-3)
k_B = 1.380662e-23; % J/K
T = 300; % K
q = 1.602192e-19;
V_T = k_B*T/q;

for i = 1 : num_Ncc
        Ncc = power(10, 16 + 8*(i - 1)/(num_Ncc - 1)); %From 10e16 m^(-3) to 10e24 m^(-3)
      
        phi = ini_phi; %The initialization of phi
        
    for j = 1 : rep

        Equation = Ncc + n_int * exp(-phi/V_T) - n_int * exp(phi/V_T);
        Taylor_first_order = -n_int/V_T * exp(-phi/V_T) - n_int/V_T * exp(phi/V_T) ;

        del_phi = -Equation/Taylor_first_order;
        phi = phi + del_phi;

    end
    
    solution(i, 1) = Ncc;
    solution(i, 2) = phi;
    solution(i, 3) = asinh(Ncc/(2*n_int))*V_T;
end
 