%% clear paramenter %%
clear all
clc
close all

%% parameter %%

h_bar = 6.626176e-34/(2*pi);
m0 = 9.109534e-31;
q = 1.602192e-19;
k_B = 1.380662e-23;
T = 300;
L = [100e-9;100e-9;5e-9];
m = [0.19;0.19;0.91];
n_max = 10;
coef = 2*L(1,1)*L(2,1)/(2*pi)*sqrt(m(1,1)*m(2,1))*m0/(h_bar^2)*(k_B*T);

%% Calculation %%
total_electron = 0;
number_Vg = 51;
Vg = -0.1 + 0.2*transpose([0:number_Vg - 1])/(number_Vg - 1);  % set the range of gate voltage (-0.1 V to 0.1 V)

for n = 1:n_max
    
    E = -Vg*q + h_bar^2/(2*m(3,1)*m0)*(pi*n/L(3,1))^2;  % the energe of system 
    total_electron = total_electron + coef*log(1+exp(-E/(k_B*T))); % the calcuation of electron number
    
end

integrated_electron = total_electron/(L(1,1)*L(2,1)); % normalize electron number