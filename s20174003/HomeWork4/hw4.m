close all;
clear all;

Voltage_sweep = 1;
Voltage_range = 0:0.1:1;%0~1V 11point step.
Charge_q = 1.602192e-19;
Boltz = 1.380662e-23; %Boltzmann constant (J/K)
Temp = 300.0; %Temperature(K)
N_acc = 1e24; %Dopant density of Acc(1e18/cm^3)
Ni = 1.075e16; %Instrinsic carrier of Si (1/cm^3)



E_0 = 8.8541878176e-12;%Permittivity(F/m)
E_r1 = 11.7;%Relative Permittivity of Si
E_r2 = 3.9; %Relative Permitivity of SiO2

E_1 = E_0 * E_r1;
E_2 = E_0 * E_r2;


Delta_X = 0.1e-9; %0.1nm
Number_N = 61; %For total 6nm thick(0~6nm)
Number_SiO2_Si = 6; %Location 0.5nm
Number_Si_SiO2 = 56; %Location 5.5nm
X = Delta_X*[0:Number_N-1]';



%%Make the Equation AX=B%%%
Equation_A = zeros(Number_N,Number_N);
Equation_A(1,1) = 1.0;
Equation_A(Number_N,Number_N) = 1.0;
for (iter =2:1:Number_N-1)
   if(iter < Number_SiO2_Si)
       Equation_A(iter,iter-1) = E_r2;
       Equation_A(iter,iter) = -2*E_r2;
       Equation_A(iter,iter+1) = E_r2;
   elseif(iter ==Number_SiO2_Si)
       Equation_A(iter,iter-1) = E_r2;
       Equation_A(iter,iter) = -E_r2 - E_r1;
       Equation_A(iter,iter+1) = E_r1;
   elseif(iter < Number_Si_SiO2)
       Equation_A(iter,iter-1) = E_r1;
       Equation_A(iter,iter) = -2*E_r1;
       Equation_A(iter,iter+1) = E_r1;
   elseif(iter==Number_Si_SiO2)
       Equation_A(iter,iter-1) = E_r1;
       Equation_A(iter,iter) = -E_r1 - E_r2;
       Equation_A(iter,iter+1) = E_r2;
   elseif(iter > Number_Si_SiO2)
       Equation_A(iter,iter-1) = E_r2;
       Equation_A(iter,iter) = -2*E_r2;
       Equation_A(iter,iter+1) = E_r2;
   end
end

Equation_B = zeros(Number_N,1);
if(Voltage_sweep==0)
    Equation_B(1,1) = 0.33374;
    Equation_B(Number_N,1) = 0.33374;
else
    Equation_B(1,1) = 0;
    Equation_B(Number_N,1) = 0;
end


for(iter=Number_SiO2_Si : 1 : Number_Si_SiO2)
   if(iter==Number_SiO2_Si)
       Equation_B(iter,1) = (Delta_X) * (Delta_X) * (Charge_q) * (N_acc) * (0.5) / (E_0);
   elseif(iter==Number_Si_SiO2)
       Equation_B(iter,1) = (Delta_X) * (Delta_X) * (Charge_q) * (N_acc) * (0.5) / (E_0);
   else
       Equation_B(iter,1) = (Delta_X) * (Delta_X) * (Charge_q) * (N_acc) * (1) / (E_0);
   end
end

if(Voltage_sweep==0)
    Potential_phi = Equation_A \ Equation_B;
else
    Potential_phi=cell(11,1);
    for(iter=1:1:11)
        Equation_B(1,1) = Voltage_range(1,iter);
        Equation_B(Number_N,1) = Voltage_range(1,iter);
        Potential_phi{iter,1} = Equation_A \ Equation_B;
    end
end

if(Voltage_sweep==0)
    Electron_density = zeros(Number_N,1);

    for(iter=Number_SiO2_Si:1:Number_Si_SiO2)
        Electron_density(iter,1) = Ni * exp(Charge_q*Potential_phi(iter,1)/(Boltz*Temp));
    end
else
    Electron_density = cell(11,1);
    for(iter=1:1:11)%%Initialize
        Electron_density{iter,1} = zeros(Number_N,1);
    end
    for(iteriter=1:1:11)
        for(iter=Number_SiO2_Si:1:Number_Si_SiO2)
            Electron_density{iteriter,1}(iter,1) = Ni * exp(Charge_q*Potential_phi{iteriter,1}(iter,1)/(Boltz*Temp));
        end
    end
end


if(Voltage_sweep==0)
    figure(1)
    plot(X/1e-9,Potential_phi,'LineWidth',2);
    title('SW HW4 Potential');
    xlabel('Distance [nm]');
    ylabel('Potential [V]');


    figure(2)
    plot(X/1e-9,Electron_density*1e-6,'LineWidth',2);
    title('SW HW4 Electron Density');
    xlabel('Position [nm]');
    ylabel('Electron density [cm^-^3]');
else
    figure(1)
    for(iter=1:1:11)
        plot(X/1e-9,Potential_phi{iter,1},'LineWidth',2);hold on;
    end
    title('SW HW4 Potential');
    xlabel('Distance [nm]');
    ylabel('Potential [V]');
    legend('0V','0.1V','0.2V','0.3V','0.4V','0.5V','0.6V','0.7V','0.8V','0.9V','1V','Location','best');

    figure(2)
    for(iter=1:1:11)
        plot(X/1e-9,Electron_density{iter,1}*1e-6,'LineWidth',2);hold on;
    end
    title('SW HW4 Electron Density');
    xlabel('Position [nm]');
    ylabel('Electron density [cm^-^3]');
    legend('0V','0.1V','0.2V','0.3V','0.4V','0.5V','0.6V','0.7V','0.8V','0.9V','1V','Location','best');
end
