close all;
clear all;

Voltage_sweep = 1;
Voltage_range = 0:0.1:1;%0~1V 11point step.
External_Voltage = 0.0;
Charge_q = 1.602192e-19; %charge Q
Boltz = 1.380662e-23; %Boltzmann constant (J/K)
Temp = 300.0; %Temperature(K)
thermal = Boltz*Temp/Charge_q;

N_acc = 1e24; %Dopant density of Acc(1e18/cm^3)
Ni = 1.075e16; %Instrinsic carrier of Si (1/cm^3)



E_0 = 8.8541878176e-12;%Permittivity(F/m)
E_r1 = 11.7;%Relative Permittivity of Si
E_r2 = 3.9; %Relative Permitivity of SiO2

E_1 = E_0 * E_r1;%Si
E_2 = E_0 * E_r2;%SiO2


Delta_X = 0.1e-9; %0.1nm
Number_N = 61; %For total 6nm thick(0~6nm)
Number_SiO2_Si = 6; %Location 0.5nm
Number_Si_SiO2 = 56; %Location 5.5nm
X = Delta_X*[0:Number_N-1]';
coef = Delta_X*Delta_X*Charge_q*E_0; %Coefficient for simple code

Phi_sweep =cell(11,1);

Volt_iter = 0;
for External_Voltage= 0 : 0.1 :1
    Phi = zeros(Number_N,1);
    Phi(:,1) = 0.33374 + External_Voltage;
    Residue = zeros(Number_N,1);
    Jaco = sparse(Number_N,Number_N);
    Res(1,1) = Phi(1,1) - (0.33374 + External_Voltage);
    Jaco(1,1) = 1.0;
    Res(Number_N,1) = Phi(Number_N,1) - (0.33374 + External_Voltage);
    Jaco(Number_N,Number_N) = 1.0;

    for newton=1:1:10
        for iter=2:1:Number_N-1%%%%Laplacian Part
           if(iter< Number_SiO2_Si || iter > Number_Si_SiO2)
               Res(iter,1) = E_2*Phi(iter+1,1) - 2*E_2*Phi(iter,1) + E_2*Phi(iter-1,1);
               Jaco(iter,iter-1) = E_2;
               Jaco(iter,iter) = -2*E_2;
               Jaco(iter,iter+1) = E_2;
           elseif(iter==Number_SiO2_Si)
               Res(iter,1) = E_1*Phi(iter+1,1) - (E_1+E_2)*Phi(iter,1) + E_2*Phi(iter-1,1);
               Jaco(iter,iter-1) = E_2;
               Jaco(iter,iter) = -(E_1+E_2);
               Jaco(iter,iter+1) = E_1;
           elseif(iter==Number_Si_SiO2)
               Res(iter,1) = E_2*Phi(iter+1,1) - (E_1+E_2)*Phi(iter,1) + E_1*Phi(iter-1,1);
               Jaco(iter,iter-1) = E_1;
               Jaco(iter,iter) = -(E_1+E_2);
               Jaco(iter,iter+1) = E_2;
           else
               Res(iter,1) = E_1*Phi(iter+1,1) - 2*E_1*Phi(iter,1) + E_1*Phi(iter-1,1);
               Jaco(iter,iter-1) = E_1;
               Jaco(iter,iter) = -2*E_1;
               Jaco(iter,iter+1) = E_1;
           end
        end


        for iter=Number_SiO2_Si:Number_Si_SiO2
            if(iter==Number_SiO2_Si)
                Res(iter,1) = Res(iter,1) - coef*(N_acc + Ni*exp(Phi(iter,1)/thermal))*0.5;
                Jaco(iter,iter) = Jaco(iter,iter) - coef*Ni*exp(Phi(iter,1)/thermal)/thermal*0.5;
            elseif(iter==Number_Si_SiO2)
                Res(iter,1) = Res(iter,1) - coef*(N_acc + Ni*exp(Phi(iter,1)/thermal))*0.5;
                Jaco(iter,iter) = Jaco(iter,iter) - coef*Ni*exp(Phi(iter,1)/thermal)/thermal*0.5;
            else
                Res(iter,1) = Res(iter,1) - coef*(N_acc + Ni*exp(Phi(iter,1)/thermal));
                Jaco(iter,iter) = Jaco(iter,iter) - coef*Ni*exp(Phi(iter,1)/thermal)/thermal;
            end
        end

        Update = Jaco \ (-Res);
        Phi = Phi + Update;
    end
    Volt_iter = Volt_iter + 1;
    Phi_sweep{Volt_iter,1} = Phi;
end
    %plot(X,Phi);
    
    %plot(X,Phi_sweep{1,1});
    

Integrated_Electron_density = zeros(11,1);
for External_Voltage= 1 : 1 :11
    Electron_density = zeros(Number_N,1);
    for iter=Number_SiO2_Si:Number_Si_SiO2
       Electron_density(iter,1) = Ni * exp(Charge_q*Phi_sweep{External_Voltage,1}(iter,1)/(Boltz*Temp));
    end

    Integrated_Electron_density_temp = zeros(Number_N,1);
    for iter=Number_SiO2_Si:Number_Si_SiO2
        Integrated_Electron_density_temp(iter+1,1) = Electron_density(iter,1) + Integrated_Electron_density_temp(iter,1);
    end

    Integrated_Electron_density(External_Voltage,1) = max(Integrated_Electron_density_temp) * (X(Number_Si_SiO2,1)-X(Number_SiO2_Si)) * (1e100);%%1e100 used 1/cm^2

end


figure(1)

plot(X/1e-9,Phi_sweep{1,1},'k','LineWidth',2); hold on;
plot(X/1e-9,Phi_sweep{2,1},'b','LineWidth',2); hold on;
plot(X/1e-9,Phi_sweep{3,1},'g','LineWidth',2); hold on;
plot(X/1e-9,Phi_sweep{4,1},'r','LineWidth',2); hold on;
plot(X/1e-9,Phi_sweep{5,1},'c','LineWidth',2); hold on;
plot(X/1e-9,Phi_sweep{6,1},'m','LineWidth',2); hold on;
plot(X/1e-9,Phi_sweep{7,1},'y','LineWidth',2); hold on;
plot(X/1e-9,Phi_sweep{8,1},'--','LineWidth',2); hold on;
plot(X/1e-9,Phi_sweep{9,1},':','LineWidth',2); hold on;
plot(X/1e-9,Phi_sweep{10,1},'-.','LineWidth',2); hold on;
plot(X/1e-9,Phi_sweep{11,1},'x','LineWidth',2); hold on;
xlabel('Distance [nm]');
ylabel('Potential [V]');
grid on;
legend('0V','0.1V','0.2V','0.3V','0.4V','0.5V','0.6V','0.7V','0.8V','0.9V','1V','Location','best');

figure(2)
plot(Voltage_range ,Integrated_Electron_density,'-o','LineWidth',2);
xlabel('Voltage [V]');
ylabel('Integrated Electron Density [1/cm^2]');
grid on;



figure(3)
plot(X/1e-9,Phi_sweep{1,1},'b','LineWidth',2); 
xlabel('Distance [nm]');
ylabel('Potential [V]');
grid on;

