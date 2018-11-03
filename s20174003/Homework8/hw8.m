close all;
clear all;



h = 6.626176e-34; %Planck constant (Js)
h_bar = h / (2*pi); %(Js)
Charge_q = 1.602192e-19; %charge Q
m0 = 9.109534e-31; %Electron reset mass (Kg)
Boltz = 1.380662e-23; %Boltzmann constant (J/K)
E_0 = 8.8541878176e-12;%Permittivity(F/m)
E_r1 = 11.7;%Relative Permittivity of Si
E_r2 = 3.9; %Relative Permitivity of SiO2
E_Si = E_0 * E_r1;%Si
E_SiO2 = E_0 * E_r2;%SiO2
Temp = 300.0; %Temperature(K)
thermal = Boltz*Temp/Charge_q; %Thermal voltage(V)

 Lx = 100e-9; %Length of x-direction (m)
 Ly = 100e-9;  %Length of y-direction (m)
% Lz = 5e-9;   %Length of z-direction (m)
mxx = 0.19;  %Mass
myy = 0.19;  %Mass
mzz = 0.91;  %Mass
Delta_Z = 0.1e-9; %0.1nm spacing
Number_N = 61; %For total 6nm thick(0~6nm)
Number_SiO2_Si = 6; %Location 0.5nm
Number_Si_SiO2 = 56; %Location 5.5nm
Z = Delta_Z*[0:Number_N-1]';
N_acc = 1e24; %Dopant density of Acc(1e18/cm^3)
Ni = 1.075e16; %Instrinsic carrier of Si (1/cm^3)
coef_Poi = Delta_Z*Delta_Z*Charge_q*E_0; %coefficient for Poi
coef_Sch = 2*Lx*Ly/(2*pi)*sqrt(mxx*myy)*m0/(h_bar^2)*(Boltz*Temp);%coefficient for Schi
Ec_Ei = 0.561004; %%E_c - E_i, (eV)
%External_Voltage = 0;

Phi_poi_sweep =cell(11,1);
Phi_sch_sweep =cell(11,1);
Elec_poi_sweep = cell(11,1);
Elec_sch_sweep = cell(11,1);

Volt_iter = 1;
elec_poi = zeros(Number_N,1);
Poi_Integrated_elec_density = zeros(11,1);
Sch_Integrated_elec_density = zeros(11,1);

for External_Voltage= 0 : 0.1 :1
    Phi = zeros(Number_N,1);
    Phi(:,1) = 0.33374 + External_Voltage;
    for newton=1:1:100
        
    Residue = zeros(Number_N,1);
    Jaco = sparse(Number_N,Number_N);
    Res(1,1) = Phi(1,1) - (0.33374 + External_Voltage);
    Jaco(1,1) = 1.0;
    Res(Number_N,1) = Phi(Number_N,1) - (0.33374 + External_Voltage);
    Jaco(Number_N,Number_N) = 1.0;

    
        for iter=2:1:Number_N-1%%%%Laplacian Part
           if(iter< Number_SiO2_Si || iter > Number_Si_SiO2)
               Res(iter,1) = E_SiO2*Phi(iter+1,1) - 2*E_SiO2*Phi(iter,1) + E_SiO2*Phi(iter-1,1);
               Jaco(iter,iter-1) = E_SiO2;
               Jaco(iter,iter) = -2*E_SiO2;
               Jaco(iter,iter+1) = E_SiO2;
           elseif(iter==Number_SiO2_Si)
               Res(iter,1) = E_Si*Phi(iter+1,1) - (E_Si+E_SiO2)*Phi(iter,1) + E_SiO2*Phi(iter-1,1);
               Jaco(iter,iter-1) = E_SiO2;
               Jaco(iter,iter) = -(E_Si+E_SiO2);
               Jaco(iter,iter+1) = E_Si;
           elseif(iter==Number_Si_SiO2)
               Res(iter,1) = E_SiO2*Phi(iter+1,1) - (E_Si+E_SiO2)*Phi(iter,1) + E_Si*Phi(iter-1,1);
               Jaco(iter,iter-1) = E_Si;
               Jaco(iter,iter) = -(E_Si+E_SiO2);
               Jaco(iter,iter+1) = E_SiO2;
           else
               Res(iter,1) = E_Si*Phi(iter+1,1) - 2*E_Si*Phi(iter,1) + E_Si*Phi(iter-1,1);
               Jaco(iter,iter-1) = E_Si;
               Jaco(iter,iter) = -2*E_Si;
               Jaco(iter,iter+1) = E_Si;
           end
        end


        for iter=Number_SiO2_Si:Number_Si_SiO2
            if(iter==Number_SiO2_Si)
                Res(iter,1) = Res(iter,1) - coef_Poi*(N_acc + Ni*exp(Phi(iter,1)/thermal))*0.5;
                Jaco(iter,iter) = Jaco(iter,iter) - coef_Poi*Ni*exp(Phi(iter,1)/thermal)/thermal*0.5;
            elseif(iter==Number_Si_SiO2)
                Res(iter,1) = Res(iter,1) - coef_Poi*(N_acc + Ni*exp(Phi(iter,1)/thermal))*0.5;
                Jaco(iter,iter) = Jaco(iter,iter) - coef_Poi*Ni*exp(Phi(iter,1)/thermal)/thermal*0.5;
            else
                Res(iter,1) = Res(iter,1) - coef_Poi*(N_acc + Ni*exp(Phi(iter,1)/thermal));
                Jaco(iter,iter) = Jaco(iter,iter) - coef_Poi*Ni*exp(Phi(iter,1)/thermal)/thermal;
            end
        end
       
        Update = Jaco \ (-Res);
        Phi = Phi + Update;
        
        for iter=Number_SiO2_Si:Number_Si_SiO2
            elec_poi(iter,1) = Ni*exp(Phi(iter,1)/thermal);
        end
        
    end
    
    Phi_poi_sweep{Volt_iter,1} = Phi;
    Elec_poi_sweep{Volt_iter,1} = elec_poi;
    for iter=1:1:size(elec_poi,1)
            Poi_Integrated_elec_density(Volt_iter,1) = sum(elec_poi(:,1)) * (Delta_Z); 
    end
    
   nSubband = 10;
   totalNumber=0;
   elec = zeros(Number_N);
    
   for iNewton=1:20
      
      Residue = zeros(Number_N,1);
      Jaco = sparse(Number_N,Number_N);
      Res(1,1) = Phi(1,1);
      Jaco(1,1) = 1.0;
      Res(Number_N,1) = Phi(Number_N,1);
      Jaco(Number_N,Number_N) = 1.0;
      
      for iValley=1:3
          mass = ones(3)*0.19;
          mass(iValley) = 0.91;
          coef_Sch = 2*Lx*Ly/(2*pi)*sqrt(mass(1)*mass(2))*m0/(h_bar^2)*(Boltz*Temp);
          
          V = Charge_q *Ec_Ei - Charge_q*Phi; %Potential energy (J)
          Nbulk = Number_Si_SiO2 - Number_SiO2_Si -1;
          Hamil = zeros(Nbulk,Nbulk);
          Hamil(1,1) = -2;
          Hamil(1,2) = 1;

          for iter=2:Nbulk-1
             Hamil(iter,iter+1) = 1;
             Hamil(iter,iter) = -2;
             Hamil(iter,iter-1) =1;
          end
          Hamil(Nbulk,Nbulk) = -2;
          Hamil(Nbulk,Nbulk-1) = 1;

          for iter=1:Nbulk
              Hamil(iter,iter) = Hamil(iter,iter) -2*mass(3)*m0*(Delta_Z/h_bar)^2*V(iter+Number_SiO2_Si,1);
          end
          
          [eigenvectors,eigenvalues] = eig(Hamil);
          scaledEz = diag(eigenvalues)/(-2*mass(3)*m0*(Delta_Z/h_bar)^2); %Eigen_energy(J)
          [sortedEz,sortedIndex] = sort(scaledEz);
          
          for iter=1:nSubband
            Ez = sortedEz(iter,1);
            wavefunction2 = eigenvectors(:,sortedIndex(iter)).^2;
            normalization = sum(wavefunction2)*Delta_Z;
            wavefunction2 = wavefunction2 / normalization;
            subbandNumber = coef_Sch*log(1+exp(-Ez/(Boltz*Temp)));
            totalNumber = totalNumber + subbandNumber;
            elec(Number_SiO2_Si+1:Number_Si_SiO2-1,1) = elec(Number_SiO2_Si+1:Number_Si_SiO2-1,1) + 1/(Lx*Ly)*wavefunction2*subbandNumber;
          end
      end
      
      for iter=2:1:Number_N-1%%%%Laplacian Part
           if(iter< Number_SiO2_Si || iter > Number_Si_SiO2)
               Res(iter,1) = E_SiO2*Phi(iter+1,1) - 2*E_SiO2*Phi(iter,1) + E_SiO2*Phi(iter-1,1);
               Jaco(iter,iter-1) = E_SiO2;
               Jaco(iter,iter) = -2*E_SiO2;
               Jaco(iter,iter+1) = E_SiO2;
           elseif(iter==Number_SiO2_Si)
               Res(iter,1) = E_Si*Phi(iter+1,1) - (E_Si+E_SiO2)*Phi(iter,1) + E_SiO2*Phi(iter-1,1);
               Jaco(iter,iter-1) = E_SiO2;
               Jaco(iter,iter) = - (E_Si+E_SiO2);
               Jaco(iter,iter+1) = E_Si;
           elseif(iter==Number_Si_SiO2)
               Res(iter,1) = E_SiO2*Phi(iter+1,1) - (E_Si+E_SiO2)*Phi(iter,1) + E_Si*Phi(iter-1,1);
               Jaco(iter,iter-1) = E_Si;
               Jaco(iter,iter) = -(E_Si+E_SiO2);
               Jaco(iter,iter+1) = E_SiO2;
           else
               Res(iter,1) = E_Si*Phi(iter+1,1) - 2*E_Si*Phi(iter,1) + E_Si*Phi(iter-1,1);
               Jaco(iter,iter-1) = E_Si;
               Jaco(iter,iter) = -2*E_Si;
               Jaco(iter,iter+1) = E_Si;
           end
      end
      
      for iter=Number_SiO2_Si:Number_Si_SiO2
            if(iter==Number_SiO2_Si)
                Res(iter,1) = Res(iter,1) - coef_Poi*(elec(iter,1))*0.5;
                Jaco(iter,iter) = Jaco(iter,iter);
            elseif(iter==Number_Si_SiO2)
                Res(iter,1) = Res(iter,1) - coef_Poi*(elec(iter,1))*0.5;
                Jaco(iter,iter) = Jaco(iter,iter);
            else
                Res(iter,1) = Res(iter,1) - coef_Poi*(elec(iter,1));
                Jaco(iter,iter) = Jaco(iter,iter);
            end
      end
      
      Update = Jaco \ (-Res);
      Phi = Phi + Update;    
   end

   Phi_sch_sweep{Volt_iter,1} = Phi;
   Elec_sch_sweep{Volt_iter,1} = elec;
   for iter=1:1:size(elec,1)
           Sch_Integrated_elec_density(Volt_iter,1) = sum(elec(:,1)) * (Delta_Z); 
   end
   
   Volt_iter = Volt_iter + 1;

   
end
 



Voltage_range = 0:0.1:1;%0~1V 11point step.
figure(1)
semilogy(Voltage_range,Poi_Integrated_elec_density/1e6,'b');
hold on;
semilogy(Voltage_range,Sch_Integrated_elec_density/1e6,'r');
xlabel('Gate Voltage[V]');
ylabel('Integrated Electron density [cm^-2]');
grid on;




figure(2)

semilogy(Z/1e-9,Elec_poi_sweep{1,1}/1e6,'k','LineWidth',2); hold on;
semilogy(Z/1e-9,Elec_poi_sweep{2,1}/1e6,'b','LineWidth',2); hold on;
semilogy(Z/1e-9,Elec_poi_sweep{3,1}/1e6,'g','LineWidth',2); hold on;
semilogy(Z/1e-9,Elec_poi_sweep{4,1}/1e6,'r','LineWidth',2); hold on;
semilogy(Z/1e-9,Elec_poi_sweep{5,1}/1e6,'c','LineWidth',2); hold on;
semilogy(Z/1e-9,Elec_poi_sweep{6,1}/1e6,'m','LineWidth',2); hold on;
semilogy(Z/1e-9,Elec_poi_sweep{7,1}/1e6,'y','LineWidth',2); hold on;
semilogy(Z/1e-9,Elec_poi_sweep{8,1}/1e6,'--','LineWidth',2); hold on;
semilogy(Z/1e-9,Elec_poi_sweep{9,1}/1e6,':','LineWidth',2); hold on;
semilogy(Z/1e-9,Elec_poi_sweep{10,1}/1e6,'-.','LineWidth',2); hold on;
semilogy(Z/1e-9,Elec_poi_sweep{11,1}/1e6,'x','LineWidth',2); 
xlabel('Distance [nm]');
ylabel('Electron density [cm^-3]');
grid on;
legend('0V','0.1V','0.2V','0.3V','0.4V','0.5V','0.6V','0.7V','0.8V','0.9V','1V','Location','best');



figure(3)
semilogy(Z/1e-9,Elec_sch_sweep{1,1}/1e6,'k','LineWidth',2); hold on;
semilogy(Z/1e-9,Elec_sch_sweep{2,1}/1e6,'b','LineWidth',2); hold on;
semilogy(Z/1e-9,Elec_sch_sweep{3,1}/1e6,'g','LineWidth',2); hold on;
semilogy(Z/1e-9,Elec_sch_sweep{4,1}/1e6,'r','LineWidth',2); hold on;
semilogy(Z/1e-9,Elec_sch_sweep{5,1}/1e6,'c','LineWidth',2); hold on;
semilogy(Z/1e-9,Elec_sch_sweep{6,1}/1e6,'m','LineWidth',2); hold on;
semilogy(Z/1e-9,Elec_sch_sweep{7,1}/1e6,'y','LineWidth',2); hold on;
semilogy(Z/1e-9,Elec_sch_sweep{8,1}/1e6,'--','LineWidth',2); hold on;
semilogy(Z/1e-9,Elec_sch_sweep{9,1}/1e6,':','LineWidth',2); hold on;
semilogy(Z/1e-9,Elec_sch_sweep{10,1}/1e6,'-.','LineWidth',2); hold on;
semilogy(Z/1e-9,Elec_sch_sweep{11,1}/1e6,'x','LineWidth',2); 
xlabel('Distance [nm]');
ylabel('Electron density[cm^-3]');
grid on;
legend('0V','0.1V','0.2V','0.3V','0.4V','0.5V','0.6V','0.7V','0.8V','0.9V','1V');
