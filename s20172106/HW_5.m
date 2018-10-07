clear all;
clc

q_ele = 1; % elementary charge, [electron coulomb]
k_B= 8.6173303e-5; % Boltzmann constant, [eV/K]
T=300.0; % temperature, % Temp [K]
Volt_Thermal = (k_B*T)/q_ele; % V_Temp [V]
n_int=1.0e10; % Intrinsic carrier density for Si on T=300K, [/cm^3]

Newton_method_Number = 1000;

N = 1000;

Phi_Init = 10; % Initial electrostatic Potential [V]
Final_data = zeros(N,7);

N_acc = logspace(10,18,N);
Final_data(:,1) = N_acc';
% 
%% Numerical solution using Newoton method

for ii = 1:N;
    phi = Phi_Init; % The estimated potential as the statrting point [V]
    ftn_plus =  N_acc(ii) + n_int*(exp(-(phi/Volt_Thermal))) - n_int*(exp((phi/Volt_Thermal)));
    Jaco_ftn_plus = (n_int/Volt_Thermal)*(-exp(-(phi/Volt_Thermal)) - exp((phi/Volt_Thermal)));
    
    Delta_phi =  Jaco_ftn_plus\(-ftn_plus);
    for jj = 1:Newton_method_Number;
        phi = phi + Delta_phi;
        ftn_plus =  N_acc(ii) + n_int*(exp(-(phi/Volt_Thermal))) - n_int*(exp((phi/Volt_Thermal)));
        Jaco_ftn_plus = (n_int/Volt_Thermal)*(-exp(-(phi/Volt_Thermal)) - exp((phi/Volt_Thermal)));
        Delta_phi =  Jaco_ftn_plus\(-ftn_plus);
    end
    
    Final_data(ii,2) = phi;
    Final_data(ii,4) = Delta_phi;
    
    phi = -Phi_Init; % The estimated potential as the statrting point [V]
    ftn_minus =  - N_acc(ii) + n_int*(exp(-(phi/Volt_Thermal))) - n_int*(exp((phi/Volt_Thermal)));
    Jaco_ftn_minus = (n_int/Volt_Thermal)*(-exp(-(phi/Volt_Thermal)) - exp((phi/Volt_Thermal)));
    
    Delta_phi =  Jaco_ftn_minus\(-ftn_minus);
    for kk = 1:Newton_method_Number;
        phi = phi + Delta_phi;
        ftn_minus =  -N_acc(ii) + n_int*(exp(-(phi/Volt_Thermal))) - n_int*(exp((phi/Volt_Thermal)));
        Jaco_ftn_minus = (n_int/Volt_Thermal)*(-exp(-(phi/Volt_Thermal)) - exp((phi/Volt_Thermal)));
        Delta_phi =  Jaco_ftn_minus\(-ftn_minus);
    end
    
    Final_data(ii,3) = phi;
    Final_data(ii,5) = Delta_phi;
end

%% Analytic solution using Newoton method

Anal_phi = zeros(N,3); %
Anal_phi(:,1) = N_acc';

for kk = 1:N;
    Anal_phi(kk,2) = Volt_Thermal*asinh(Anal_phi(kk,1)/(2*n_int));
    Anal_phi(kk,3) = Volt_Thermal*asinh(-Anal_phi(kk,1)/(2*n_int));
end

Final_data(:,4) = Anal_phi(:,2);
Final_data(:,5) = Anal_phi(:,3);

%% Plot

figure(1);
semilogx(Final_data(:,1),Final_data(:,2),'r',Final_data(:,1),Final_data(:,3),'b'); hold on;
semilogx(Final_data(:,1),Final_data(:,4),'o','MarkerEdgeColor',[1 0 0]); hold on;
semilogx(Final_data(:,1),Final_data(:,5),'o','MarkerEdgeColor',[0 0 1]); hold on;
legend('Numerical N^+','Numerical N^-','Analytic N^+','Analytic N^-') 
xlabel('N^+ or N^- (cm^-^3)');
ylabel('Potential (V)');

figure(2)
subplot(2,1,1);
plot(Final_data(:,1),Final_data(:,2)-Final_data(:,4),'o','markersize',6,'markerEdgecolor','r');
xlabel('N^+ (cm^-^3)');
ylabel('Error for potential (V)');

subplot(2,1,2);
plot(Final_data(:,1),Final_data(:,3)-Final_data(:,5),'o','markersize',6,'markerEdgecolor','b');
xlabel('N^- (cm^-^3)');
ylabel('Error for potential (V)');













