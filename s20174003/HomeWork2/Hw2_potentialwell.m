close all;
clear all;

N = 5; %%N = 5,50,500
%N = 50;
%N = 500;
M_0 = 9.10938356e-31; % Electron rest mass (Kg)
Planck_constant = 1.054571800e-34; %Planck constant divided 2pi(Js)

M = 0.19 * M_0; %%Effective electron mass
a = 5e-9; %%5nm 

Analytic_K_square = (pi^2) / (a^2);
Analytic_Energy = ((Planck_constant^2) * Analytic_K_square) / (2*M); %%[J]
Analytic_Energy_electron = Analytic_Energy / (1.6e-19);%%[eV]

%%%%Initialize for N=5,50,500%%%%
Energy = zeros(3,1);
Energy_electron = zeros(3,1);
%%%%%%%%%%%%%%%%%%

for(iter=1:1:3)
    if(iter==1)
        N=5;
    elseif(iter==2)
        N=50;
    else
        N=500;
    end
    delta_x = a/(N-1);
    
    Equation_A = zeros(N-2);
    Equation_A = Equation_A + (-2*eye(N-2));
    Equation_A = Equation_A + diag(ones(1,N-3),1);
    Equation_A = Equation_A + diag(ones(1,N-3),-1);

    [Eigen_V, Eigen_D] = eig(Equation_A);
    Find_Min_K_num = find(Eigen_D);
    Find_Min_K = (-1) * max(Eigen_D(Find_Min_K_num));

    K_square = Find_Min_K / (delta_x^2);

    Energy(iter,1) = ((Planck_constant^2) * K_square) / (2*M);

    Energy_electron(iter,1) = Energy(iter,1) / (1.6e-19)
end


figure
stem(Energy_electron(:,1));hold on;
hy = graph2d.constantline(Analytic_Energy_electron, 'LineStyle',':','Color',[0 .0 .0]);
title('SW HW Potential Well Problem in N=5,50,and 500 [eV]');
legend('N=5,50,500','Analytic solution','Location','best');
ylim([0.075 0.0795])





