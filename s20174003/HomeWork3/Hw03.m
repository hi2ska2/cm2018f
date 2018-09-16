close all;
clear all;

N = 5; %%N = 5,50,500

E_0 = 8.8541878176e-12; %Permittivity(F/m)
E_r1 = 11.7; %Relative Permittivity of Si
E_r2 = 3.9; %Relative Permittivity of SiO2
E_1 = E_0 * E_r1;
E_2 = E_0 * E_r2;

T = 5e-9; %%5nm Thickness (Space divided 2.5nm/2.nm)


Equation_A = [1 0 0 0 0; E_1 (-2*E_1) E_1 0 0; 0 E_1 -(E_1+E_2) E_2 0; 0 0 E_2 (-2*E_2) E_2; 0 0 0 0 1];
Equation_B = [0;0;0;0;1];
X = Equation_A \ Equation_B;%Potential value

C_1 = E_1/(T/2);%Capacitance of Si(F/cm^2)
C_2 = E_2/(T-T/2);%Capacitance of Si(F/cm^2)
C_Total = (1/C_1 + 1/C_2)^-1;

sprintf('%0.10f',C_1)%%[F/cm^2]
sprintf('%0.10f',C_2)%%[F/cm^2]
sprintf('%0.10f',C_Total)%%[F/cm^2]


Equation_Analytic_Total = ones(1,5);
iter = 1;
for distance_x_2=0 : 5/4 : 5;
    if(distance_x_2 <= 2.5)
        Equation_Analytic_Total(1,iter) = (X(3)/(T/2))*distance_x_2*(1e-9);
    else
        Equation_Analytic_Total(1,iter) = ((X(5)-X(3))/(T-T/2))*(distance_x_2-2.5)*(1e-9) + X(3);
    end
    iter = iter + 1;
end
%%%%%%%%%%%%%%%%%%%%

distance_x_2 = 0:5/4:5;

 figure(1)
 distance_x = 0:5/4:5;%distance(5/4nm space)
 stem(distance_x,X,'r','LineStyle',':');
 hold on;
 plot(distance_x_2,Equation_Analytic_Total);
 title('SW HW3 Potential');
 xlabel('Distance [nm]');
 ylabel('Potential [V]');
 legend('N=5 Discrete','Analytic','Location','best');
