%% Define the geometry
Thickness = [20, 30, 20, 30, 20, 30, 20, 30, 5, 10]*1e-9; % [m]
reps = [127, 300, 127, 300, 127, 300, 127, 300, 3.9, 11.7]; % relative permittivity
Voltage = [0, 100]; % [Voltage of left, and right]
Lstep = [2, 3, 2, 3, 2, 3, 2, 3, 0.5 0.5]*1e-9; % [m]






DoAnalytical = 0; % [5nm, 5nm], [11.7, 3.9], [0, 1], [0.1nm 0.1nm]에 대한 analytical solution 생성.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eps = reps * 8.85e-12;


%% Matrix of the poission equation
LastX = zeros(length(Thickness), 1);
for k = 1 : length(Thickness)
    X{k, 1} = 0+LastX(k, 1) : Lstep(k) : Thickness(k)+LastX(k, 1);
    LastX(k+1, 1) = LastX(1, 1) + max(X{k, 1}(1, :));
end
x = unique(horzcat(X{:, :}))';
dx = diff(x);

for k = 1 : length(x)
    Junction(k, :) = logical(x(k) == LastX(2:end, 1));
end


Matrix = zeros(length(x), length(x));
Matrix(1, 1) = 1;
Matrix(end, end) = 1;
Matrix(2, 1) = eps(1)/dx(1);
Matrix(end-1, end) = eps(end)/dx(end);
iindex = 1 : length(x);
for k = 1 : length(Junction(1, :)) - 1
    Junctioniindex(k, 1) = iindex(Junction(:, k));
end
sectioni = vertcat(0, Junctioniindex, iindex(end));


for i = 2 : length(x) - 1
    for j = 2 : length(x) - 1
        for k = 1 : length(Junction(1, :)) 
            if Junction(i, k) == 1 && Junction(j, k) == 1 % Junction에 해당하는 i, j에서
                Matrix(i, j) = -eps(k+1)/dx(i) -eps(k)/dx(i-1);
            elseif Junction(i, k) == 1 && Junction(j+1, k) == 1 % Junction 주변
                Matrix(i, j) = eps(k)/dx(i-1);
            elseif Junction(i, k) == 1 && Junction(j-1, k) == 1 % Junction 주변
                Matrix(i, j) = eps(k+1)/dx(i);
            elseif i > sectioni(k, 1) && i < sectioni(k+1, 1) % Junction이 아닌곳
                if i == j
                    Matrix(i, j) = -2*eps(k)/dx(i-1);
                elseif j == i - 1
                    Matrix(i, j) = eps(k)/dx(i-1);
                elseif j == i + 1
                    Matrix(i, j) = eps(k)/dx(i-1);
                end
            end
        end
    end
end


%% Take the solution of potential
Boundary = zeros(length(x), 1);
Boundary(1, 1) = Voltage(1);
Boundary(end, 1) = Voltage(2);

Nphi = Matrix \ Boundary;

%% Do analytical solution and compare with numerical solution
if DoAnalytical == 1
    syms a y
    Aphi1 = y/(2*a);
    Aphi2 = (3*y)/(2*a) - (1/2);
    
    Aphi1N = double(subs(subs(Aphi1, a, x(end)), y, x(1:Junctioniindex, 1)));
    Aphi2N = double(subs(subs(Aphi2, a, x(end)), y, x(Junctioniindex+1:end, 1)));
    AphiN = vertcat(Aphi1N, Aphi2N);
    
    plot(x*1e+9, Nphi, 'ob', x*1e+9, AphiN, '-r')
    title('Comparison with numerical and analytical solutions')
    xlabel('Position [nm]')
    ylabel('Voltage [V]')
    legend('Numerical solution', 'Analytical solution')
else
    plot(x*1e+9, Nphi, 'ob')
    title('Numerical solution of a heterostructure')
    xlabel('Position [nm]')
    ylabel('Voltage [V]')
    legend('Numerical solution')
end