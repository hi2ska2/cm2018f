%% System parameters
meff = 0.19; % effective of mass for the electrons in a periodic potential.
L = 5e-9; % square well width.
Npoints = 500; % number of discretized points.








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Elementary constants (SI units)
hbar = 6.626e-34;
m0 = 9.11e-31;

%% Analytical solutions
syms N X
FAES = sqrt(2/L)*sin((N*pi/ L)*X); % Function of analytical eigenstates.
FAEV = (N*pi*hbar)^2/(2*(meff*m0)*L^2); % Function of nalytical eigenvalues.

%% Numerical solutions
xtot = linspace(0, L, Npoints);
x = xtot(2 : end-1); % x points without boundaries.

% Simplified Hamiltonian
SH = zeros(length(x), length(x));
for i = 1 : length(x)
    for j = 1 : length(x)
        if i == j
            SH(i, j) = -2;
        elseif i == j+1
            SH(i, j) = 1;
        elseif j == i+1
            SH(i, j) = 1;
        end
    end
end

% Numerical eigenstates and values of the hamiltonian.
[ES, EV] = eig(-SH); % the minus makes positive results, just trick.
NES = sqrt((length(xtot)-1)/L)*ES; % normalizing coefficient and boundary.-
NEV = diag((hbar^2/(2*meff*m0))*(EV/(x(2)-x(1))^2))';

%% Comparison between analytical and numerical solutions
% generated quantum numbers
n = 1 : length(x);

% generate analytical eigenstates and values.
AES = zeros(length(x), length(n));
for i = 1 : length(x)
    tic
    fprintf('Iteration of sampling from analytical solution : (%d/%d)\n', i, length(x))
    AES(i, :) = double(subs(subs(FAES, N, n), X, x(i)));
    toc
end
AEV = double(subs(FAEV, N, n));

% electron density.
BDR = zeros(1, length(x));
DAES = [BDR; AES.*AES; BDR];
DNES = [BDR; NES.*NES; BDR];

% estimation of the error between eigenvalues.
RatioEVerr = sqrt((NEV - AEV).*(NEV - AEV))./sqrt(AEV.*AEV);
RatioGRerr = RatioEVerr(1); % Ratio error of ground state.

%% Visualizing first three states only
guideposition = linspace(0, L, 1000);
figure(1)
subplot(3, 1, 1)
guideline1 = double(subs(subs(FAES, N, 1), X, guideposition));
plot(xtot*1e+9, DAES(:, 1), '*r', xtot*1e+9, DNES(:, 1), 'ob-', guideposition*1e+9, guideline1.^2, '-r')
title(sprintf('Visualizing first three states, Npoints = %d', Npoints))
legend(sprintf('Analytic, Energy : %2.2f [eV]', AEV(1)/1.6e-19), sprintf('Numeric, Energy : %2.2f [eV]', NEV(1)/1.6e-19))

subplot(3, 1, 2)
guideline2 = double(subs(subs(FAES, N, 2), X, guideposition));
plot(xtot*1e+9, DAES(:, 2), '*r', xtot*1e+9, DNES(:, 2), 'ob-', guideposition*1e+9, guideline2.^2, '-r')
legend(sprintf('Analytic, Energy : %2.2f [eV]', AEV(2)/1.6e-19), sprintf('Numeric, Energy : %2.2f [eV]', NEV(2)/1.6e-19))
ylabel('Squared eigenstates')

subplot(3, 1, 3)
guideline3 = double(subs(subs(FAES, N, 3), X, guideposition));
plot(xtot*1e+9, DAES(:, 3), '*r', xtot*1e+9, DNES(:, 3), 'ob-', guideposition*1e+9, guideline3.^2, '-r')
legend(sprintf('Analytic, Energy : %2.2f [eV]', AEV(3)/1.6e-19), sprintf('Numeric, Energy : %2.2f [eV]', NEV(3)/1.6e-19))
xlabel('Position [nm]')