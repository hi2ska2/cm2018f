%% Defining variables
% Basic constants
q = 1.602e-19; % Elementary charge [C]
eps0 = 8.85e-12; % Vacuum permittivity [F/m]
kB = 1.38e-23; % Boltzmann constant [J/K]
T = 300; % Temperature [K]

ThermalV = kB*T/q; % Thermal voltage [V]


% Geometry
DeltaX = 0.1e-9; % 0.1 nm spacing [m]
N = 61; % 
x = DeltaX * transpose(0:N-1); % real space [m]
interface1 = 6; % At x = 0.5 nm
interface2 = 56; % At x = 5.5 nm

coef = DeltaX * DeltaX * q/eps0;


% System
eps_si = 11.7; eps_ox = 3.9; % Relative permittivity
N_acc = 1e24; % 1e18 / cm^3
ni = 1.075e16; % 1.075e10 / cm^3
GateVoltages = linspace(0, 0.2, 101);

for jj = 1 : length(GateVoltages)
    %% Nonlinear Poission Equation
    for newton = 1 : 10
        % Initialization
        Residue = zeros(N, 1);
        Jacobian = sparse(N, N);

        phi = zeros(N, 1);
        phi(:, 1) = 0.33374 + GateVoltages(jj);
        

        % Boundary conditions
        Residue(1, 1) = phi(1, 1) - 0.33374;
        Jacobian(1, 1) = 1;
        Residue(N, 1) = phi(N, 1) - 0.33374;
        Jacobian(N, N) = 1;

        % Laplacian part
        for ii = 2: N-1
            if ii < interface1 || ii > interface2
                Residue(ii, 1) = eps_ox * phi(ii+1, 1) - 2*eps_ox * phi(ii, 1) + eps_ox * phi(ii-1, 1);
                Jacobian(ii, ii-1) = eps_ox; Jacobian(ii, ii) = -2*eps_ox; Jacobian(ii, ii+1) = eps_ox;
            elseif ii == interface1
                Residue(ii, 1) = eps_si * phi(ii+1, 1) - (eps_si+eps_ox) * phi(ii, 1) + eps_ox * phi(ii-1, 1);
                Jacobian(ii, ii-1) = eps_ox; Jacobian(ii, ii) = -(eps_si+eps_ox); Jacobian(ii, ii+1) = eps_si;
            elseif ii == interface2
                Residue(ii, 1) = eps_ox * phi(ii+1, 1) - (eps_ox+eps_si) * phi(ii, 1) + eps_si * phi(ii-1, 1);
                Jacobian(ii, ii-1) = eps_si; Jacobian(ii, ii) = -(eps_ox+eps_si); Jacobian(ii, ii+1) = eps_ox;
            else
                Residue(ii, 1) = eps_si * phi(ii+1, 1) - 2*eps_si * phi(ii, 1) + eps_si * phi(ii-1, 1);
                Jacobian(ii, ii-1) = eps_si; Jacobian(ii, ii) = -2*eps_si; Jacobian(ii, ii+1) = eps_si;
            end
        end

        % Charge part
        for ii = interface1 : interface2
            if ii == interface1
                Residue(ii, 1) = Residue(ii, 1) - coef * (N_acc + ni*exp(phi(ii, 1)/ThermalV)) * 0.5;
                Jacobian(ii, ii) = Jacobian(ii, ii) - coef * ni*exp(phi(ii, 1)/ThermalV)/ThermalV * 0.5;
            elseif ii == interface2
                Residue(ii, 1) = Residue(ii, 1) - coef * (N_acc + ni*exp(phi(ii, 1)/ThermalV)) * 0.5;
                Jacobian(ii, ii) = Jacobian(ii, ii) - coef * ni*exp(phi(ii, 1)/ThermalV)/ThermalV * 0.5;
            else
                Residue(ii, 1) = Residue(ii, 1) - coef * (N_acc + ni*exp(phi(ii, 1)/ThermalV));
                Jacobian(ii, ii) = Jacobian(ii, ii) - coef * ni*exp(phi(ii, 1)/ThermalV)/ThermalV;
            end
        end

        % Updating potential
        Delphi = Jacobian \ (-Residue);
        phi = phi + Delphi;
    end

    %% 3D and 2D electron density
    n3D = ni*exp(q*phi/(kB*T)); % [#/m^3]
    n2D(jj) = trapz(x, n3D); % [#/m^2]
end
plot(GateVoltages, (n2D/1e+4)/1e+12)
xlabel('Gate Voltage [V]')
ylabel('n_{2D} [10^{12}/cm^2]')