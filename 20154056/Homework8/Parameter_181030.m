function [SysParam] = Parameter_181030()

%% GEOMETRICAL PARAMETERS (only for double gate : oxide | meterial | oxide)
SysParam.L      = [100 100] * 1e-9; % Lx, Ly [m]
SysParam.Deltaz = 0.1e-9; % 0.1 nm spacing [m]
SysParam.Nz     = 61; % 6 nm thick [#]
SysParam.Intf   = [6 56]; % Index of the two interfaces [#]

%% SAMPLE PROPERTIES
SysParam.T         = 300; % Temperature [K]
SysParam.m         = [0.19 0.19 0.91]; % Effective mass : mxx, myy, mzz [#]
SysParam.eps       = [3.9 11.7 3.9]; % Relative permittivity [#]
SysParam.Nacc      = 1e24; % Accepter concentration [#/m^3]
SysParam.ni        = 1.075e16; % Intrinsic carrier concentration [#/m^3]
SysParam.Ei        = 0.33374;
SysParam.EcminusEi = 0.561004; % Ec - Ei [eV]

%% ELEMENTARY PARAMETERS
SysParam.h       = 6.626176e-34; % Plank constant [J-s]
SysParam.hbar    = SysParam.h/(2*pi); % Reduced Plank constant [J-s]
SysParam.q       = 1.602192e-19; % Elementary charge [C]
SysParam.m0      = 9.109534e-31; % Electron rest mass [kg]
SysParam.k_B     = 1.380662e-23; % Boltzmann constant [J/K]
SysParam.eps0    = 8.854187817e-12; % Vacuum permittivity [F/m]
end