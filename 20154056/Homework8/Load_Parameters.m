%% LOAD PARAMETERS
h         = SysParam.h;
hbar      = SysParam.hbar;
q         = SysParam.q;
m0        = SysParam.m0;
k_B       = SysParam.k_B;
eps0      = SysParam.eps0;

L         = SysParam.L;
Deltaz    = SysParam.Deltaz;
Nz        = SysParam.Nz;
Intf      = SysParam.Intf;

T         = SysParam.T;
m         = SysParam.m;
eps       = SysParam.eps;
Nacc      = SysParam.Nacc;
ni        = SysParam.ni;
Ei        = SysParam.Ei;
EcminusEi = SysParam.EcminusEi;

%% DEFINITION OF FREQUENTLY USED PARAMETERS
z        = Deltaz * transpose([0:Nz-1]); % real space [m]
thermal  = k_B * T / q; % Thermal voltage [V]
Coef_Poi = Deltaz^2 * q / eps0; % Coefficient for fixted source Poisson equation.
Coef_Sch = 2 * L(1)*L(2)/(2*pi) * sqrt(m(1)*m(2))*m0/hbar^2 * k_B*T; % Coefficient for when it calulates total number of electrons.