
//constants
const double q = 1.602192e-19;
const double ep0 = 8.854187817e-12;
const double k_B = 1.380662e-23;
const double T = 300.0;
const double thermal = k_B*T/q;
const double eps = 11.7, epo = 3.9;

//settings
const int Nz = 61, Ny = 31;
const double Lz = 6.0e-9, Ly = 120.0e-9;
const double delz = Lz/(double)(Nz-1), dely = Ly/(double)(Ly-1);
const int interfacez1 = 5, interfacez2 = 55;
const int interfacey1 = 10, interfacey2 = 20;
const double Hdop = 5.0e25, Ldop = 2.0e23;
const double ni = 1.075e16;
