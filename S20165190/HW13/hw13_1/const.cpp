/*
 * const.cpp
 *
 *  Created on: 2018. 11. 18.
 *      Author: han
 */

#include "LA.h"

// Constants
const double q = 1.602192e-19; // Elementary charge, C
const double eps0 = 8.854187817e-12; // Vacuum permittivity, F/m
const double k_B = 1.380662e-23; // Boltzmann constant, J/K
const double pi = 3.141592; // PI
const double ni = 1.075e16; // Intrinsic electron density, /m^3
const double eSi = 11.7; const double eOx = 3.9; // Relative permittivity
const double Dn = 0.01;

// Settings
double Length = 600e-9; // Length of system, m
const double T = 300; // Temperature, K
const double thermal = k_B*T/q; // Thermal voltage, V
double Hdop = 5e23;
double Ldop = 2e21;

// Discretization
int N;
const int VN = 10;
int interface1;
int interface2;
double termx;
double coef;

