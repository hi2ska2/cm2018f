//fucntions for learning DFT
#include <cmath>
#include <iomanip>

#ifndef E_function
#define E_function

#define PI 3.14159265359
#define h_bar 1.05457180e-34
#define h_bar_ev (h_bar/q_e)
#define m_e 9.10938356e-31
#define q_e 1.6021766208e-19
#define es 8.854187817e-12
#define k_b 1.380662e-23

double E_1D(int n, double a){return n*n*PI*PI*h_bar*h_bar/2/m_e/(a*a);}
double psi_1D(int n, double x, double a){return sqrt(2/a)*sin(n*PI*x/a);}
double psi_1D_d2(int n, double x, double a){return -(double)n*n*PI*PI/(a*a)*sqrt(2/a)*sin(n*PI*x/a);}
double N_x(int n, double x, double a){return (n+1)/a-sin((n+1)*(PI*x/a))/(a*sin(PI*x/a))*cos(n*PI*x/a);}
double test_P(double x, double e0, double e1){
	if(x <= 2.5e-9) return e1/(e0+e1)/(2.5e-9)*x;
	else return e0/(e0+e1)*(x-5e-9)/(2.5e-9)+1;
}
double es_r(double r){return es*r;}
#endif
