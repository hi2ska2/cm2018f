#include <iostream>
#include <cmath>
#include <cstdio>
#include <iomanip>

#define n_int 10e16
#define V_T 0.026

using namespace std;
double A(double N){return V_T*log(N/(2*n_int)+sqrt((N/(2*n_int))*(N/(2*n_int))+1));}
double F(double phi, double N){return N+(n_int*(exp(-(phi/V_T))))-(n_int*(exp(phi/V_T)));}
double dF(double phi){return -(n_int/V_T)*exp(-(phi/V_T))-(n_int/V_T)*exp(phi/V_T);}

int main(int argc, char *argv[]){

	int i, count = 0, points = 65;
	double phi_const=-0.5, phi, dphi, N = -1e16;

	for(i = 0; i < points; i++){
		
		phi = phi_const;
		dphi = phi_const;

		while( abs(dphi/phi) > 1e-13){
	
			dphi = -F(phi, N)/dF(phi);
			phi += dphi; //Newton-Raphson method for 1 equation.

		}

		cout  << setprecision(13) << -N << ", " << A(N) << endl;

		N *= pow(10, 0.125);

	}
	

	return 0;
}

