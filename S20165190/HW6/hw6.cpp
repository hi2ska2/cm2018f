#include "Constants.h"
#include <iostream>
using namespace std;

void N_Rmethod(double *J, double *r, double *phi);
int N(int points, double *pois, double *n);

int main(int argc, char *argv[]){
	int i, j, points = 71;
	double *J, *r, *phi, *n, Vi, V = atof(argv[1]), term = (Lo*2+Ls)/(points-1), sum = 0;

	for(j = 0; j < 11; j++){
	sum = 0;
	Vi = j*0.1;

	J = new double[points*points];
	r = new double[points];
	n = new double[points];
	phi = new double[points];

	for(i = 0; i < points; i++) phi[i] = Vi+0.33374;
	for(i = 0; i < 30; i++) N_Rmethod(J, r, phi);

	N(points, phi, n);

	for(i = 0; i < points; i++) {
		if(i == 10 || i == 60) sum += n[i]/2;
		sum += n[i];
	}

	cout << Vi << ", " << sum << endl;
	}

}

void N_Rmethod(double *J, double *r, double *phi){

	int points=71, i, j;
	double term = (Lo*2+Ls)/(points-1);

	for(i = 0; i < points; i++){
	       	for(j = 0; j < points; j++) J[points*i+j] = 0;

		if(i == 0 || i == points-1) {
			r[i] = 0;
			J[points*i+i] = 1.0;
		}

		else if(i == 10){
			J[points*i+i] = -(eSi+eOx)/term - term*q*ni/VT*exp((*(phi+i))/VT)*0.5;
			J[points*i+i-1] = eOx/term;
			J[points*i+i+1] = eSi/term;
			r[i] = 1/term*(phi[i+1]*eSi-(eOx+eSi)*phi[i]+phi[i-1]*eOx) - term*q*(Nacc+ni*exp(phi[i]/VT))*0.5;
		} // Interface1

		else if(i == 60){
			J[points*i+i] = -(eSi+eOx)/term - term*q*ni/VT*exp((*(phi+i))/VT)*0.5;
			J[points*i+i-1] = eSi/term;
			J[points*i+i+1] = eOx/term;
			r[i] = 1/term*(phi[i+1]*eOx-(eSi+eOx)*phi[i]+phi[i-1]*eSi) - term*q*(Nacc+ni*exp(phi[i]/VT))*0.5;
		} // Interface2

		else if(i > 10 && i < 60){
			J[points*i+i] = (-2.0)*(eSi/term) - term*q*ni/VT*exp((*(phi+i))/VT);
			J[points*i+i-1] = eSi/term;
			J[points*i+i+1] = eSi/term;
			r[i] = eSi/term*(phi[i+1]-2*phi[i]+phi[i-1]) - term*q*(Nacc+ni*exp(phi[i]/VT));
		} // silicon layer

		else {
			J[points*i+i] = (-2.0)*(eOx/term);
			J[points*i+i-1] = eOx/term;
			J[points*i+i+1] = eOx/term;
			r[i] = eOx/term*(phi[i+1]-2*phi[i]+phi[i-1]);
			
		} // dioxide layer
	}

	for(i = 0; i < points; i++) *(r+i) = -(*(r+i));

       	Solver(points, J, r);	

	for(i = 0; i < points; i++) *(phi+i) += *(r+i);

}

int N(int points, double *pois, double *n){
        int i;
        double term = (Lo*2+Ls)/(points-1);
        for(i = 0; i < points; i++){
                if(term*i < Lo || term*i > Lo+Ls) *(n+i) = 0;
                else *(n+i) = ni*exp(q*(*(pois+i))/(k_B*T));
        }
        return i==(points-1);
} // Calculation of Electron density
 

