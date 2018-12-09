#include "LA.h"
#include <fstream>
#include "const.cpp"
using namespace std;
int Solver(int, double*, double*);
double *transRC(double freq);
int main(int argc, char *argv[]){

	double freqi=atof(argv[1]);
	transRC(freqi);

	return 0;

}

double *transRC(double freq){
	
	LA A(5);
	double *result = new double[5], *resultOld = new double[5], t=0, *b = new double[5], deltat, w, I;
	result[0] = 0; result[1] = 0; result[2] = 0; result[3] = 1.0; result[4] = 0;
	for(int i = 0; i < 5; i++) b[i] = 0;
	
	deltat = 1.0/freq/100.0;
	w = 2.0*pi*freq;

	A(0,3) = 1.0;
	A(1,1) = 1.0; A(1,3) = -C/deltat; A(1,4) = C/deltat;
	A(2,2) = 1.0; A(2,4) = -1.0/R;
	A(3,0) = 1.0; A(3,1) = 1.0;
	A(4,1) = -1.0; A(4,2) = 1.0;


	for(int i = 1; i < 1000; i++){
		t = (double)i*deltat;
		for(int j = 0; j < 5; j++) resultOld[j] = result[j];
		b[0] = cos(2.0*pi*freq*t);
		b[1] = -C/deltat*(resultOld[3]-resultOld[4]);
		result = A / b;
		I = w*w*R*C*C/(1+w*w*R*R*C*C)*cos(w*t) - w*C/(1+w*w*R*R*C*C)*sin(w*t);	
		cout << t << '\t' << result[4] << '\t' << I*R << endl;
	}

	A.~LA();

	return result;
}


