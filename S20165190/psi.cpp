#include <iostream>
#include <iomanip>
#include "psi_f.h"

using namespace std;

int main(int argc, char *argv[]){
	
	if(argc != 3) {cout << "There should be 2 arguments only." << endl; return 1;}
	int N = atoi(argv[1]), DN = atoi(argv[2])-2, T;

	double *E = new double[DN], *psi = new double[DN*DN], *E_d = new double[DN], term = 5e-9/(DN+1);
/*
	for(int i = 1; i <= DN; i++) {
		T = psix(N, i, psi, *(E_d+i));
		cout << setprecision(9) << i << ", " << *(E_d+i) << ", " << E_1D(N, 5e-9) << ", " << T << endl;
	}
*/	//code to print eigen energy


	
	T = psix(N, DN, psi, E); // This instruction decides what to print among N(x) and psi(x)
	
	for(int i = 0; i < N; i++) {
		cout << "Eigen energy is " << *(E+i) << "J" << endl;
		cout << "Psis are " << endl;
		cout << "0, 0" << endl;
		for(int j = 0; j < DN; j++) cout << term*(j+1) << ", "<< *(psi+DN*i+j) << endl;
		cout << "5e-9, 0" << endl;
	}

	//Code to print N(x) or psi(x)


	return 0;
}	
