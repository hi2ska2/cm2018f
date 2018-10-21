#include "Constant.h"

using namespace std;

int main(int argc, char *argv[]){

	int Nz = 51;
	double subbandNumber, totalNumber = 0, term = Lz/(Nz-1);
	double Ef=-0.1*q, Ez, *elec = new double[Nz], sum = 0;

	for(int k=0; k < 11; k++){ 
		for(int i=1; i<=nmax; i++){
			for(int j=0; j < Nz; j++){ 
				Ez = (hbar*hbar)/(2*mzz*m0)*(PI*i/Lz)*(PI*i/Lz);
				subbandNumber = coef*log(1+exp((Ef-Ez)/(k_B*T)));
				totalNumber += subbandNumber;
				*(elec+j) += 2/(Lx*Ly*Lz)*(sin(i*PI*j*term/Lz)*sin(i*PI*j*term/Lz))*subbandNumber;

			}
		}

		for(int j=0; j < Nz; j++) {
			sum += *(elec+j)*term;
	      	 	*(elec+j) = 0;
		}
		cout << Ef/q << ", " << sum << endl;
		Ef += 0.02*q;
		sum = 0;
		
	}
	return 0;
}
	



