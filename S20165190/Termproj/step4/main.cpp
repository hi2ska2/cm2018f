/*
 * main.cpp
 *
 *  Created on: 2018. 11. 18.
 *      Author: han
 */

#include "LA.h"
#include "const.cpp"

double* F_0(double*);

int main(int argc, char *argv[]){
	double *f0,  Vi = 0, *psi = new double[31];
	ifstream file_f0;

	file_f0.open("05_01_4.dat");
	for(int j = 0 ; j < 61; j++) {
		
		for(int i = 0 ; i < 31; i++) file_f0 >> psi[i];
		f0 = F_0(psi);
		for(int i=0; i<N; i++) cout << 0.1e-9*j << ", " << 4e-9*i << ", " << f0[i] << endl;
	}

		file_f0.close();

	return 0;

}

double* F_0(double *psi){
	LA s(N);
	double *V = new double[N];
	double *b = new double[N];

	double c1, c2;

	for(int i=0; i < N; i++) V[i] = Ec_Ei - *(psi+i);

	const double fs = sqrt(2*pi)/(1+exp(q*H/(k_B*T)));
	const double fd = sqrt(2*pi)/(1+exp(q*(H+V[N-1])/(k_B*T)));

	for(int i=0; i < N; i++)
		if(i == 0 || i == N-1) s(i,i) = 1.0;
		else
		{
			c1 = H + 0.5*(V[i]+V[i-1]);
			c2 = H + 0.5*(V[i]+V[i+1]);
			s(i,i-1) = c1; s(i,i) = -c1-c2; s(i,i+1) = c2;
		}

	for(int i=0; i < N; i++)
		if(i == 0) b[i] = fs;
		else if(i == N-1) b[i] = fd;
		else b[i] = 0;

	return s/b;
}
