#include <iostream>
#include <cstdio>
#include <cmath>
#include "E_function.h"

#ifndef MOSfunctions
#define MOSfunctions

using namespace std;

extern "C" {
	void dgetrf_(int* M, int* N, double* A, int* LDA, int* IPIV, int* INFO);
	void dgetrs_(char* TRANS, int* N, int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, int* INFO);
	void dgetri_(int* N, double* A, int* LDA, int* IPIV, double* WORK, int* LWORK, int* INFO);
} // refer to lapack library


int Solver(int N, double *A, double *b){

	int i, j, *ipiv = new int[N], nrhs = 1, info, lwork = 4*N;
        char trans = 'N';
	double *work = new double[lwork], *temp = new double[N];

	dgetrf_(&N, &N, A, &N, ipiv, &info);
	dgetri_(&N, A, &N, ipiv, work, &lwork, &info);

	for(i = 0; i < N; i++) {*(temp+i) = *(b+i); *(b+i) = 0;}
        for(i = 0; i < N; i++) for(j = 0; j < N; j++) *(b+i) += A[N*i+j]*temp[j];

	delete ipiv, work, temp;
        return info;
} // Linear Solver


int Pois(int points, double *pois, double *n, double cof1, double cof2, double iv, int Switch){
	int i, j, info;
	double *temp1 = new double[points*points], a = 5e-9, d1 = 1e-9, d2 = 1e-9, term = (a+d1+d2)/(points-1);

	double e0 = cof1*es, e1 = cof2*es, ep, N_acc = 1e24, n_i = 1.075e16, T = 300;

	for(i = 0; i < points; i++){
		if (Switch == 1) N_acc = 1e24 - *(n+i);
		if(i == points-1 || i == 0) *(pois+i) = iv + 0.33374;
		else if(term * i == d1 || term * i == d1+a) *(pois+i) = term*term*q_e*N_acc/2;
		else if(term * i > d1 && term * i < d1+a) *(pois+i) = term*term*q_e*N_acc;
		else *(pois+i) = 0;
	}

	for(i = 0; i < points; i++){
		for(j = 0; j < points; j++) *(temp1+i*points+j) = 0;

		if(i == 0 || i == points-1) {*(temp1+i*points+i) = 1; continue;}
			
		if(i*term == d1) {
			*(temp1+i*points+i) = -e0-e1;
			*(temp1+i*points+i-1) = e0;
			*(temp1+i*points+i+1) = e1;
		}

		else if(i*term == a+d1){
			*(temp1+i*points+i) = -e0-e1;
			*(temp1+i*points+i-1) = e1;
			*(temp1+i*points+i+1) = e0;
		}
		else {
			if(i*term < d1 || i*term > a+d1) ep = e0;
			else ep = e1;

			*(temp1+i*points+i) = -2*ep;
			*(temp1+i*points+i-1) = ep;
			*(temp1+i*points+i+1) = ep;
		}
	}

	info = Solver(points, temp1, pois);

	delete temp1;

	return (info == 0)? true:false;
} // Solve Poisson's equation to find out electrostatic potential

int N_acc(int points, double *pois, double *n){
	int i;
	double n_i = 1.075e16, T = 300, a = 5e-9, d1 = 1e-9, d2 = 1e-9, term = (a+d1+d2)/(points-1);
	for(i = 0; i < points; i++){
		if(term*i < d1 || term*i > d1+a) *(n+i) = 0;
		else *(n+i) = n_i*exp(q_e*(*(pois+i))/(k_b*T));
	}
	return i==(points-1);

}

#endif
