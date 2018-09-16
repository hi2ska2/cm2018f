#include <iostream>
#include <cstdio>
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


int Pois(int points, double *pois, double cof1, double cof2){

	int i, j, info;
	double *temp1 = new double[points*points], d1 = 25e-10, d2 = 25e-10, term = (d1+d2)/(points-1);
	double e0 = cof1*es, e1 = cof2*es, ep;

	for(i = 0; i < points; i++){
		if(i == points-1) *(pois+i) = 1;
		else *(pois+i) = 0;
	}

	for(i = 0; i < points; i++){
		for(j = 0; j < points; j++) *(temp1+i*points+j) = 0;
		if(i == 0 || i == points-1) {*(temp1+i*points+i) = 1; continue;}
		if(i*term != d1) {
			if(i*term < d1) ep = e0;
			else if(i*term > d1) ep = e1;

			*(temp1+i*points+i) = -2*ep;
			*(temp1+i*points+i-1) = ep;
			*(temp1+i*points+i+1) = ep;
		}
		else{
			*(temp1+i*points+i) = -e0-e1;
			*(temp1+i*points+i-1) = e0;
			*(temp1+i*points+i+1) = e1;

		}
	}
	
	info = Solver(points, temp1, pois);

	delete temp1;

	return (info == 0)? true:false;
} // Solve Poisson's equation to find out electrostatic potential

#endif
