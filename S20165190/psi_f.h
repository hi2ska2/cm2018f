#include <iostream>
#include <cstdio>
#include "E_function.h"

#ifndef Nfunctions
#define Nfunctions

using namespace std;

extern "C" void dsyevx_(char* JOBZ, char* RANGE, char* UPLO, int* N, double* A, int* LDA, double* VL, double* VU, int* IL, int* IU, double* ABSTOL, int* M, double* W, double* Z, int* LDZ, double* WORK, int* LWORK, int* IWORK, int* IFAIL, int *INFO);
extern "C" void dgesv_(int* N, int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, int* INFO);

void setMatA1(int n, double *A); //Set matrix for 1D infinity potential well
int psix(int N, int DN, double *psi, double *E); // N x DN matrix for psi(x)
int Nx(int N, int DN, double *n); // DN matrix for N(x)


void setMatA1(int n, double *A){

	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) A[n*i+j] = 0;

		A[n*i+i] += 2;

		if(!i)  A[n*i+1] += -1;
		else if(i==(n-1)) A[n*i+n-2] += -1;
		else{
			A[n*i+i-1] += -1;
			A[n*i+i+1] += -1;
		}
	}
}

int psix(int N, int DN, double *psi, double *E){
	
	double a = 5e-9, sum = 0, term = a/(DN+1);
	int judge = 0;

	double *A;
        A = new double[DN*DN];

	setMatA1(DN, A);
	
	char jobz = 'V', range = 'I', uplo = 'U';
	int n, lda, info, ldz, ldvr, lwork, il=1, iu=N, m;
	double vl, vu;
	n=lda=ldz=DN;

	lwork=8*DN;
	int *iwork = new int[DN*5];
	int *ifail = new int[DN];
	double *w = new double[DN];
	double *work = new double[lwork];

	dsyevx_(&jobz, &range, &uplo, &n, A, &lda, &vl, &vu, &il, &iu, &vl, &m, w, psi, &ldz, work, &lwork, iwork, ifail, &info);
	for(int j = 0; j < N; j++){
  		for(int i = 0; i < DN; i++) sum += psi[DN*j+i]*psi[DN*j+i]*term;
		for(int i = 0; i < DN; i++) psi[DN*j+i] /= sqrt(sum);
		sum = 0;
	}
	for(int i = 0; i < N; i++){ judge += (ifail[i] || (E <= 0)); *(E+i) = w[i]*(h_bar/(2*m_e*0.19*term*term)*h_bar);}

	
	delete A, iwork, ifail, w, work;
	return judge;

}

int Nx(int N, int DN, double *n){
	int T, i, j;
	double *temp1 = new double[DN], *temp2 = new double[DN*DN], temp3;

	T = psix(N, DN, temp2, temp1);
	for(i = 0; i < DN; i++) for(j = 0; j < i; j++){
		temp3 = *(temp2+DN*i+j);
		*(temp2+DN*i+j) = *(temp2+DN*j+i);
		*(temp2+DN*j+i) = temp3;
	}

	for(i = 0; i < DN; i++) *(n+i) = 0;
	for(i = 0; i < DN; i++) for(j = 0; j < N; j++) *(n+i) += (*(temp2+DN*i+j))*(*(temp2+DN*i+j));

	delete temp1, temp2;
	return T;
}



#endif
