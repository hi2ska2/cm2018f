#ifndef Rfunctions
#define Rfunctions

#include "step3.h"


extern "C" void dsyevx_(char* JOBZ, char* RANGE, char* UPLO, int* N, double* A, int* LDA, double* VL, double* VU, int* IL, int* IU, double* ABSTOL, int* M, double* W, double* Z, int* LDZ, double* WORK, int* LWORK, int* IWORK, int* IFAIL, int *INFO);
// Reffer to lapack library

int R_Sch(double *sch, double *pois, double *E, int valley){

	int i, j, judge = 0, SP = points;
	double *mass = new double[3];
	for(i = 0; i < 3; i++) mass[i] = 0.19*m0;
	mass[valley] = 0.91*m0;

	double coef_Sch=(2*Lx*Ly/(2*PI)*sqrt(mass[0]*mass[1])/(hbar*hbar)*(k_B*T));
	double *temp1 = new double[SP*SP], *V_sch = new double[SP], sum = 0, term = 0.1e-9;

	for(i = 0; i < SP; i++) *(V_sch+i) = q*Ec_Ei - q*(*(pois+i+5));
	for(i = 0; i < SP; i++){
		for(j = 0; j < SP; j++) *(temp1+i*SP+j) = 0;
		*(temp1+SP*i+i) = 2 + *(V_sch+i)/hbar*2*mass[2]/hbar*term*term;
		if(i > 0) *(temp1+SP*i+i-1) = -1;
		if(i < (SP-1))*(temp1+SP*i+i+1) = -1;
	}

	char jobz = 'V', range = 'I', uplo = 'U';
	int n, lda, info, ldz, ldvr, lwork, il=1, iu=1, m;
	double vl, vu;
	n=lda=ldz=SP;

	lwork=8*SP;
	int *iwork = new int[SP*5];
	int *ifail = new int[SP];
	double *w = new double[SP];
	double *work = new double[lwork];

	dsyevx_(&jobz, &range, &uplo, &n, temp1, &lda, &vl, &vu, &il, &iu, &vl, &m, w, sch, &ldz, work, &lwork, iwork, ifail, &info);
	for(j = 0; j < SP;j++){
		for(i = 0; i < SP; i++) sum += sch[SP*j+i]*sch[SP*j+i]*term;
		for(i = 0; i < SP; i++) sch[SP*j+i] /= sqrt(sum);
		sum = 0;
	}


	for(i = 0; i < iu; i++){ judge += (ifail[i] || (E <= 0)); *(E+i) = w[i]*(hbar/(2*mass[2]*term*term)*hbar);}
	
	delete temp1, iwork, ifail, w, work;
	return judge;
}


int R_Nx(double *n, double *sch, double *E, int valley){
	int i, j;
	int SP = points;       
       	double *mass = new double[3];
        for(i = 0; i < 3; i++) mass[i] = 0.19*m0;
        mass[valley] = 0.91*m0;
	double coef_Sch=(2*Lx*Ly/(2*PI)*sqrt(mass[0]*mass[1])/(hbar*hbar)*(k_B*T));
	for(i = 0; i < SP; i++) *(n+i) = 0;
	for(i = 0; i < 1; i++) for(j = 0; j < SP; j++)
		*(n+j) += ((*(sch+SP*i+j))*(*(sch+SP*i+j)))/(Lx*Ly)*coef_Sch*log(1+exp(-(*(E+i))/(k_B*T)));
	return j;
}
#endif
