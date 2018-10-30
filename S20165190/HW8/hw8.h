#include "Constants.h"

#ifndef Nfunctions
#define Nfunctions

using namespace std;

extern "C" void dsyevx_(char* JOBZ, char* RANGE, char* UPLO, int* N, double* A, int* LDA, double* VL, double* VU, int* IL, int* IU, double* ABSTOL, int* M, double* W, double* Z, int* LDZ, double* WORK, int* LWORK, int* IWORK, int* IFAIL, int *INFO);
// Reffer to lapack library

int R_Sch(double *sch, double *pois, double *E){

	int i, j, judge = 0, SP = interface2 - interface1 - 1;

	double *temp1 = new double[SP*SP], *V_sch = new double[SP], sum = 0, term = (2*Lo+Ls)/(points-1);

	for(i = 0; i < SP; i++) *(V_sch+i) = q*Ec_Ei - q*(*(pois+i+interface1+1));
	for(i = 0; i < SP; i++){
		for(j = 0; j < SP; j++) *(temp1+i*SP+j) = 0;
		*(temp1+SP*i+i) = 2 + *(V_sch+i)/hbar*2*mzz/hbar*term*term;
		if(i > 0) *(temp1+SP*i+i-1) = -1;
		if(i < (SP-1))*(temp1+SP*i+i+1) = -1;
	}

	char jobz = 'V', range = 'I', uplo = 'U';
	int n, lda, info, ldz, ldvr, lwork, il=1, iu=En, m;
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


	for(i = 0; i < iu; i++){ judge += (ifail[i] || (E <= 0)); *(E+i) = w[i]*(hbar/(2*mzz*term*term)*hbar);}
	
	delete temp1, iwork, ifail, w, work;
	return judge;
}


int R_Nx(double *n, double *sch, double *E){
	int i, j;
	int SP = interface2 - interface1 -1;
	for(i = 0; i < SP; i++) *(n+i) = 0;
	for(i = 0; i < En; i++) for(j = 0; j < SP; j++)
		*(n+j) += ((*(sch+SP*i+j))*(*(sch+SP*i+j)))/(Lx*Ly)*coef_Sch*log(1+exp(-(*(E+i))/(k_B*T)));
	return j;
}

int R_Pois(double *pois, double *n, double iv){
	int i, j, info;
	double *temp1 = new double[points*points], ep, term = (2*Lo+Ls)/(points-1);

	for(i = 0; i < points; i++){
		if(i == points-1 || i == 0) *(pois+i) = iv + 0.33374;
		else if(i == interface1 || i == interface2) *(pois+i) = term*term*q*Nacc/2;
		else if(i > interface1 && i < interface2) *(pois+i) = term*term*q*(Nacc+n[i-interface1-1]);
		else *(pois+i) = 0;
	}

	for(i = 0; i < points; i++){
		for(j = 0; j < points; j++) *(temp1+i*points+j) = 0;

		if(i == 0 || i == points-1) {*(temp1+i*points+i) = 1; continue;}

		if(i == interface1) {
			*(temp1+i*points+i) = -eOx-eSi;
			*(temp1+i*points+i-1) = eOx;
			*(temp1+i*points+i+1) = eSi;
		}

		else if(i == interface2){
			*(temp1+i*points+i) = -eOx-eSi;
			*(temp1+i*points+i-1) = eSi;
			*(temp1+i*points+i+1) = eOx;
		}
		else {
			if(i < interface1 || i > interface2) ep = eOx;
			else ep = eSi;

			*(temp1+i*points+i) = -2*ep;
			*(temp1+i*points+i-1) = ep;
			*(temp1+i*points+i+1) = ep;
		}
	}

	info = Solver(points, temp1, pois);

	delete temp1;

	return (info == 0)? true:false;
} // Solve Poisson's equation to find out electrostatic potential

void N_Rmethod(double *J, double *r, double *phi){

	int i, j;
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

int N_classic(double *pois, double *n){
        int i, SP = interface2 - interface1 - 1;
        double term = (Lo*2+Ls)/(points-1);
        for(i = 0; i < SP; i++){
                 *(n+i) = ni*exp(q*(*(pois+i+interface1+1))/(k_B*T));
        }
        return i==(points-1);
} // Calculation of Electron density




// This code is used for self-consistent

#endif
