/*
 * hw12.cpp
 *
 *  Created on: 2018. 11. 25.
 *      Author: han
 */

#include "LA.h"
#include "const.cpp"

double **SCeq(double Vapp);
double Ber(double x);

int main(int argc, char* argv[]){
	double **result, *phi, *elec, *J, Vapp=0.0;
	N = atoi(argv[2]);
	J = new double[N-1];
	if(*argv[1] == 's' || *argv[1] == 'S') {
		Hdop = 5e25;
		Ldop = 2e23;
		Length = 120e-9;
		interface1 = (N-1)/3;
		interface2 = (N-1)*2/3;
	}
	else{
		interface1 = (N-1)/6;
		interface2 = (N-1)*5/6;
	}
	termx = Length/(double)(N-1);
	coef = termx*termx*q/eps0;

	for(int j = 0; j < 11; j++)
	{
		result = SCeq(Vapp);
		phi = *result;
		elec = *(result+1);
	//	for(int i = 0; i < N; i++) cout << elec[i] << '\t' << phi[i] << endl;
		for(int i = 0; i < N-1; i++) J[i] = q*Dn/termx*(elec[i+1]*Ber((phi[i+1]-phi[i])/thermal)-elec[i]*Ber((phi[i]-phi[i+1])/thermal));

		cout << Vapp << " ";
		for(int i = 0; i < N-1; i++) cout << J[i] << " ";
		cout << endl;
		Vapp += 0.05;
		delete(phi, elec, result);
	}
	return 0;
}


double **SCeq(double Vapp)
{

	LA jaco(2*N);
	double **result = new double*[2];

	double *res = new double[2*N], *update, *elec = new double[N], *phi = new double[N];
	double *Ndon = new double[N], *Cvector = new double[2*N], *Rvector = new double[2*N];
	double n_av, dphidx, delecdx, Jn, sum;
	result[0] = phi; result[1] = elec;


	for(int i = 0; i < N; i++)
	{
		if(i <= interface1 || i >= interface2) Ndon[i] = Hdop;
		else Ndon[i] = Ldop;

		phi[i] = thermal*log(Ndon[i]/ni);
		elec[i] = ni*exp(phi[i]/thermal);
	}

	for(int r = 0; r < VN; r++)
	{
		for(int i = 0; i < 2*N; i++) {
			res[i] = 0;
			for(int j = 0; j < 2*N; j++) jaco(i,j) = 0;
		}

		res[0] = phi[0] - thermal*log(Ndon[0]/ni);
		jaco(0,0) = 1.0;

		for(int i = 1; i < N-1; i++)
		{
			res[2*i] = eSi*(phi[i+1] - 2*phi[i] + phi[i-1]) + coef*(Ndon[i]-elec[i]);
			jaco(2*i,2*i-2) = eSi;
			jaco(2*i,2*i) = -2.0*eSi;
			jaco(2*i,2*i+2) = eSi;
			jaco(2*i,2*i+1) = -coef;
		}
		res[2*N-2] = phi[N-1] - thermal*log(Ndon[N-1]/ni) - Vapp;
		jaco(2*N-2, 2*N-2) = 1.0;

		res[1] = elec[0] - Ndon[0];
		jaco(1,1) = 1.0;
		for(int i = 1; i < N-1; i++)
		{
			n_av = 0.5*(elec[i+1]+elec[i]);
			dphidx = (phi[i+1]-phi[i])/termx;
			delecdx = (elec[i+1]-elec[i])/termx;
			Jn = n_av * dphidx - thermal * delecdx;
			res[2*i+1] += Jn;
			jaco(2*i+1, 2*i+3) += 0.5*dphidx - thermal / termx;
			jaco(2*i+1, 2*i+1) += 0.5*dphidx + thermal / termx;
			jaco(2*i+1, 2*i+2) += n_av / termx;
			jaco(2*i+1, 2*i) += - n_av / termx;
			res[2*i+3] += -Jn;
			jaco(2*i+3, 2*i+3) += -0.5*dphidx + thermal / termx;
			jaco(2*i+3, 2*i+1) += -0.5*dphidx - thermal / termx;
			jaco(2*i+3, 2*i+2) += - n_av / termx;
			jaco(2*i+3, 2*i) += n_av / termx;
		}
		res[2*N-1] = elec[N-1] - Ndon[N-1];
		jaco(2*N-1,2*N-1) = 1.0;

		for(int i = 0; i < N; i++)
		{
			Cvector[2*i] = thermal;
			Cvector[2*i+1] = sqrt((double)2*interface1*(Hdop*Hdop)+(double)(N-interface1-interface2)*(Ldop*Ldop));
		}
		jaco.Rdiag(Cvector);
		sum = jaco.abssum(1);
		for(int i = 0; i < 2*N; i++) Rvector[i] = 1.0/sum;
		jaco.Ldiag(Rvector);
		for(int i = 0; i < 2*N; i++) res[i] *= Rvector[i];
		update = -jaco / (res);
		for(int i = 0; i < 2*N; i++) update[i] *= Cvector[i];

		for(int i = 0; i < N; i++) phi[i] += update[2*i];
		for(int i = 0; i < N; i++) elec[i] += update[2*i+1];
	}

	jaco.~LA();
	delete(res, Ndon, Cvector, Rvector);
	return result;
}

double Ber(double x){ return x/(exp(x)-1.0); }
