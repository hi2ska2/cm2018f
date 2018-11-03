#include "hw8.h"
using namespace std;


int main(int argc, char *argv[]){
	int i, j, Swit = atoi(argv[1]), SP =  interface2 - interface1 - 1;
	double *J, *r, *phi, *n, *sch, *E, Vi = 0, term = (Lo*2+Ls)/(points-1), sum = 0;
	char filename[100];

	FILE *fp;

	J = new double[points*points];
	r = new double[points];
	n = new double[SP];
	phi = new double[points];
	sch = new double[SP*SP];
	E = new double[En];



	sprintf(filename, "%d_int.txt", Swit);
	fp = fopen(filename, "w");
	for(j = 0; j < 11; j++){

	for(i = 0; i < points; i++) phi[i] = Vi+0.33374;
	for(i = 0; i < 30; i++) N_Rmethod(J, r, phi);

	N_classic(phi, n);

	if(Swit == 1){
		R_Sch(sch, phi, E);
		R_Nx(n, sch, E);
		for(i = 0 ; i < 20; i++){
			R_Pois(phi, sch, Vi);
			R_Sch(sch, phi, E);
			R_Nx(n, sch, E);
		}
	}

	for(i = 0; i < SP; i++ ) {
		sum += n[i]*term;

	}
	fprintf(fp, "%lf, %lf\n",  Vi, sum);
	sum = 0;
	Vi += 0.1;


	}
	fclose(fp);
	return 0;


}
