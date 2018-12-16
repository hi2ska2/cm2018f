#include "step3.h"
#include "step3function.cpp"
#include <string>

int main(int argc, char *argv[]){

	string num;	
	double *sch, *phi, *E, *n;
	int valley = atoi(argv[1]);

	sch = new double[points*points];
	phi = new double[points+12];
	n = new double[points];
	E = new double[points];

	ifstream input("05_05.dat");

	for(int j = 0; j < 31; j++){
		for(int i = 0; i < points+12; i++) {input >> num; phi[i] = stof(num);}
		R_Sch(sch,phi,E,valley);
		R_Nx(n,sch,E,valley);
		for(int i = 0; i < points; i++) cout << j*4e-9 << ", " << 0.6e-9+i*0.1e-9 <<", " <<  n[i] << endl;
	}

	return 0;
}
