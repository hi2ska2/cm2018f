#include <iostream>
#include "MOS_function.h"

int main(int argc, char *argv[]){
	int points = atoi(argv[1]), repeat = atoi(argv[4]),  Switch = 0;
	double cof1 = atof(argv[2]), cof2 = atof(argv[3]);

	double *pois = new double[points], *n = new double[points], *p2 = new double[points];
        double d1 = 1e-9, d2 = 1e-9, a = 5e-9, term = (d1+d2+a)/(points-1), iv = 0.1;
	for(int j = 0; j < repeat; j++){
		iv = 0.1*j;
		Switch = 0;
		Pois(points, pois, n, cof1, cof2, iv, Switch);
		N_acc(points, pois, n);
		Switch = 1;
		Pois(points, p2, n, cof1, cof2, iv, Switch);
	
		cout << setprecision(9) << iv << ", "  << *(pois+points/2) - *(p2+points/2)  << endl;
	}

	return 0;

}
