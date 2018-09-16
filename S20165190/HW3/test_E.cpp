#include <iostream>
#include "MOS_function.h"

int main(int argc, char *argv[]){
	int points = atoi(argv[1]);
	double cof1 = atof(argv[2]), cof2 = atof(argv[3]);

	double *pois = new double[points], d1 = 25e-10, d2 = 25e-10, term = (d1+d2)/(points-1);
	
	for(int i = 0; i < points; i++) cout << setprecision(9) <<i*term << ", " << test_P(i*term, cof1, cof2) << endl;

	return 0;

}
