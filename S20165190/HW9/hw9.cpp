#include <iostream>
using namespace std;
#include <iomanip>

#define numY 9
#define numZ 5

extern "C" {
        void dgetrf_(int* M, int* N, double* A, int* LDA, int* IPIV, int* INFO);
        void dgetrs_(char* TRANS, int* N, int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, int* INFO);
        void dgetri_(int* N, double* A, int* LDA, int* IPIV, double* WORK, int* LWORK, int* INFO);
}
// refer to lapack library 

int Solver(int, double*, double*);

// comments for toy problem
int main(int argc, char* argv[]){

	int i, j, k, points = numZ*numY, Case = atoi(argv[1]);

	double *toy = new double[points*points], *vb = new double[points], C1, C2, C3;
	
	switch (Case){
		case 1:
			C1 = 1;
			C2 = 0;
			C3 = 0;
			break;

		case 2:
			C1 = 0;
			C2 = 1;
			C3 = 0;
			break;

		case 3:
			C1 = 0;
			C2 = 0;
			C3 = 1;
			break;

		case 4:
			C1 = 1;
			C2 = 1;
			C3 = 1;
			break;

		default:
			cout << "error!\n" << endl;
			return 1;
	}
		

	for(i = 0; i < numZ; i++){

		for(j = 0; j < numY; j++){

			for(k = 0; k < points; k++) toy[(i*numY+j)*points+k] = 0;

			if(i == 0){
				if(j < 2) toy[(i*numY+j)*points+i*numY+j] = 1; // for first blue points
				else if(j >= numY-2) toy[(i*numY+j)*points+i*numY+j] = 1; // for red points
				else {
					toy[(i*numY+j)*points+i*numY+j] = -2;
					toy[(i*numY+j)*points+(i+1)*numY+j] = 1;
					toy[(i*numY+j)*points+i*numY+j-1] = 0.5;
					toy[(i*numY+j)*points+i*numY+j+1] = 0.5;
				}// for black points;
			}
			else if(i == numZ - 1) toy[(i*numY+j)*points+i*numY+j] = 1; // for blue points in last line

			else if(j == 0) {
				toy[(i*numY+j)*points+i*numY+j] = -2;
				toy[(i*numY+j)*points+(i+1)*numY+j] = 0.5;
				toy[(i*numY+j)*points+(i-1)*numY+j] = 0.5;
				toy[(i*numY+j)*points+i*numY+j+1] = 1;
			} // for black points in left line

			else if(j == numY - 1) {
				toy[(i*numY+j)*points+i*numY+j] = -2;
				toy[(i*numY+j)*points+(i+1)*numY+j] = 0.5;
				toy[(i*numY+j)*points+(i-1)*numY+j] = 0.5;
				toy[(i*numY+j)*points+i*numY+j-1] = 1;
			} // for black points in right line

			else {
				toy[(i*numY+j)*points+i*numY+j] = -4;
				toy[(i*numY+j)*points+(i-1)*numY+j] = 1;
				toy[(i*numY+j)*points+(i+1)*numY+j] = 1;
				toy[(i*numY+j)*points+i*numY+(j-1)] = 1;
				toy[(i*numY+j)*points+i*numY+(j+1)] = 1;
			} // for other points;

		}
	}

	for(i = 0; i < numZ; i++) for(j = 0; j < numY; j++){
		if(i == 0){
			if(j < 2) vb[i*numY+j] = C2; //for first blue points
			else if(j >= numY - 2) vb[i*numY+j] = C1; //for first red points
			else vb[i*numY+j] = 0; 	
		}

		else if(i == numZ-1) vb[i*numY+j] = C3;

		else if(j == 0) vb[i*numY+j] = 0;

		else if(j == numY - 1) vb[i*numY+j] = 0;

		else vb[i*numY+j] = 0;

	}

	Solver(points, toy, vb);

	for(i = 0; i < numZ; i++) { for(j = 0; j < numY; j++) cout << setprecision(16) << i+1 << ", " << j+1 << ", " << vb[i*numY+j] << endl;}

	return 0;

}


int Solver(int N, double *A, double *b){

        int i, j, *ipiv = new int[N], info, lwork = 4*N;
        double *work = new double[lwork], *temp = new double[N];

        dgetrf_(&N, &N, A, &N, ipiv, &info);
        dgetri_(&N, A, &N, ipiv, work, &lwork, &info);

        for(i = 0; i < N; i++) {*(temp+i) = *(b+i); *(b+i) = 0;}
        for(i = 0; i < N; i++) for(j = 0; j < N; j++) *(b+i) += A[N*i+j]*temp[j];

       // delete ipiv, work, temp;
        return info;
}
// Linear Solver


