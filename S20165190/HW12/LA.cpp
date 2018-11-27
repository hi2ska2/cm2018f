/*
 * LA.cpp
 *
 *  Created on: 2018. 11. 18.
 *      Author: han
 */
#include "LA.h"

extern "C" {
        void dgetrf_(int* M, int* N, double* A, int* LDA, int* IPIV, int* INFO);
        void dgetrs_(char* TRANS, int* N, int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, int* INFO);
        void dgetri_(int* N, double* A, int* LDA, int* IPIV, double* WORK, int* LWORK, int* INFO);
}
// refer to lapack library


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

		void LA::setMat(int x, int y, double phi){
			*(arr+size*x+y) = phi;
		}
		void LA::Ldiag(double *vector){
			for (int i = 0; i < size; i++)
				for(int j = 0; j < size; j++) *(arr+size*i+j) *= vector[i];
		}
		void LA::Rdiag(double *vector){
			for (int i = 0; i < size; i++)
				for(int j = 0; j < size; j++) *(arr+size*j+i) *= vector[i];
		}

		double* LA::getP(void){
			return arr;
		}

		double LA::getMat(int x, int y){
			return *(arr+size*x+y);
		}

		double LA::getSize(void){
			return size;
		}

		double LA::abssum(int n){
			double sum = 0;
			for(int i = 0; i < size; i++) sum += *(arr+n*size+i);
			return sum;
		}

		LA LA::operator- (){
			LA temp(size);
			for(int i = 0; i < size; i++) for(int j = 0; j < size; j++) temp(i,j) = -(*this)(i,j);
			return temp;
		}

		LA LA::operator- (LA& B){
				if(size != B.getSize()) { cout << "Error (different size)" << endl; return 3;}
				LA temp(size);
				for(int i = 0; i < size; i++) for(int j = 0; j < size; j++) temp(i,j) = (*this)(i,j) - B(i,j);
				return temp;
		}

		LA LA::operator+ (LA& B){
			if(size != B.getSize()) { cout << "Error (different size)" << endl; return 3;}
			LA temp(size);
			for(int i = 0; i < size; i++) for(int j = 0; j < size; j++) temp(i,j) = (*this)(i,j) + B(i,j);
			return temp;
		}

		double* LA::operator/ (double* vector){

			double *temp = new double[size];
			for(int i=0;i<size;i++)
				temp[i] = vector[i];

			Solver(size, (*this).getP(), temp);
			return temp;
		}

		double& LA::operator() (int x, int y) {
			return *(arr+size*x+y);
		}


		int LA::operator= (LA B){
			if(size != B.getSize()) { cout << "Error (different size)" << endl; return 3;}
			for(int i = 0; i < size; i++) for(int j = 0; j < size; j++) (*this)(i,j) = B(i,j);
			return 0;
		}

		int LA::operator= (double* B){
			for(int i = 0; i < size; i++) for(int j = 0; j < size; j++) (*this)(i,j) = B[size*i+j];
			return 0;
		}

