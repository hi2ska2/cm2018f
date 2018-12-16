/*
 * LA.h
 *
 *  Created on: 2018. 11. 18.
 *      Author: han
 */

#ifndef LA_H_
#define LA_H_


#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>

using namespace std;


class LA3D {

	private:
		int ysize, zsize, points;
		double *arr;

	public:
		LA3D(int zN, int yN){
			points = zN*yN;
			arr = new double[points*points];
			for (int i = 0; i < points*points; i++) arr[i] = 0;
			ysize = yN; zsize = zN;
		}

		~LA3D(){
		}

		int operator= (LA3D B){
			for(int i = 0; i < points*points; i++) *((this->getP())+i) = *(B.getP()+i);
		}

		void setMat(int s1, int s2, int z, int y, double phi){
			arr[(s1*ysize+s2)*points+z*ysize+y] = phi;
		}

		double* getP(void){return arr;}
		double getMat(int s1, int s2, int x, int y){return arr[(s1*ysize+s2)*points+x*ysize+y];}
		double getySize(void){return ysize;}
		double getzSize(void){return zsize;}

		double& operator() (int s1, int s2, int x, int y){
			return arr[(s1*ysize+s2)*points+x*ysize+y];
		}

		double* operator/ (double *vector);

};

class LA {
	private:
		int size;
		double *arr;

	public:
		LA(int N){
			arr = new double[N*N];
			for (int i = 0; i < N; i++) for(int j = 0; j < N; j++) *(arr+N*i+j) = 0;
			size = N;
		}

		LA(int N, double *vector){
			arr = new double[N*N];
			for (int i = 0; i < N; i++) {
				for(int j = 0; j < N; j++) *(arr+N*i+j) = 0;
				*(arr+N*i+i) = vector[i];
			}
			size = N;
		}

		~LA(){
		}

		int operator= (LA B);
		int operator= (double* B);
		void setMat(int x, int y, double phi);
		void Ldiag(double *vector);
		void Rdiag(double *vector);
		double* getP(void);
		double getMat(int x, int y);
		double getSize(void);
		double abssum(int n);
		LA operator- ();
		LA operator- (LA &B);
		LA operator+ (LA &B);
		double& operator() (int x, int y);
		double* operator/ (double *vector);

};


#endif /* LA_H_ */
