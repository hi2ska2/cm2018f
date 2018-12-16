#include "LA.h"
#include "const.cpp"
double *NP2D(double, double, double);

int main(void){
	double Vg = 0.5, Vs = 0.0, Vd = 0.5, *phi;

	phi = NP2D(Vg, Vs, Vd);

	for(int i= 0; i < Nz; i++) for(int j = 0; j < Ny; j++)  cout << phi[Ny*i+j] << endl;

	return 0;
	
}

double *NP2D(double Vg, double Vs, double Vd){
	LA3D Jaco(Nz, Ny);
	double *res = new double[Nz*Ny], *phi = new double[Nz*Ny], *update = new double[Nz*Ny], *Ndon = new double[Nz*Ny];


	for(int i = 0; i < Nz; i++){
		for(int j = 0; j < Ny; j++){
			if(i < interfacez1 || i > interfacez2) Ndon[i*Ny+j] = ni;
			else if(i == interfacez1 || i == interfacez2) 
				if(j < interfacey1 || j > interfacey2) Ndon[i*Ny+j] = Hdop/2.0;
				else Ndon[i*Ny+j] = Ldop/2.0;
			else if(j < interfacey1 || j > interfacey2) Ndon[i*Ny+j] = Hdop;
			else Ndon[i*Ny+j] = Ldop;
		}
	}

	for(int i = 0; i < Nz; i++)
		for( int j = 0; j < Ny; j++){ phi[i*Ny+j] = thermal*log(Ndon[i*Ny+j]/ni); }

	for(int NR = 0; NR < 5; NR++){

	for(int i = 0; i < Nz; i++){
		for(int j = 0; j < Ny; j++){

			if(i < interfacez1 || i > interfacez2){
				if(i == 0 && j == 0){
	                                Jaco(i,j,i,j) = -0.5*epo*(delz/dely+dely/delz);
	                                Jaco(i,j,i,j+1) = 0.5*epo*delz/dely;
	                                Jaco(i,j,i+1,j) = 0.5*epo*dely/delz;
	                                res[i*Ny+j] = Jaco(i,j,i,j+1)*phi[i*Ny+j+1]+Jaco(i,j,i,j)*phi[i*Ny+j]+Jaco(i,j,i+1,j)*phi[(i+1)*Ny+j];
				}
				else if(i == 0 && j == Ny-1){
	                                Jaco(i,j,i,j) = -0.5*epo*(delz/dely+dely/delz);
	                                Jaco(i,j,i,j-1) = 0.5*epo*delz/dely;
	                                Jaco(i,j,i+1,j) = 0.5*epo*dely/delz;
	                                res[i*Ny+j] = Jaco(i,j,i,j-1)*phi[i*Ny+j-1]+Jaco(i,j,i,j)*phi[i*Ny+j]+Jaco(i,j,i+1,j)*phi[(i+1)*Ny+j];
				}
				else if(i == Nz-1 && j == 0){
	                                Jaco(i,j,i,j) = -0.5*epo*(delz/dely+dely/delz);
	                                Jaco(i,j,i,j+1) = 0.5*epo*delz/dely;
	                                Jaco(i,j,i-1,j) = 0.5*epo*dely/delz;
	                                res[i*Ny+j] = Jaco(i,j,i,j+1)*phi[i*Ny+j+1]+Jaco(i,j,i,j)*phi[i*Ny+j]+Jaco(i,j,i-1,j)*phi[(i-1)*Ny+j];
				}
				else if(i == Nz-1 && j == Ny-1){
	                                Jaco(i,j,i,j) = -0.5*epo*(delz/dely+dely/delz);
	                                Jaco(i,j,i,j-1) = 0.5*epo*delz/dely;
	                                Jaco(i,j,i-1,j) = 0.5*epo*dely/delz;
	                                res[i*Ny+j] = Jaco(i,j,i,j-1)*phi[i*Ny+j-1]+Jaco(i,j,i,j)*phi[i*Ny+j]+Jaco(i,j,i-1,j)*phi[(i-1)*Ny+j];
				}
				else if(i == Nz-1 && j >= interfacey1 && j <= interfacey2){
					Jaco(i,j,i,j) = 1.0;
					res[i*Ny+j] = phi[i*Ny+j] - thermal*log(Ndon[i*Ny+j]/ni) - (Vg + 0.33374);
				}
				else if(i == Nz-1){
					Jaco(i,j,i,j) = -1.0*epo*(delz/dely+dely/delz);
					Jaco(i,j,i,j-1) = 0.5*epo*delz/dely;
					Jaco(i,j,i,j+1) = 0.5*epo*delz/dely;
					Jaco(i,j,i-1,j) = epo*dely/delz;
					res[i*Ny+j] = Jaco(i,j,i,j-1)*phi[i*Ny+j-1]+Jaco(i,j,i,j+1)*phi[i*Ny+j+1]+Jaco(i,j,i,j)*phi[i*Ny+j]+Jaco(i,j,i-1,j)*phi[(i-1)*Ny+j];
				}
				else if(i == 0 && j >= interfacey1 && j <= interfacey2){
					Jaco(i,j,i,j) = 1.0;
					res[i*Ny+j] = phi[i*Ny+j] - thermal*log(Ndon[i*Ny+j]/ni) - (Vg + 0.33374);
				}
				else if(i == 0){
					Jaco(i,j,i,j) = -1.0*epo*(delz/dely+dely/delz);
					Jaco(i,j,i,j-1) = 0.5*epo*delz/dely;
					Jaco(i,j,i,j+1) = 0.5*epo*delz/dely;
					Jaco(i,j,i+1,j) = epo*dely/delz;
					res[i*Ny+j] = Jaco(i,j,i,j-1)*phi[i*Ny+j-1]+Jaco(i,j,i,j+1)*phi[i*Ny+j+1]+Jaco(i,j,i,j)*phi[i*Ny+j]+Jaco(i,j,i+1,j)*phi[(i+1)*Ny+j];
				}
				else if(j == Ny-1){
					Jaco(i,j,i,j) = -1.0*epo*(delz/dely+dely/delz);
					Jaco(i,j,i,j-1) = epo*delz/dely;
					Jaco(i,j,i-1,j) = 0.5*epo*dely/delz;
					Jaco(i,j,i+1,j) = 0.5*epo*dely/delz;
					res[i*Ny+j] = Jaco(i,j,i,j-1)*phi[i*Ny+j-1]+Jaco(i,j,i,j)*phi[i*Ny+j]+Jaco(i,j,i-1,j)*phi[(i-1)*Ny+j]+Jaco(i,j,i+1,j)*phi[(i+1)*Ny+j];
	
				}
				else if(j == 0){
					Jaco(i,j,i,j) = -1.0*epo*(delz/dely+dely/delz);
					Jaco(i,j,i,j+1) = epo*delz/dely;
					Jaco(i,j,i-1,j) = 0.5*epo*dely/delz;
					Jaco(i,j,i+1,j) = 0.5*epo*dely/delz;
					res[i*Ny+j] = Jaco(i,j,i,j+1)*phi[i*Ny+j+1]+Jaco(i,j,i,j)*phi[i*Ny+j]+Jaco(i,j,i-1,j)*phi[(i-1)*Ny+j]+Jaco(i,j,i+1,j)*phi[(i+1)*Ny+j];
				}
				else {
					Jaco(i,j,i,j) = -2.0*epo*(delz/dely+dely/delz);
					Jaco(i,j,i,j-1) = epo*delz/dely;
					Jaco(i,j,i,j+1) = epo*delz/dely;
					Jaco(i,j,i-1,j) = epo*dely/delz;
					Jaco(i,j,i+1,j) = epo*dely/delz;
					res[i*Ny+j] = Jaco(i,j,i,j-1)*phi[i*Ny+j-1]+Jaco(i,j,i,j+1)*phi[i*Ny+j+1]+Jaco(i,j,i,j)*phi[i*Ny+j]+Jaco(i,j,i-1,j)*phi[(i-1)*Ny+j]+Jaco(i,j,i+1,j)*phi[(i+1)*Ny+j];
				}
			}
			else if(i == interfacez1){
				Jaco(i,j,i,j) = -(epo+eps);
				Jaco(i,j,i-1,j) = epo;
				Jaco(i,j,i+1,j) = eps;
				res[i*Ny+j] = Jaco(i,j,i-1,j)*phi[(i-1)*Ny+j]+Jaco(i,j,i+1,j)*phi[(i+1)*Ny+j]+Jaco(i,j,i,j)*phi[i*Ny+j];
			}
			else if(i == interfacez2){
				Jaco(i,j,i,j) = -(epo+eps);
				Jaco(i,j,i-1,j) = eps;
				Jaco(i,j,i+1,j) = epo;
				res[i*Ny+j] = Jaco(i,j,i-1,j)*phi[(i-1)*Ny+j]+Jaco(i,j,i+1,j)*phi[(i+1)*Ny+j]+Jaco(i,j,i,j)*phi[i*Ny+j];
			}

			else if(i > interfacez1 && i < interfacez2 && j==0 ){
				Jaco(i,j,i,j) = 1.0;
				res[i*Ny+j] = phi[i*Ny+j] - thermal*log(Ndon[i*Ny+j]/ni) - 0.33374 - Vs;
			}
			else if(i > interfacez1 && i < interfacez2 && j==Ny-1){
				Jaco(i,j,i,j) = 1.0;
				res[i*Ny+j] = phi[i*Ny+j] - thermal*log(Ndon[i*Ny+j]/ni) - 0.33374 - Vd;
			}
			
			else {
				Jaco(i,j,i,j) = -2.0*eps*(delz/dely+dely/delz);
				Jaco(i,j,i,j-1) = eps*delz/dely;
				Jaco(i,j,i,j+1) = eps*delz/dely;
				Jaco(i,j,i-1,j) = eps*dely/delz;
				Jaco(i,j,i+1,j) = eps*dely/delz;
				res[i*Ny+j] = Jaco(i,j,i,j-1)*phi[i*Ny+j-1]+Jaco(i,j,i,j+1)*phi[i*Ny+j+1]+Jaco(i,j,i,j)*phi[i*Ny+j]+Jaco(i,j,i-1,j)*phi[(i-1)*Ny+j]+Jaco(i,j,i+1,j)*phi[(i+1)*Ny+j];

				res[i*Ny+j] -= q/ep0*dely*delz*(ni*exp(phi[i*Ny+j]/thermal)-Ndon[i*Ny+j]);
				Jaco(i,j,i,j) -= q/ep0*dely*delz*ni/thermal*exp(phi[i*Ny+j]/thermal);
			}

//			cout <<i << " " <<j <<" "<<  res[i*Ny+j] << endl;	

		}
	}


	for( int i = 0; i < Ny*Nz; i++) res[i] = -res[i];
	update = Jaco / res;
	for( int i = 0; i < Ny*Nz; i++) phi[i] += update[i];

	}
	Jaco.~LA3D();
	delete(update, res, Ndon);
	return phi;
	

}
