#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

#include "trsv.h"

using namespace std;

double errCalc(int n, double *exact, double *calc);

#define PREC 1e-5
//#define RAND_MAX 1000//

int main() {

	//System size:
	int n=10;
	
	//Build A and x, then compute b:
	double *A_U = new double[n*n];
	double *A_L = new double[n*n];
	double *A_U_D = new double[n*n];	//=A_U with 1s on diagonal
	double *A_L_D = new double[n*n];	//=A_L with 1s on diagonal
	double *x = new double[n];
	double *b_U = new double[n];
	double *b_L = new double[n];
	double *b_U_T = new double[n];		//=b_L
	double *b_L_T = new double[n];		//=b_U
	double *b_U_D = new double[n];
	double *b_L_D = new double[n];
	double *b_U_T_D = new double[n];
	double *b_L_T_D = new double[n];
	
	//Fill A, x with random numbers:
	srand(5);
    for (int i=0; i<n; ++i) {
        for (int j=i; j<n; ++j) {
        
        	A_U[j*n +i] = 0;
        	A_L[i*n +j] = 0;
        	
        	double tmp=rand()/(double)RAND_MAX;
            A_U[i*n +j] = tmp;
            A_L[j*n +i] = tmp;
            
            A_U_D[i*n +j] = A_U[i*n +j];
            A_L_D[j*n +i] = A_L[j*n +i];
            
        }
        
        A_U_D[i*n +i] = 1;
        A_L_D[i*n +i] = 1;
        
        x[i] = rand()/(double)RAND_MAX;
        
    }
    
    //b = Ax
    //cout << "b:" << endl;
    for (int i=0; i<n; ++i) {
    
    	b_U[i] = 0;
    	b_L[i] = 0;
    	
    	b_U_D[i] = 0;
    	b_L_D[i] = 0;
    	
    	for (int j=0; j<n; ++j) {
    	
    		b_U[i] += A_U[i*n +j] * x[j];
    		b_L[i] += A_L[i*n +j] * x[j];
    		
    		b_U_D[i] += A_U_D[i*n +j] * x[j];
    		b_L_D[i] += A_L_D[i*n +j] * x[j];
    		
    	}
    	
    	b_U_T[i] = b_L[i];
    	b_L_T[i] = b_U[i];
    	
    	b_U_T_D[i] = b_L_D[i];
    	b_L_T_D[i] = b_U_D[i];
    	
    	//cout << b_U_T[i] << endl;
    	
    }
	
	//Solve using the solver:
	int lda = n;
	int incX = 1;

	//Upper:
	trsv('U', 'N', 'N', n, A_U, lda, b_U, incX);
	trsv('U', 'T', 'N', n, A_U, lda, b_U_T, incX);
	trsv('U', 'N', 'U', n, A_U_D, lda, b_U_D, incX);
	trsv('U', 'T', 'U', n, A_U_D, lda, b_U_T_D, incX);

	//Lower:
	trsv('L', 'N', 'N', n, A_L, lda, b_L, incX);
	trsv('L', 'T', 'N', n, A_L, lda, b_L_T, incX);
	trsv('L', 'N', 'U', n, A_L_D, lda, b_L_D, incX);
	trsv('L', 'T', 'U', n, A_L_D, lda, b_L_T_D, incX);

	//Output results:
	cout << "\nResults:\n" << endl;
	
	double err_U = errCalc(n,x,b_U);
	double err_U_T = errCalc(n,x,b_U_T);
	double err_U_D = errCalc(n,x,b_U_D);
	double err_U_T_D = errCalc(n,x,b_U_T_D);
	
	double err_L = errCalc(n,x,b_L);
	double err_L_T = errCalc(n,x,b_L_T);
	double err_L_D = errCalc(n,x,b_L_D);
	double err_L_T_D = errCalc(n,x,b_L_T_D);
	
	cout << fixed << setprecision(4) << "Upper: " << err_U << "\t\t\t ----> " << ((err_U<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Upper transpose: " << err_U_T << "\t\t ----> " << ((err_U_T<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Upper unit: " << err_U_D << "\t\t ----> " << ((err_U_D<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Upper unit transpose: " << err_U_T_D << "\t ----> " << ((err_U_T_D<PREC)?"PASS":"FAIL") << endl;
		
	cout << fixed << setprecision(4) << "Lower: " << err_L << "\t\t\t ----> " << ((err_L<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Lower transpose: " << err_L_T << "\t\t ----> " << ((err_L_T<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Lower unit: " << err_L_D << "\t\t ----> " << ((err_L_D<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Lower unit transpose: " << err_L_T_D << "\t ----> " << ((err_L_T_D<PREC)?"PASS":"FAIL") << endl;

	return 0;
	
}

double errCalc(int n, double *exact, double *calc) {

	double error = 0;
	
	for (int i=0; i<n; ++i) {
	
		error += fabs(exact[i] - calc[i]);
		
	}
	
	return error;
		
}