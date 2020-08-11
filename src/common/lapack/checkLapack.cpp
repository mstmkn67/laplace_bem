#include "LapackFunctions.h"
#include <iostream>
using namespace std;


int main(){
	cout << "check of solve_liner_equations_general" << endl;
	int n=3;
	double a[]={1.0,2.5,3.5,
	            2.0,2.0,2.0,
	            3.0,2.0,4.0};//transpose 
	double b[]={1.0,2.0,3.0};
	int info=lap::solve_linear_equations_general(n,a,b);
	cout << "info: "  << info << endl;
	cout << "solution: " <<  b[0] << " " << b[1] << " " << b[2] << endl;
	cout << endl << endl;
	/////////////////////////////////////////////////////////////////
	cout << "check of get_inverse_matrix_symmeteric" << endl;
	n=3;
	double c[]={1.0,1.0,1.0,
		              2.0,1.0,
		                  4.0};
	double d[]={1.0,1.0,1.0,
		              2.0,1.0,
		                  4.0};
	info=lap::get_inverse_matrix_symmetric(n,d);
	cout << "info: "  << info << endl;
	cout << "matrix times inverse matirx : " << endl;
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			double s=0.0;
			for(int k=0;k<n;k++){
				s+=c[lap::sindex(i,k,n)]*d[lap::sindex(k,j,n)];
			}
			cout << s << " ";
		}
		cout << endl;
	}
	cout << endl << endl;
	////////////////////////////////////////////////////////////////////
	cout << "check of get_inverse_matrix_general " << endl;
	n=3;
	double e[]={3.0,2.0,3.0,
	            2.5,3.0,2.0,
	            3.5,2.0,4.0};//transpose
	double f[]={3.0,2.0,3.0,
	            2.5,3.0,2.0,
	            3.5,2.0,4.0};//transpose
	info=lap::get_inverse_matrix_general(n,e);
	cout << "info: "  << info << endl;
	cout << "matrix times inverse matirx : " << endl;
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			double s=0.0;
			for(int k=0;k<n;k++){
				s+=e[lap::index(i,k,n)]*f[lap::index(k,j,n)];
			}
			cout << s << " ";
		}
		cout << endl;
	}
	cout << endl << endl;
	////////////////////////////////////////////////////////////////////
	cout << "check of get_Cholesky_decomposition " << endl;
	n=3;
	double g[]={1.0,1.0,1.0,
		              2.0,1.0,
		                  4.0};
	double h[]={1.0,1.0,1.0,
		              2.0,1.0,
		                  4.0};
	info=lap::get_Cholesky_decomposition(n,g);
	cout << "info: "  << info << endl;
	cout << "decomposed matrix times transposed matrix- original matrix: " << endl;
	for(int i=0;i<n;i++){
		for(int j=i;j<n;j++){
			double s=0.0;
			for(int k=0;k<n;k++){
				if(i>=k && j>=k){
					s+=g[lap::sindex(i,k,n)]*g[lap::sindex(j,k,n)];
				}
			}
			cout << s-h[lap::sindex(i,j,n)] << " ";
		}
		cout << endl;
	}
	cout << endl << endl;
	////////////////////////////////////////////////////////////////////
	cout << "get_eigenvalues_and_vectors" << endl;
	double w[n];
	double o[]={1.0,1.0,1.0,
		          9.9,2.0,1.0,
		          9.9,9.9,4.0};//elements of 9.9 are not used
	double p[]={1.0,1.0,1.0,
		          1.0,2.0,1.0,
		          1.0,1.0,4.0};
	info=lap::get_eigenvalues_and_vectors_symmetric(n,o,w);
	cout << "info: "  << info << endl;
	cout << "eigen value 0: " << w[0]
	     << " eigen vector 0: " << o[0] << " " << o[1] << " " << o[2] << endl;
	cout << "eigen value 1: " << w[1]
	     << " eigen vector 1: " << o[3] << " " << o[4] << " " << o[5] << endl;
	cout << "eigen value 2: " << w[2]
	     << " eigen vector 2: " << o[6] << " " << o[7] << " " << o[8] << endl;
	cout << "matrix reproduced from spectrum decomposition - original matrix: " << endl;
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			double s=0.0;
			for(int k=0;k<n;k++){ 
				s+=w[k]*o[n*k+i]*o[n*k+j];
			}
			cout << s-p[lap::index(i,j,n)] << " ";
		}
		cout << endl;
	}
	cout << endl << endl;
	////////////////////////////////////////////////////////////////////
	cout << "get_eigenvalues" << endl;
	double q[n];
	double r[]={1.0,1.0,1.0,
		          9.9,2.0,1.0,
		          9.9,9.9,4.0};//elements of 9.9 are not used
	info=lap::get_eigenvalues_symmetric(n,r,q);
	cout << "info: "  << info << endl;
	cout << "eigen values: " << q[0] << " " << q[1] << " " << q[2] << endl;
	cout << endl << endl;
	////////////////////////////////////////////////////////////////////
	cout << "check of solve_linear_equations_symmetric " << endl;
	n=3;
	double s[]={1.0,1.0,1.0,
		          1.0,2.0,1.0,
		          1.0,1.0,4.0};
	double t[]={1.0,2.0,3.0};
	info=lap::solve_linear_equations_symmetric(n,s,t);
	cout << "info: "  << info << endl;
	cout << "solution: " <<  t[0] << " " << t[1] << " " << t[2] << endl;
	cout << endl << endl;
	return 0;
}
