#ifndef _LAPACK_FUNCTIONS_H_
#define _LAPACK_FUNCTIONS_H_

#include <iostream>
using namespace std;

////////////original subroutines in LAPACK////////////////////////
extern"C"{
	//int ilaenv_(int* ISPEC,char* NAME,char* OPTS,int* N1,int* N2,int* N3,int* N4);
	
	//void dgbsv_(int* N,int* KL,int* KU,int* NRHS,double* AB,int* LDAB,int* IPIV,double* B,int* INFO);
	void dgesv_(int* N,int* NRHS,double* A,int* LDA,int* IPIV,double* B,int* LDB,int* INFO);
	//	
	void dsysv_(char* UPLO,int* N,int* NRHS,double* A,int* LDA,int* IPIV,
	            double* B,int* LDB,double* WORK,int* LWORK,int* INFO);

	void dgetrf_(int* M,int* N,double* A,int* LDA,int* IPIV,int* INFO);
	void dgetri_(int* N,double* A,int* LDA,int* IPIV,double* WORK,int* LWORK,int* INFO);

	void dpptrf_(char* UPLO,int* N,double* AP,int* INFO);
	void dpptri_(char* UPLO,int* N,double* AP,int* INFO);

	void dsyev_(char* JOBS,char* UPLO,int* N,double* A,int* LDA,double* W,double* WORK,int* LWORK,int* INFO);
}

namespace lap{
//////////////////////////////////////////////////////////////////
///////////interface for c++ ////////////////////////////////
//void solve_linear_equations_band(int* n,int* kl,int* ku,);
// strage method of matrix is fortran type 
int solve_linear_equations_general(int n,double* A,double* B);
int solve_linear_equations_general(int n,double* A,double* B,double* X);//solutions are not written in B

int solve_linear_equations_symmetric(int n,double* A,double* B);
int solve_linear_equations_symmetric(int n,double* A,double* B,double* X);

int get_inverse_matrix_symmetric(int n,double* A);//Cholesky decomposition
int get_inverse_matrix_general(int n,double* A);//LU decomposition

int get_Cholesky_decomposition(int n,double* A);

int get_eigenvalues_symmetric(int n,double* A,double* W);
int get_eigenvalues_and_vectors_symmetric(int n,double* A,double* W);


//general matrix index
// i,j (=0...n-1) index of n x n matrix 
inline int index(int i,int j,int n){
	return n*i+j;
}
inline int left(int index){
	return index-1;
}
inline int right(int index){
	return index+1;
}
inline int up(int index,int n){
	return index-n;
}
inline int down(int index,int n){
	return index+n;
}

//a,b (=0,1,2) index of  ei,ej sub matrix of 3n x 3n matrix
inline int mindex(int ei,int ej,int a,int b,int n){
	return 3*n*(3*ei+a)+3*ej+b;
}
inline int mleft(int index){
	return index-1;
}
inline int mright(int index){
	return index+1;
}
inline int mup(int index,int n){
	return index-3*n;
}
inline int mdown(int index,int n){
	return index+3*n;
}


//symmetric compressed matrix index
// upper triangle matrix in C or C++ (lower one in Fortran) when you use symmetric case
// i,j (=0...n-1) index of n x n matrix 
inline int sindex(int i,int j,int n){ 
	return i>j ? (j*(2*n-j-1)/2+i):(i*(2*n-i-1)/2+j);
}
inline int sleft(int index){
	return index-1;
}
inline int sright(int index){
	return index+1;
}
inline int sup(int index,int i,int n){
	return index-n+i;
}
inline int sdown(int index,int i,int n){
	return index+n-i-1;
}
//a,b (=0,1,2) index of  ei,ej sub matrix of 3n x 3n matrix
inline int msindex(int ei,int ej,int a,int b,int n){
	int i=3*ei+a,j=3*ej+b;
	return i>j ? (j*(6*n-j-1)/2+i):(i*(6*n-i-1)/2+j);
}
inline int msleft(int index){
	return index-1;
}
inline int msright(int index){
	return index+1;
}
inline int msup(int index,int i,int a,int n){
	return index-3*n+3*i+a;
}
inline int msdown(int index,int i,int a,int n){
	return index+3*n-3*i-a-1;
}

//index classes 
/*
class Indexl{
public:
	Indexl(int n_=0):n(n_){index=0;};
	Indexl(int i,int j,int n_):n(n_){index=n*i+j;};
	~Indexl(){};
	inline int seek(int i,int j){return index=n*i+j;};
	inline int reset(int n_){n=n_;return index=0;};
	
	inline int left(){return index-1;};
	inline int right(){return index+1;};
	inline int up(){return index-n;};
	inline int down(){return index+n;};
	
	inline int move_left(){return index-=1;};
	inline int move_right(){return index+=1;};
	inline int move_up(){return index-=n;};
	inline int move_down(){return index+=n;};
	
	inline int operator()(){return index;};
	inline int operator++(){return index+=1;};
	inline int operator--(){return index-=1;};
	
	int n;//dimension of matrix
	int index;//current index
};

class Mindexl{
public:
	Mindexl(int n_=0):n(n_){index=0;};
	Mindexl(int ei,int a,int ej,int b,int n_):n(n_){index=3*n*(3*ei+a)+3*ej+b;};
	~Mindexl(){};
	inline int seek(int ei,int ej,int a,int b){return index=3*n*(3*ei+a)+3*ej+b;};
	inline int reset(int n_){n=n_;return index=0;};
	
	inline int left(){return index-1;};
	inline int right(){return index+1;};
	inline int up(){return index-3*n;};
	inline int down(){return index+3*n;};
	
	inline int move_left(){return index-=1;};
	inline int move_right(){return index+=1;};
	inline int move_up(){return index-=3*n;};
	inline int move_down(){return index+=3*n;};
	
	inline int sub_begin(int ei,int ej){return index=3*n*(3*ei)+3*ej;};
	inline int sub_get(int a,int b){return index+3*n*a+b;}
	
	inline int operator()(){return index;};
	inline int operator++(){return index+=1;};
	inline int operator--(){return index-=1;};
	
	int n;//dimension of matrix
	int index;//current index
};

class Sindexl{
public:
	Sindexl(int n_=0):n(n_){index=0;};
	Sindexl(int i,int j,int n_):n(n_){index=i>j ? (j*(2*n-j-1)/2+i):(i*(2*n-i-1)/2+j);};
	~Sindexl(){};
	inline int seek(int i,int j){return index=i>j ? (j*(2*n-j-1)/2+i):(i*(2*n-i-1)/2+j);};
	inline int reset(int n_){n=n_;return index=0;};
	
	inline int left(){return index-1;};
	inline int right(){return index+1;};
	inline int up(){return index-n+1;};
	inline int down(){return index+n-i-1;};
	
	inline int move_left(){return index+=-1;};
	inline int move_right(){return index+=+1;};
	inline int move_up(){return index+=-n+1;};
	inline int move_down(){return index+=n-i-1;};
	
	inline int operator()(){return index};
	inline int operator++(){return index+=1;};
	inline int operator--(){return index-=1;};
	
	int n;//dimension of matrix
	int index;//current index
};
class Msindexl{
public:
	Msindexl(int n_=0):n(n_){index=0;};
	Msindexl(int ei,int a,int ej,int b,int n_):n(n_){
		int i=3*ei+a,j=3*ej+b;
		index=i>j ? (j*(6*n-j-1)/2+i):(i*(6*n-i-1)/2+j);
	}
	~Msindexl(){};
	inline int seek(int ei,int ej,int a,int b){
		int i=3*ei+a,j=3*ej+b;
		return index=i>j ? (j*(6*n-j-1)/2+i):(i*(6*n-i-1)/2+j);
	};
	inline int reset(int n_){n=n_;index=0;};
	
	inline int left(){return index-1;};
	inline int right(){return index+1;};
	inline int up(){return index-3*n+3*i+a;};
	inline int down(){return index+3*n-3*i-a-1;};
	
	inline int move_left(){return index+=-1;};
	inline int move_right(){return index+=+1;};
	inline int move_up(){return index+=-3*n+3*i+a;};
	inline int move_down(){return index+=3*n-3*i-a-1;};

	inline int sub_begin(int ei,int ej){//ei <= ej
		int i=3*ei,j=3*ej;
		return index=(i*(6*n-i-1)/2+j);
	};
	inline int sub_get(int a,int b){
		
	};
	
	inline int operator()(){return index;};
	inline int operator++()return index+=1;};
	inline int operator--(){return index-=1;};
	
	int n;//dimension of matrix
	int index;//current index
};*/

}
/////////////////////////////////////////

#endif // _LAPACK_FUNCTIONS_H_
