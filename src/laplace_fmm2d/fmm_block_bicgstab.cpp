#include "fmm_block_bicgstab.h"


FMMBlockBiCGSTAB::FMMBlockBiCGSTAB(int n,FMMTree_c* tree,Mesh_c* mesh,double* x,double* b,
                    int moments_number,int local_number,int max_ite,double tol,int ipn,const string& greenFunc,
                    double phi,const complex<double>& gphi)
:FMMNonBiCGSTAB(n,tree,mesh,x,b,moments_number,local_number,max_ite,tol,ipn,greenFunc,phi,gphi){
}

FMMBlockBiCGSTAB::~FMMBlockBiCGSTAB(){}

void FMMBlockBiCGSTAB::psolve(double* xx,double* bb){
	vector<Cell_c*>::iterator c=tree->leaf.begin();
	for(;c!=tree->leaf.end();c++){
		int sub_n=(*c)->edge.size();
		double* sub_b=new double[sub_n];
		double* sub_A=new double[sub_n*sub_n];
		vector<Edge_c*>::iterator ei=(*c)->edge.begin();
		for(int i=0;ei!=(*c)->edge.end();ei++,i++){
			sub_b[i]=bb[(*ei)->id];
			complex<double>& ri=(*ei)->position;
			complex<double>& ri0=(*ei)->vertex[0]->position;
			complex<double>& ri1=(*ei)->vertex[1]->position;
			vector<Edge_c*>::iterator ej=(*c)->edge.begin();
			for(int j=0;ej!=(*c)->edge.end();ej++,j++){
				complex<double>& rj0=(*ej)->vertex[0]->position;
				complex<double>& rj1=(*ej)->vertex[1]->position;
				complex<double>& nj=(*ej)->normal;
				double length_j=(*ej)->length;
				if((*ei)!=(*ej)){
					if((*ej)->type=="phi"){
						sub_A[j*sub_n+i]= calc_Aab(ri,rj0,rj1,length_j);
					}else if((*ej)->type=="q"){
						sub_A[j*sub_n+i]=-calc_Bab(ri,rj0,rj1,nj,length_j);
					}
				}else{
					if((*ei)->type=="phi"){
						sub_A[i*sub_n+i]=calc_Aaa(ri,ri0,ri1,(*ei)->length);
					}else if((*ei)->type=="q"){
						sub_A[i*sub_n+i]=0.5;
					}
				}
			}
		}
		lap::solve_linear_equations_general(sub_n,sub_A,sub_b);
		ei=(*c)->edge.begin();
		for(int i=0;ei!=(*c)->edge.end();ei++,i++){
			xx[(*ei)->id]=sub_b[i];
		}
		delete[] sub_b;delete[] sub_A;
	}
}
