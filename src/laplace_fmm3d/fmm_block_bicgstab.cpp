#include "fmm_block_bicgstab.h"


FMMBlockBiCGSTAB::FMMBlockBiCGSTAB(int n,FMMTree3d* tree,Mesh3d* mesh,double* x,double* b,
                    int moments_number,int local_number,int max_ite,double tol,int ipn,const string& g,
                    double p,const Vector3d& gp)
:FMMNonBiCGSTAB(n,tree,mesh,x,b,moments_number,local_number,max_ite,tol,ipn,g,p,gp){
}

FMMBlockBiCGSTAB::~FMMBlockBiCGSTAB(){}

void FMMBlockBiCGSTAB::psolve(double* xx,double* bb){
	vector<Cell3d*>::iterator c=tree->leaf.begin();
	for(;c!=tree->leaf.end();c++){
		int sub_n=(*c)->face.size();
		double* sub_b=new double[sub_n];
		double* sub_A=new double[sub_n*sub_n];
		vector<Face3d*>::iterator fi=(*c)->face.begin();
		for(int i=0;fi!=(*c)->face.end();fi++,i++){
			sub_b[i]=bb[(*fi)->id];
			Vector3d& ri=(*fi)->position;
			vector<Face3d*>::iterator fj=(*c)->face.begin();
			for(int j=0;fj!=(*c)->face.end();fj++,j++){
				Vector3d& rj0=(*fj)->vertex[0]->position;
				Vector3d& rj1=(*fj)->vertex[1]->position;
				Vector3d& rj2=(*fj)->vertex[2]->position;
				Vector3d& nj=(*fj)->normal;
				double area_j=(*fj)->area;
				if((*fi)!=(*fj)){
					if((*fj)->type=="phi"){
						sub_A[j*sub_n+i]= calc_Aab(ri,rj0,rj1,rj2,area_j);
					}else if((*fj)->type=="q"){
						sub_A[j*sub_n+i]=-calc_Bab(ri,rj0,rj1,rj2,nj,area_j);
					}
				}else{
					if((*fi)->type=="phi"){
						sub_A[i*sub_n+i]=calc_Aaa(ri,(*fi)->vertex[0]->position,
						                             (*fi)->vertex[1]->position,
						                             (*fi)->vertex[2]->position);
					}else if((*fi)->type=="q"){
						sub_A[i*sub_n+i]=0.5;
					}
				}
			}
		}
		lap::solve_linear_equations_general(sub_n,sub_A,sub_b);
		fi=(*c)->face.begin();
		for(int i=0;fi!=(*c)->face.end();fi++,i++){
			xx[(*fi)->id]=sub_b[i];
		}
		delete[] sub_b;delete[] sub_A;
	}
}
