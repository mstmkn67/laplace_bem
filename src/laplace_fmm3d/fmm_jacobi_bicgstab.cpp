#include "fmm_jacobi_bicgstab.h"


FMMJacobiBiCGSTAB::FMMJacobiBiCGSTAB(int n,FMMTree3d* cell,Mesh3d* mesh,double* x,double* b,
                    int moments_number,int local_number,int max_ite,double tol,int ipn,const string& g,
                    double p,const Vector3d& gp)
:FMMNonBiCGSTAB(n,cell,mesh,x,b,moments_number,local_number,max_ite,tol,ipn,g,p,gp){
}

FMMJacobiBiCGSTAB::~FMMJacobiBiCGSTAB(){}

void FMMJacobiBiCGSTAB::psolve(double* xx,double* bb){
	vector<Face3d*>::iterator fi=mesh->face.begin();
	for(int i=0;fi!=mesh->face.end();fi++,i++){
		Vector3d& ri=(*fi)->position;
		string& type_i=(*fi)->type;
		double s;
		if(type_i=="phi"){
			s=0.5;
		}else if(type_i=="q"){
			s=calc_Aaa((*fi)->position,(*fi)->vertex[0]->position,
			           (*fi)->vertex[1]->position,(*fi)->vertex[2]->position);
		}
		xx[i]=bb[i]/s;
	}
}
