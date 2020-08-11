#include "fmm_jacobi_bicgstab.h"


FMMJacobiBiCGSTAB::FMMJacobiBiCGSTAB(int n,FMMTree_c* cell,Mesh_c* mesh,double* x,double* b,
                    int moments_number,int local_number,int max_ite,double tol,int ipn,const string& greenFunc,
                    double phi,const complex<double>& gphi)
:FMMNonBiCGSTAB(n,cell,mesh,x,b,moments_number,local_number,max_ite,tol,ipn,greenFunc,phi,gphi){
}

FMMJacobiBiCGSTAB::~FMMJacobiBiCGSTAB(){}

void FMMJacobiBiCGSTAB::psolve(double* xx,double* bb){
	vector<Edge_c*>::iterator ei=mesh->edge.begin();
	for(int i=0;ei!=mesh->edge.end();ei++,i++){
		complex<double>& ri=(*ei)->position;
		complex<double>& ri0=(*ei)->vertex[0]->position;
		complex<double>& ri1=(*ei)->vertex[1]->position;
		string& type_i=(*ei)->type;
		double s;
		if(type_i=="phi"){
			//s=-0.5*(*ei)->phi;
			s=0.5;
		}else if(type_i=="q"){
			//s=-(*ei)->q*calc_Aaa(ri,ri0,ri1,(*ei)->length);
			s=calc_Aaa(ri,ri0,ri1,(*ei)->length);
		}
		xx[i]=bb[i]/s;
	}
}
