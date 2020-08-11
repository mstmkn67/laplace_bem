#include "LaplaceBEM_FMM3d.h"

LaplaceBEM_FMM3d::LaplaceBEM_FMM3d(
	FMMTree3d* t,Mesh3d* m,int mn,int ln,
	int ite,double tol,int ipn,const string& pre,const string& green,double p,const Vector3d& gp)
:tree(t),mesh(m),moments_number(mn),local_number(ln){
	face_number=mesh->face.size();
	vector_b=new double[face_number];
	vector_x=new double[face_number];
	if(pre=="point_Jacobi"){
		bicgstab=new FMMJacobiBiCGSTAB(face_number,tree,mesh,vector_x,vector_b,mn,ln,ite,tol,ipn,green,p,gp);
	}else if(pre=="block_solution_in_leaf"){
		bicgstab=new FMMBlockBiCGSTAB(face_number,tree,mesh,vector_x,vector_b,mn,ln,ite,tol,ipn,green,p,gp);
	}else{
		bicgstab=new FMMNonBiCGSTAB(face_number,tree,mesh,vector_x,vector_b,mn,ln,ite,tol,ipn,green,p,gp);
	}
}

LaplaceBEM_FMM3d::~LaplaceBEM_FMM3d(){
	delete bicgstab;
	delete vector_x;
	delete vector_b;
}

void LaplaceBEM_FMM3d::update(){
	bicgstab->update();
	vector<Face3d*>::iterator i=mesh->face.begin();
	for(int index=0;i!=mesh->face.end();i++,index++){
		if((*i)->type=="phi"){
			((LaplaceFace3d*)(*i))->q=vector_x[index];
		}else{
			((LaplaceFace3d*)(*i))->phi=vector_x[index];
		}
	}
}

void LaplaceBEM_FMM3d::evaluate(const Vector3d& r,double& phi,Cell3d*& cell){
	cell=tree->get_cell_at_evaluation_position(r);
	phi=bicgstab->calc_evaluation_point(cell,r);
}

