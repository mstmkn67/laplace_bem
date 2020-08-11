#include "LaplaceBEM_FMM2d.h"

LaplaceBEM_FMM2d::LaplaceBEM_FMM2d(
	FMMTree_c* t,Mesh_c* m,int mn,int ln,
	int ite,double tol,int ipn,const string& pre,const string& greenFunc,
	double phi,const complex<double>& gphi)
:tree(t),mesh(m),moments_number(mn),local_number(ln){
	potential_green_func_c::calc_factorial(local_number+moments_number);
	edge_number=mesh->edge.size();
	vector_b=new double[edge_number];
	vector_x=new double[edge_number];
	if(pre=="point_Jacobi"){
		bicgstab=new FMMJacobiBiCGSTAB(edge_number,tree,mesh,vector_x,vector_b,mn,ln,ite,tol,ipn,greenFunc,phi,gphi);
	}else if(pre=="block_solution_in_leaf"){
		bicgstab=new FMMBlockBiCGSTAB(edge_number,tree,mesh,vector_x,vector_b,mn,ln,ite,tol,ipn,greenFunc,phi,gphi);
	}else{
		bicgstab=new FMMNonBiCGSTAB(edge_number,tree,mesh,vector_x,vector_b,mn,ln,ite,tol,ipn,greenFunc,phi,gphi);
	}
}

LaplaceBEM_FMM2d::~LaplaceBEM_FMM2d(){
	delete bicgstab;
	delete vector_x;
	delete vector_b;
}

void LaplaceBEM_FMM2d::update(){
	bicgstab->update();
	vector<Edge_c*>::iterator i=mesh->edge.begin();
	for(int index=0;i!=mesh->edge.end();i++,index++){
		if((*i)->type=="phi"){
			((LaplaceEdge_c*)(*i))->q=vector_x[index];
		}else{
			((LaplaceEdge_c*)(*i))->phi=vector_x[index];
		}
	}
}

void LaplaceBEM_FMM2d::evaluate(
	const complex<double>& r,double& phi,Cell_c*& cell){
	cell=tree->get_cell_at_evaluation_position(r);
	phi=bicgstab->calc_evaluation_point(cell,r);
}

