#include "LaplaceBEM3d.h"


LaplaceBEM3d::LaplaceBEM3d(Mesh3d* me,int m,double t,int ipn,const string& p,const string& g,double ph,
	const Vector3d& gp)
:mesh(me),phi(ph),gphi(gp){
	face_number=mesh->face.size();
	matrix_A=new double[face_number*face_number];
	vector_b=new double[face_number];vector_x=new double[face_number];
	integral= new ParameterOfLegendreGaussFormulaTiangle(ipn);
	if(p=="LU"){
		cout << "The solver of LU_decomposition is used." << endl;
		bicgstab=0;
	}else if(p=="point_Jacobi"){
		bicgstab=new JacobiBiCGSTAB(face_number,matrix_A,vector_x,vector_b,m,t);
	}else if(p=="ILU"){
		bicgstab=new ILUBiCGSTAB(face_number,matrix_A,vector_x,vector_b,m,t);
	}else{
		bicgstab=new NonBiCGSTAB(face_number,matrix_A,vector_x,vector_b,m,t);
	}
	if(g=="zero_potential_z0"){
		greenFunc=new potential_green_func::GreenFuncZeroPotentialOnWall;
	}else if(g=="no_flux_z0"){
		greenFunc=new potential_green_func::GreenFuncNoFluxOnWall;
	}else{
		greenFunc=new potential_green_func::GreenFunc;
	}
}

LaplaceBEM3d::~LaplaceBEM3d(){
	delete greenFunc;
	if(bicgstab!=0)delete bicgstab;
	delete integral;delete[] vector_x;
	delete[] vector_b;delete[] matrix_A;
}

void LaplaceBEM3d::calc_bicgstab(){
	vector<Face3d*>::iterator i=mesh->face.begin();
	for(int a=0;i!=mesh->face.end();i++,a++){
		Vector3d& ri=(*i)->position;
		Vector3d& ri0=(*i)->vertex[0]->position;
		Vector3d& ri1=(*i)->vertex[1]->position;
		Vector3d& ri2=(*i)->vertex[2]->position;
		string& type_i=(*i)->type;
		//vector_b
		if(type_i=="phi"){
			vector_b[a]=-0.5*((LaplaceFace3d*)(*i))->phi +phi+gphi*ri;
		}else if(type_i=="q"){
			vector_b[a]=-((LaplaceFace3d*)(*i))->q*calc_Aaa(ri,ri0,ri1,ri2)+phi+gphi*ri;
		}
		//matrix_A
		vector<Face3d*>::iterator j=mesh->face.begin();
		for(int b=0;j!=mesh->face.end();j++,b++){
			Vector3d& rj0=(*j)->vertex[0]->position;
			Vector3d& rj1=(*j)->vertex[1]->position;
			Vector3d& rj2=(*j)->vertex[2]->position;
			Vector3d& nj=(*j)->normal;
			double area_j=(*j)->area;
			string& type_j=(*j)->type;
			if(a==b){
				if(type_j=="phi"){
					matrix_A[a*face_number+b]=calc_Aaa(ri,ri0,ri1,ri2);
				}else if(type_j=="q"){
					matrix_A[a*face_number+b]=0.5;
				}
			}else{
				if(type_j=="phi"){
					matrix_A[a*face_number+b]=calc_Aab(ri,rj0,rj1,rj2,area_j);
					vector_b[a]+=((LaplaceFace3d*)(*j))->phi*calc_Bab(ri,rj0,rj1,rj2,nj,area_j);
				}else if(type_j=="q"){
					matrix_A[a*face_number+b]=-calc_Bab(ri,rj0,rj1,rj2,nj,area_j);
					vector_b[a]+=-((LaplaceFace3d*)(*j))->q*calc_Aab(ri,rj0,rj1,rj2,area_j);
				}
			}
		}
	}
}

void LaplaceBEM3d::calc_lapack(){
	vector<Face3d*>::iterator i=mesh->face.begin();
	for(int a=0;i!=mesh->face.end();i++,a++){
		Vector3d& ri=(*i)->position;
		Vector3d& ri0=(*i)->vertex[0]->position;
		Vector3d& ri1=(*i)->vertex[1]->position;
		Vector3d& ri2=(*i)->vertex[2]->position;
		string& type_i=(*i)->type;
		//vector_b
		if(type_i=="phi"){
			vector_b[a]=-0.5*((LaplaceFace3d*)(*i))->phi+ phi+gphi*ri;
		}else if(type_i=="q"){
			vector_b[a]=-((LaplaceFace3d*)(*i))->q*calc_Aaa(ri,ri0,ri1,ri2)+ phi+gphi*ri;
		}
		//matrix_A
		vector<Face3d*>::iterator j=mesh->face.begin();
		for(int b=0;j!=mesh->face.end();j++,b++){
			Vector3d& rj0=(*j)->vertex[0]->position;
			Vector3d& rj1=(*j)->vertex[1]->position;
			Vector3d& rj2=(*j)->vertex[2]->position;
			Vector3d& nj=(*j)->normal;
			double area_j=(*j)->area;
			string& type_j=(*j)->type;
			if(a==b){
				if(type_j=="phi"){
					matrix_A[b*face_number+a]=calc_Aaa(ri,ri0,ri1,ri2);
				}else if(type_j=="q"){
					matrix_A[b*face_number+a]=0.5;
				}
			}else{
				if(type_j=="phi"){
					matrix_A[b*face_number+a]=calc_Aab(ri,rj0,rj1,rj2,area_j);
					vector_b[a]+=((LaplaceFace3d*)(*j))->phi*calc_Bab(ri,rj0,rj1,rj2,nj,area_j);
				}else if(type_j=="q"){
					matrix_A[b*face_number+a]=-calc_Bab(ri,rj0,rj1,rj2,nj,area_j);
					vector_b[a]+=-((LaplaceFace3d*)(*j))->q*calc_Aab(ri,rj0,rj1,rj2,area_j);
				}
			}
		}
	}
}

void LaplaceBEM3d::update(){
	if(bicgstab==0){
		calc_lapack();
	}else{
		calc_bicgstab();
	}
	for(int i=0;i<face_number;i++)vector_x[i]=0.0;
	if(bicgstab==0){
		lap::solve_linear_equations_general(face_number,matrix_A,vector_b,vector_x);
	}else{
		int info=bicgstab->update();
	}
	for(int a=0;a<face_number;a++){
		if(mesh->face[a]->type=="phi"){
			((LaplaceFace3d*)(mesh->face[a]))->q=vector_x[a];
		}else{
			((LaplaceFace3d*)(mesh->face[a]))->phi=vector_x[a];
		}
	}
}

void LaplaceBEM3d::evaluate(const Vector3d& r,double& phi){
	phi=0.0;//grad_phi.real()=0.0;grad_phi.imag()=0.0;
	vector<Face3d*>::iterator i=mesh->face.begin();
	for(;i!=mesh->face.end();i++){
		const Vector3d& ri0=(*i)->vertex[0]->position;
		const Vector3d& ri1=(*i)->vertex[1]->position;
		const Vector3d& ri2=(*i)->vertex[2]->position;
		const Vector3d& ni=(*i)->normal;
		double area=(*i)->area;
		phi+=-((LaplaceFace3d*)(*i))->q*calc_Aab(r,ri0,ri1,ri2,area)
		     +((LaplaceFace3d*)(*i))->phi*calc_Bab(r,ri0,ri1,ri2,ni,area);
	}
	phi+=this->phi+gphi*r;
}

double LaplaceBEM3d::calc_Aab(
	const Vector3d& r0,const Vector3d& ra,const Vector3d& rb,const Vector3d& rc,double area)const{
	double s=0.0;
	for(int l=0;l<integral->point;l++){
		Vector3d r=integral->xi[l].x*(rb-ra)+integral->xi[l].y*(rc-ra)+ra;
		s+=integral->w[l]*greenFunc->get_G(r0,r);
		//s+=integral->w[l]*potential_green_func::get_free_G(r0,r);
	}
	return area*s;
}

double LaplaceBEM3d::calc_Aaa(const Vector3d& r0,const Vector3d& ra,const Vector3d& rb,const Vector3d& rc)const{
	double s0=0.5*((ra-r0)^(rb-r0)).length();double s1=0.5*((rb-r0)^(rc-r0)).length();
	double s2=0.5*((rc-r0)^(ra-r0)).length();
	double a=calc_Aab(r0,r0,ra,rb,s0);
	a+=calc_Aab(r0,r0,rb,rc,s1);a+=calc_Aab(r0,r0,rc,ra,s2);
	return a;
}

double LaplaceBEM3d::calc_Bab(
	const Vector3d& r0,const Vector3d& ra,const Vector3d& rb,const Vector3d& rc,const Vector3d& n,double area)const{
	double s=0.0;
	for(int l=0;l<integral->point;l++){
		Vector3d r=integral->xi[l].x*(rb-ra)+integral->xi[l].y*(rc-ra)+ra;
		s+=integral->w[l]*greenFunc->get_Tn(r0,r,n);
		//s+=integral->w[l]*n*potential_green_func::get_free_T(r,r0);
	}
	return area*s;
}

