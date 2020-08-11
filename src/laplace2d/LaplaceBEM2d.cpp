#include "LaplaceBEM2d.h"


LaplaceBEM2d::LaplaceBEM2d(Mesh_c* me,int m,double t,int ipn,const string& p,const string& g,
	double ph,const complex<double>& gp)
:mesh(me),phi(ph),gphi(gp){
	edge_number=mesh->edge.size();
	matrix_A=new double[edge_number*edge_number];
	vector_b=new double[edge_number];vector_x=new double[edge_number];
	integral= new ParameterOfLegendreGaussFormula1d(ipn);
	if(p=="LU"){
		cout << "The solver of LU_decomposition is used." << endl;
		bicgstab=0;
	}else if(p=="point_Jacobi"){
		bicgstab=new JacobiBiCGSTAB(edge_number,matrix_A,vector_x,vector_b,m,t);
	}else if(p=="ILU"){
		bicgstab=new ILUBiCGSTAB(edge_number,matrix_A,vector_x,vector_b,m,t);
	}else{
		bicgstab=new NonBiCGSTAB(edge_number,matrix_A,vector_x,vector_b,m,t);
	}
	if(g=="free_analytic_integral"){
		greenFunc=0;
	}else if(g=="free"){
		greenFunc=new potential_green_func_c::GreenFunc;
	}else if(g=="zero_potential_y0"){
		greenFunc=new potential_green_func_c::GreenFuncZeroPotentialOnWall;
	}else if(g=="no_flux_y0"){
		greenFunc=new potential_green_func_c::GreenFuncNoFluxOnWall;
	}
}

LaplaceBEM2d::~LaplaceBEM2d(){
	if(greenFunc!=0)delete greenFunc;
	if(bicgstab!=0)delete bicgstab;
	delete integral;delete[] vector_x;
	delete[] vector_b;delete[] matrix_A;
}

void LaplaceBEM2d::calc_bicgstab(){
	vector<Edge_c*>::iterator i=mesh->edge.begin();
	for(int a=0;i!=mesh->edge.end();i++,a++){
		complex<double>& ri=(*i)->position;
		complex<double>& ri0=(*i)->vertex[0]->position;
		complex<double>& ri1=(*i)->vertex[1]->position;
		string& type_i=(*i)->type;
		//vector_b
		if(type_i=="phi"){
			vector_b[a]=-0.5*((LaplaceEdge_c*)(*i))->phi +phi+ri.real()*gphi.real()+ri.imag()*gphi.imag();
		}else if(type_i=="q"){
			vector_b[a]=-((LaplaceEdge_c*)(*i))->q*calc_phi_Aaa(ri,ri0,ri1,(*i)->length)+phi+ri.real()*gphi.real()+ri.imag()*gphi.imag();
		}
		//matrix_A
		vector<Edge_c*>::iterator j=mesh->edge.begin();
		for(int b=0;j!=mesh->edge.end();j++,b++){
			complex<double>& rj0=(*j)->vertex[0]->position;
			complex<double>& rj1=(*j)->vertex[1]->position;
			complex<double>& nj=(*j)->normal;
			double length_j=(*j)->length;
			string& type_j=(*j)->type;
			if(a==b){
				if(type_j=="phi"){
					matrix_A[a*edge_number+b]=calc_phi_Aaa(ri,rj0,rj1,length_j);
				}else if(type_j=="q"){
					matrix_A[a*edge_number+b]=0.5;
				}
			}else{
				if(type_j=="phi"){
					matrix_A[a*edge_number+b]=calc_phi_Aab(ri,rj0,rj1,length_j);
					vector_b[a]+=((LaplaceEdge_c*)(*j))->phi*calc_q_Aab(ri,rj0,rj1,nj,length_j);
				}else if(type_j=="q"){
					matrix_A[a*edge_number+b]=-calc_q_Aab(ri,rj0,rj1,nj,length_j);
					vector_b[a]+=-((LaplaceEdge_c*)(*j))->q*calc_phi_Aab(ri,rj0,rj1,length_j);
				}
			}
		}
	}
}

void LaplaceBEM2d::calc_lapack(){
	vector<Edge_c*>::iterator i=mesh->edge.begin();
	for(int a=0;i!=mesh->edge.end();i++,a++){
		complex<double>& ri=(*i)->position;
		complex<double>& ri0=(*i)->vertex[0]->position;
		complex<double>& ri1=(*i)->vertex[1]->position;
		string& type_i=(*i)->type;
		//vector_b
		if(type_i=="phi"){
			vector_b[a]=-0.5*((LaplaceEdge_c*)(*i))->phi;
		}else if(type_i=="q"){
			vector_b[a]=-((LaplaceEdge_c*)(*i))->q*calc_phi_Aaa(ri,ri0,ri1,(*i)->length);
		}
		//matrix_A
		vector<Edge_c*>::iterator j=mesh->edge.begin();
		for(int b=0;j!=mesh->edge.end();j++,b++){
			complex<double>& rj0=(*j)->vertex[0]->position;
			complex<double>& rj1=(*j)->vertex[1]->position;
			complex<double>& nj=(*j)->normal;
			double length_j=(*j)->length;
			string& type_j=(*j)->type;
			if(a==b){
				if(type_j=="phi"){
					matrix_A[b*edge_number+a]=calc_phi_Aaa(ri,rj0,rj1,length_j);
				}else if(type_j=="q"){
					matrix_A[b*edge_number+a]=0.5;
				}
			}else{
				if(type_j=="phi"){
					matrix_A[b*edge_number+a]=calc_phi_Aab(ri,rj0,rj1,length_j);
					vector_b[a]+=((LaplaceEdge_c*)(*j))->phi*calc_q_Aab(ri,rj0,rj1,nj,length_j);
				}else if(type_j=="q"){
					matrix_A[b*edge_number+a]=-calc_q_Aab(ri,rj0,rj1,nj,length_j);
					vector_b[a]+=-((LaplaceEdge_c*)(*j))->q*calc_phi_Aab(ri,rj0,rj1,length_j);
				}
			}
		}
	}
}

void LaplaceBEM2d::update(){
	if(bicgstab==0){
		calc_lapack();
	}else{
		calc_bicgstab();
	}
	for(int i=0;i<edge_number;i++)vector_x[i]=0.0;
	if(bicgstab==0){
		lap::solve_linear_equations_general(edge_number,matrix_A,vector_b,vector_x);
	}else{
		int info=bicgstab->update();
	}
	for(int a=0;a<edge_number;a++){
		if(mesh->edge[a]->type=="phi"){
			((LaplaceEdge_c*)mesh->edge[a])->q=vector_x[a];
		}else{
			((LaplaceEdge_c*)mesh->edge[a])->phi=vector_x[a];
		}
	}
}

void LaplaceBEM2d::evaluate(const complex<double>& r,double& phi){
	phi=0.0;//grad_phi.real()=0.0;grad_phi.imag()=0.0;
	vector<Edge_c*>::iterator i=mesh->edge.begin();
	for(;i!=mesh->edge.end();i++){
		complex<double>& ri0=(*i)->vertex[0]->position;
		complex<double>& ri1=(*i)->vertex[1]->position;
		complex<double>& ni=(*i)->normal;
		double length_i=(*i)->length;
		phi+=-((LaplaceEdge_c*)(*i))->q*calc_phi_Aab(r,ri0,ri1,length_i)
		     +((LaplaceEdge_c*)(*i))->phi*calc_q_Aab(r,ri0,ri1,ni,length_i);
	}
	phi+=this->phi+r.real()*gphi.real()+r.imag()*gphi.imag();
}

double LaplaceBEM2d::calc_q_Aab(
	const complex<double>& r0,const complex<double>& ra,const complex<double>& rb,
	const complex<double>& n,double length)const{
	if(greenFunc==0){
		complex<double> dr1=ra-r0;
		complex<double> dr2=rb-r0;
		double dl1=sqrt(norm(dr1));
		double c1=dr1.real()/dl1;double s1=dr1.imag()/dl1;
		double dx2r=dr2.real()*c1+dr2.imag()*s1;
		double dy2r=-dr2.real()*s1+dr2.imag()*c1;
		double da=atan2(dy2r,dx2r);
		return da*0.5/3.14159265358979;
	}else{
		double s=0.0;
		for(int l=0;l<integral->point;l++){
			complex<double> r=0.5*((1.0-integral->xi[l])*ra+(integral->xi[l]+1.0)*rb);
			s+=integral->w[l]*greenFunc->get_Tn(r0,r,n);
		}
		return 0.5*length*s;
	}
}

double LaplaceBEM2d::calc_phi_Aab(
	const complex<double>& r0,const complex<double>& ra,const complex<double>& rb,double length)const{
	if(greenFunc==0){
		complex<double> dr1=ra-r0;
		complex<double> dr2=rb-r0;
		double dl1=sqrt(norm(dr1));double dl2=sqrt(norm(dr2));
		double c1=dr1.real()/dl1;double s1=dr1.imag()/dl1;
		double dx2r=dr2.real()*c1+dr2.imag()*s1;
		double dy2r=-dr2.real()*s1+dr2.imag()*c1;
		double theta=atan2(dy2r,dx2r);
		complex<double> tt=(rb-ra)/length;
		return 0.5/3.14159265358979*((-dr2.real()*tt.real()-dr2.imag()*tt.imag())*log(dl2)
			                 +( dr1.real()*tt.real()+dr1.imag()*tt.imag())*log(dl1)
			                 +length
			                 +(-dr1.real()*tt.imag()+dr1.imag()*tt.real())*theta);
	}else{
		double s=0.0;
		for(int l=0;l<integral->point;l++){
			complex<double> r=0.5*((1.0-integral->xi[l])*ra+(integral->xi[l]+1.0)*rb);
			s+=integral->w[l]*greenFunc->get_G(r0,r);
		}
		return 0.5*length*s;
	}
}

double LaplaceBEM2d::calc_phi_Aaa(const complex<double>& r0,const complex<double>& ra,
	                                const complex<double>& rb,double length)const{
	if(greenFunc==0){
		double l=0.5*length;
		return l*(-log(l)+1.0)/3.14159265358979;
	}else{
		double la=sqrt(norm(r0-ra));
		return calc_phi_Aab(r0,ra,r0,la)+calc_phi_Aab(r0,r0,rb,length-la);
	}
}

