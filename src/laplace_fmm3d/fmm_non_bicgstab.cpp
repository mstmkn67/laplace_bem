#include "fmm_non_bicgstab.h"

#define PI 3.14159265358979

FMMNonBiCGSTAB::FMMNonBiCGSTAB(int n,FMMTree3d* t,Mesh3d* m,double* vx,double* vb,
	                             int mn,int ln,int max_ite,double tol,int ipn,const string& g,
                               double p,const Vector3d& gp)
:BiCGSTAB(n,vx,vb,max_ite,tol),tree(t),mesh(m),moments_number(mn),local_number(ln),phi(p),gphi(gp){
	integral= new ParameterOfLegendreGaussFormulaTiangle(ipn);
	//for(int i=0;i<n;i++)x[i]=0.0;
	//input the initial guess about x 
	vector<Face3d*>::iterator i=mesh->face.begin();
	for(int index=0;i!=mesh->face.end();i++,index++){
		if((*i)->type=="phi"){
			x[index]=((LaplaceFace3d*)(*i))->q;
		}else{
			x[index]=((LaplaceFace3d*)(*i))->phi;
		}
	}
	//
	if(tree->generation_table.size()<=2){
		cout << "This simulation is same calculation of the conventional boundary element method." << endl;
	}
	if(g=="zero_potential_z0"){
		greenFunc=new potential_green_func::GreenFuncZeroPotentialOnWall;
	}else if(g=="no_flux_z0"){
		greenFunc=new potential_green_func::GreenFuncNoFluxOnWall;
	}else{
		greenFunc=new potential_green_func::GreenFunc;
	}
}

FMMNonBiCGSTAB::~FMMNonBiCGSTAB(){
	delete greenFunc;
	delete integral;
}

void FMMNonBiCGSTAB::psolve(double* xx,double* bb){
	for(int i=0;i<N;i++){
		xx[i]=bb[i];
	}
}

int FMMNonBiCGSTAB::update(){
	calc_vector_b();
	int info=BiCGSTAB::update();
	return info;
}

void FMMNonBiCGSTAB::calc_vector_b(){
	rnm.resize(local_number);
	calc_moment_leaf_b();
	calc_moment_M2M();
	calc_local_expansion_2nd_generation();
	calc_local_expansion();
	vector<Face3d*>::iterator fi=mesh->face.begin();
	for(int i=0;fi!=mesh->face.end();fi++,i++){
		Vector3d& ri=(*fi)->position;
		Vector3d& ri0=(*fi)->vertex[0]->position;
		Vector3d& ri1=(*fi)->vertex[1]->position;
		Vector3d& ri2=(*fi)->vertex[2]->position;
		string& type_i=(*fi)->type;
		//from self effect
		if(type_i=="phi"){
			b[i]=-0.5*((LaplaceFace3d*)(*fi))->phi +phi+gphi*ri;
		}else if(type_i=="q"){
			b[i]=-((LaplaceFace3d*)(*fi))->q*calc_Aaa(ri,ri0,ri1,ri2) +phi+gphi*ri;
		}
		//from faces including in same cell 
		b[i]+=calc_face_b(*fi,(*fi)->cell);
		//from edges including in adjacent cells
		vector<Cell3d*>::iterator c=(*fi)->cell->adjacent.begin();
		for(;c!=(*fi)->cell->adjacent.end();c++){
			b[i]+=calc_face_b(*fi,*c);
		}
		//from far cells using local expansion
		Vector3d r=(*fi)->position-(*fi)->cell->position;
		rnm.update(r);
		Mnm& local_L=((LaplaceMoment3d*)((*fi)->cell->local))->M;
		//Mnm& local_K=((LaplaceMoment3d*)((*fi)->cell->local))->N;
		for(int n=0;n<=local_number;n++){
			for(int m=-n;m<=n;m++){
				b[i]+=(0.25/PI*rnm.get(n,m)*local_L.get(n,m)).real();
				//b[i]+=(-0.25/PI*rnm.get(n,m)*local_L.get(n,m)).real();
				//b[i]+=( 0.25/PI*rnm.get(n,m)*local_K.get(n,m)).real();
			}
		}
	}
	cout << "\n\tcalc_vector_b function is done" << endl;
}

double FMMNonBiCGSTAB::calc_face_b(Face3d* fi,Cell3d* cell){
	double s=0.0;
	if(cell->face.size()==0){
		vector<Cell3d*>::iterator c=cell->child.begin();
		for(;c!=cell->child.end();c++){
			s+=calc_face_b(fi,*c);
		}
		return s;
	}
	Vector3d& ri=fi->position;
	vector<Face3d*>::iterator fj=cell->face.begin();
	for(;fj!=cell->face.end();fj++){
		Vector3d& rj0=(*fj)->vertex[0]->position;
		Vector3d& rj1=(*fj)->vertex[1]->position;
		Vector3d& rj2=(*fj)->vertex[2]->position;
		Vector3d& nj=(*fj)->normal;
		double area_j=(*fj)->area;
		string& type_j=(*fj)->type;
		if(fi!=(*fj)){
			if(type_j=="phi"){
				s+= ((LaplaceFace3d*)(*fj))->phi*calc_Bab(ri,rj0,rj1,rj2,nj,area_j);
			}else if(type_j=="q"){
				s+=-((LaplaceFace3d*)(*fj))->q*calc_Aab(ri,rj0,rj1,rj2,area_j);
			}
		}
	}
	return s;
}

double FMMNonBiCGSTAB::calc_face_x(Face3d* fi,Cell3d* cell){
	double s=0.0;
	if(cell->face.size()==0){
		vector<Cell3d*>::iterator c=cell->child.begin();
		for(;c!=cell->child.end();c++){
			s+=calc_face_x(fi,*c);
		}
		return s;
	}
	Vector3d& ri=fi->position;
	vector<Face3d*>::iterator fj=cell->face.begin();
	for(;fj!=cell->face.end();fj++){
		Vector3d& rj0=(*fj)->vertex[0]->position;
		Vector3d& rj1=(*fj)->vertex[1]->position;
		Vector3d& rj2=(*fj)->vertex[2]->position;
		Vector3d& nj=(*fj)->normal;
		double area_j=(*fj)->area;
		string& type_j=(*fj)->type;
		if(fi!=(*fj)){
			if(type_j=="q"){
				s+=-((LaplaceFace3d*)(*fj))->phi*calc_Bab(ri,rj0,rj1,rj2,nj,area_j);
			}else if(type_j=="phi"){
				s+= ((LaplaceFace3d*)(*fj))->q*calc_Aab(ri,rj0,rj1,rj2,area_j);
			}
		}
	}
	return s;
}

double FMMNonBiCGSTAB::calc_face_evaluation_point(const Vector3d& ri,Cell3d* cell){
	double s=0.0;
	if(cell->face.size()==0){
		vector<Cell3d*>::iterator c=cell->child.begin();
		for(;c!=cell->child.end();c++){
			s+=calc_face_evaluation_point(ri,*c);
		}
		return s;
	}
	vector<Face3d*>::iterator fj=cell->face.begin();
	for(;fj!=cell->face.end();fj++){
		Vector3d& rj0=(*fj)->vertex[0]->position;
		Vector3d& rj1=(*fj)->vertex[1]->position;
		Vector3d& rj2=(*fj)->vertex[2]->position;
		Vector3d& nj=(*fj)->normal;
		double area_j=(*fj)->area;
		s+=-((LaplaceFace3d*)(*fj))->q*calc_Aab(ri,rj0,rj1,rj2,area_j);
		s+= ((LaplaceFace3d*)(*fj))->phi*calc_Bab(ri,rj0,rj1,rj2,nj,area_j);
	}
	s+=phi+gphi*ri;
	return s;
}

void FMMNonBiCGSTAB::matvec(double alpha,double *xx,double beta,double* yy){
	//potential_green_func::Rnm rnm(local_number);
	rnm.resize(local_number);
	vector<Face3d*>::iterator i=mesh->face.begin();
	for(int index=0;i!=mesh->face.end();i++,index++){
		if((*i)->type=="phi"){
			((LaplaceFace3d*)(*i))->q=xx[index];
		}else{
			((LaplaceFace3d*)(*i))->phi=xx[index];
		}
	}
	calc_moment_leaf_x();
	calc_moment_M2M();
	calc_local_expansion_2nd_generation();
	calc_local_expansion();
	vector<Face3d*>::iterator fi=mesh->face.begin();
	for(int i=0;fi!=mesh->face.end();fi++,i++){
		yy[i]=beta*yy[i];
		Vector3d& ri=(*fi)->position;
		Vector3d& ri0=(*fi)->vertex[0]->position;
		Vector3d& ri1=(*fi)->vertex[1]->position;
		Vector3d& ri2=(*fi)->vertex[2]->position;
		string& type_i=(*fi)->type;
		//from self effect
		if(type_i=="q"){
			yy[i]+=alpha*0.5*((LaplaceFace3d*)(*fi))->phi;
		}else if(type_i=="phi"){
			yy[i]+=alpha*((LaplaceFace3d*)(*fi))->q*calc_Aaa(ri,ri0,ri1,ri2);
		}
		//from faces including in same cell 
		yy[i]+=alpha*calc_face_x(*fi,(*fi)->cell);
		//from faces including in adjacent cells
		vector<Cell3d*>::iterator c=(*fi)->cell->adjacent.begin();
		for(;c!=(*fi)->cell->adjacent.end();c++){
			yy[i]+=alpha*calc_face_x(*fi,*c);
		}
		//from far cells using local expansion
		Vector3d r=(*fi)->position-(*fi)->cell->position;
		rnm.update(r);
		Mnm& local_L=((LaplaceMoment3d*)((*fi)->cell->local))->M;
		//Mnm& local_K=((LaplaceMoment3d*)((*fi)->cell->local))->N;
		for(int n=0;n<=local_number;n++){
			for(int m=-n;m<=n;m++){
				yy[i]-=alpha*( 0.25/PI*rnm.get(n,m)*local_L.get(n,m)).real();
				//yy[i]+=alpha*( 0.25/PI*rnm.get(n,m)*local_L.get(n,m)).real();
				//yy[i]+=alpha*(-0.25/PI*rnm.get(n,m)*local_K.get(n,m)).real();
			}
		}
	}
	cout << "\tmatvec function is done (" << resid << ")" << endl;
}

void FMMNonBiCGSTAB::init_evaluation_point(){
	calc_moment_leaf_evaluation_point();
	calc_moment_M2M();
	calc_local_expansion_2nd_generation();
	calc_local_expansion();
}

double FMMNonBiCGSTAB::calc_evaluation_point(Cell3d* ce,const Vector3d& ri){
	//potential_green_func::Rnm rnm(local_number);
	rnm.resize(local_number);
	double phi=0.0;
	//from faces including in same cell 
	phi+=calc_face_evaluation_point(ri,ce);
	//from faces including in adjacent cells
	vector<Cell3d*>::iterator c=ce->adjacent.begin();
	for(;c!=ce->adjacent.end();c++){
		phi+=calc_face_evaluation_point(ri,*c);
	}
	//from far cells using local expansion
	Vector3d r=ri - ce->position;
	Mnm& local_L=((LaplaceMoment3d*)(ce->local))->M;
	//Mnm& local_K=((LaplaceMoment3d*)(ce->local))->N;
	rnm.update(r);
	for(int n=0;n<=local_number;n++){
		for(int m=-n;m<=n;m++){
			phi+=(0.25/PI*rnm.get(n,m)*local_L.get(n,m)).real();
			//phi+=(-0.25/PI*rnm.get(n,m)*local_L.get(n,m)).real();
			//phi+=( 0.25/PI*rnm.get(n,m)*local_K.get(n,m)).real();
		}
	}
	return phi;
}

void FMMNonBiCGSTAB::calc_local_expansion(){
	if(tree->generation_table.size()<=2)return;
	//potential_green_func::Rnm rnm(local_number);
	//potential_green_func::Snm snm(moments_number+local_number);
	rnm.resize(local_number);
	snm.resize(moments_number+local_number);
	vector<vector<Cell3d*> >::iterator g=tree->generation_table.begin()+3;
	for(;g!=tree->generation_table.end();g++){
		vector<Cell3d*>::iterator i=g->begin();// child
		for(;i!=g->end();i++){
			((LaplaceMoment3d*)((*i)->local))->zero_reset();
			//Local expansion from L2L being parent
			Vector3d r=(*i)->position-(*i)->parent->position;
			rnm.update(r);
			for(int n=0;n<=local_number;n++){
				for(int m=0;m<=n;m++){
					for(int np=0;np<=local_number;np++){
						for(int mp=-np;mp<=np;mp++){
							((LaplaceMoment3d*)((*i)->local))->M.moment[n][m]+=
								rnm.get(np-n,mp-m)*((LaplaceMoment3d*)((*i)->parent->local))->M.get(np,mp);
							//((LaplaceMoment3d*)((*i)->local))->N.moment[n][m]+=
							//	rnm.get(np-n,mp-m)*((LaplaceMoment3d*)((*i)->parent->local))->N.get(np,mp);
						}
					}
				}
			}
			
			//M2L from interaction cells
			vector<Cell3d*>::iterator j=(*i)->interaction.begin();
			for(;j!=(*i)->interaction.end();j++){
				Vector3d r=(*i)->position-(*j)->position;
				snm.update(r);
				for(int n=0;n<=local_number;n++){
					for(int m=0;m<=n;m++){
						for(int np=0;np<=moments_number;np++){
							for(int mp=-np;mp<=np;mp++){
								double c=1.0;
								if(n%2==1)c=-1;
								//Local expansion is evaluated by M2L only.
								((LaplaceMoment3d*)((*i)->local))->M.moment[n][m]+=
									c*conj(snm.get(n+np,m+mp))*((LaplaceMoment3d*)((*j)->moment))->M.get(np,mp);
								//((LaplaceMoment3d*)((*i)->local))->N.moment[n][m]+=
								//	c*conj(snm.get(n+np,m+mp))*((LaplaceMoment3d*)((*j)->moment))->N.get(np,mp);
							}
						}
					}
				}
			}
		}
	}
}

void FMMNonBiCGSTAB::calc_local_expansion_2nd_generation(){
	//potential_green_func::Snm snm(moments_number+local_number);
	snm.resize(moments_number+local_number);
	if(tree->generation_table.size()<=2)return;
	vector<Cell3d*>::iterator i=tree->generation_table[2].begin();
	for(;i!=tree->generation_table[2].end();i++){
		((LaplaceMoment3d*)((*i)->local))->zero_reset();
		vector<Cell3d*>::iterator j=(*i)->interaction.begin();
		for(;j!=(*i)->interaction.end();j++){
			Vector3d r=(*i)->position-(*j)->position;
			snm.update(r);
			for(int n=0;n<=local_number;n++){
				for(int m=0;m<=n;m++){
					for(int np=0;np<=moments_number;np++){
						for(int mp=-np;mp<=np;mp++){
							double c=1.0;
							if(n%2==1)c=-1;
							//Local expansion is evaluated by M2L only.
							((LaplaceMoment3d*)((*i)->local))->M.moment[n][m]+=
								c*conj(snm.get(n+np,m+mp))*((LaplaceMoment3d*)((*j)->moment))->M.get(np,mp);
							//((LaplaceMoment3d*)((*i)->local))->N.moment[n][m]+=
							//	c*conj(snm.get(n+np,m+mp))*((LaplaceMoment3d*)((*j)->moment))->N.get(np,mp);
						}
					}
				}
			}
		}
	}
}

void FMMNonBiCGSTAB::calc_moment_M2M(){
	if(tree->generation_table.size()<=2)return;
	rnm.resize(moments_number);
	//from generation_max-1 (moments in leaves are calculated in calc_moment_leaf)
	vector<vector<Cell3d*> >::iterator i=tree->generation_table.end()-2;
	for(;i!=tree->generation_table.begin()+1;i--){//from max_generation-1 to generation 2 
		vector<Cell3d*>::iterator j=i->begin();// parent
		for(;j!=i->end();j++){
			if((*j)->face.size()!=0)continue;//if you are leaf
			((LaplaceMoment3d*)(*j)->moment)->zero_reset();
			//get the moments from his children
			vector<Cell3d*>::iterator c=(*j)->child.begin();// child
			for(;c!=(*j)->child.end();c++){
				Vector3d r=(*c)->position-(*j)->position;
				rnm.update(r);
				for(int n=0;n<=moments_number;n++){
					for(int m=0;m<=n;m++){
						for(int np=0;np<=n;np++){
							for(int mp=-np;mp<=np;mp++){
								((LaplaceMoment3d*)((*j)->moment))->M.moment[n][m]+=
									rnm.get(np,mp)*((LaplaceMoment3d*)((*c)->moment))->M.get(n-np,m-mp);
								//((LaplaceMoment3d*)((*j)->moment))->N.moment[n][m]+=
								//	rnm.get(np,mp)*((LaplaceMoment3d*)((*c)->moment))->N.get(n-np,m-mp);
							}
						}
					}
				}
			}
		}
	}
}

void FMMNonBiCGSTAB::calc_moment_leaf_b(){
	rnm.resize(moments_number);
	//moments in leaves are calculated
	vector<Cell3d*>::iterator i=tree->leaf.begin();
	for(;i!=tree->leaf.end();i++){
		((LaplaceMoment3d*)((*i)->moment))->zero_reset();
		vector<Face3d*>::iterator j=(*i)->face.begin();
		for(;j!=(*i)->face.end();j++){
			Vector3d& rj0=(*j)->vertex[0]->position;
			Vector3d& rj1=(*j)->vertex[1]->position;
			Vector3d& rj2=(*j)->vertex[2]->position;
			double area_j=(*j)->area;
			for(int a=0;a<integral->point;a++){
				Vector3d r=integral->xi[a].x*(rj1-rj0)+integral->xi[a].y*(rj2-rj0)+rj0
				           -(*i)->position;
				rnm.update(r);
				for(int n=0;n<=moments_number;n++){
					for(int m=0;m<=n;m++){
						if((*j)->type=="q"){
							((LaplaceMoment3d*)((*i)->moment))->M.moment[n][m]-=
								integral->w[a]*area_j*rnm.get(n,m)*((LaplaceFace3d*)(*j))->q;
							//((LaplaceMoment3d*)((*i)->moment))->M.moment[n][m]+=
							//	integral->w[a]*area_j*rnm.get(n,m)*((LaplaceFace3d*)(*j))->q;
						}else{
							((LaplaceMoment3d*)((*i)->moment))->M.moment[n][m]+=
								integral->w[a]*area_j*rnm.get_dRdn(n,m,(*j)->normal)*((LaplaceFace3d*)(*j))->phi;
							//((LaplaceMoment3d*)((*i)->moment))->N.moment[n][m]+=
							//	integral->w[a]*area_j*rnm.get_dRdn(n,m,(*j)->normal)*((LaplaceFace3d*)(*j))->phi;
						}
					}
				}
			}
		}
	}
}

void FMMNonBiCGSTAB::calc_moment_leaf_x(){
	rnm.resize(moments_number);
	//moments in leaves are calculated
	vector<Cell3d*>::iterator i=tree->leaf.begin();
	for(;i!=tree->leaf.end();i++){
		((LaplaceMoment3d*)(*i)->moment)->zero_reset();
		vector<Face3d*>::iterator j=(*i)->face.begin();
		for(;j!=(*i)->face.end();j++){
			Vector3d& rj0=(*j)->vertex[0]->position;
			Vector3d& rj1=(*j)->vertex[1]->position;
			Vector3d& rj2=(*j)->vertex[2]->position;
			double area_j=(*j)->area;
			for(int a=0;a<integral->point;a++){
				Vector3d r=integral->xi[a].x*(rj1-rj0)+integral->xi[a].y*(rj2-rj0)+rj0
				           -(*i)->position;
				rnm.update(r);
				for(int n=0;n<=moments_number;n++){
					for(int m=0;m<=n;m++){
						if((*j)->type=="phi"){
							((LaplaceMoment3d*)((*i)->moment))->M.moment[n][m]-=
								integral->w[a]*area_j*rnm.get(n,m)*((LaplaceFace3d*)(*j))->q;
							//((LaplaceMoment3d*)((*i)->moment))->M.moment[n][m]+=
							//	integral->w[a]*area_j*rnm.get(n,m)*((LaplaceFace3d*)(*j))->q;
						}else{
							((LaplaceMoment3d*)((*i)->moment))->M.moment[n][m]+=
								integral->w[a]*area_j*rnm.get_dRdn(n,m,(*j)->normal)*((LaplaceFace3d*)(*j))->phi;
							//((LaplaceMoment3d*)((*i)->moment))->N.moment[n][m]+=
							//	integral->w[a]*area_j*rnm.get_dRdn(n,m,(*j)->normal)*((LaplaceFace3d*)(*j))->phi;
						}
					}
				}
			}
		}
	}
}

void FMMNonBiCGSTAB::calc_moment_leaf_evaluation_point(){
	rnm.resize(moments_number);
	//moments in leaves are calculated
	vector<Cell3d*>::iterator i=tree->leaf.begin();
	for(;i!=tree->leaf.end();i++){
		((LaplaceMoment3d*)((*i)->moment))->zero_reset();
		vector<Face3d*>::iterator j=(*i)->face.begin();
		for(;j!=(*i)->face.end();j++){
			Vector3d& rj0=(*j)->vertex[0]->position;
			Vector3d& rj1=(*j)->vertex[1]->position;
			Vector3d& rj2=(*j)->vertex[2]->position;
			double area_j=(*j)->area;
			for(int a=0;a<integral->point;a++){
				Vector3d r=integral->xi[a].x*(rj1-rj0)+integral->xi[a].y*(rj2-rj0)+rj0
				           -(*i)->position;
				rnm.update(r);
				for(int n=0;n<=moments_number;n++){
					for(int m=0;m<=n;m++){
						((LaplaceMoment3d*)((*i)->moment))->M.moment[n][m]+=
							-integral->w[a]*area_j*rnm.get(n,m)*((LaplaceFace3d*)(*j))->q;
							+integral->w[a]*area_j*rnm.get_dRdn(n,m,(*j)->normal)*((LaplaceFace3d*)(*j))->phi;
						//((LaplaceMoment3d*)((*i)->moment))->M.moment[n][m]+=
						//	integral->w[a]*area_j*rnm.get(n,m)*((LaplaceFace3d*)(*j))->q;
						//((LaplaceMoment3d*)((*i)->moment))->N.moment[n][m]+=
						//	integral->w[a]*area_j*rnm.get_dRdn(n,m,(*j)->normal)*((LaplaceFace3d*)(*j))->phi;
					}
				}
			}
		}
	}
}

double FMMNonBiCGSTAB::calc_Aab(
	const Vector3d& r0,const Vector3d& ra,const Vector3d& rb,const Vector3d& rc,double area)const{
	double s=0.0;
	for(int l=0;l<integral->point;l++){
		Vector3d r=integral->xi[l].x*(rb-ra)+integral->xi[l].y*(rc-ra)+ra;
		s+=integral->w[l]*greenFunc->get_G(r0,r);
	}
	return area*s;
}

double FMMNonBiCGSTAB::calc_Bab(
	const Vector3d& r0,const Vector3d& ra,const Vector3d& rb,const Vector3d& rc,const Vector3d& n,double area)const{
	double s=0.0;
	for(int l=0;l<integral->point;l++){
		Vector3d r=integral->xi[l].x*(rb-ra)+integral->xi[l].y*(rc-ra)+ra;
		s+=integral->w[l]*greenFunc->get_Tn(r0,r,n);
	}
	return area*s;
}

double FMMNonBiCGSTAB::calc_Aaa(const Vector3d& r0,const Vector3d& ra,const Vector3d& rb,const Vector3d& rc)const{
	double s0=0.5*((ra-r0)^(rb-r0)).length();double s1=0.5*((rb-r0)^(rc-r0)).length();
	double s2=0.5*((rc-r0)^(ra-r0)).length();
	double a=calc_Aab(r0,r0,ra,rb,s0);
	a+=calc_Aab(r0,r0,rb,rc,s1);a+=calc_Aab(r0,r0,rc,ra,s2);
	return a;
}

#undef PI
