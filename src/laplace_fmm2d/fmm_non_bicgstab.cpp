#include "fmm_non_bicgstab.h"

#define PI 3.14159265358979

FMMNonBiCGSTAB::FMMNonBiCGSTAB(int n,FMMTree_c* t,Mesh_c* m,double* vx,double* vb,
	                             int mn,int ln,int max_ite,double tol,int ipn,const string& gFunc,
	                             double ph,const complex<double>& gp)
:BiCGSTAB(n,vx,vb,max_ite,tol),tree(t),mesh(m),moments_number(mn),local_number(ln),phi(ph),gphi(gp){
	integral= new ParameterOfLegendreGaussFormula1d(ipn);
	//for(int i=0;i<n;i++)x[i]=0.0;
	//input the initial guess about x 
	vector<Edge_c*>::iterator i=mesh->edge.begin();
	for(int index=0;i!=mesh->edge.end();i++,index++){
		if((*i)->type=="phi"){
			x[index]=((LaplaceEdge_c*)(*i))->q;
		}else{
			x[index]=((LaplaceEdge_c*)(*i))->phi;
		}
	}
	if(tree->generation_table.size()<=2){
		cout << "This simulation is same calculation of the conventional boundary element method." << endl;
	}
	if(gFunc=="free_analytic_integral"){
		greenFunc=0;
	}else if(gFunc=="free"){
		greenFunc=new potential_green_func_c::GreenFunc;
	}else if(gFunc=="zero_potential_y0"){
		greenFunc=new potential_green_func_c::GreenFuncZeroPotentialOnWall;
	}else if(gFunc=="no_flux_y0"){
		greenFunc=new potential_green_func_c::GreenFuncNoFluxOnWall;
	}
}

FMMNonBiCGSTAB::~FMMNonBiCGSTAB(){
	delete integral;
	if(greenFunc!=0)delete greenFunc;
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

double FMMNonBiCGSTAB::calc_edge_b(Edge_c* ei,Cell_c* cell){
	double s=0.0;
	if(cell->edge.size()==0){
		vector<Cell_c*>::iterator c=cell->child.begin();
		for(;c!=cell->child.end();c++){
			s+=calc_edge_b(ei,*c);
		}
		return s;
	}
	complex<double>& ri=ei->position;
	vector<Edge_c*>::iterator ej=cell->edge.begin();
	for(;ej!=cell->edge.end();ej++){
		complex<double>& rj0=(*ej)->vertex[0]->position;
		complex<double>& rj1=(*ej)->vertex[1]->position;
		complex<double>& nj=(*ej)->normal;
		double length_j=(*ej)->length;
		string& type_j=(*ej)->type;
		if(ei!=(*ej)){
			if(type_j=="phi"){
				s+= ((LaplaceEdge_c*)(*ej))->phi*calc_Bab(ri,rj0,rj1,nj,length_j);
			}else if(type_j=="q"){
				s+=-((LaplaceEdge_c*)(*ej))->q*calc_Aab(ri,rj0,rj1,length_j);
			}
		}
	}
	return s;
}

double FMMNonBiCGSTAB::calc_edge_x(Edge_c* ei,Cell_c* cell){
	double s=0.0;
	if(cell->edge.size()==0){
		vector<Cell_c*>::iterator c=cell->child.begin();
		for(;c!=cell->child.end();c++){
			s+=calc_edge_x(ei,*c);
		}
		return s;
	}
	complex<double>& ri=ei->position;
	vector<Edge_c*>::iterator ej=cell->edge.begin();
	for(;ej!=cell->edge.end();ej++){
		complex<double>& rj0=(*ej)->vertex[0]->position;
		complex<double>& rj1=(*ej)->vertex[1]->position;
		complex<double>& nj=(*ej)->normal;
		double length_j=(*ej)->length;
		string& type_j=(*ej)->type;
		if(ei!=(*ej)){
			if(type_j=="q"){
				s+=-((LaplaceEdge_c*)(*ej))->phi*calc_Bab(ri,rj0,rj1,nj,length_j);
			}else if(type_j=="phi"){
				s+= ((LaplaceEdge_c*)(*ej))->q*calc_Aab(ri,rj0,rj1,length_j);
			}
		}
	}
	return s;
}

double FMMNonBiCGSTAB::calc_edge_evaluation_point(const complex<double>& ri,Cell_c* cell){
	double s=0.0;
	if(cell->edge.size()==0){
		vector<Cell_c*>::iterator c=cell->child.begin();
		for(;c!=cell->child.end();c++){
			s+=calc_edge_evaluation_point(ri,*c);
		}
		return s;
	}
	vector<Edge_c*>::iterator ej=cell->edge.begin();
	for(;ej!=cell->edge.end();ej++){
		complex<double>& rj0=(*ej)->vertex[0]->position;
		complex<double>& rj1=(*ej)->vertex[1]->position;
		complex<double>& nj=(*ej)->normal;
		double length_j=(*ej)->length;
		s+=-((LaplaceEdge_c*)(*ej))->q*calc_Aab(ri,rj0,rj1,length_j);
		s+= ((LaplaceEdge_c*)(*ej))->phi*calc_Bab(ri,rj0,rj1,nj,length_j);
	}
	return s;
}

void FMMNonBiCGSTAB::calc_vector_b(){
	calc_moment_leaf_b();
	calc_moment_M2M();
	calc_local_expansion_2nd_generation();
	calc_local_expansion();
	vector<Edge_c*>::iterator ei=mesh->edge.begin();
	for(int i=0;ei!=mesh->edge.end();ei++,i++){
		complex<double>& ri=(*ei)->position;
		complex<double>& ri0=(*ei)->vertex[0]->position;
		complex<double>& ri1=(*ei)->vertex[1]->position;
		string& type_i=(*ei)->type;
		//from self effect
		if(type_i=="phi"){
			b[i]=-0.5*((LaplaceEdge_c*)(*ei))->phi +phi+gphi.real()*ri.real()+gphi.imag()*ri.imag();
		}else if(type_i=="q"){
			b[i]=-((LaplaceEdge_c*)(*ei))->q*calc_Aaa(ri,ri0,ri1,(*ei)->length)+phi+gphi.real()*ri.real()+gphi.imag()*ri.imag();
		}
		//from edges including in same cell 
		b[i]+=calc_edge_b(*ei,(*ei)->cell);
		//from edges including in adjacent cells
		vector<Cell_c*>::iterator c=(*ei)->cell->adjacent.begin();
		for(;c!=(*ei)->cell->adjacent.end();c++){
			b[i]+=calc_edge_b(*ei,*c);
		}
		//from far cells using local expansion
		complex<double> z=(*ei)->position-(*ei)->cell->position;
		vector<complex<double> >& local_L=((LaplaceMoment_c*) ((*ei)->cell->local))->M;//local_L
		//vector<complex<double> >& local_K=((LaplaceMoment_c*) ((*ei)->cell->local))->N;//local_K
		//for(int ii=0;ii<=local_L.size();ii++){
		for(int ii=0;ii<=local_number;ii++){
			b[i]+=(local_L[ii]*potential_green_func_c::get_Ik(z,ii)).real();
			//b[i]+=-(local_L[ii]*potential_green_func_c::get_Ik(z,ii)).real();
			//b[i]+= (local_K[ii]*potential_green_func_c::get_Ik(z,ii)).real();
		}
	}
	cout << "\n\tcalc_vector_b function is done" << endl;
}

void FMMNonBiCGSTAB::matvec(double alpha,double *xx,double beta,double* yy){
	vector<Edge_c*>::iterator i=mesh->edge.begin();
	for(int index=0;i!=mesh->edge.end();i++,index++){
		if((*i)->type=="phi"){
			((LaplaceEdge_c*)(*i))->q=xx[index];
		}else{
			((LaplaceEdge_c*)(*i))->phi=xx[index];
		}
	}
	calc_moment_leaf_x();
	calc_moment_M2M();
	calc_local_expansion_2nd_generation();
	calc_local_expansion();
	vector<Edge_c*>::iterator ei=mesh->edge.begin();
	for(int i=0;ei!=mesh->edge.end();ei++,i++){
		yy[i]=beta*yy[i];
		complex<double>& ri=(*ei)->position;
		complex<double>& ri0=(*ei)->vertex[0]->position;
		complex<double>& ri1=(*ei)->vertex[1]->position;
		string& type_i=(*ei)->type;
		//from self effect
		if(type_i=="q"){
			yy[i]+=alpha*0.5*((LaplaceEdge_c*)(*ei))->phi;
		}else if(type_i=="phi"){
			yy[i]+=alpha*((LaplaceEdge_c*)(*ei))->q*calc_Aaa(ri,ri0,ri1,(*ei)->length);
		}
		//from edges including in same cell 
		yy[i]+=alpha*calc_edge_x(*ei,(*ei)->cell);
		//from edges including in adjacent cells
		vector<Cell_c*>::iterator c=(*ei)->cell->adjacent.begin();
		for(;c!=(*ei)->cell->adjacent.end();c++){
			yy[i]+=alpha*calc_edge_x(*ei,*c);
		}
		//from far cells using local expansion
		complex<double> z=(*ei)->position-(*ei)->cell->position;
		vector<complex<double> >& local_L=((LaplaceMoment_c*) ((*ei)->cell->local)) ->M; //L;
		//vector<complex<double> >& local_K=((LaplaceMoment_c*) ((*ei)->cell->local)) ->N; //K;
		//for(int ii=0;ii<=local_L.size();ii++){
		for(int ii=0;ii<=local_number;ii++){
			yy[i]+=-alpha*(local_L[ii]*potential_green_func_c::get_Ik(z,ii)).real();
			//yy[i]+=-alpha*(local_K[ii]*potential_green_func_c::get_Ik(z,ii)).real();
			//yy[i]+= alpha*(local_L[ii]*potential_green_func_c::get_Ik(z,ii)).real();
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

double FMMNonBiCGSTAB::calc_evaluation_point(Cell_c* ei,const complex<double>& ri){
	double phi=this->phi+gphi.real()*ri.real()+gphi.imag()*ri.imag();
	//from edges including in same cell 
	phi+=calc_edge_evaluation_point(ri,ei);
	//from edges including in adjacent cells
	vector<Cell_c*>::iterator c=ei->adjacent.begin();
	for(;c!=ei->adjacent.end();c++){
		phi+=calc_edge_evaluation_point(ri,*c);
	}
	//from far cells using local expansion
	complex<double> z=ri - ei->position;
	vector<complex<double> >& local_L=((LaplaceMoment_c*) (ei->local))->M;  //L;
	//vector<complex<double> >& local_K=((LaplaceMoment_c*) (ei->local))->N;  //K;
	//for(int ii=0;ii<=local_L.size();ii++){
	for(int ii=0;ii<=local_number;ii++){
		phi+=(local_L[ii]*potential_green_func_c::get_Ik(z,ii)).real();
		//phi+=-(local_L[ii]*potential_green_func_c::get_Ik(z,ii)).real();
		//phi+= (local_K[ii]*potential_green_func_c::get_Ik(z,ii)).real();
	}
	return phi;
}

void FMMNonBiCGSTAB::calc_moment_leaf_b(){
	//moments in leaves are calculated
	vector<Cell_c*>::iterator i=tree->leaf.begin();
	for(;i!=tree->leaf.end();i++){
		(*i)->moment->zero_reset();
		vector<Edge_c*>::iterator j=(*i)->edge.begin();
		for(;j!=(*i)->edge.end();j++){
			//complex<double> z=(*j)->position-(*i)->position;
			complex<double>& rj0=(*j)->vertex[0]->position;
			complex<double>& rj1=(*j)->vertex[1]->position;
			for(int k=0;k<=moments_number;k++){
				for(int a=0;a<integral->point;a++){
					complex<double> r=0.5*((1.0-integral->xi[a])*rj0+(integral->xi[a]+1.0)*rj1);
					complex<double> z=r-(*i)->position;
					if((*j)->type=="q"){
						//((StokesMoment_c*) ((*i)->moment)) ->M[k]-=0.5*integral->w[a]*(*j)->length*potential_green_func_c::get_Ik(z,k)*((LaplaceEdge_c*)(*j))->q;
						//((LaplaceMoment_c*) ((*i)->moment)) ->M[k]+=0.5*integral->w[a]*(*j)->length*potential_green_func_c::get_Ik(z,k)*((LaplaceEdge_c*)(*j))->q;
					}else{
						//((StokesMoment_c*) ((*i)->moment)) ->M[k]+=0.5*integral->w[a]*(*j)->length*((LaplaceEdge_c*)(*j))->phi*(*j)->normal*potential_green_func_c::get_Ik(z,k-1);
						//((LaplaceMoment_c*) ((*i)->moment)) ->N[k]+=0.5*integral->w[a]*(*j)->length*((LaplaceEdge_c*)(*j))->phi*(*j)->normal*potential_green_func_c::get_Ik(z,k-1);
					}
				}
			}
		}
	}
}

void FMMNonBiCGSTAB::calc_moment_leaf_x(){
	//moments in leaves are calculated
	vector<Cell_c*>::iterator i=tree->leaf.begin();
	for(;i!=tree->leaf.end();i++){
		(*i)->moment->zero_reset();
		vector<Edge_c*>::iterator j=(*i)->edge.begin();
		for(;j!=(*i)->edge.end();j++){
			complex<double>& rj0=(*j)->vertex[0]->position;
			complex<double>& rj1=(*j)->vertex[1]->position;
			//complex<double> z=(*j)->position-(*i)->position;
			for(int k=0;k<=moments_number;k++){
				for(int a=0;a<integral->point;a++){
					complex<double> r=0.5*((1.0-integral->xi[a])*rj0+(integral->xi[a]+1.0)*rj1);
					complex<double> z=r-(*i)->position;
					if((*j)->type=="phi"){
						((LaplaceMoment_c*) ((*i)->moment)) ->M[k]-=0.5*integral->w[a]*(*j)->length*potential_green_func_c::get_Ik(z,k)*((LaplaceEdge_c*)(*j))->q;
						//((LaplaceMoment_c*) ((*i)->moment)) ->M[k]+=0.5*integral->w[a]*(*j)->length*potential_green_func_c::get_Ik(z,k)*((LaplaceEdge_c*)(*j))->q;
					}else{
						((LaplaceMoment_c*) ((*i)->moment)) ->M[k]+=0.5*integral->w[a]*(*j)->length*((LaplaceEdge_c*)(*j))->phi*(*j)->normal*potential_green_func_c::get_Ik(z,k-1);
						//((LaplaceMoment_c*) ((*i)->moment)) ->N[k]+=0.5*integral->w[a]*(*j)->length*((LaplaceEdge_c*)(*j))->phi*(*j)->normal*potential_green_func_c::get_Ik(z,k-1);
					}
				}
			}
		}
	}
}

void FMMNonBiCGSTAB::calc_moment_leaf_evaluation_point(){
	//moments in leaves are calculated
	vector<Cell_c*>::iterator i=tree->leaf.begin();
	for(;i!=tree->leaf.end();i++){
		(*i)->moment->zero_reset();
		vector<Edge_c*>::iterator j=(*i)->edge.begin();
		for(;j!=(*i)->edge.end();j++){
			complex<double>& rj0=(*j)->vertex[0]->position;
			complex<double>& rj1=(*j)->vertex[1]->position;
			//complex<double> z=(*j)->position-(*i)->position;
			for(int k=0;k<=moments_number;k++){
				for(int a=0;a<integral->point;a++){
					complex<double> r=0.5*((1.0-integral->xi[a])*rj0+(integral->xi[a]+1.0)*rj1);
					complex<double> z=r-(*i)->position;
					((LaplaceMoment_c*) ((*i)->moment)) ->M[k]+=-0.5*integral->w[a]*(*j)->length*potential_green_func_c::get_Ik(z,k)*((LaplaceEdge_c*)(*j))->q;
					                                            +0.5*integral->w[a]*(*j)->length*((LaplaceEdge_c*)(*j))->phi*(*j)->normal*potential_green_func_c::get_Ik(z,k-1);
					//((LaplaceMoment_c*) ((*i)->moment)) ->N[k]+=0.5*integral->w[a]*(*j)->length*((LaplaceEdge_c*)(*j))->phi*(*j)->normal*potential_green_func_c::get_Ik(z,k-1);
				}
			}
		}
	}
}

void FMMNonBiCGSTAB::calc_moment_M2M(){
	if(tree->generation_table.size()<=2)return;
	//from generation_max-1 (moments in leaves are calculated in calc_moment_leaf)
	vector<vector<Cell_c*> >::iterator i=tree->generation_table.end()-2;
	for(;i!=tree->generation_table.begin()+1;i--){//from max_generation-1 to generation 2 
		vector<Cell_c*>::iterator j=i->begin();// parent
		for(;j!=i->end();j++){
			if((*j)->edge.size()!=0)continue;//if you are leaf
			(*j)->moment->zero_reset();
			//get the moments from his children
			vector<Cell_c*>::iterator c=(*j)->child.begin();// child
			for(;c!=(*j)->child.end();c++){
				complex<double> z=(*c)->position-(*j)->position;
				for(int k=0;k<=moments_number;k++){
					for(int l=0;l<=k;l++){
						((LaplaceMoment_c*) ((*j)->moment))->M[k]+=
								potential_green_func_c::get_Ik(z,k-l)*((LaplaceMoment_c*)((*c)->moment))->M[l];
						//((LaplaceMoment_c*) ((*j)->moment))->N[k]+=
						//		potential_green_func_c::get_Ik(z,k-l)*((LaplaceMoment_c*)((*c)->moment))->N[l];
					}
				}
			}
		}
	}
}

void FMMNonBiCGSTAB::calc_local_expansion_2nd_generation(){
	if(tree->generation_table.size()<=2)return;
	vector<Cell_c*>::iterator i=tree->generation_table[2].begin();
	for(;i!=tree->generation_table[2].end();i++){
		(*i)->local->zero_reset();
		vector<Cell_c*>::iterator j=(*i)->interaction.begin();
		for(;j!=(*i)->interaction.end();j++){
			complex<double> z=(*i)->position-(*j)->position;
			for(int k=0;k<=moments_number;k++){
				for(int l=0;l<=local_number;l++){
					double c=pow(-1.0,l)*0.5/PI;
					//Local expansion is evaluated by M2L only.
					//potential_green_func_c::get_Ok(z,l+k) << endl;
					((LaplaceMoment_c*) ((*i)->local))->M[l]+=
							c*potential_green_func_c::get_Ok(z,l+k)*((LaplaceMoment_c*) ((*j)->moment))->M[k];
					//((LaplaceMoment_c*) ((*i)->local))->N[l]+=
					//		c*potential_green_func_c::get_Ok(z,l+k)*((LaplaceMoment_c*) ((*j)->moment))->N[k];
				}
			}
		}
	}
}

void FMMNonBiCGSTAB::calc_local_expansion(){
	if(tree->generation_table.size()<=2)return;
	vector<vector<Cell_c*> >::iterator g=tree->generation_table.begin()+3;
	for(;g!=tree->generation_table.end();g++){
		vector<Cell_c*>::iterator i=g->begin();// child
		for(;i!=g->end();i++){
			(*i)->local->zero_reset();
			//Local expansion from L2L being parent
			complex<double> z=(*i)->position-(*i)->parent->position;
			for(int l=0;l<=local_number;l++){
				for(int m=l;m<=local_number;m++){
					((LaplaceMoment_c*) ((*i)->local))->M[l]+=potential_green_func_c::get_Ik(z,m-l)*((LaplaceMoment_c*)((*i)->parent)->local)->M[m];
					//((LaplaceMoment_c*) ((*i)->local))->N[l]+=potential_green_func_c::get_Ik(z,m-l)*((LaplaceMoment_c*)((*i)->parent)->local)->N[m];
				}
			}
			//M2L from interaction cells
			vector<Cell_c*>::iterator j=(*i)->interaction.begin();
			for(;j!=(*i)->interaction.end();j++){
				complex<double> z=(*i)->position-(*j)->position;
				for(int k=0;k<=moments_number;k++){
					for(int l=0;l<=local_number;l++){
						double c=pow(-1.0,l)*0.5/PI;
						((LaplaceMoment_c*) ((*i)->local))->M[l]+=c*potential_green_func_c::get_Ok(z,l+k)*((LaplaceMoment_c*) ((*j)->moment))->M[k];
						//((LaplaceMoment_c*) ((*i)->local))->N[l]+=c*potential_green_func_c::get_Ok(z,l+k)*((LaplaceMoment_c*) ((*j)->moment))->N[k];
					}
				}
			}
		}
	}
}

double FMMNonBiCGSTAB::calc_Aab(
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
		return 0.5/PI*((-dr2.real()*tt.real()-dr2.imag()*tt.imag())*log(dl2)
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

double FMMNonBiCGSTAB::calc_Bab(
	const complex<double>& r0,const complex<double>& ra,
	const complex<double>& rb,const complex<double>& n,double length)const{
	if(greenFunc==0){
		complex<double> dr1=ra-r0;
		complex<double> dr2=rb-r0;
		double dl1=sqrt(norm(dr1));
		double c1=dr1.real()/dl1;double s1=dr1.imag()/dl1;
		double dx2r=dr2.real()*c1+dr2.imag()*s1;
		double dy2r=-dr2.real()*s1+dr2.imag()*c1;
		double da=atan2(dy2r,dx2r);
		return da*0.5/PI;
	}else{
		double s=0.0;
		for(int l=0;l<integral->point;l++){
			complex<double> r=0.5*((1.0-integral->xi[l])*ra+(integral->xi[l]+1.0)*rb);
			s+=integral->w[l]*greenFunc->get_Tn(r0,r,n);
		}
		return 0.5*length*s;
	}
}

double FMMNonBiCGSTAB::calc_Aaa(const complex<double>& r0,const complex<double>& ra,
	                              const complex<double>& rb,double length)const{
	if(greenFunc==0){
		double l=0.5*length;
		return l*(-log(l)+1.0)/PI;
	}else{
		double la=sqrt(norm(r0-ra));
		return calc_Aab(r0,ra,r0,la)+calc_Aab(r0,r0,rb,length-la);
	}
}

#undef PI
