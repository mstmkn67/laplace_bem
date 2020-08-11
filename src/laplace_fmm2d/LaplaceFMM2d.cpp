#include "LaplaceFMM2d.h"

LaplaceFMM2d::LaplaceFMM2d(UDFManager* i,UDFManager* o)
:in(i),out(o){
	cout << "Simulation system is setting ....  " << flush;
	timer.start();
	assign_mesh();
	assign_region_condition();
	assign_initial_guess();
	create_cell();
	assign_bem();
	cout << " done (" << timer << " seconds )" << endl;
}

LaplaceFMM2d::~LaplaceFMM2d(){
	out->put("output.bicgstab.iteration",bem->bicgstab->max_iter);
	out->put("output.bicgstab.tolerance",bem->bicgstab->tol);
	delete bem;
	//vector<vector<Cell_c*> >::iterator c=tree->generation_table.begin();
	//for(;c!=tree->generation_table.end();c++){
	//	vector<Cell_c*>::iterator d=c->begin();
	//	for(;d!=c->end();d++){
	//		delete (*d)->moment;
	//		delete (*d)->local;
	//	}
	//}
	delete tree;
	vector<Edge_c*>::iterator i=mesh->edge.begin();
	for(;i!=mesh->edge.end();i++)delete *i;
	vector<Vertex_c*>::iterator j=mesh->vertex.begin();
	for(;j!=mesh->vertex.end();j++)delete *j;
	double t=timer.get();
	out->put("output.time",t);
	out->write();
}

void LaplaceFMM2d::run(){
	cout << "Calculation of BEM-FMM .... "  << flush;
	bem->update();
	cout << " done (" << timer << " seconds )" << endl;
	cout << "Reporting the imformation of edges and cells .... " << flush;
	report_edge();
	report_cell();
	cout << " done (" << timer << " seconds )" << endl;
	cout << "Evaluating the values at evaluation_point .... " << flush;
	report_evaluation_point();
	cout << " done (" << timer << " seconds )" << endl;
	
}

void LaplaceFMM2d::assign_mesh(){
	mesh=new Mesh_c;
	//assign vertex
	Location vloc("input.vertex[]");
	for(int i=0;i<in->size(vloc);i++){
		vloc.next();
		complex<double> r(in->d(vloc.sub("position.x")),in->d(vloc.sub("position.y")));
		mesh->vertex.push_back(new Vertex_c(r));
	}
	cout << "\n\tnumber of vertex: " << mesh->vertex.size() << endl;
	//assign face
	Location floc("input.edge[]");
	for(int j=0;j<in->size(floc);j++){
		floc.next();
		vector<int>& v=in->iarray(floc.sub("vertex[]"));
		int id0=in->getLocation("Vertex",v[0]).getIndex().get()[0];
		int id1=in->getLocation("Vertex",v[1]).getIndex().get()[0];
		//int id=in->i(floc.sub("id"));
		mesh->edge.push_back(new LaplaceEdge_c(mesh->vertex[id0],mesh->vertex[id1],j));
		//mesh->edge.push_back(new Edge_c(mesh->vertex[id0],mesh->vertex[id1],id));
	}
	cout << "\tnumber of edge: " << mesh->edge.size() << endl;
	mesh->calc_position_length_normal();
}

//void LaplaceFMM2d::assign_evaluation_point(){
//	Location eloc("input.evaluation_point[]");
//	for(int m=0;m<in->size(eloc);m++){
//		eloc.next();
//		evaluation_point.push_back(complex<double>(in->d(eloc.sub("x")),in->d(eloc.sub("y"))));
//	}
//}

void LaplaceFMM2d::assign_region_condition(){
	Location rloc("input.region_condition[]");
	for(int i=0;i<in->size(rloc);i++){
		rloc.next();
		string type=in->s(rloc.sub("type"));
		vector<int>& e=in->iarray(rloc.sub("edge[]"));
		vector<int>::iterator j=e.begin();
		for(;j!=e.end();j++){
			int id=in->getLocation("Edge",(*j)).getIndex().get()[0];
			if(type=="phi"){
				mesh->edge[id]->type="phi";
				((LaplaceEdge_c*)(mesh->edge[id]))->phi=in->d(rloc.sub("phi"));
				((LaplaceEdge_c*)(mesh->edge[id]))->q=0.0;
			}else if(type=="q"){
				mesh->edge[id]->type="q";
				((LaplaceEdge_c*)(mesh->edge[id]))->phi=0.0;
				((LaplaceEdge_c*)(mesh->edge[id]))->q=in->d(rloc.sub("q"));
			}
		}
	}
}

void LaplaceFMM2d::assign_initial_guess(){
	if(mesh->edge.size()!=in->size("output.edge[]"))return;
	Location loc("output.edge[]");
	vector<Edge_c*>::iterator e=mesh->edge.begin();
	for(;e!=mesh->edge.end();e++){
		loc.next();
		if((*e)->type=="phi"){
			((LaplaceEdge_c*)(*e))->q=in->d(loc.sub("q"));
		}else if((*e)->type=="q"){
			((LaplaceEdge_c*)(*e))->phi=in->d(loc.sub("phi"));
		}
	}
}

void LaplaceFMM2d::create_cell(){
	double size=in->d("input.system.size");
	complex<double> center(in->d("input.system.center.x"),in->d("input.system.center.y"));
	int pop=in->i("input.fmm.number_of_edges_per_leaf");
	int moment_number=in->i("input.fmm.number_of_terms_moments");
	int local_number=in->i("input.fmm.number_of_terms_local_expansions");
	tree=new LaplaceFMMTree_c(mesh->edge,size,center,pop,moment_number,local_number);
	//tree->update();
	//vector<vector<Cell_c*> >::iterator c=tree->generation_table.begin();
	//for(;c!=tree->generation_table.end();c++){
	//	vector<Cell_c*>::iterator d=c->begin();
	//	for(;d!=c->end();d++){
	//		(*d)->moment=new LaplaceMoment_c;
	//		(*d)->moment->resize(moment_number+1);
	//		(*d)->local=new LaplaceMoment_c;
	//		(*d)->local->resize(local_number+1);
	//	}
	//}
}

void LaplaceFMM2d::assign_bem(){
	int moment_number=in->i("input.fmm.number_of_terms_moments");
	int local_number=in->i("input.fmm.number_of_terms_local_expansions");
	int max_ite=in->i("input.bicgstab.max_iteration");
	double tolerance=in->d("input.bicgstab.tolerance");
	string temp=in->s("input.Gauss_Legendre_integral_points_number");
	int ipn=atoi(temp.c_str());
	string pre_condition=in->s("input.bicgstab.pre_condition");
	//string greenFunc=in->s("input.GreenFunction");
	double phi=in->d("input.external_field.phi");
	complex<double> gphi(in->d("input.external_field.gradient_phi.x"),
	                     in->d("input.external_field.gradient_phi.y"));
	bem=new LaplaceBEM_FMM2d(tree,mesh,moment_number,local_number,
	                         max_ite,tolerance,ipn,pre_condition,"free_analytic_integral",
	                         phi,gphi);
}

void LaplaceFMM2d::report_edge(){
	Location loc("output.edge[]");
	vector<Edge_c*>::iterator i=mesh->edge.begin();
	for(;i!=mesh->edge.end();i++){
		loc.next();
		out->put(loc.sub("phi"),((LaplaceEdge_c*)(*i))->phi);
		out->put(loc.sub("q"),((LaplaceEdge_c*)(*i))->q);
	}
}

void LaplaceFMM2d::report_evaluation_point(){
	Location iloc("input.evaluation_point[]");
	Location oloc("output.evaluation_point[]");
	int n=in->size(iloc);
	if(n==0)return;
	bem->bicgstab->init_evaluation_point();//calc local expansions
	for(int i=0;i<n;i++){
		iloc.next();oloc.next();
		const complex<double> r(in->d(iloc.sub("position.x")),in->d(iloc.sub("position.y")));
		double phi;/*complex<double> grad_phi;*/Cell_c* cell;
		bem->evaluate(r,phi,cell);
		out->put(oloc.sub("phi"),phi);
		//out->put(oloc.sub("grad_phi.x"),grad_phi.real());
		//out->put(oloc.sub("grad_phi.y"),grad_phi.imag());
		out->put(oloc.sub("cell"),cell->id);
	}
}

void  LaplaceFMM2d::report_cell(){
	vector<Cell_c*> cell_list;
	create_cell_vector(cell_list,tree->root);
	//report cell
	Location loc("output.cell[]");
	vector<Cell_c*>::iterator k=cell_list.begin();
	for(;k!=cell_list.end();k++){
		loc.next();
		//out->put(loc.sub("cell.id"),(*k)->id);
		out->put(loc.sub("generation"),(*k)->generation);
		out->put(loc.sub("size"),(*k)->size);
		complex<double> r=(*k)->position;
		out->put(loc.sub("position.x"),r.real());
		out->put(loc.sub("position.y"),r.imag());
		//children
		vector<int> ids;
		vector<Cell_c*>::iterator i=(*k)->child.begin();
		for(;i!=(*k)->child.end();i++){
			ids.push_back((*i)->id);
		}
		out->putArray(loc.sub("child[]"),ids);
		//adjacent
		ids.clear();
		vector<Cell_c*>::iterator a=(*k)->adjacent.begin();
		for(;a!=(*k)->adjacent.end();a++){
			ids.push_back((*a)->id);
		}
		out->putArray(loc.sub("adjacent[]"),ids);
		//interaction
		ids.clear();
		a=(*k)->interaction.begin();
		for(;a!=(*k)->interaction.end();a++){
			ids.push_back((*a)->id);
		}
		out->putArray(loc.sub("interaction[]"),ids);
		//edge
		ids.clear();
		vector<Edge_c*>::iterator j=(*k)->edge.begin();
		for(;j!=(*k)->edge.end();j++){
			ids.push_back((*j)->id);
		}
		out->putArray(loc.sub("edge[]"),ids);
	}
}

void LaplaceFMM2d::create_cell_vector(vector<Cell_c*>& cell_list,Cell_c* cell){
	cell_list.push_back(cell);
	vector<Cell_c*>::const_iterator i=cell->child.begin();
	for(;i!=cell->child.end();i++){
		create_cell_vector(cell_list,*i);
	}
}

