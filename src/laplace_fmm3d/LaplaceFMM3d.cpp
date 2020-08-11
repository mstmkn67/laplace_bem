#include "LaplaceFMM3d.h"

LaplaceFMM3d::LaplaceFMM3d(UDFManager* i,UDFManager* o):in(i),out(o){
	cout << "Simulation system is setting ....  " << flush;
	timer.start();
	assign_mesh();
	assign_region_condition();
	assign_initial_guess();
	create_cell();
	assign_bem();
	cout << " done (" << timer << " seconds )" << endl;
}

LaplaceFMM3d::~LaplaceFMM3d(){
	out->put("output.bicgstab.iteration",bem->bicgstab->max_iter);
	out->put("output.bicgstab.tolerance",bem->bicgstab->tol);
	delete bem;
	delete tree;
	vector<Face3d*>::iterator i=mesh->face.begin();
	for(;i!=mesh->face.end();i++)delete *i;
	vector<Vertex3d*>::iterator j=mesh->vertex.begin();
	for(;j!=mesh->vertex.end();j++)delete *j;
	double t=timer.get();
	out->put("output.time",t);
	out->write();
}

void LaplaceFMM3d::run(){
	cout << "Calculation of BEM-FMM .... "  << flush;
	bem->update();
	cout << " done (" << timer << " seconds )" << endl;
	cout << "Reporting the imformation of faces and cells .... " << flush;
	report_face();
	report_cell();
	cout << " done (" << timer << " seconds )" << endl;
	cout << "Evaluating the values at evaluation_point .... " << flush;
	report_evaluation_point();
	cout << " done (" << timer << " seconds )" << endl;
	
}

void LaplaceFMM3d::assign_mesh(){
	mesh=new Mesh3d;
	//assign vertex
	Location vloc("input.vertex[]");
	for(int i=0;i<in->size(vloc);i++){
		vloc.next();
		Vector3d r(in->d(vloc.sub("position.x")),in->d(vloc.sub("position.y")),in->d(vloc.sub("position.z")));
		mesh->vertex.push_back(new Vertex3d(r));
	}
	cout << "\n\tnumber of vertex: " << mesh->vertex.size() << endl;
	//assign face
	Location floc("input.face[]");
	for(int j=0;j<in->size(floc);j++){
		floc.next();
		vector<int>& v=in->iarray(floc.sub("vertex[]"));
		int id0=in->getLocation("Vertex",v[0]).getIndex().get()[0];
		int id1=in->getLocation("Vertex",v[1]).getIndex().get()[0];
		int id2=in->getLocation("Vertex",v[2]).getIndex().get()[0];
		//int id=in->i(floc.sub("id"));
		mesh->face.push_back(new LaplaceFace3d(mesh->vertex[id0],mesh->vertex[id1],mesh->vertex[id2],j));
	}
	cout << "\tnumber of face: " << mesh->face.size() << endl;
	mesh->calc_position_area_normal();
}

void LaplaceFMM3d::assign_region_condition(){
	Location rloc("input.region_condition[]");
	for(int i=0;i<in->size(rloc);i++){
		rloc.next();
		string type=in->s(rloc.sub("type"));
		vector<int>& e=in->iarray(rloc.sub("face[]"));
		vector<int>::iterator j=e.begin();
		for(;j!=e.end();j++){
			int id=in->getLocation("Face",(*j)).getIndex().get()[0];
			if(type=="phi"){
				mesh->face[id]->type="phi";
				((LaplaceFace3d*)(mesh->face[id]))->phi=in->d(rloc.sub("phi"));
				((LaplaceFace3d*)(mesh->face[id]))->q=0.0;
			}else if(type=="q"){
				mesh->face[id]->type="q";
				((LaplaceFace3d*)(mesh->face[id]))->phi=0.0;
				((LaplaceFace3d*)(mesh->face[id]))->q=in->d(rloc.sub("q"));
			}
		}
	}
}

void LaplaceFMM3d::assign_initial_guess(){
	if(mesh->face.size()!=in->size("output.face[]"))return;
	Location loc("output.face[]");
	vector<Face3d*>::iterator f=mesh->face.begin();
	for(;f!=mesh->face.end();f++){
		loc.next();
		if((*f)->type=="phi"){
			((LaplaceFace3d*)(*f))->q=in->d(loc.sub("q"));
		}else if((*f)->type=="q"){
			((LaplaceFace3d*)(*f))->phi=in->d(loc.sub("phi"));
		}
	}
}

void LaplaceFMM3d::create_cell(){
	double size=in->d("input.system.size");
	Vector3d center(in->d("input.system.center.x"),in->d("input.system.center.y"),in->d("input.system.center.z"));
	int pop=in->i("input.fmm.number_of_faces_per_leaf");
	int moment_number=in->i("input.fmm.number_of_terms_moments");
	int local_number=in->i("input.fmm.number_of_terms_local_expansions");
	tree=new LaplaceFMMTree3d(mesh->face,size,center,pop,moment_number,local_number);
	//tree->update();
	//vector<vector<Cell3d*> >::iterator c=tree->generation_table.begin();
	//for(;c!=tree->generation_table.end();c++){
	//	vector<Cell3d*>::iterator d=c->begin();
	//	for(;d!=c->end();d++){
	//		(*d)->moment=new LaplaceMoment3d;
	//		(*d)->moment->resize(moment_number);
	//		(*d)->local=new LaplaceMoment3d;
	//		(*d)->local->resize(local_number);
	//	}
	//}
}

void LaplaceFMM3d::assign_bem(){
	int moment_number=in->i("input.fmm.number_of_terms_moments");
	int local_number=in->i("input.fmm.number_of_terms_local_expansions");
	int max_ite=in->i("input.bicgstab.max_iteration");
	double tolerance=in->d("input.bicgstab.tolerance");
	string temp=in->s("input.Gauss_Legendre_integral_points_number");
	int ipn=atoi(temp.c_str());
	string pre_condition=in->s("input.bicgstab.pre_condition");
	//string green=in->s("input.GreenFunction");
	double phi=in->d("input.external_field.phi");
	Vector3d gphi(in->d("input.external_field.gradient_phi.x"),
	              in->d("input.external_field.gradient_phi.y"),
	              in->d("input.external_field.gradient_phi.z"));
	bem=new LaplaceBEM_FMM3d(tree,mesh,moment_number,local_number,
	                         max_ite,tolerance,ipn,pre_condition,"free",phi,gphi);
}

void LaplaceFMM3d::report_face(){
	Location loc("output.face[]");
	vector<Face3d*>::iterator i=mesh->face.begin();
	for(;i!=mesh->face.end();i++){
		loc.next();
		out->put(loc.sub("phi"),((LaplaceFace3d*)(*i))->phi);
		out->put(loc.sub("q"),((LaplaceFace3d*)(*i))->q);
	}
}

void LaplaceFMM3d::report_evaluation_point(){
	Location iloc("input.evaluation_point[]");
	Location oloc("output.evaluation_point[]");
	int n=in->size(iloc);
	if(n==0)return;
	bem->bicgstab->init_evaluation_point();//calc local expansions
	for(int i=0;i<n;i++){
		iloc.next();oloc.next();
		const Vector3d r(in->d(iloc.sub("position.x")),in->d(iloc.sub("position.y")),in->d(iloc.sub("position.z")));
		double phi;/*complex<double> grad_phi;*/Cell3d* cell;
		bem->evaluate(r,phi,cell);
		out->put(oloc.sub("phi"),phi);
		//out->put(oloc.sub("grad_phi.x"),grad_phi.real());
		//out->put(oloc.sub("grad_phi.y"),grad_phi.imag());
		out->put(oloc.sub("cell"),cell->id);
	}
}

void  LaplaceFMM3d::report_cell(){
	vector<Cell3d*> cell_list;
	create_cell_vector(cell_list,tree->root);
	//report cell
	Location loc("output.cell[]");
	vector<Cell3d*>::iterator k=cell_list.begin();
	for(;k!=cell_list.end();k++){
		loc.next();
		//out->put(loc.sub("cell.id"),(*k)->id);
		out->put(loc.sub("generation"),(*k)->generation);
		out->put(loc.sub("size"),(*k)->size);
		const Vector3d& r=(*k)->position;
		out->put(loc.sub("position.x"),r.x);
		out->put(loc.sub("position.y"),r.y);
		out->put(loc.sub("position.z"),r.z);
		//children
		vector<int> ids;
		vector<Cell3d*>::iterator i=(*k)->child.begin();
		for(;i!=(*k)->child.end();i++){
			ids.push_back((*i)->id);
		}
		out->putArray(loc.sub("child[]"),ids);
		//adjacent
		ids.clear();
		vector<Cell3d*>::iterator a=(*k)->adjacent.begin();
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
		vector<Face3d*>::iterator j=(*k)->face.begin();
		for(;j!=(*k)->face.end();j++){
			ids.push_back((*j)->id);
		}
		out->putArray(loc.sub("face[]"),ids);
	}
}

void LaplaceFMM3d::create_cell_vector(vector<Cell3d*>& cell_list,Cell3d* cell){
	cell_list.push_back(cell);
	vector<Cell3d*>::const_iterator i=cell->child.begin();
	for(;i!=cell->child.end();i++){
		create_cell_vector(cell_list,*i);
	}
}

