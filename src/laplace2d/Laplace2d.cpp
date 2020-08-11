#include "Laplace2d.h"

Laplace2d::Laplace2d(UDFManager* i,UDFManager* o)
:in(i),out(o){
	cout << "Simulation system is setting ....  ";
	timer.start();
	assign_mesh();
	assign_region_condition();
	assign_bem();
	cout << " done " << endl;
}

Laplace2d::~Laplace2d(){
	if(bem->bicgstab!=0){
		out->put("output.bicgstab.iteration",bem->bicgstab->max_iter);
		out->put("output.bicgstab.tolerance",bem->bicgstab->tol);
	}
	delete bem;
	vector<Edge_c*>::iterator i=mesh->edge.begin();
	for(;i!=mesh->edge.end();i++)delete *i;
	vector<Vertex_c*>::iterator j=mesh->vertex.begin();
	for(;j!=mesh->vertex.end();j++)delete *j;
	delete mesh;
	double t=timer.get();
	out->put("output.time",t);
	out->write();
}

void Laplace2d::run(){
	cout << "Calculation of BEM .... ";
	bem->update();
	cout << " done " << endl;
	cout << "Reporting the imformation of edges .... ";
	report_edge();
	cout << " done " << endl;
	cout << "Evaluating the values at evaluation_point .... ";
	report_evaluation_point();
	cout << " done " << endl;
}

void Laplace2d::assign_mesh(){
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
		int id=in->i(floc.sub("id"));
		mesh->edge.push_back(new LaplaceEdge_c(mesh->vertex[id0],mesh->vertex[id1],id));
	}
	cout << "\tnumber of edge: " << mesh->edge.size() << endl;
	mesh->calc_position_length_normal();
}

void Laplace2d::assign_region_condition(){
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

void Laplace2d::assign_bem(){
	int max_ite=in->i("input.bicgstab.max_iteration");
	double tolerance=in->d("input.bicgstab.tolerance");
	string pre_condition=in->s("input.bicgstab.pre_condition");
	string temp=in->s("input.Gauss_Legendre_integral_points_number");
	//string green=in->s("input.GreenFunction");
	int ipn=atoi(temp.c_str());
	//
	Location loc("input.additional_information[]");
	for(int i=0;i<in->size(loc);i++){
		loc.next();
		string item=in->s(loc.sub("item"));
		if(item=="SOLVER"){
			if(in->s(loc.sub("value"))=="LU"){
				pre_condition="LU";
			}
		}
	}
	double phi=in->d("input.external_field.phi");
	complex<double> gphi(in->d("input.external_field.gradient_phi.x"),
	                     in->d("input.external_field.gradient_phi.y"));
	//
	bem=new LaplaceBEM2d(mesh,max_ite,tolerance,ipn,pre_condition,"free_analytic_integral",phi,gphi);
}

void Laplace2d::report_edge(){
	Location loc("output.edge[]");
	vector<Edge_c*>::iterator i=mesh->edge.begin();
	for(;i!=mesh->edge.end();i++){
		loc.next();
		out->put(loc.sub("phi"),((LaplaceEdge_c*)(*i))->phi);
		out->put(loc.sub("q"),((LaplaceEdge_c*)(*i))->q);
	}
}

void Laplace2d::report_evaluation_point(){
	Location iloc("input.evaluation_point[]");
	Location oloc("output.evaluation_point[]");
	for(int i=0;i<in->size(iloc);i++){
		iloc.next();oloc.next();
		const complex<double> r(in->d(iloc.sub("position.x")),in->d(iloc.sub("position.y")));
		double phi;//complex<double> grad_phi;
		bem->evaluate(r,phi);
		out->put(oloc.sub("phi"),phi);
		//out->put(oloc.sub("grad_phi.x"),grad_phi.real());
		//out->put(oloc.sub("grad_phi.y"),grad_phi.imag());
	}
}
