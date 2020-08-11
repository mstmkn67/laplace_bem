#include "Laplace3d.h"

Laplace3d::Laplace3d(UDFManager* i,UDFManager* o)
:in(i),out(o){
	cout << "Simulation system is setting ....  ";
	timer.start();
	assign_mesh();
	assign_region_condition();
	assign_bem();
	cout << " done " << endl;
}

Laplace3d::~Laplace3d(){
	if(bem->bicgstab!=0){
		out->put("output.bicgstab.iteration",bem->bicgstab->max_iter);
		out->put("output.bicgstab.tolerance",bem->bicgstab->tol);
	}
	delete bem;
	vector<Face3d*>::iterator i=mesh->face.begin();
	for(;i!=mesh->face.end();i++)delete *i;
	vector<Vertex3d*>::iterator j=mesh->vertex.begin();
	for(;j!=mesh->vertex.end();j++)delete *j;
	delete mesh;
	double t=timer.get();
	out->put("output.time",t);
	out->write();
}

void Laplace3d::run(){
	cout << "Calculation of BEM .... ";
	bem->update();
	cout << " done " << endl;
	cout << "Reporting the imformation of faces .... ";
	report_face();
	cout << " done " << endl;
	cout << "Evaluating the values at evaluation_point .... ";
	report_evaluation_point();
	cout << " done " << endl;
}

void Laplace3d::assign_mesh(){
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
		int id=in->i(floc.sub("id"));
		mesh->face.push_back(new LaplaceFace3d(mesh->vertex[id0],mesh->vertex[id1],mesh->vertex[id2],id));
	}
	cout << "\tnumber of face: " << mesh->face.size() << endl;
	mesh->calc_position_area_normal();
}

void Laplace3d::assign_region_condition(){
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

void Laplace3d::assign_bem(){
	int max_ite=in->i("input.bicgstab.max_iteration");
	double tolerance=in->d("input.bicgstab.tolerance");
	string pre_condition=in->s("input.bicgstab.pre_condition");
	string temp=in->s("input.Gauss_Legendre_integral_points_number");
	int ipn=atoi(temp.c_str());
	//string green=in->s("input.GreenFunction");
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
	Vector3d gphi(in->d("input.external_field.gradient_phi.x"),
	              in->d("input.external_field.gradient_phi.y"),
	              in->d("input.external_field.gradient_phi.z"));
	//
	bem=new LaplaceBEM3d(mesh,max_ite,tolerance,ipn,pre_condition,"free",phi,gphi);
}

void Laplace3d::report_face(){
	Location loc("output.face[]");
	vector<Face3d*>::iterator i=mesh->face.begin();
	for(;i!=mesh->face.end();i++){
		loc.next();
		out->put(loc.sub("phi"),((LaplaceFace3d*)(*i))->phi);
		out->put(loc.sub("q"),((LaplaceFace3d*)(*i))->q);
	}
}

void Laplace3d::report_evaluation_point(){
	Location iloc("input.evaluation_point[]");
	Location oloc("output.evaluation_point[]");
	for(int i=0;i<in->size(iloc);i++){
		iloc.next();oloc.next();
		const Vector3d r(in->d(iloc.sub("position.x")),in->d(iloc.sub("position.y")),in->d(iloc.sub("position.z")));
		double phi;//complex<double> grad_phi;
		bem->evaluate(r,phi);
		out->put(oloc.sub("phi"),phi);
	}
}
