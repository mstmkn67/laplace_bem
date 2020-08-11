#ifndef _MESHC_H_
#define _MESHC_H_

#include <complex>
#include <vector>
#include <cmath>
#include <string>
using namespace std;

class Cell_c;

class Vertex_c{
public:
	Vertex_c(){}
	Vertex_c(const complex<double>& r){position=r;}
	virtual ~Vertex_c(){}
	complex<double> position;
};

class Edge_c{
public:
	Edge_c();
	Edge_c(Vertex_c* v0,Vertex_c* v1,int i=-1,const string& t="")
	:id(i),type(t){vertex[0]=v0;vertex[1]=v1;}
	virtual ~Edge_c(){}
	
	int id;
	Vertex_c* vertex[2];
	complex<double> position;
	complex<double> normal;
	double length;
	string type;
	
	Cell_c* cell;
	
	virtual void calc_position_length_normal();
private:
	static int total_id;
	
};

class Mesh_c{
public:
	Mesh_c(){};
	virtual ~Mesh_c(){};
	virtual void calc_position_length_normal();
	vector<Vertex_c*> vertex;
	vector<Edge_c*> edge;
};

#endif // _MESH_C_H_
