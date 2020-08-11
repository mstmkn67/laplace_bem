#ifndef _MESH_3D_H_
#define _MESH_3D_H_

#include "../Vector3d.h"
#include <vector>
#include <cmath>
#include <string>
using namespace std;

class Cell3d;

class Vertex3d{
public:
	Vertex3d(){}
	Vertex3d(const Vector3d& r){position=r;}
	virtual ~Vertex3d(){}
	Vector3d position;
};

class Face3d{
public:
	Face3d(){};
	Face3d(Vertex3d* v0,Vertex3d* v1,Vertex3d* v2,int i=-1,const string& t="")
	:id(i),type(t){vertex[0]=v0;vertex[1]=v1;vertex[2]=v2;}
	virtual ~Face3d(){}
	
	int id;
	Vertex3d* vertex[3];
	Vector3d position;
	Vector3d normal;
	double area;
	string type;

	
	Cell3d* cell;
	
	virtual void calc_position_area_normal();
private:
	static int total_id;
};

class Mesh3d{
public:
	Mesh3d(){};
	virtual ~Mesh3d(){};
	virtual void calc_position_area_normal();
	vector<Vertex3d*> vertex;
	vector<Face3d*> face;
};

#endif // _MESH_3D_H_
