#ifndef _LAPLACE_MESH_3D_H_
#define _LAPLACE_MESH_3D_H_

#include "Mesh3d.h"

class LaplaceFace3d:public Face3d{
public:
	LaplaceFace3d(){};
	LaplaceFace3d(Vertex3d* v0,Vertex3d* v1,Vertex3d* v2,int i=-1)
	:Face3d(v0,v1,v2,i,"phi"),q(0.0),phi(0.0){}
	virtual ~LaplaceFace3d(){};
	double q;
	double phi;

};

#endif // _LAPLACE_MESH_3D_H_
