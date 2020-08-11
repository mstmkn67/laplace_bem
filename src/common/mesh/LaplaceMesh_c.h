#ifndef _LAPLACE_MESH_C_H_
#define _LAPLACE_MESH_C_H_

#include "Mesh_c.h"

class LaplaceEdge_c:public Edge_c{
public:
	LaplaceEdge_c():q(0.0),phi(0.0){};
	LaplaceEdge_c(Vertex_c* v0,Vertex_c* v1,int i=-1)
	:Edge_c(v0,v1,i,"phi"),q(0.0),phi(0.0){}
	virtual ~LaplaceEdge_c(){};
	double q;
	double phi;

};

#endif // _LAPLACE_MESH_C_H_
