#include "Mesh_c.h"

void Edge_c::calc_position_length_normal(){
	const complex<double>& r0=vertex[0]->position;
	const complex<double>& r1=vertex[1]->position;
	position=(r0+r1)/2.0;
	normal=r1-r0;
	length=sqrt(norm(normal));
	normal/=length;
	normal=complex<double>(-normal.imag(),normal.real());
}

void Mesh_c::calc_position_length_normal(){
	vector<Edge_c*>::iterator i=edge.begin();
	for(;i!=edge.end();i++){
		(*i)->calc_position_length_normal();
	}
}
