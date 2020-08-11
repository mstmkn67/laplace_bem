#include "Mesh3d.h"

void Face3d::calc_position_area_normal(){
	const Vector3d& r0=vertex[0]->position;
	const Vector3d& r1=vertex[1]->position;
	const Vector3d& r2=vertex[2]->position;
	position=(r0+r1+r2)/3.0;
	normal=((r1-r0)^(r2-r0));
	area=0.5*normal.length();
	normal/=normal.length();
}

void Mesh3d::calc_position_area_normal(){
	vector<Face3d*>::iterator i=face.begin();
	for(;i!=face.end();i++){
		(*i)->calc_position_area_normal();
	}
}
