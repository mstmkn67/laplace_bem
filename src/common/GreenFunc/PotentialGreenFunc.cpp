#include "PotentialGreenFunc.h"

#define GPI (0.25/3.14159285358979)

namespace potential_green_func{

double get_free_G(const Vector3d& r,const Vector3d& rp){
	Vector3d dr=r-rp;
	double l=1./dr.length();
	return GPI*l;
}

Vector3d get_free_T(const Vector3d& r,const Vector3d& rp){
	Vector3d dr=r-rp;
	double l2=1./dr.length2();
	return GPI*l2*sqrt(l2)*dr;
}

double get_free_Tn(const Vector3d& r,const Vector3d& rp,const Vector3d& np){
	Vector3d dr=r-rp;
	double l2=1./dr.length2();
	return GPI*l2*sqrt(l2)*dr*np;
}

/////////////////////
double get_n_grad_free_G(const Vector3d& r0,const Vector3d& n0,const Vector3d& rp){
	Vector3d dr=r0-rp;
	double l2=1./dr.length2();
	return -GPI*l2*sqrt(l2)*(n0*dr);
}

Vector3d get_grad_free_G(const Vector3d& r0,const Vector3d& rp){
	Vector3d dr=r0-rp;
	double l2=1./dr.length2();
	return -GPI*l2*sqrt(l2)*dr;
}
/////////////////////
double get_free_F(const Vector3d& r0,const Vector3d& rp,const Vector3d& np){
	Vector3d dr=r0-rp;
	double l2=1./dr.length2();
	return GPI*l2*sqrt(l2)*dr*np;
}

double get_free_K(const Vector3d& r0,const Vector3d& rp,const Vector3d& n0){
	Vector3d dr=rp-r0;
	double l2=1./dr.length2();
	return GPI*l2*sqrt(l2)*dr*n0;
}

double get_free_H(const Vector3d& r0,const Vector3d& rp,const Vector3d& n0,const Vector3d& np){
	Vector3d dr=r0-rp;
	double l2=1./dr.length2();
	return GPI*l2*sqrt(l2)*(3.*(n0*dr)*(np*dr)*l2-n0*np);
}
//////////////////////
//
double get_zero_potential_on_wall_G(const Vector3d& r,const Vector3d& rp){
	Vector3d rpi=Vector3d(rp.x,rp.y,-rp.z);
	return get_free_G(r,rp)-get_free_G(r,rpi);
}

Vector3d get_zero_potential_on_wall_T(const Vector3d& r,const Vector3d& rp){
	Vector3d rpi=Vector3d(rp.x,rp.y,-rp.z);
	return get_free_T(r,rp)-get_free_T(r,rpi);
}

double get_zero_potential_on_wall_Tn(const Vector3d& r,const Vector3d& rp,const Vector3d& np){
	Vector3d rpi=Vector3d(rp.x,rp.y,-rp.z);
	return get_free_Tn(r,rp,np)-get_free_Tn(r,rpi,np);
}

double get_no_flux_on_wall_G(const Vector3d& r,const Vector3d& rp){
	Vector3d rpi=Vector3d(rp.x,rp.y,-rp.z);
	return get_free_G(r,rp)+get_free_G(r,rpi);
}

Vector3d get_no_flux_on_wall_T(const Vector3d& r,const Vector3d& rp){
	Vector3d rpi=Vector3d(rp.x,rp.y,-rp.z);
	return get_free_T(r,rp)+get_free_T(r,rpi);
}

double get_no_flux_on_wall_Tn(const Vector3d& r,const Vector3d& rp,const Vector3d& np){
	Vector3d rpi=Vector3d(rp.x,rp.y,-rp.z);
	return get_free_Tn(r,rp,np)+get_free_Tn(r,rpi,np);
}
//
double GreenFunc::get_G(const Vector3d& r0,const Vector3d& rp){
	return get_free_G(r0,rp);
}

Vector3d GreenFunc::get_T(const Vector3d& r0,const Vector3d& rp){
	return get_free_T(r0,rp);
}

double GreenFunc::get_Tn(const Vector3d& r0,const Vector3d& rp,const Vector3d& np){
	return get_free_Tn(r0,rp,np);
}
//
double GreenFuncZeroPotentialOnWall::get_G(const Vector3d& r0,const Vector3d& rp){
	return get_zero_potential_on_wall_G(r0,rp);
}

Vector3d GreenFuncZeroPotentialOnWall::get_T(const Vector3d& r0,const Vector3d& rp){
	return get_zero_potential_on_wall_T(r0,rp);
}

double GreenFuncZeroPotentialOnWall::get_Tn(const Vector3d& r0,const Vector3d& rp,const Vector3d& np){
	return get_zero_potential_on_wall_Tn(r0,rp,np);
}
//
double GreenFuncNoFluxOnWall::get_G(const Vector3d& r0,const Vector3d& rp){
	return get_no_flux_on_wall_G(r0,rp);
}

Vector3d GreenFuncNoFluxOnWall::get_T(const Vector3d& r0,const Vector3d& rp){
	return get_no_flux_on_wall_T(r0,rp);
}

double GreenFuncNoFluxOnWall::get_Tn(const Vector3d& r0,const Vector3d& rp,const Vector3d& np){
	return get_no_flux_on_wall_Tn(r0,rp,np);
}



}

