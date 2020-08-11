#include "StokesGreenFunc.h"

#define GPI (0.125/3.14159285358979)
#define TPI (0.75 /3.14159285358979)

namespace stokes_green_func{

Tensor3x3 delta(1.0,0.0,0.0,
                0.0,1.0,0.0,
                0.0,0.0,1.0);

/////////////////////////////////////////////////////////////
Tensor3x3 get_free_G(const Vector3d& x,const Vector3d& y){
	Vector3d dr=y-x;
	double l2=1./dr.length2();
	return GPI*(delta+dyad(dr,dr)*l2)*sqrt(l2);
}

Tensor3x3x3 get_free_T(const Vector3d& x,const Vector3d& y){
	Vector3d dr=y-x;
	double l2=1./dr.length2();
	return -TPI*dyad(dr,dyad(dr,dr))*l2*l2*sqrt(l2);
}

Tensor3x3 get_free_T(const Vector3d& x,const Vector3d& y,const Vector3d& ny){
	Vector3d dr=y-x;
	double l2=1./dr.length2();
	return -TPI*(dr*ny)*dyad(dr,dr)*l2*l2*sqrt(l2);
}

Vector3d get_free_G(const Vector3d& x,const Vector3d& y,const Vector3d& fy){
	Vector3d dr=y-x;
	double l2=1./dr.length2();
	return GPI*sqrt(l2)*(fy+l2*(dr*fy)*dr);
}

Vector3d get_free_T(const Vector3d& x,const Vector3d& y,const Vector3d& ny,const Vector3d& fy){
	Vector3d dr=y-x;
	double l2=1./dr.length2();
	return -TPI*l2*l2*sqrt(l2)*(dr*ny)*(dr*fy)*dr;
}


}

#undef GPI
#undef TPI
