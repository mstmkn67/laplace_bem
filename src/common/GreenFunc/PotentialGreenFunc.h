#ifndef _POTENTIAL_GREEN_FUNC_H_
#define _POTENTIAL_GREEN_FUNC_H_

#include "../Vector3d.h"
#include <complex>
#include <vector>
using namespace std;

namespace potential_green_func{
//r0 is the evaluation point, rp is the soruce position
double get_free_G(const Vector3d& r0,const Vector3d& rp);
Vector3d get_free_T(const Vector3d& r0,const Vector3d& rp);
double get_free_Tn(const Vector3d& r0,const Vector3d& rp,const Vector3d& np);
//////////////////////////////
//for indirect boundary element method
//diffrential operator acts on r0
double get_n_grad_free_G(const Vector3d& r0,const Vector3d& n0,const Vector3d& rp);
Vector3d get_grad_free_G(const Vector3d& r0,const Vector3d& rp);
///////////////////////////
double get_free_F(const Vector3d& r0,const Vector3d& rp,const Vector3d& np);//this is same to get_free_Tn
double get_free_K(const Vector3d& r0,const Vector3d& rp,const Vector3d& n0);
double get_free_H(const Vector3d& r0,const Vector3d& rp,const Vector3d& n0,const Vector3d& np);
////////////////////////

//The wall is located on z=0.
double get_zero_potential_on_wall_G(const Vector3d& r,const Vector3d& rp);
Vector3d get_zero_potential_on_wall_T(const Vector3d& r,const Vector3d& rp);
double get_zero_potential_on_wall_Tn(const Vector3d& r,const Vector3d& rp,const Vector3d& np);

double get_no_flux_on_wall_G(const Vector3d& r,const Vector3d& rp);
Vector3d get_no_flux_on_wall_T(const Vector3d& r,const Vector3d& rp);
double get_no_flux_on_wall_Tn(const Vector3d& r,const Vector3d& rp,const Vector3d& np);

//

class GreenFunc{
public:
	virtual double get_G(const Vector3d& r0,const Vector3d& rp);
	virtual Vector3d get_T(const Vector3d& r0,const Vector3d& rp);
	virtual double get_Tn(const Vector3d& r0,const Vector3d& rp,const Vector3d& np);
};

class GreenFuncZeroPotentialOnWall:public GreenFunc{
public:
	virtual double get_G(const Vector3d& r0,const Vector3d& rp);
	virtual Vector3d get_T(const Vector3d& r0,const Vector3d& rp);
	virtual double get_Tn(const Vector3d& r0,const Vector3d& rp,const Vector3d& np);
};

class GreenFuncNoFluxOnWall:public GreenFunc{
public:
	virtual double get_G(const Vector3d& r0,const Vector3d& rp);
	virtual Vector3d get_T(const Vector3d& r0,const Vector3d& rp);
	virtual double get_Tn(const Vector3d& r0,const Vector3d& rp,const Vector3d& np);
};

}


#endif // _POTENTIAL_GREEN_FUNC_H_
