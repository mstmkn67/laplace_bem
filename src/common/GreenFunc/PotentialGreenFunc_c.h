#ifndef _POTENTIAL_GREEN_FUNC_C_H_
#define _POTENTIAL_GREEN_FUNC_C_H_

#include <complex>
#include <vector>
using namespace std;

namespace potential_green_func_c{

double get_free_G(const complex<double>& r0,const complex<double>& rp);
complex<double> get_free_T(const complex<double>& r0,const complex<double>& rp);
double get_free_Tn(const complex<double>& r0,const complex<double>& rp,const complex<double>& np);
//////////////////////////////
//for indirect boundary element method
//diffrential operator acts on r0
double get_n_grad_free_G(const complex<double>& r0,const complex<double>& n0,const complex<double>& rp);
complex<double> get_grad_free_G(const complex<double>& r0,const complex<double>& rp);

//The wall is located on y=0.
double get_zero_potential_on_wall_G(const complex<double>& r,const complex<double>& rp);
complex<double> get_zero_potential_on_wall_T(const complex<double>& r,const complex<double>& rp);
double get_zero_potential_on_wall_Tn(const complex<double>& r,const complex<double>& rp,const complex<double>& np);

double get_no_flux_on_wall_G(const complex<double>& r,const complex<double>& rp);
complex<double> get_no_flux_on_wall_T(const complex<double>& r,const complex<double>& rp);
double get_no_flux_on_wall_T(const complex<double>& r,const complex<double>& rp,const complex<double>& np);

//for multipole moment

complex<double> get_Ik(const complex<double>& z,int k);
complex<double> get_Ok(const complex<double>& z,int k);

void calc_factorial(int n);

//
class GreenFunc{
public:
	virtual double get_G(const complex<double>& r0,const complex<double>& rp);
	virtual complex<double> get_T(const complex<double>& r0,const complex<double>& rp);
	virtual double get_Tn(const complex<double>& r0,const complex<double>& rp,const complex<double>& np);
};

class GreenFuncZeroPotentialOnWall:public GreenFunc{
public:
	virtual double get_G(const complex<double>& r0,const complex<double>& rp);
	virtual complex<double> get_T(const complex<double>& r0,const complex<double>& rp);
	virtual double get_Tn(const complex<double>& r0,const complex<double>& rp,const complex<double>& np);
};

class GreenFuncNoFluxOnWall:public GreenFunc{
public:
	virtual double get_G(const complex<double>& r0,const complex<double>& rp);
	virtual complex<double> get_T(const complex<double>& r0,const complex<double>& rp);
	virtual double get_Tn(const complex<double>& r0,const complex<double>& rp,const complex<double>& np);
};

}


#endif // _POTENTIAL_GREEN_FUNC_C_H_
