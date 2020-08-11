#ifndef _STOKES_GREEN_FUNC_C_H_
#define _STOKES_GREEN_FUNC_C_H_

#include <complex>
using namespace std;
#include "../Tensor2x2x2_c.h"

namespace stokes_green_func_c{

Tensor2x2_c get_free_G(const complex<double>& r0,const complex<double>& rp);
Tensor2x2x2_c get_free_T(const complex<double>& r0,const complex<double>& rp);
Tensor2x2_c get_free_Tn(const complex<double>& r0,const complex<double>& rp,
                        const complex<double>& np);
//
Tensor2x2_c get_free_Kn(const complex<double>& r0,const complex<double>& rp,
	                      const complex<double>& n0);
Tensor2x2_c get_free_Hnn(const complex<double>& r0,const complex<double>& rp,
	                       const complex<double>& n0,const complex<double>& np);
//


}

#endif // _STOKES_GREEN_FUNC_C_H_
