#include "StokesGreenFunc_c.h"

#define GPI (0.25/3.14159285358979)
#define TPI (1.0/3.14159285358979)

namespace stokes_green_func_c{

Tensor2x2_c delta(1.0,0.0,
                  0.0,1.0);

Tensor2x2_c get_free_G(const complex<double>& r0,const complex<double>& rp){
	complex<double> dr=rp-r0;
	double l2=norm(dr);
	return GPI*(-(log(sqrt(l2))+0.5)*delta+dyad(dr,dr)/l2);
}

Tensor2x2x2_c get_free_T(const complex<double>& r0,const complex<double>& rp){
	complex<double> dr=rp-r0;
	double l2=norm(dr);
	return -TPI*dyad(dr,dyad(dr,dr))/l2/l2;
}
Tensor2x2_c get_free_Tn(const complex<double>& r0,const complex<double>& rp,
                        const complex<double>& np){
	complex<double> dr=rp-r0;
	double l2=norm(dr);
	return -TPI*(dr.real()*np.real()+dr.imag()*np.imag())/l2/l2*dyad(dr,dr);
}

//
Tensor2x2_c get_free_Kn(const complex<double>& r0,const complex<double>& rp,
	                      const complex<double>& n0){
	complex<double> dr=r0-rp;
	double l2=norm(dr);
	return TPI*(dr.real()*n0.real()+dr.imag()*n0.imag())/l2/l2*dyad(dr,dr);
}

Tensor2x2_c get_free_Hnn(const complex<double>& r0,const complex<double>& rp,
	                       const complex<double>& n0,const complex<double>& np){
	complex<double> dr=r0-rp;
	double l2=norm(dr);
	double rn0=dr.real()*n0.real()+dr.imag()*n0.imag();
	double rnp=dr.real()*np.real()+dr.imag()*np.imag();
	double n0np=n0.real()*np.real()+n0.imag()*np.imag();
	Tensor2x2_c a=((rn0*rnp)/l2)*(delta-8*dyad(dr,dr)/l2)+(rnp/l2)*dyad(dr,n0);
	a+=(rn0/l2)*dyad(np,dr)+(n0np/l2)*dyad(dr,dr)+dyad(n0,np);
	return (TPI/l2)*a;
}

}

#undef GPI
#undef TPI
