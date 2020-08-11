#include "PotentialGreenFunc_c.h"

#define GPI ((0.5)/(3.14159285358979))
vector<double> factorial;
#include <iostream>
namespace potential_green_func_c{

double get_free_G(const complex<double>& r0,const complex<double>& rp){
	return -GPI*log(sqrt(norm(r0-rp)));
}

complex<double> get_free_T(const complex<double>& r0,const complex<double>& rp){
	return conj(GPI/(r0-rp));
}

double get_free_Tn(const complex<double>& r0,const complex<double>& rp,const complex<double>& np){
	return ((GPI/(r0-rp))*np).real();
	//return (GPI/norm(r0-rp)*(r0-rp)*conj(np)).real();
}
//
double get_n_grad_free_G(const complex<double>& r0,const complex<double>& n0,const complex<double>& rp){
	return ((-GPI/(r0-rp))*n0).real();
}

complex<double> get_grad_free_G(const complex<double>& r0,const complex<double>& rp){
	return conj(-GPI/(r0-rp));
}

//
double get_zero_potential_on_wall_G(const complex<double>& r,const complex<double>& rp){
	return get_free_G(r,rp)-get_free_G(r,conj(rp));
}

complex<double> get_zero_potential_on_wall_T(const complex<double>& r,const complex<double>& rp){
	return get_free_T(r,rp)-get_free_T(r,conj(rp));
}

double get_zero_potential_on_wall_Tn(const complex<double>& r,const complex<double>& rp,const complex<double>& np){
	return get_free_Tn(r,rp,np)-get_free_Tn(r,conj(rp),np);
}
//

double get_no_flux_on_wall_G(const complex<double>& r,const complex<double>& rp){
	return get_free_G(r,rp)+get_free_G(r,conj(rp));
}

complex<double> get_no_flux_on_wall_T(const complex<double>& r,const complex<double>& rp){
	return get_free_T(r,rp)+get_free_T(r,conj(rp));
}

double get_no_flux_on_wall_Tn(const complex<double>& r,const complex<double>& rp,const complex<double>& np){
	return get_free_Tn(r,rp,np)+get_free_Tn(r,conj(rp),np);
}

//

complex<double> get_Ik(const complex<double>& z,int k){
	if(k>=0)return pow(z,k)/factorial[k];
	else return complex<double>(0.0,0.0);
}

complex<double> get_Ok(const complex<double>& z,int k){
	if(k>0){
		return factorial[k-1]/pow(z,k);
	}else if(k==0){
		return -log(z);
	}else{
		return complex<double>(0.0,0.0);
	}
}

void calc_factorial(int n){
	factorial.clear();
	factorial.resize(n+1);
	factorial[0]=1.0;
	for(int i=1;i<=n;i++){
		factorial[i]=i*factorial[i-1];
	}
}
//

double GreenFunc::get_G(const complex<double>& r0,const complex<double>& rp){
	return get_free_G(r0,rp);
}

complex<double> GreenFunc::get_T(const complex<double>& r0,const complex<double>& rp){
	return get_free_T(r0,rp);
}

double GreenFunc::get_Tn(const complex<double>& r0,const complex<double>& rp,const complex<double>& np){
	return get_free_Tn(r0,rp,np);
}
//
double GreenFuncZeroPotentialOnWall::get_G(const complex<double>& r0,const complex<double>& rp){
	return get_zero_potential_on_wall_G(r0,rp);
}

complex<double> GreenFuncZeroPotentialOnWall::get_T(const complex<double>& r0,const complex<double>& rp){
	return get_zero_potential_on_wall_T(r0,rp);
}

double GreenFuncZeroPotentialOnWall::get_Tn(const complex<double>& r0,const complex<double>& rp,const complex<double>& np){
	return get_zero_potential_on_wall_Tn(r0,rp,np);
}
//
double GreenFuncNoFluxOnWall::get_G(const complex<double>& r0,const complex<double>& rp){
	return get_no_flux_on_wall_G(r0,rp);
}

complex<double> GreenFuncNoFluxOnWall::get_T(const complex<double>& r0,const complex<double>& rp){
	return get_no_flux_on_wall_T(r0,rp);
}

double GreenFuncNoFluxOnWall::get_Tn(const complex<double>& r0,const complex<double>& rp,const complex<double>& np){
	return get_no_flux_on_wall_Tn(r0,rp,np);
}


}

#undef GPI
