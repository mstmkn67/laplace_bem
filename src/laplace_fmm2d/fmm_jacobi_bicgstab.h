#ifndef _FMM_JACOBI_BICGSTAB_H_
#define _FMM_JACOBI_BICGSTAB_H_

#include "fmm_non_bicgstab.h"

class FMMJacobiBiCGSTAB:public FMMNonBiCGSTAB{
public:
	FMMJacobiBiCGSTAB(int n,FMMTree_c* cell,Mesh_c* mesh,double* x,double* b,
                    int moments_number,int local_number,int max_ite,double tol,
                    int integral_point_num,const string& greenFunc,
                    double phi,const complex<double>& gphi);
	virtual ~FMMJacobiBiCGSTAB();
	virtual void psolve(double* xx,double* bb);
private:
};

#endif // _FMM_JACOBI_BICGSTAB_H_
