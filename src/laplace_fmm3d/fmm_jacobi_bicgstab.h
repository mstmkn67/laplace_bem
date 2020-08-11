#ifndef _FMM_JACOBI_BICGSTAB_H_
#define _FMM_JACOBI_BICGSTAB_H_

#include "fmm_non_bicgstab.h"

class FMMJacobiBiCGSTAB:public FMMNonBiCGSTAB{
public:
	FMMJacobiBiCGSTAB(int n,FMMTree3d* cell,Mesh3d* mesh,double* x,double* b,
                    int moments_number,int local_number,int max_ite,double tol,
                    int integral_point_num,const string& green,
                    double phi,const Vector3d& gphi);
	virtual ~FMMJacobiBiCGSTAB();
	virtual void psolve(double* xx,double* bb);
private:
};

#endif // _FMM_JACOBI_BICGSTAB_H_
