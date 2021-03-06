#ifndef _LAPLACE_BEM_FMM_2D_H_
#define _LAPLACE_BEM_FMM_2D_H_

#include "../common/mesh/FMMTree3d.h"
#include "fmm_non_bicgstab.h"
#include "fmm_jacobi_bicgstab.h"
#include "fmm_block_bicgstab.h"
#include <set>
#include <iostream>
using namespace std;

class LaplaceBEM_FMM3d{
public:
	LaplaceBEM_FMM3d(FMMTree3d* tree,Mesh3d* mesh,
                   int moments_number,int local_number,
                   int max_ite,double tolerance,int integral_point_number,
                   const string& preconditioner,const string& greenFunc,
                   double phi,const Vector3d& gphi);
	virtual ~LaplaceBEM_FMM3d();
	virtual void update();
	virtual void evaluate(const Vector3d& r,double& phi,Cell3d*& cell);
	FMMNonBiCGSTAB* bicgstab;
protected:
private:
	FMMTree3d* tree;
	Mesh3d* mesh;
	int face_number;
	double* vector_b;
	double* vector_x;
	int moments_number;
	int local_number;
};

#endif // _LAPLACE_BEM_FMM_3D_H_
