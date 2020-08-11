#ifndef _LAPLACE_BEM_FMM_2D_H_
#define _LAPLACE_BEM_FMM_2D_H_

#include "../common/mesh/FMMTree_c.h"
#include "../common/mesh/Moment_c.h"
#include "fmm_non_bicgstab.h"
#include "fmm_jacobi_bicgstab.h"
#include "fmm_block_bicgstab.h"
#include <set>
#include <iostream>
using namespace std;

class LaplaceBEM_FMM2d{
public:
	LaplaceBEM_FMM2d(FMMTree_c* tree,Mesh_c* mesh,
                   int moments_number,int local_number,
                   int max_ite,double tolerance,int integral_point_number,
                   const string& preconditioner,const string& greenFunc,
                   double phi,const complex<double>& gphi);
	virtual ~LaplaceBEM_FMM2d();
	virtual void update();
	virtual void evaluate(const complex<double>& r,
	                      double& phi,
	                      Cell_c*& cell);
	FMMNonBiCGSTAB* bicgstab;
protected:
private:
	FMMTree_c* tree;
	Mesh_c* mesh;
	int edge_number;
	double* vector_b;
	double* vector_x;
	int moments_number;
	int local_number;
};

#endif // _LAPLACE_BEM_FMM_2D_H_
