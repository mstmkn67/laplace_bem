#ifndef _LAPLACE_BEM_2D_H_
#define _LAPLACE_BEM_2D_H_

#include "../common/mesh/LaplaceMesh_c.h"
#include "../common/GreenFunc/PotentialGreenFunc_c.h"
#include "../common/bicgstab/ilu_bicgstab.h"
#include "../common/bicgstab/jacobi_bicgstab.h"
#include "../common/bicgstab/non_bicgstab.h"
#include "../common/integration/LegendreGaussFormula1d.h"
#include "../common/lapack/LapackFunctions.h"
#include <iostream>
using namespace std;

class LaplaceBEM2d{
public:
	LaplaceBEM2d(Mesh_c* mesh,int max_ite,double tolerance,
               int integral_point_number,const string& preconditioner,const string& greenFunc,
               double phi,const complex<double>& gphi);
	virtual ~LaplaceBEM2d();
	virtual void update();
	virtual void evaluate(const complex<double>& r,double& phi);
	BiCGSTAB* bicgstab;
protected:
	virtual void calc_bicgstab();
	virtual void calc_lapack();
	//
	virtual double calc_q_Aab(const complex<double>& r0,const complex<double>& ra,
	                          const complex<double>& rb,const complex<double>& n,
	                          double length)const;
	virtual double calc_phi_Aab(const complex<double>& r0,const complex<double>& ra,
	                            const complex<double>& rb,double length)const;
	virtual double calc_phi_Aaa(const complex<double>& r0,const complex<double>& ra,
	                            const complex<double>& rb,double length)const;
	//
private:
	Mesh_c* mesh;
	ParameterOfLegendreGaussFormula1d* integral;
	potential_green_func_c::GreenFunc* greenFunc;
	
	double phi;
	complex<double> gphi;
	
	int edge_number;
	double *matrix_A;
	double *vector_b;
	double *vector_x;
};

#endif // _LAPLACE_BEM_2D_H_
