#ifndef _LAPLACE_BEM_3D_H_
#define _LAPLACE_BEM_3D_H_

#include "../common/mesh/LaplaceMesh3d.h"
#include "../common/GreenFunc/PotentialGreenFunc.h"
#include "../common/bicgstab/ilu_bicgstab.h"
#include "../common/bicgstab/jacobi_bicgstab.h"
#include "../common/bicgstab/non_bicgstab.h"
#include "../common/integration/LegendreGaussFormulaTriangle.h"
#include "../common/lapack/LapackFunctions.h"
#include <iostream>
using namespace std;

class LaplaceBEM3d{
public:
	LaplaceBEM3d(Mesh3d* mesh,int max_ite,double tolerance,
               int integral_point_number,const string& preconditioner,const string& greenFunc,
               double phi,const Vector3d& gphi);
	virtual ~LaplaceBEM3d();
	virtual void update();
	virtual void evaluate(const Vector3d& r,double& phi);
	BiCGSTAB* bicgstab;
protected:
	virtual void calc_bicgstab();
	virtual void calc_lapack();
	//
	virtual double calc_Aab(const Vector3d& r0,
	                        const Vector3d& ra,const Vector3d& rb,const Vector3d& rc,
	                        double area)const;
	virtual double calc_Aaa(const Vector3d& r0,const Vector3d& ra,const Vector3d& rb,const Vector3d& rc)const;
	virtual double calc_Bab(const Vector3d& r0,
	                        const Vector3d& ra,const Vector3d& rb,const Vector3d& rc,
	                        const Vector3d& n,double area)const;
	//
private:
	Mesh3d* mesh;
	ParameterOfLegendreGaussFormulaTiangle* integral;
	potential_green_func::GreenFunc* greenFunc;
	
	double phi;
	Vector3d gphi;
	
	int face_number;
	double *matrix_A;
	double *vector_b;
	double *vector_x;
};

#endif // _LAPLACE_BEM_3D_H_
