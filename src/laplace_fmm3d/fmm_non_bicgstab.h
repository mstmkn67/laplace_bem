#ifndef _FMM_NON_BICGSTAB_H_
#define _FMM_NON_BICGSTAB_H_

#include "../common/bicgstab/bicgstab.h"
#include "../common/mesh/FMMTree3d.h"
#include "../common/mesh/Moment3d.h"
#include "../common/mesh/LaplaceMesh3d.h"
#include "../common/GreenFunc/PotentialGreenFunc.h"
#include "../common/integration/LegendreGaussFormulaTriangle.h"

class FMMNonBiCGSTAB:public BiCGSTAB{
public:
	FMMNonBiCGSTAB(int n,FMMTree3d* tree,Mesh3d* mesh,double* x,double* b,
                 int moments_number,int local_number,int max_ite,double tol,
                 int integral_point_num,const string& green,
                 double phi,const Vector3d& gphi);
	virtual ~FMMNonBiCGSTAB();
	
	virtual int update();
	
	virtual void matvec(double alpha,double* xx,double beta,double* yy);
	virtual void psolve(double* xx,double* bb);
	virtual void init_evaluation_point();
	virtual double calc_evaluation_point(Cell3d* cell,const Vector3d& r);
protected:
	virtual void calc_vector_b();
	//virtual double get_matrix(int i,int j);
	virtual double calc_Aab(const Vector3d& r0,const Vector3d& ra,const Vector3d& rb,const Vector3d& rc,
	                        double area)const;
	virtual double calc_Bab(const Vector3d& r0,const Vector3d& ra,const Vector3d& rb,const Vector3d& rc,
	                        const Vector3d& n,double area)const;
	virtual double calc_Aaa(const Vector3d& r0,const Vector3d& ra,const Vector3d& rb,const Vector3d& rc)const;
	
	virtual double calc_face_b(Face3d* fi,Cell3d* cell);
	virtual double calc_face_x(Face3d* fi,Cell3d* cell);
	virtual double calc_face_evaluation_point(const Vector3d& ri,Cell3d* cell);
	virtual void calc_moment_leaf_b();
	virtual void calc_moment_leaf_x();
	virtual void calc_moment_leaf_evaluation_point();
	virtual void calc_moment_M2M();//upward procedure
	virtual void calc_local_expansion_2nd_generation();
	virtual void calc_local_expansion();//downward procedure

	FMMTree3d* tree;
	Mesh3d* mesh;
	ParameterOfLegendreGaussFormulaTiangle* integral;
	potential_green_func::GreenFunc* greenFunc;
	Rnm rnm;
	Snm snm;
	
	double phi;
	Vector3d gphi;
private:
	int moments_number;
	int local_number;
};

#endif // _FMM_NON_BICGSTAB_H_
