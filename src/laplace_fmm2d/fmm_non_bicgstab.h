#ifndef _FMM_NON_BICGSTAB_H_
#define _FMM_NON_BICGSTAB_H_

#include "../common/bicgstab/bicgstab.h"
#include "../common/mesh/LaplaceMesh_c.h"
#include "../common/mesh/Cell_c.h"
#include "../common/mesh/FMMTree_c.h"
#include "../common/GreenFunc/PotentialGreenFunc_c.h"
#include "../common/integration/LegendreGaussFormula1d.h"

class FMMNonBiCGSTAB:public BiCGSTAB{
public:
	FMMNonBiCGSTAB(int n,FMMTree_c* tree,Mesh_c* mesh,double* x,double* b,
                 int moments_number,int local_number,int max_ite,double tol,
                 int integral_point_num,const string& greenFunc,
                 double phi,const complex<double>& gphi);
	virtual ~FMMNonBiCGSTAB();
	
	virtual int update();
	
	virtual void matvec(double alpha,double* xx,double beta,double* yy);
	virtual void psolve(double* xx,double* bb);
	virtual void init_evaluation_point();
	virtual double calc_evaluation_point(Cell_c* cell,const complex<double>& r);
protected:
	virtual void calc_vector_b();
	//virtual double get_matrix(int i,int j);
	virtual double calc_Aab(const complex<double>& r0,const complex<double>& ra,
	                            const complex<double>& rb,double length)const;
	virtual double calc_Bab(const complex<double>& r0,const complex<double>& ra,
	                          const complex<double>& rb,const complex<double>& n,
	                          double length)const;
	virtual double calc_Aaa(const complex<double>& r0,const complex<double>& ra,
	                        const complex<double>& rb,double length)const;
	//virtual double calc_Aaa(double length)const;
	
	virtual double calc_edge_b(Edge_c* ei,Cell_c* cell);
	virtual double calc_edge_x(Edge_c* ei,Cell_c* cell);
	virtual double calc_edge_evaluation_point(const complex<double>& ri,Cell_c* cell);
	virtual void calc_moment_leaf_b();
	virtual void calc_moment_leaf_x();
	virtual void calc_moment_leaf_evaluation_point();
	virtual void calc_moment_M2M();//upward procedure
	virtual void calc_local_expansion_2nd_generation();
	virtual void calc_local_expansion();//downward procedure

	FMMTree_c* tree;
	Mesh_c* mesh;
	ParameterOfLegendreGaussFormula1d* integral;
	potential_green_func_c::GreenFunc* greenFunc;
	
	double phi;
	complex<double> gphi;
private:
	int moments_number;
	int local_number;
};

#endif // _FMM_NON_BICGSTAB_H_
