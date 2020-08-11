#ifndef _LAPLACE_FMM_2D_H_
#define _LAPLACE_FMM_2D_H_

#include "../common/mesh/LaplaceFMMTree_c.h"
#include "LaplaceBEM_FMM2d.h"
#include "../common/Timer.h"
#include "udfmanager.h"

class LaplaceFMM2d{
public:
	LaplaceFMM2d(UDFManager* in,UDFManager* out);
	virtual ~LaplaceFMM2d();
	
	virtual void run();
	
protected:
	virtual void assign_mesh();
	virtual void assign_region_condition();
	virtual void assign_initial_guess();
	//virtual void assign_evaluation_point();
	virtual void create_cell();
	virtual void assign_bem();
	virtual void create_cell_vector(vector<Cell_c*>& cell_list,Cell_c* cell);
	virtual void report_edge();
	virtual void report_evaluation_point();
	virtual void report_cell();
private:
	UDFManager* in;
	UDFManager* out;
	FMMTree_c* tree;
	Mesh_c* mesh;
	//vector<complex<double> > evaluation_point;
	LaplaceBEM_FMM2d* bem;
	Timer timer;
};

#endif // _LAPLACE_FMM_2D_H_
