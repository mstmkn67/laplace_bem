#ifndef _LAPLACE_FMM_3D_H_
#define _LAPLACE_FMM_3D_H_

#include "../common/mesh/LaplaceFMMTree3d.h"
#include "LaplaceBEM_FMM3d.h"
#include "../common/Timer.h"
#include "udfmanager.h"

class LaplaceFMM3d{
public:
	LaplaceFMM3d(UDFManager* in,UDFManager* out);
	virtual ~LaplaceFMM3d();
	
	virtual void run();
	
protected:
	virtual void assign_mesh();
	virtual void assign_region_condition();
	virtual void assign_initial_guess();
	//virtual void assign_evaluation_point();
	virtual void create_cell();
	virtual void assign_bem();
	virtual void create_cell_vector(vector<Cell3d*>& cell_list,Cell3d* cell);
	virtual void report_face();
	virtual void report_evaluation_point();
	virtual void report_cell();
private:
	UDFManager* in;
	UDFManager* out;
	FMMTree3d* tree;
	Mesh3d* mesh;
	//vector<complex<double> > evaluation_point;
	LaplaceBEM_FMM3d* bem;
	Timer timer;
};

#endif // _LAPLACE_FMM_3D_H_
