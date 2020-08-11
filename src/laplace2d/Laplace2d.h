#ifndef _LAPLACE_2D_H_
#define _LAPLACE_2D_H_

#include "LaplaceBEM2d.h"
#include "../common/Timer.h"
#include "udfmanager.h"

class Laplace2d{
public:
	Laplace2d(UDFManager* in,UDFManager* out);
	virtual ~Laplace2d();
	
	virtual void run();
	
protected:
	virtual void assign_mesh();
	virtual void assign_region_condition();
	//virtual void assign_evaluation_point();
	virtual void assign_bem();
	virtual void report_edge();
	virtual void report_evaluation_point();
private:
	UDFManager* in;
	UDFManager* out;
	Mesh_c* mesh;
	//vector<complex<double> > evaluation_point;
	LaplaceBEM2d* bem;
	Timer timer;
};

#endif // _LAPLACE_2D_H_
