#ifndef _LAPLACE_3D_H_
#define _LAPLACE_3D_H_

#include "LaplaceBEM3d.h"
#include "../common/Timer.h"
#include "udfmanager.h"

class Laplace3d{
public:
	Laplace3d(UDFManager* in,UDFManager* out);
	virtual ~Laplace3d();
	
	virtual void run();
	
protected:
	virtual void assign_mesh();
	virtual void assign_region_condition();
	//virtual void assign_evaluation_point();
	virtual void assign_bem();
	virtual void report_face();
	virtual void report_evaluation_point();
private:
	UDFManager* in;
	UDFManager* out;
	Mesh3d* mesh;
	//vector<complex<double> > evaluation_point;
	LaplaceBEM3d* bem;
	Timer timer;
};

#endif // _Laplace_3D_H_
