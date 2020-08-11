#ifndef _LAPLACE_FMM_TREE_3D_H_
#define _LAPLACE_FMM_TREE_3D_H_

#include "FMMTree3d.h"

class LaplaceFMMTree3d:public FMMTree3d{
public:
	LaplaceFMMTree3d(vector<Face3d*>& face,
	                 double size,
	                 const Vector3d& position,
                   int population_per_leaf,
                   int moment_number,
                   int local_number);
	virtual ~LaplaceFMMTree3d();
	virtual void update();
	virtual void update(double size,const Vector3d& position);
private:
	int moment_number;
	int local_number;
};

#endif // _LAPLACE_FMM_TREE_3D_H_
