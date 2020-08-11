#ifndef _LAPLACE_FMM_TREE_C_H_
#define _LAPLACE_FMM_TREE_C_H_

#include "FMMTree_c.h"

class LaplaceFMMTree_c:public FMMTree_c{
public:
	LaplaceFMMTree_c(vector<Edge_c*>& edge,
	                 double size,
	                 const complex<double>& position,
                   int population_per_leaf,
                   int moment_number,
                   int local_number);
	virtual ~LaplaceFMMTree_c();
	virtual void update();
	virtual void update(double size,const complex<double>& position);
private:
	int moment_number;
	int local_number;
};

#endif // _LAPLACE_FMM_TREE_C_H_
