#ifndef _FMM_TREE_C_H_
#define _FMM_TREE_C_H_

#include "Mesh_c.h"
#include "Cell_c.h"
#include <vector>
#include <iostream>
using namespace std;

class Cell_c;

//Cell manager of system
//Square system is considered
//quadra tree is produced
// each cell has maximum four children.
//eldest is upper right, second upper left, next lower left, lower right  
//+----+
//|1  0|
//|2  3|
//+----+
class FMMTree_c{ 
public:
	FMMTree_c(vector<Edge_c*>& edge,
	          double size,
	          const complex<double>& position,
            int population_per_leaf);
	virtual ~FMMTree_c();
	virtual void update();
	virtual void update(double size,const complex<double>& position);
	
	virtual Cell_c* get_cell_at_evaluation_position(const complex<double>& r);
	
	double size;
	complex<double> position;
	vector<Edge_c*>& edge;
	Cell_c* root;
	vector<Cell_c*> leaf;
	//generation table
	vector<vector<Cell_c*> > generation_table;
protected:
	//virtual void create_children(vector<Edge_c*>& edge);
private:
	int generation_max;
	virtual void getLeaf(Cell_c* cell);
	virtual void creatGenerationTable();
	virtual void registerGenerationTable(Cell_c* cell);
};



#endif // _FMM_TREE_C_H_
