#ifndef _FMM_TREE_3D_H_
#define _FMM_TREE_3D_H_

#include "Mesh3d.h"
#include "Cell3d.h"
#include <vector>
#include <iostream>
using namespace std;

class Cell3d;

class FMMTree3d{
public:
	FMMTree3d(vector<Face3d*>& face,
	          double size,
	          const Vector3d& position,
            int population_per_leaf);
	virtual ~FMMTree3d();
	virtual void update();
	virtual void update(double size,const Vector3d& position);
	virtual Cell3d* get_cell_at_evaluation_position(const Vector3d& r);
	double size;
	Vector3d position;
	Cell3d* root;
	vector<Cell3d*> leaf;
	//generation table
	vector<vector<Cell3d*> > generation_table;

	vector<Face3d*>& face;
protected:
	//virtual void create_children(vector<Edge_c*>& edge);
private:
	int generation_max;
	virtual void getLeaf(Cell3d* cell);
	virtual void creatGenerationTable();
	virtual void registerGenerationTable(Cell3d* cell);
};

#endif // _FMM_TREE_3D_H_
