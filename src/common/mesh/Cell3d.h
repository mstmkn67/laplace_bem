#ifndef _CELL_3D_H_
#define _CELL_3D_H_

#include "Mesh3d.h"
#include "Moment3d.h"
#include <vector>
#include <complex>
#include <iostream>
using namespace std;


class Cell3d{
public:
	Cell3d(vector<Face3d*> face,Cell3d* parent,
	       double size, const Vector3d& position,int generation);
	virtual ~Cell3d();
	
	const int id; //cell id
	const int generation;
	//
	const double size;      //linear size of cell
	const Vector3d position;//center of cell
	//tree
	Cell3d* parent;
	vector<Cell3d*> child;
	//
	vector<Cell3d*> adjacent;
	vector<Cell3d*> interaction;
	//used when cell is leaf
	vector<Face3d*> face;
	//moment
	Moment3d* moment;
	//local_expansion
	Moment3d* local;
	//
	static int population_per_leaf;
	static int total_cell_id;

	virtual void teach_children_adjacent_interaction_cells();//parent has the obligation to educate his children
	virtual Cell3d* get_cell_at_evaluation_position(const Vector3d& r);
protected:
	virtual void create_children(vector<Face3d*> face);
private:
};

#endif // _CELL_3D_H_
