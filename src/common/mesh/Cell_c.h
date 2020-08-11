#ifndef _CELL_C_H_
#define _CELL_C_H_

#include "Mesh_c.h"
#include "Moment_c.h"
#include <vector>
#include <iostream>
using namespace std;

class Cell_c{
public:
	Cell_c(vector<Edge_c*> edge,Cell_c* parent,
	       double size, const complex<double>& position,int generation);
	virtual ~Cell_c();
	
	const int id; //cell id
	const int generation;
	//
	const double size;      //linear size of cell
	const complex<double> position;//center of cell
	//tree
	Cell_c* parent;
	vector<Cell_c*> child;
	//
	vector<Cell_c*> adjacent;
	vector<Cell_c*> interaction;
	//used when cell is leaf
	vector<Edge_c*> edge;
	//moment
	Moment_c* moment;
	//local_expansion
	Moment_c* local;
	//
	static int population_per_leaf;
	static int total_cell_id;

	virtual void teach_children_adjacent_interaction_cells();//parent has the obligation to educate his children
	virtual Cell_c* get_cell_at_evaluation_position(const complex<double>& r);
	virtual void create_children(vector<Edge_c*> edge);
private:
};

#endif // _CELL_C_H_
