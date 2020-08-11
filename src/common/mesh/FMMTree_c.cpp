#include "FMMTree_c.h"

FMMTree_c::FMMTree_c(vector<Edge_c*>& e,double s,const complex<double>& p,int n)
:root(0),edge(e),size(s),position(p),generation_max(0){
	Cell_c::population_per_leaf=n;
	update();
}

FMMTree_c::~FMMTree_c(){
	delete root;
}

void FMMTree_c::update(){
	if(root!=0){
		delete root;
		leaf.clear();
		generation_table.clear();
	}
	Cell_c::total_cell_id=0;
	root=new Cell_c(edge,0,size,position,0);
	root->teach_children_adjacent_interaction_cells();
	getLeaf(root);
	creatGenerationTable();
}

void FMMTree_c::update(double s,const complex<double>& p){
	size=s;
	position=p;
	update();
}

Cell_c* FMMTree_c::get_cell_at_evaluation_position(const complex<double>& r){
	return root->get_cell_at_evaluation_position(r);
}

void FMMTree_c::getLeaf(Cell_c* cell){
	if(cell->child.size()==0){
		leaf.push_back(cell);
		if(cell->generation > generation_max){
			generation_max=cell->generation;
		}
		return;
	}
	vector<Cell_c*>::const_iterator i=cell->child.begin();
	for(;i!=cell->child.end();i++){
		getLeaf(*i);
	}
}

void FMMTree_c::creatGenerationTable(){
	generation_table.resize(generation_max+1);
	registerGenerationTable(root);
}

void FMMTree_c::registerGenerationTable(Cell_c* cell){
	generation_table[cell->generation].push_back(cell);
	vector<Cell_c*>::iterator i=cell->child.begin();
	for(;i!=cell->child.end();i++){
		registerGenerationTable(*i);
	}
}

