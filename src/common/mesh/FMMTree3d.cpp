#include "FMMTree3d.h"

FMMTree3d::FMMTree3d(vector<Face3d*>& f,double s,const Vector3d& p,int n)
:root(0),face(f),size(s),position(p),generation_max(0){
	Cell3d::population_per_leaf=n;
	update();
}

FMMTree3d::~FMMTree3d(){
	delete root;
}

void FMMTree3d::update(){
	if(root!=0){
		delete root;
		leaf.clear();
		generation_table.clear();
	}
	Cell3d::total_cell_id=0;
	root=new Cell3d(face,0,size,position,0);
	root->teach_children_adjacent_interaction_cells();
	getLeaf(root);
	creatGenerationTable();
}

void FMMTree3d::update(double s,const Vector3d& p){
	size=s;
	position=p;
	update();
}

Cell3d* FMMTree3d::get_cell_at_evaluation_position(const Vector3d& r){
	return root->get_cell_at_evaluation_position(r);
}

void FMMTree3d::getLeaf(Cell3d* cell){
	if(cell->child.size()==0){
		leaf.push_back(cell);
		if(cell->generation > generation_max){
			generation_max=cell->generation;
		}
		return;
	}
	vector<Cell3d*>::const_iterator i=cell->child.begin();
	for(;i!=cell->child.end();i++){
		getLeaf(*i);
	}
}

void FMMTree3d::creatGenerationTable(){
	generation_table.resize(generation_max+1);
	registerGenerationTable(root);
}

void FMMTree3d::registerGenerationTable(Cell3d* cell){
	generation_table[cell->generation].push_back(cell);
	vector<Cell3d*>::iterator i=cell->child.begin();
	for(;i!=cell->child.end();i++){
		registerGenerationTable(*i);
	}
}

