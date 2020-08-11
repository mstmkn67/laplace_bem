#include "Cell3d.h"

int Cell3d::population_per_leaf=1;
int Cell3d::total_cell_id=0;

	
Cell3d::Cell3d(vector<Face3d*> f,Cell3d* p,double s,const Vector3d& pos,int g)
:parent(p),size(s),position(pos),generation(g),id(total_cell_id++){
	create_children(f);
}

Cell3d::~Cell3d(){
	//parent kills his childeren before his death.
	vector<Cell3d*>::iterator i=child.begin();
	for(;i!=child.end();i++)delete *i;
}

Cell3d* Cell3d::get_cell_at_evaluation_position(const Vector3d& er){
	//search the child room 
	int en=-1;
	if(er.z<position.z){
		if(er.y>position.y){ 
			if(er.x>position.x)en=0;
			else en=1;
		}else{
			if(er.x>position.x)en=3;
			else en=2;
		}
	}else{
		if(er.y>position.y){ 
			if(er.x>position.x)en=4;
			else en=5;
		}else{
			if(er.x>position.x)en=7;
			else en=6;
		}
	}
	//if child exists in room, put the evaluation point
	//if not, this point is treated by this cell
	vector<Cell3d*>::iterator i=child.begin();
	for(;i!=child.end();i++){
		const Vector3d& cr=(*i)->position;
		if(cr.z<position.z){
			if(cr.y>position.y){ 
				if(cr.x>position.x){
					if(en==0)return (*i)->get_cell_at_evaluation_position(er);
				}else{
					if(en==1)return (*i)->get_cell_at_evaluation_position(er);
				}
			}else{
				if(cr.x>position.x){
					if(en==3)return (*i)->get_cell_at_evaluation_position(er);
				}else{
					if(en==2)return (*i)->get_cell_at_evaluation_position(er);
				}
			}
		}else{
			if(cr.y>position.y){ 
				if(cr.x>position.x){
					if(en==4)return (*i)->get_cell_at_evaluation_position(er);
				}else{
					if(en==5)return (*i)->get_cell_at_evaluation_position(er);
				}
			}else{
				if(cr.x>position.x){
					if(en==7)return (*i)->get_cell_at_evaluation_position(er);
				}else{
					if(en==6)return (*i)->get_cell_at_evaluation_position(er);
				}
			}
		}
	}
	return this;
}

void Cell3d::create_children(vector<Face3d*> f){
	if(f.size()<=population_per_leaf){//this cell is leaf
		face=f;
		vector<Face3d*>::iterator i=face.begin();
		for(;i!=face.end();i++){
			(*i)->cell=this;
		}
		return;
	}
	vector<Face3d*> ne[8];// edges for next generation
	double ccs=0.25*size;
	Vector3d nc[8];//
	nc[0]=position+ccs*Vector3d( 1.0, 1.0,-1.0);nc[1]=position+ccs*Vector3d(-1.0, 1.0,-1.0);
	nc[2]=position+ccs*Vector3d(-1.0,-1.0,-1.0);nc[3]=position+ccs*Vector3d( 1.0,-1.0,-1.0);
	nc[4]=position+ccs*Vector3d( 1.0, 1.0, 1.0);nc[5]=position+ccs*Vector3d(-1.0, 1.0, 1.0);
	nc[6]=position+ccs*Vector3d(-1.0,-1.0, 1.0);nc[7]=position+ccs*Vector3d( 1.0,-1.0, 1.0);
	//
	vector<Face3d*>::iterator i=f.begin();
	for(;i!=f.end();i++){
		const Vector3d& r=(*i)->position;
		if(r.z<position.z){
			if(r.y>position.y){ 
				if(r.x>position.x)ne[0].push_back((*i));//child 0
				else ne[1].push_back((*i));//child 1
			}else{
				if(r.x>position.x)ne[3].push_back((*i));//child 3
				else ne[2].push_back((*i)); //child 2
			}
		}else{
			if(r.y>position.y){ 
				if(r.x>position.x)ne[4].push_back((*i));//child 4
				else ne[5].push_back((*i));//child 5
			}else{
				if(r.x>position.x)ne[7].push_back((*i));//child 7
				else ne[6].push_back((*i)); //child 6
			}
		}
	}
	for(int i=0;i<8;i++){
		if(ne[i].size()>=1){
			child.push_back(new Cell3d(ne[i],this,0.5*size,nc[i],generation+1));
		}
	}
}

void Cell3d::teach_children_adjacent_interaction_cells(){
	if(child.size()==0){//if you don't have children, you don't have any obligations.
		return;
	}else if(parent==0){//if you don't have parent, children are adjacent each other.
		vector<Cell3d*>::iterator i=child.begin();
		for(;i!=child.end();i++){
			vector<Cell3d*>::iterator j=i+1;
			for(;j!=child.end();j++){
				(*i)->adjacent.push_back(*j);
				(*j)->adjacent.push_back(*i);
			}
		}
	}else{//your negiborhood's information is thaught to your children
		vector<Cell3d*>::iterator c=child.begin();
		for(;c!=child.end();c++){
			//introduction to brothers 
			vector<Cell3d*>::iterator k=c+1;
			for(;k!=child.end();k++){
				(*c)->adjacent.push_back(*k);
				(*k)->adjacent.push_back(*c);
			}
			//introduction to friend from your adjacent
			const Vector3d& pc=(*c)->position;
			vector<Cell3d*>::iterator i=adjacent.begin();
			for(;i!=adjacent.end();i++){
				//if adjacent cell is leaf, it is adjacent of children
				if((*i)->child.size()==0){
					(*c)->adjacent.push_back(*i);
					continue;
				}
				//introduce the friend of adjacent's children
				vector<Cell3d*>::iterator j=(*i)->child.begin();
				for(;j!=(*i)->child.end();j++){
					//check for same generation of children
					const Vector3d& pj=(*j)->position;
					double l2=(pj-pc).length2();
					if(l2<0.8*size*size){//adjacent, (l/2)^2 or (l*0.5*sqrt(2))^2 or 3/4 l^2
						(*c)->adjacent.push_back(*j);
					}else{//interaction
						(*c)->interaction.push_back(*j);
					}
				}
			}
		}
	}
	vector<Cell3d*>::iterator c=child.begin();//parent teaches teaching
	for(;c!=child.end();c++){
		(*c)->teach_children_adjacent_interaction_cells();
	}
}

