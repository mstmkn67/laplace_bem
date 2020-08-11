#include "Cell_c.h"

int Cell_c::population_per_leaf=1;
int Cell_c::total_cell_id=0;

	
Cell_c::Cell_c(vector<Edge_c*> e,Cell_c* p,double s,const complex<double>& pos,int g)
:parent(p),size(s),position(pos),generation(g),id(total_cell_id++){
	create_children(e);
}

Cell_c::~Cell_c(){
	//parent kills his childeren before his death.
	vector<Cell_c*>::iterator i=child.begin();
	for(;i!=child.end();i++)delete *i;
}

Cell_c* Cell_c::get_cell_at_evaluation_position(const complex<double>& er){
	//search the child room 
	int en=-1;
	if(er.imag()>position.imag()){ 
		if(er.real()>position.real())en=0;
		else en=1;
	}else{
		if(er.real()>position.real())en=3;
		else en=2;
	}
	//if child exists in room, put the evaluation point
	//if not, this point is treated by this cell
	vector<Cell_c*>::iterator i=child.begin();
	for(;i!=child.end();i++){
		const complex<double>& cr=(*i)->position;
		if(cr.imag()>position.imag()){ 
			if(cr.real()>position.real()){
				if(en==0)return (*i)->get_cell_at_evaluation_position(er);
			}else{
				if(en==1)return (*i)->get_cell_at_evaluation_position(er);
			}
		}else{
			if(cr.real()>position.real()){
				if(en==3)return (*i)->get_cell_at_evaluation_position(er);
			}
			else{
				if(en==2)return (*i)->get_cell_at_evaluation_position(er);
			}
		}
	}
	return this;
}

void Cell_c::create_children(vector<Edge_c*> e){
	if(e.size()<=population_per_leaf){//this cell is leaf
		edge=e;
		vector<Edge_c*>::iterator i=edge.begin();
		for(;i!=edge.end();i++){
			(*i)->cell=this;
		}
		return;
	}
	vector<Edge_c*> ne[4];// edges for next generation
	double ccs=0.25*size;
	complex<double> nc[4];//
	nc[0]=position+ccs*complex<double>( 1.0, 1.0);
	nc[1]=position+ccs*complex<double>(-1.0, 1.0);
	nc[2]=position+ccs*complex<double>(-1.0,-1.0);
	nc[3]=position+ccs*complex<double>( 1.0,-1.0);
	//
	vector<Edge_c*>::iterator i=e.begin();
	for(;i!=e.end();i++){
		complex<double>& r=(*i)->position;
		if(r.imag()>position.imag()){ 
			if(r.real()>position.real()){//child 0
				ne[0].push_back((*i));
			}else{             //child 1
				ne[1].push_back((*i));
			}
		}else{
			if(r.real()>position.real()){//child 3 
				ne[3].push_back((*i));
			}else{             //child 2
				ne[2].push_back((*i));
			}
		}
	}
	for(int i=0;i<4;i++){
		if(ne[i].size()>=1){
			child.push_back(new Cell_c(ne[i],this,0.5*size,nc[i],generation+1));
		}
	}
}

void Cell_c::teach_children_adjacent_interaction_cells(){
	if(child.size()==0){//if you don't have children, you don't have any obligations.
		return;
	}else if(parent==0){//if you don't have parent, children are adjacent each other.
		vector<Cell_c*>::iterator i=child.begin();
		for(;i!=child.end();i++){
			vector<Cell_c*>::iterator j=i+1;
			for(;j!=child.end();j++){
				(*i)->adjacent.push_back(*j);
				(*j)->adjacent.push_back(*i);
			}
		}
	}else{//your negiborhood's information is thaught to your children
		vector<Cell_c*>::iterator c=child.begin();
		for(;c!=child.end();c++){
			//introduction to brothers 
			vector<Cell_c*>::iterator k=c+1;
			for(;k!=child.end();k++){
				(*c)->adjacent.push_back(*k);
				(*k)->adjacent.push_back(*c);
			}
			//introduction to friend from your adjacent
			const complex<double>& pc=(*c)->position;
			vector<Cell_c*>::iterator i=adjacent.begin();
			for(;i!=adjacent.end();i++){
				//if adjacent cell is leaf, it is adjacent of children
				if((*i)->child.size()==0){
					(*c)->adjacent.push_back(*i); //I believe this is better.
					//(*c)->interaction.push_back(*i);
					continue;
				}
				//introduce the friend of adjacent's children
				vector<Cell_c*>::iterator j=(*i)->child.begin();
				for(;j!=(*i)->child.end();j++){
					//check for same generation of children
					const complex<double>& pj=(*j)->position;
					double l2=norm(pj-pc);
					if(l2<0.6*size*size){//adjacent, exactly, (l/2)^2 or (l*0.5*sqrt(2))^2
						(*c)->adjacent.push_back(*j);
					}else{//interaction
						(*c)->interaction.push_back(*j);
					}
				}
			}
		}
	}
	vector<Cell_c*>::iterator c=child.begin();//parent teaches teaching
	for(;c!=child.end();c++){
		(*c)->teach_children_adjacent_interaction_cells();
	}
}
