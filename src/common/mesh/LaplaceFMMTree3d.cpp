#include "LaplaceFMMTree3d.h"

LaplaceFMMTree3d::LaplaceFMMTree3d(vector<Face3d*>& face,double size,const Vector3d& position,
	int population_per_leaf,int mn,int ln)
:FMMTree3d(face,size,position,population_per_leaf),moment_number(mn),local_number(ln){
	vector<vector<Cell3d*> >::iterator c=generation_table.begin();
	for(;c!=generation_table.end();c++){
		vector<Cell3d*>::iterator d=c->begin();
		for(;d!=c->end();d++){
			(*d)->moment=new LaplaceMoment3d;
			(*d)->moment->resize(moment_number);
			(*d)->local=new LaplaceMoment3d;
			(*d)->local->resize(local_number);
		}
	}
}

LaplaceFMMTree3d::~LaplaceFMMTree3d(){
	vector<vector<Cell3d*> >::iterator c=generation_table.begin();
	for(;c!=generation_table.end();c++){
		vector<Cell3d*>::iterator d=c->begin();
		for(;d!=c->end();d++){
			delete (*d)->moment;
			delete (*d)->local;
		}
	}
}

void LaplaceFMMTree3d::update(){
	vector<vector<Cell3d*> >::iterator c=generation_table.begin();
	for(;c!=generation_table.end();c++){
		vector<Cell3d*>::iterator d=c->begin();
		for(;d!=c->end();d++){
			delete (*d)->moment;
			delete (*d)->local;
		}
	}
	FMMTree3d::update();
	c=generation_table.begin();
	for(;c!=generation_table.end();c++){
		vector<Cell3d*>::iterator d=c->begin();
		for(;d!=c->end();d++){
			(*d)->moment=new LaplaceMoment3d;
			(*d)->moment->resize(moment_number);
			(*d)->local=new LaplaceMoment3d;
			(*d)->local->resize(local_number);
		}
	}
}

void LaplaceFMMTree3d::update(double s,const Vector3d& p){
	vector<vector<Cell3d*> >::iterator c=generation_table.begin();
	for(;c!=generation_table.end();c++){
		vector<Cell3d*>::iterator d=c->begin();
		for(;d!=c->end();d++){
			delete (*d)->moment;
			delete (*d)->local;
		}
	}
	update(size,position);
	c=generation_table.begin();
	for(;c!=generation_table.end();c++){
		vector<Cell3d*>::iterator d=c->begin();
		for(;d!=c->end();d++){
			(*d)->moment=new LaplaceMoment3d;
			(*d)->moment->resize(moment_number);
			(*d)->local=new LaplaceMoment3d;
			(*d)->local->resize(local_number);
		}
	}
}
