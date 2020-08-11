#include "LaplaceFMMTree_c.h"

LaplaceFMMTree_c::LaplaceFMMTree_c(vector<Edge_c*>& edge,double size,const complex<double>& position,
	int population_per_leaf,int mn,int ln)
:FMMTree_c(edge,size,position,population_per_leaf),moment_number(mn),local_number(ln){
	vector<vector<Cell_c*> >::iterator c=generation_table.begin();
	for(;c!=generation_table.end();c++){
		vector<Cell_c*>::iterator d=c->begin();
		for(;d!=c->end();d++){
			(*d)->moment=new LaplaceMoment_c;
			(*d)->moment->resize(moment_number+1);
			(*d)->local=new LaplaceMoment_c;
			(*d)->local->resize(local_number+1);
		}
	}
}

LaplaceFMMTree_c::~LaplaceFMMTree_c(){
	vector<vector<Cell_c*> >::iterator c=generation_table.begin();
	for(;c!=generation_table.end();c++){
		vector<Cell_c*>::iterator d=c->begin();
		for(;d!=c->end();d++){
			delete (*d)->moment;
			delete (*d)->local;
		}
	}
}

void LaplaceFMMTree_c::update(){
	vector<vector<Cell_c*> >::iterator c=generation_table.begin();
	for(;c!=generation_table.end();c++){
		vector<Cell_c*>::iterator d=c->begin();
		for(;d!=c->end();d++){
			delete (*d)->moment;
			delete (*d)->local;
		}
	}
	FMMTree_c::update();
	c=generation_table.begin();
	for(;c!=generation_table.end();c++){
		vector<Cell_c*>::iterator d=c->begin();
		for(;d!=c->end();d++){
			(*d)->moment=new LaplaceMoment_c;
			(*d)->moment->resize(moment_number+1);
			(*d)->local=new LaplaceMoment_c;
			(*d)->local->resize(local_number+1);
		}
	}
}

void LaplaceFMMTree_c::update(double size,const complex<double>& position){
	vector<vector<Cell_c*> >::iterator c=generation_table.begin();
	for(;c!=generation_table.end();c++){
		vector<Cell_c*>::iterator d=c->begin();
		for(;d!=c->end();d++){
			delete (*d)->moment;
			delete (*d)->local;
		}
	}
	FMMTree_c::update(size,position);
	c=generation_table.begin();
	for(;c!=generation_table.end();c++){
		vector<Cell_c*>::iterator d=c->begin();
		for(;d!=c->end();d++){
			(*d)->moment=new LaplaceMoment_c;
			(*d)->moment->resize(moment_number+1);
			(*d)->local=new LaplaceMoment_c;
			(*d)->local->resize(local_number+1);
		}
	}
}
