#include "Moment_c.h"

void LaplaceMoment_c::resize(int n){
	M.resize(n);
	//N.resize(n);
}

void LaplaceMoment_c::zero_reset(){
	vector<complex<double> >::iterator i=M.begin();
	for(;i!=M.end();i++){
		(*i)=complex<double>(0.0,0.0);
	}
	//vector<complex<double> >::iterator j=N.begin();
	//for(;j!=N.end();j++){
	//	j->real()=0.0;j->imag()=0.0;
	//}
}


void StokesMoment_c::resize(int n){
	M.resize(n);
	N.resize(n);
}

void StokesMoment_c::zero_reset(){
	vector<complex<double> >::iterator i=M.begin();
	for(;i!=M.end();i++){
		(*i)=complex<double>(0.0,0.0);
	}
	vector<complex<double> >::iterator j=N.begin();
	for(;j!=N.end();j++){
		(*j)=complex<double>(0.0,0.0);
	}
}

