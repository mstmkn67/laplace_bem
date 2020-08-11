#include "Moment3d.h"

Mnm::Mnm(int ns){resize(ns);}
Mnm::~Mnm(){}
void Mnm::resize(int ns){
	moment.resize(ns+1);
	for(int i=0;i<=ns;i++){
	 moment[i].resize(i+1,complex<double>(0.0,0.0));
	}
}
complex<double> Mnm::get(int n,int m){
	if(n<0 || n<m || -n>m)return 0.0;
	if(m>=0){
		return moment[n][m];
	}else{
		if((-m)%2==0){
			return conj(moment[n][-m]);
		}else{
			return -conj(moment[n][-m]);
		}
	}
}
void Mnm::zero_reset(){
	vector<vector<complex<double> > >::iterator i=moment.begin();
	for(;i!=moment.end();i++){
		vector<complex<double> >::iterator j=i->begin();
		for(;j!=i->end();j++){
			(*j)=complex<double>(0.0,0.0);
		}
	}
}
//
Rnm::Rnm(int p):Mnm(p){}
Rnm::~Rnm(){}
void Rnm::update(const Vector3d& x){
	int p=moment.size()-1;
	moment[0][0]=complex<double>(1.0,0.0);
	for(int n=1;n<=p;n++){
		moment[n][n]=0.5/n*complex<double>(x.x,x.y)*moment[n-1][n-1];
	}
	for(int n=0;n<p;n++){
		moment[n+1][n]=x.z*moment[n][n];
	}
	for(int m=0;m<p-1;m++){
		for(int n=m+2;n<=p;n++){
			double a=1./(n+m)/(n-m);
			moment[n][m]=a*((2*n-1)*x.z*moment[n-1][m]-x.length2()*moment[n-2][m]);
		}
	}
}
complex<double> Rnm::get_dRdn(int n,int m,const Vector3d& nn){
	return 0.5*nn.x*(get(n-1,m-1)-get(n-1,m+1))
		    +0.5*nn.y*complex<double>(0.0,1.0)*(get(n-1,m-1)+get(n-1,m+1))
		    +nn.z*get(n-1,m);
}
complex<double> Rnm::get_dRdxi(int n,int m,int i){
	if(i==0)     return 0.5*(get(n-1,m-1)-get(n-1,m+1));
	else if(i==1)return 0.5*complex<double>(0.0,1.0)*(get(n-1,m-1)+get(n-1,m+1));
	else if(i==2)return get(n-1,m);
}
//
Snm::Snm(int p):Mnm(p){}
Snm::~Snm(){}
void Snm::update(const Vector3d& x){
	int p=moment.size()-1;
	moment[0][0]=complex<double>(1./x.length(),0.0);
	for(int n=1;n<=p;n++){
		moment[n][n]=(2*n-1.)/x.length2()*complex<double>(x.x,x.y)*moment[n-1][n-1];
	}
	for(int n=0;n<p;n++){
		moment[n+1][n]=(2*n+1.)/x.length2()*x.z*moment[n][n];
	}
	for(int m=0;m<p-1;m++){
		for(int n=m+2;n<=p;n++){
			moment[n][m]=1./x.length2()*((2*n-1.)*x.z*moment[n-1][m]-(n+m-1.)*(n-m-1.)*moment[n-2][m]);
		}
	}
}
complex<double> Snm::get_dSdn(int n,int m,const Vector3d& nn){
	return 0.5*nn.x*(get(n+1,m-1)-get(n+1,m+1))
		    +0.5*nn.y*complex<double>(0.0,1.0)*(get(n+1,m+1)+get(n+1,m-1))
		    -nn.z*get(n+1,m);
}
complex<double> Snm::get_dSdxi(int n,int m,int i){
	if(i==0)return 0.5*(get(n+1,m-1)-get(n+1,m+1));
	else if(i==1)return 0.5*complex<double>(0.0,1.0)*(get(n+1,m+1)+get(n+1,m-1));
	else if(i==2)return -get(n+1,m);
}
//
//Pnm::Pnm(int p):rnm(p){}
//Pnm::~Pnm(){}
//void Pnm::update(const Vector3d& x){
//	rnm.update(x);
//	rr=x;
//}
//complex<double> Pnm::get(int i,int j,int n,int m){
//	complex<double> s;
//	if(i==j)s=rnm.get(n,m);
//	
//}
//complex<double> Pnm::get(int i,int n,int m){
//	return rnm.get_dRdxi(n,m,i);
//}
	
//
void LaplaceMoment3d::resize(int n){
	M.resize(n);
	//N.resize(n);
}

void LaplaceMoment3d::zero_reset(){
	M.zero_reset();
	//N.zero_reset();
}

void StokesMoment3d::resize(int n){
	M1x.resize(n);
	M1y.resize(n);
	M1z.resize(n);
	M2.resize(n);
}

void StokesMoment3d::zero_reset(){
	M1x.zero_reset();
	M1y.zero_reset();
	M1z.zero_reset();
	M2.zero_reset();
}
