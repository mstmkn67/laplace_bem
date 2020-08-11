#ifndef _LEGENDRE_GAUSS_FORMULA_1D_H_
#define _LEGENDRE_GAUSS_FORMULA_1D_H_

#include <iostream>
using namespace std;

//base class of 1 dimensional Legendre-Gauss integration formula 
template<class returnType,class argumentType>
class LegendreGaussFormula1d{
public:
	LegendreGaussFormula1d(int point);
	virtual ~LegendreGaussFormula1d();
	virtual returnType integrand(const argumentType& r)const=0;
	virtual returnType integral(const argumentType& rmin,const argumentType& rmax,double length)const;
	double* xi;
	double* w;
	int point;
protected:
};

struct ParameterOfLegendreGaussFormula1d:public LegendreGaussFormula1d<double,double>{
public:
	ParameterOfLegendreGaussFormula1d(int point):LegendreGaussFormula1d<double,double>(point){};
	virtual double integrand(const double& r)const{return 0;};
private:
};

template<class returnType,class argumentType>
LegendreGaussFormula1d<returnType,argumentType>::LegendreGaussFormula1d(int p):point(p){
	xi=new double[point];w=new double[point];
	if(point==1){
		xi[0]=0.0;w[0]=2.0;
	}else if(point==2){
		xi[0]=1.0/sqrt(3.);xi[1]=-xi[0];w[0]=w[1]=1.0;
	}else if(point==3){
		xi[0]=0.0;w[0]=8.0/9.0;
		xi[1]=sqrt(0.6);xi[2]=-xi[1];w[1]=w[2]=5.0/9.0;
	}else if(point==4){
		xi[0]=0.861136311594053;xi[1]=-xi[0];w[0]=w[1]=0.347854845137454;
		xi[2]=0.339981043584856;xi[3]=-xi[2];w[2]=w[3]=0.652145154862546;
	}else if(point==5){
		xi[0]=0.0;w[0]=0.568888888888889;
		xi[1]=0.906179845938664;xi[2]=-xi[1];w[1]=w[2]=0.236926885056189;
		xi[3]=0.538469310105683;xi[4]=-xi[3];w[3]=w[4]=0.478628670499366;
	}else if(point==6){
		xi[0]=0.932469514203152;xi[1]=-xi[0];w[0]=w[1]=0.171324492379170;
		xi[2]=0.661209386466265;xi[3]=-xi[2];w[2]=w[3]=0.360761573048139;
		xi[4]=0.238619186083197;xi[5]=-xi[4];w[4]=w[5]=0.467913934572691;
	}else if(point==7){
		xi[0]=0.0;w[0]=0.417959183673469;
		xi[1]=0.405845151377397;xi[2]=-xi[1];w[1]=w[2]=0.381830050505119;
		xi[3]=0.741531185599394;xi[4]=-xi[3];w[3]=w[4]=0.279705391489277;
		xi[5]=0.949107912342759;xi[6]=-xi[5];w[5]=w[6]=0.129484966168870;
	}else if(point==8){
		xi[0]=0.183435652495650;xi[1]=-xi[0];w[0]=w[1]=0.362683783378362;
		xi[2]=0.535532409916329;xi[3]=-xi[2];w[2]=w[3]=0.313706645877887;
		xi[4]=0.796666477413627;xi[5]=-xi[4];w[4]=w[5]=0.222381034453374;
		xi[6]=0.960289856497536;xi[7]=-xi[6];w[6]=w[7]=0.101228536290376;
	}else if(point==10){
		xi[0]=0.148874338981631;xi[1]=-xi[0];w[0]=w[1]=0.295524224714753;
		xi[2]=0.433395394129247;xi[3]=-xi[2];w[2]=w[3]=0.269266719309996;
		xi[4]=0.679409568299024;xi[5]=-xi[4];w[4]=w[5]=0.219086362515982;
		xi[6]=0.865063366688985;xi[7]=-xi[6];w[6]=w[7]=0.149451349150581;
		xi[7]=0.973906528517172;xi[8]=-xi[7];w[8]=w[9]=0.066671344208688;
	}else if(point==12){
		xi[0] =0.125233408511469;xi[1]= -xi[0]; w[0]= w[1]= 0.249147045813403;
		xi[2] =0.367831498998180;xi[3]= -xi[2]; w[2]= w[3]= 0.233492536538355;
		xi[4] =0.587317954286617;xi[5]= -xi[4]; w[4]= w[5]= 0.203167426723066;
		xi[6] =0.769902674194305;xi[7]= -xi[6]; w[6]= w[7]= 0.160078328543346;
		xi[7] =0.904117256370475;xi[8]= -xi[9]; w[8]= w[9]= 0.106939325995318;
		xi[10]=0.981560634246719;xi[11]=-xi[10];w[10]=w[11]=0.047175336386512;
	}else if(point==16){
		xi[0] =0.095012509837637;xi[1]= -xi[0]; w[0]= w[1]= 0.189450610455068;
		xi[2] =0.281603550779258;xi[3]= -xi[2]; w[2]= w[3]= 0.182603415044923;
		xi[4] =0.458016777657227;xi[5]= -xi[4]; w[4]= w[5]= 0.169156519395002;
		xi[6] =0.617876244402643;xi[7]= -xi[6]; w[6]= w[7]= 0.149595988816576;
		xi[7] =0.755404408355003;xi[8]= -xi[9]; w[8]= w[9]= 0.124628971255533;
		xi[10]=0.865631202387831;xi[11]=-xi[10];w[10]=w[11]=0.095158511682492;
		xi[12]=0.944575023073232;xi[13]=-xi[10];w[12]=w[13]=0.062253523938647;
		xi[14]=0.989400934991649;xi[15]=-xi[10];w[12]=w[13]=0.027152459411754;
	}else{
		cout << point << " points algorithm is ";
		cout << "not implemented in class LegendreGaussFormula1d" << endl;
		exit(1);
	}
}
template<class returnType,class argumentType>
LegendreGaussFormula1d<returnType,argumentType>::~LegendreGaussFormula1d(){
	delete[] xi;delete[] w;
}

template<class returnType,class argumentType>
returnType LegendreGaussFormula1d<returnType,argumentType>::
integral(const argumentType& rmin,const argumentType& rmax,double length)const{
	returnType s(0);
	for(int i=0;i<point;i++){
		argumentType r=0.5*((xi[i]+1.0)*rmax+(1.0-xi[i])*rmin);
		s+=integrand(r)*w[i];
	}
	return 0.5*length*s;
}

#endif // _LEGENDRE_GAUSS_FORMULA_1D_H_
