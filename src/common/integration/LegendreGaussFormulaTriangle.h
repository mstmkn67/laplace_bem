#ifndef _LEGENDRE_GAUSS_FORMULA_TRIANGLE_H_
#define _LEGENDRE_GAUSS_FORMULA_TRIANGLE_H_

#include "../Vector3d.h"

template<class returnType,class argumentType>
class LegendreGaussFormulaTriangle{
public:
	LegendreGaussFormulaTriangle(int point);
	virtual ~LegendreGaussFormulaTriangle();
	virtual returnType integrand(const argumentType& r)const=0;
	virtual returnType integral(const argumentType& r0,const argumentType& r1,const argumentType& r2,double area)const;
	Vector3d* xi;
	double* w;
	int point;
protected:
};

struct ParameterOfLegendreGaussFormulaTiangle:public LegendreGaussFormulaTriangle<double,double>{
public:
	ParameterOfLegendreGaussFormulaTiangle(int point):LegendreGaussFormulaTriangle<double,double>(point){};
	virtual double integrand(const double& r)const{return 0;};
private:
};

template<class returnType,class argumentType>
LegendreGaussFormulaTriangle<returnType,argumentType>::LegendreGaussFormulaTriangle(int p):point(p){
	xi=new Vector3d[point];w=new double[point];
	if(point==1){
		xi[0].set(1./3,1./3,1./3);
		w[0]=1.0;
	}else if(point==3){
		xi[0].set(0.5,0.5,0.0);xi[1].set(0.0,0.5,0.5);xi[2].set(0.5,0.0,0.5);
		w[0]=w[1]=w[2]=1./3;
	}else if(point==4){
		xi[0].set(1./3,1./3,1./3);xi[1].set(0.6,0.2,0.2);xi[2].set(0.2,0.6,0.2);xi[3].set(0.2,0.2,0.6);
		w[0]=-9./16;w[1]=w[2]=w[3]=25./48;
	}else if(point==7){
		xi[0].set(0.33333333,0.33333333,0.33333333);xi[1].set(0.79742699,0.10128651,0.10128651);
		xi[2].set(0.10128651,0.79742699,0.10128651);xi[3].set(0.10128651,0.10128651,0.79742699);
		xi[4].set(0.05971587,0.47014206,0.47014206);xi[5].set(0.47014206,0.05971587,0.47014206);
		xi[6].set(0.47014206,0.47014206,0.05971587);
		w[0]=0.225;w[1]=w[2]=w[3]=0.12593918;w[4]=w[5]=w[6]=0.13239415;
	}else{
		cout << point << " points algorithm is ";
		cout << "not implemented in class LegendreGaussFormulaTriangle" << endl;
		exit(1);
	}
}
template<class returnType,class argumentType>
LegendreGaussFormulaTriangle<returnType,argumentType>::~LegendreGaussFormulaTriangle(){
	delete[] xi;delete[] w;
}

template<class returnType,class argumentType>
returnType LegendreGaussFormulaTriangle<returnType,argumentType>::
integral(const argumentType& r0,const argumentType& r1,const argumentType& r2,double area)const{
	argumentType r01=r1-r0,r02=r2-r0;
	returnType s(0);
	for(int i=0;i<point;i++){
		argumentType r=xi[i].x*r01+xi[i].y*r02+r0;
		s+=integrand(r)*w[i];
	}
	return area*s;
}

#endif // _LEGENDRE_GAUSS_FORMULA_TRIANGLE_H_
