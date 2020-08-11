#include "LegendreGaussFormula1d.h"
#include <iostream>
using namespace std;

class Test:public LegendreGaussFormula1d<double,double>{
public:
	Test(int point):LegendreGaussFormula1d<double,double>(point){};
	virtual double integrand(const double& r)const{return r*r*r;}
};

int main(){
	double dx=0.05;
	double s[6]={0,0,0,0,0,0};
	LegendreGaussFormula1d<double,double>* test[6]; 
	test[0]=new Test(1);test[1]=new Test(2);test[2]=new Test(3);
	test[3]=new Test(4);test[4]=new Test(5);test[5]=new Test(6);
	for(double x=0.0;x<1.0;x+=dx){
		for(int i=0;i<6;i++){
			s[i]+=test[i]->integral(x,x+dx,dx);
		}
	}
	for(int i=0;i<6;i++){
		cout << i+1 << "points Legendre-Gauss integral: " <<s[i] << endl;
		delete test[i];
	}
	return 0;
}

