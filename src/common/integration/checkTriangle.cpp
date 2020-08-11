#include "LegendreGaussFormulaTriangle.h"
#include <iostream>
using namespace std;
#include "../Vector2d.h"

class Test:public LegendreGaussFormulaTriangle<double,Vector2d>{
public:
	Test(int point):LegendreGaussFormulaTriangle<double,Vector2d>(point){}
	virtual double integrand(const Vector2d& r)const{
		if(r.x*r.x+r.y*r.y>1.0)return 0.0;
		return sqrt(1.0-r.x*r.x-r.y*r.y);
	}
};


int main(){
	double ds=0.05;
	double s[4]={0,0,0,0};
	LegendreGaussFormulaTriangle<double,Vector2d>* test[4]; 
	test[0]=new Test(1);test[1]=new Test(3);test[2]=new Test(4);test[3]=new Test(7);
	double area=ds*ds*0.5;
	for(double x=-1.0;x<1.0;x+=ds){
		for(double y=-1.0;y<1.0;y+=ds){
			Vector2d r0(x,y),r1(x+ds,y),r2(x+ds,y+ds);
			for(int i=0;i<4;i++){
				s[i]+=test[i]->integral(r0,r1,r2,area);
			}
			r1.set(x,y+ds);
			for(int i=0;i<4;i++){
				s[i]+=test[i]->integral(r0,r1,r2,area);
			}
		}
	}
	for(int i=0;i<4;i++){
		cout << "Legendre-Gauss integral: " <<s[i] << endl;
		delete test[i];
	}
	return 0;
}
