#ifndef _MOMENT_3D_H_
#define _MOMENT_3D_H_

#include <vector>
#include <complex>
#include "../Vector3d.h"
using namespace std;

class Mnm{
public:
	Mnm(int n_size=0);
	virtual ~Mnm();
	virtual void resize(int n_size);
	virtual void zero_reset();
	virtual complex<double> get(int n,int m);
	//void set(int n,int m,const complex<double>& v);
	vector<vector<complex<double> > > moment;
};
class Rnm:public Mnm{
public:
	Rnm(int p=0);
	virtual ~Rnm();
	virtual void update(const Vector3d& x);
	//gradient of Rnm and vector nn, not necessary that nn is unit normal
	virtual complex<double> get_dRdn(int n,int m,const Vector3d& nn);
	//diferential of x_i
	virtual complex<double> get_dRdxi(int n,int m,int i);
private:
};

class Snm:public Mnm{
public:
	Snm(int p=0);
	virtual ~Snm();
	virtual void update(const Vector3d& x);
	virtual complex<double> get_dSdn(int n,int m,const Vector3d& nn);
	//diferential of x_i
	virtual complex<double> get_dSdxi(int n,int m,int i);
private:
};

//class Pnm{
//public:
//	Pnm(int p=0);
//	virtual ~Pnm();
//	virtual void update(const Vector3d& x);
//	virtual complex<double> get(int i,int j,int n,int m);
//	virtual complex<double> get(int i,int n,int m);
//	Vector3d rr;
//	Rnm rnm;
//};

class Moment3d{
public:
	virtual void resize(int n)=0;
	virtual void zero_reset()=0;
};

class LaplaceMoment3d:public Moment3d{
public:
	virtual void resize(int n);
	virtual void zero_reset();
	Mnm M;//or L
	//Mnm N;//or K
};

class StokesMoment3d:public Moment3d{
public:
	virtual void resize(int n);
	virtual void zero_reset();
	//Minm M1;
	Mnm M1x;
	Mnm M1y;
	Mnm M1z;
	Mnm M2;
};

#endif // _MOMENT_3D_H_

