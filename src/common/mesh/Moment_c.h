#ifndef _MOMENT_C_H_
#define _MOMENT_C_H_

#include <vector>
#include <complex>
using namespace std;

class Moment_c{
public:
	virtual void resize(int n)=0;
	virtual void zero_reset()=0;
};

class LaplaceMoment_c:public Moment_c{
public:
	virtual void resize(int n);
	virtual void zero_reset();
	vector<complex<double> > M;//or L
	//vector<complex<double> > N;//or K
};

class StokesMoment_c:public Moment_c{
public:
	virtual void resize(int n);
	virtual void zero_reset();
	vector<complex<double> > M;//or L
	vector<complex<double> > N;//or K
};

#endif // _MOMENT_C_H_
