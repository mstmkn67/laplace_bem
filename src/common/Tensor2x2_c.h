//Tensor2x2_c

#ifndef _TENSOR_2X2_C_H_
#define _TENSOR_2X2_C_H_
#include <iomanip>
#include <vector>
#include <complex>
using namespace std;

class Tensor2x2_c{
public:
	Tensor2x2_c(double mxx=0.0,double mxy=0.0,
              double myx=0.0,double myy=0.0):x(mxx,mxy),y(myx,myy){}
	Tensor2x2_c(const complex<double>& vecx,const complex<double>& vecy)
		:x(vecx),y(vecy){}
	//ÉÅÉìÉo
	complex<double> x,y;
	void set(double mxx=0.0,double mxy=0.0,
           double myx=0.0,double myy=0.0)
		{x=complex<double>(mxx,mxy);y=complex<double>(myx,myy);}
	void set(const complex<double>& vecx,const complex<double>& vecy)
		{x=vecx;y=vecy;}
	double det();
	double det()const;
	double trace(){return x.real()+y.imag();}
	double trace()const{return x.real()+y.imag();}
	void clear(){x=complex<double>(0.0,0.0);y=complex<double>(0.0,0.0);}
	//double& operator()(unsigned i,unsigned j);
	double get(unsigned i,unsigned j)const;
	Tensor2x2_c transpose();
	Tensor2x2_c invert();
	Tensor2x2_c transpose()const;
	Tensor2x2_c invert()const;

	Tensor2x2_c& operator = (const Tensor2x2_c& m);
	Tensor2x2_c& operator += (const Tensor2x2_c& m);
	Tensor2x2_c& operator -= (const Tensor2x2_c& m);
	Tensor2x2_c& operator *= (const double& m);
	Tensor2x2_c& operator /= (const double& m);

	Tensor2x2_c operator - ()const;
	
	friend Tensor2x2_c operator + (const Tensor2x2_c& m1,const Tensor2x2_c& m2);
	friend Tensor2x2_c operator - (const Tensor2x2_c& m1,const Tensor2x2_c& m2);
	friend Tensor2x2_c operator / (const Tensor2x2_c& m,const double& s);
	friend Tensor2x2_c operator * (const Tensor2x2_c& m1,const Tensor2x2_c& m2);
	friend Tensor2x2_c operator * (const double& d, const Tensor2x2_c& m1);
	friend Tensor2x2_c operator * (const Tensor2x2_c& m1,const double& d);
	friend complex<double>  operator * (const Tensor2x2_c& m, const complex<double>& v);
	friend Tensor2x2_c dyad(const complex<double>& v1,const complex<double>& v2);

	friend ostream& operator << (ostream& os, const Tensor2x2_c& m)
	{
		os	<< m.x << endl << m.y ;
		return os;
	}

	friend istream& operator >> (istream& is, Tensor2x2_c& m)
	{
		is	>> m.x	>> m.y;
		return is;
	}
};

//inline double& Tensor2x2_c::operator()(unsigned i,unsigned j)
//{
//	if(i==0){
//		if(j==0)return x.real();
//		return x.imag();
//	}else if(i==1){
//		if(j==0)return y.real();
//		return y.imag();
//	}
//}

inline double Tensor2x2_c::get(unsigned i,unsigned j)const
{
	if(i==0){
		if(j==0)return x.real();
		return x.imag();
	}else if(i==1){
		if(j==0)return y.real();
		return y.imag();
	}
}

inline Tensor2x2_c& Tensor2x2_c::operator = (const Tensor2x2_c& m)
{
	x=m.x;
	y=m.y;
	return *this;
}

inline Tensor2x2_c& Tensor2x2_c::operator += (const Tensor2x2_c& m)
{
	x+=m.x;
	y+=m.y;
	return *this;
}

inline Tensor2x2_c& Tensor2x2_c::operator -= (const Tensor2x2_c& m)
{
	x-=m.x;
	y-=m.y;
	return *this;
}

inline Tensor2x2_c& Tensor2x2_c::operator *= (const double& m)
{
	x*=m;
	y*=m;
	return *this;
}

inline Tensor2x2_c& Tensor2x2_c::operator /= (const double& m)
{
	x/=m;
	y/=m;
	return *this;
}

inline Tensor2x2_c Tensor2x2_c::operator - ()const
{
	return Tensor2x2_c(-x,-y);
}

inline Tensor2x2_c operator + (const Tensor2x2_c& m1, const Tensor2x2_c& m2)
{
	return Tensor2x2_c(m1.x+m2.x, m1.y+m2.y);
}

inline Tensor2x2_c operator - (const Tensor2x2_c& m1, const Tensor2x2_c& m2)
{
	return Tensor2x2_c(m1.x-m2.x, m1.y-m2.y);
}

inline Tensor2x2_c operator / (const Tensor2x2_c& m, const double& s)
{
	return Tensor2x2_c(m.x/s, m.y/s);
}

inline Tensor2x2_c operator * (const Tensor2x2_c& m1, const Tensor2x2_c& m2){
	return Tensor2x2_c(
		m1.x.real()*m2.x.real() + m1.x.imag()*m2.y.real(),
	  m1.x.real()*m2.x.imag() + m1.x.imag()*m2.y.imag(),
		m1.y.real()*m2.x.real() + m1.y.imag()*m2.y.real(),
	  m1.y.real()*m2.x.imag() + m1.y.imag()*m2.y.imag());
}

inline Tensor2x2_c operator * (const double& d,const Tensor2x2_c& m){
	return Tensor2x2_c(d*m.x,d*m.y);
}

inline Tensor2x2_c operator * (const Tensor2x2_c& m, const double& d){
	return Tensor2x2_c(d*m.x,d*m.y);
}

inline complex<double>  operator * (const Tensor2x2_c& m, const complex<double>& v){
	return complex<double>(m.x.real()*v.real()+m.x.imag()*v.imag(),
	                       m.y.real()*v.real()+m.y.imag()*v.imag());
}

inline Tensor2x2_c dyad(const complex<double>& v1,const complex<double>& v2){
	return Tensor2x2_c(v1.real()*v2.real(),v1.real()*v2.imag(),
	                   v1.imag()*v2.real(),v1.imag()*v2.imag());
}

inline Tensor2x2_c Tensor2x2_c::transpose()
{
	return Tensor2x2_c(x.real(),y.real(),
	                   x.imag(),y.imag());
}

inline Tensor2x2_c Tensor2x2_c::transpose()const
{
	return Tensor2x2_c(x.real(),y.real(),
	                   x.imag(),y.imag());
}

inline Tensor2x2_c Tensor2x2_c::invert(){
	Tensor2x2_c hinv( y.imag(),-x.imag(),
	                 -y.real(), x.real());
	double a=x.real()*y.imag()-x.imag()*y.real();
	return (hinv/a);
}

inline Tensor2x2_c Tensor2x2_c::invert()const{
	Tensor2x2_c hinv( y.imag(),-x.imag(),
	                 -y.real(), x.real());
	double a=x.real()*y.imag()-x.imag()*y.real();
	return (hinv/a);
}

inline double Tensor2x2_c::det()
{
	return x.real()*y.imag()-x.imag()*y.real();
}

inline double Tensor2x2_c::det()const
{
	return x.real()*y.imag()-x.imag()*y.real();
}


#endif //_TENSOR_2X2_C_H_

