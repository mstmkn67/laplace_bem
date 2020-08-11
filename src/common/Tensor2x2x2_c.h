//Tensor2x2x2

#ifndef _TENSOR_2X2X2_C_H_
#define _TENSOR_2X2X2_C_H_
#include <iomanip>
#include "Tensor2x2_c.h"

class Tensor2x2x2_c{
public:
	//コンストラクタ
	Tensor2x2x2_c(double mxxx= 0.0,double mxxy= 0.0,
                double mxyx= 0.0,double mxyy= 0.0,//
                double myxx= 0.0,double myxy= 0.0,
                double myyx= 0.0,double myyy= 0.0)
	:x(mxxx,mxxy,mxyx,mxyy),y(myxx,myxy,myyx,myyy){}
	//コンストラクタ
	Tensor2x2x2_c(Tensor2x2_c mx,Tensor2x2_c my):x(mx),y(my){}
	
	//メンバ
	Tensor2x2_c x,y;
	
	void clear(){
		x.clear();y.clear();
	}
	//double& operator()(unsigned i,unsigned j,unsigned k);
	double get(unsigned i,unsigned j,unsigned k)const;
	Tensor2x2x2_c& operator = (const Tensor2x2x2_c& m);
	Tensor2x2x2_c& operator += (const Tensor2x2x2_c& m);
	Tensor2x2x2_c& operator -= (const Tensor2x2x2_c& m);
	Tensor2x2x2_c& operator *= (const double& m);
	Tensor2x2x2_c& operator /= (const double& m);
	
	Tensor2x2x2_c operator - ()const;

	friend Tensor2x2x2_c operator + (const Tensor2x2x2_c& m1,const Tensor2x2x2_c& m2);
	friend Tensor2x2x2_c operator - (const Tensor2x2x2_c& m1,const Tensor2x2x2_c& m2);
	friend Tensor2x2x2_c operator / (const Tensor2x2x2_c& m,const double& s);
	friend Tensor2x2x2_c operator * (const double& d,const Tensor2x2x2_c& m1);
	friend Tensor2x2x2_c operator * (const Tensor2x2x2_c& m1,const double& d);
	
	//特殊な操作
	//A_{ijl}B_{l}の操作A・Bを行う
	friend Tensor2x2_c operator * (const complex<double>& v,const Tensor2x2x2_c& m);
	friend Tensor2x2_c operator * (const Tensor2x2x2_c& m,const complex<double>& v);
	//A_{ijl}B_{lk}の操作A・Bを行う
	friend Tensor2x2x2_c operator * (const Tensor2x2_c& m1,const Tensor2x2x2_c& m2);
	friend Tensor2x2x2_c operator * (const Tensor2x2x2_c& m1,const Tensor2x2_c& m2);
	
	friend Tensor2x2x2_c dyad(const complex<double>& v,const Tensor2x2_c& m);
	friend Tensor2x2x2_c dyad(const Tensor2x2_c& m,const complex<double>& v);
	
	friend ostream& operator << (ostream& os, const Tensor2x2x2_c& m)
	{
		os	<< m.x << endl << m.y;
		return os;
	}

	friend istream& operator >> (istream& is, Tensor2x2x2_c& m)
	{
		is	>> m.x >> m.y;
		return is;
	}
};

//inline double& Tensor2x2x2_c::operator()(unsigned i,unsigned j,unsigned k)
//{
//	if(i==0){
//		return x(j,k);
//	}else{
//		return y(j,k);
//	}
//}

inline double Tensor2x2x2_c::get(unsigned i,unsigned j,unsigned k)const
{
	if(i==0){
		return x.get(j,k);
	}else{
		return y.get(j,k);
	}
}

inline Tensor2x2x2_c& Tensor2x2x2_c::operator = (const Tensor2x2x2_c& m)
{
	x=m.x;y=m.y;
	return *this;
}

inline Tensor2x2x2_c& Tensor2x2x2_c::operator += (const Tensor2x2x2_c& m)
{
	x+=m.x;y+=m.y;
	return *this;
}

inline Tensor2x2x2_c& Tensor2x2x2_c::operator -= (const Tensor2x2x2_c& m)
{
	x-=m.x;y-=m.y;
	return *this;
}

inline Tensor2x2x2_c& Tensor2x2x2_c::operator *= (const double& c)
{
	x*=c;y*=c;
	return *this;
}

inline Tensor2x2x2_c& Tensor2x2x2_c::operator /= (const double& c)
{
	x/=c;y/=c;
	return *this;
}

inline Tensor2x2x2_c Tensor2x2x2_c::operator - () const
{
	return Tensor2x2x2_c(-x, -y);
}

inline Tensor2x2x2_c operator + (const Tensor2x2x2_c& m1, const Tensor2x2x2_c& m2)
{
	return Tensor2x2x2_c(m1.x+m2.x,m1.y+m2.y);
}

inline Tensor2x2x2_c operator - (const Tensor2x2x2_c& m1, const Tensor2x2x2_c& m2)
{
	return Tensor2x2x2_c(m1.x-m2.x,m1.y-m2.y);
}

inline Tensor2x2x2_c operator / (const Tensor2x2x2_c& m, const double& s)
{
	return Tensor2x2x2_c(m.x/s,m.y/s);
}

inline Tensor2x2x2_c operator * (const double& s, const Tensor2x2x2_c& m)
{
	return Tensor2x2x2_c(s*m.x,s*m.y);
}

inline Tensor2x2x2_c operator * (const Tensor2x2x2_c& m, const double& s)
{
	return Tensor2x2x2_c(m.x*s,m.y*s);
}

inline Tensor2x2_c operator * (const complex<double>& v,const Tensor2x2x2_c& m)
{
	return Tensor2x2_c(v.real()*m.x.x.real()+v.imag()*m.y.x.real(),
		                   v.real()*m.x.x.imag()+v.imag()*m.y.x.imag(),
	                   v.real()*m.x.y.real()+v.imag()*m.y.x.imag(),
	                     v.real()*m.x.y.imag()+v.imag()*m.y.y.imag());
}

inline Tensor2x2_c operator * (const Tensor2x2x2_c& m,const complex<double>& v)
{
	return Tensor2x2_c(m.x*v,m.y*v);
}

inline Tensor2x2x2_c operator * (const Tensor2x2_c& m1,const Tensor2x2x2_c& m2)
{
	return Tensor2x2x2_c(m1.x.real()*m2.x.x.real()+m1.x.imag()*m2.y.x.real(),
	                     m1.x.real()*m2.x.x.imag()+m1.x.imag()*m2.y.x.imag(),
	                     m1.x.real()*m2.x.y.real()+m1.x.imag()*m2.y.y.real(),
	                     m1.x.real()*m2.x.y.imag()+m1.x.imag()*m2.y.y.imag(),
		                   
	                     m1.y.real()*m2.x.x.real()+m1.y.imag()*m2.y.x.real(),
	                     m1.y.real()*m2.x.x.imag()+m1.y.imag()*m2.y.x.imag(),
	                     m1.y.real()*m2.x.y.real()+m1.y.imag()*m2.y.y.real(),
	                     m1.y.real()*m2.x.y.imag()+m1.y.imag()*m2.y.y.imag());
}

inline Tensor2x2x2_c operator * (const Tensor2x2x2_c& m1,const Tensor2x2_c& m2)
{
	return Tensor2x2x2_c(m1.x*m2,m1.y*m2);
}

inline Tensor2x2x2_c dyad(const complex<double>& v,const Tensor2x2_c& m){
	return Tensor2x2x2_c(v.real()*m,v.imag()*m);
}

inline Tensor2x2x2_c dyad(const Tensor2x2_c& m,const complex<double>& v){
	return Tensor2x2x2_c(m.x.real()*v.real(), m.x.real()*v.imag(),
	                     m.x.imag()*v.real(), m.x.imag()*v.imag(),
	                     m.y.real()*v.real(), m.y.real()*v.imag(),
	                     m.y.imag()*v.real(), m.y.imag()*v.imag());
}

#endif //_TENSOR_2X2X2_C_H_

