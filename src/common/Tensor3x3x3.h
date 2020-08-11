//Tensor3x3x3
//3x3x3の実数型、A_{ijk}の3階のテンソル

#ifndef _TENSOR_3X3X3_H_
#define _TENSOR_3X3X3_H_
#ifdef _WIN32
#pragma warning (disable : 4786)	// VC++ でデバッガーに関する不要の警告を押さえる
#endif
#include <iomanip>
#include "Vector3d.h"
#include "Tensor3x3.h"

class Tensor3x3x3{
public:
	//コンストラクタ
	Tensor3x3x3(double mxxx= 0.0,double mxxy= 0.0,double mxxz= 0.0,
              double mxyx= 0.0,double mxyy= 0.0,double mxyz= 0.0,
              double mxzx= 0.0,double mxzy= 0.0,double mxzz= 0.0,//
              double myxx= 0.0,double myxy= 0.0,double myxz= 0.0,
              double myyx= 0.0,double myyy= 0.0,double myyz= 0.0,
              double myzx= 0.0,double myzy= 0.0,double myzz= 0.0,//
              double mzxx= 0.0,double mzxy= 0.0,double mzxz= 0.0,
              double mzyx= 0.0,double mzyy= 0.0,double mzyz= 0.0,
              double mzzx= 0.0,double mzzy= 0.0,double mzzz= 0.0)
	:x(mxxx,mxxy,mxxz,mxyx,mxyy,mxyz,mxzx,mxzy,mxzz),
	 y(myxx,myxy,myxz,myyx,myyy,myyz,myzx,myzy,myzz),
	 z(mzxx,mzxy,mzxz,mzyx,mzyy,mzyz,mzzx,mzzy,mzzz){}
	//コンストラクタ
	Tensor3x3x3(Tensor3x3 mx,Tensor3x3 my,Tensor3x3 mz):x(mx),y(my),z(mz){}
	
	//メンバ
	Tensor3x3 x,y,z;
	
	void clear(){
		x.clear();y.clear();z.clear();
	}
	double& operator()(unsigned i,unsigned j,unsigned k);
	double get(unsigned i,unsigned j,unsigned k)const;
	vector<vector<vector<double> > > getSTLVector();
	Tensor3x3x3& operator = (const Tensor3x3x3& m);
	Tensor3x3x3& operator += (const Tensor3x3x3& m);
	Tensor3x3x3& operator -= (const Tensor3x3x3& m);
	Tensor3x3x3& operator *= (const double& m);
	Tensor3x3x3& operator /= (const double& m);
	
	Tensor3x3x3 operator - ()const;

	friend Tensor3x3x3 operator + (const Tensor3x3x3& m1,const Tensor3x3x3& m2);
	friend Tensor3x3x3 operator - (const Tensor3x3x3& m1,const Tensor3x3x3& m2);
	friend Tensor3x3x3 operator / (const Tensor3x3x3& m,const double& s);
	friend Tensor3x3x3 operator * (const double& d,const Tensor3x3x3& m1);
	friend Tensor3x3x3 operator * (const Tensor3x3x3& m1,const double& d);
	
	//特殊な操作
	//A_{ijl}B_{l}の操作A・Bを行う
	friend Tensor3x3 operator * (const Vector3d& v,const Tensor3x3x3& m);
	friend Tensor3x3 operator * (const Tensor3x3x3& m,const Vector3d& v);
	//A_{ijl}B_{lk}の操作A・Bを行う
	friend Tensor3x3x3 operator * (const Tensor3x3& m1,const Tensor3x3x3& m2);
	friend Tensor3x3x3 operator * (const Tensor3x3x3& m1,const Tensor3x3& m2);
	//A_{ijk}B_{jk}の操作A:Bを行う
	friend Vector3d  dot2(const Tensor3x3x3& m1, const Tensor3x3& m2);
	//A_{ij}B_{k}
	friend Tensor3x3x3 dyad(const Tensor3x3& m,const Vector3d& v);
	//B_{i}A_{jk}
	friend Tensor3x3x3 dyad(const Vector3d& v,const Tensor3x3& m);
	//A_{ik}B_{j}の操作を行う
	friend Tensor3x3x3 diamond(const Tensor3x3& m,const Vector3d& v);
	//B_{j}A_{ik}の操作を行う
	friend Tensor3x3x3 diamond(const Vector3d& v,const Tensor3x3& m);
	//座標系の変換のテンソルQを使って座標系を変換
	//G'_{abc}=Q_{ar}Q_{bp}Q_{cq}G_{rpq}
	friend Tensor3x3x3 convert(const Tensor3x3& Q,const Tensor3x3x3& G);
	//(G×r)_{abc}=ε_{ady}G_{bcd}r_{y}
	friend Tensor3x3x3 operator^(const Tensor3x3x3& G,const Vector3d& r);
	friend Tensor3x3x3 operator^(const Vector3d& r,const Tensor3x3x3& G);

	friend ostream& operator << (ostream& os, const Tensor3x3x3& m)
	{
		os	<< m.x << endl << m.y << endl << m.z;
		return os;
	}

	friend istream& operator >> (istream& is, Tensor3x3x3& m)
	{
		is	>> m.x >> m.y >> m.z;
		return is;
	}
};

inline double& Tensor3x3x3::operator()(unsigned i,unsigned j,unsigned k)
{
	if(i==0){
		return x(j,k);
	}else if(i==1){
		return y(j,k);
	}else{
		return z(j,k);
	}
}

inline double Tensor3x3x3::get(unsigned i,unsigned j,unsigned k)const
{
	if(i==0){
		return x.get(j,k);
	}else if(i==1){
		return y.get(j,k);
	}else{
		return z.get(j,k);
	}
}

inline vector<vector<vector<double> > > Tensor3x3x3::getSTLVector()
{
	vector<vector<vector<double> > > a;
	a.resize(3);
	a[0]=x.getSTLVector();
	a[1]=y.getSTLVector();
	a[2]=y.getSTLVector();
	return a;
}

inline Tensor3x3x3& Tensor3x3x3::operator = (const Tensor3x3x3& m)
{
	x=m.x;y=m.y;z=m.z;
	return *this;
}

inline Tensor3x3x3& Tensor3x3x3::operator += (const Tensor3x3x3& m)
{
	x+=m.x;y+=m.y;z+=m.z;
	return *this;
}

inline Tensor3x3x3& Tensor3x3x3::operator -= (const Tensor3x3x3& m)
{
	x-=m.x;y-=m.y;z-=m.z;
	return *this;
}

inline Tensor3x3x3& Tensor3x3x3::operator *= (const double& c)
{
	x*=c;y*=c;z*=c;
	return *this;
}

inline Tensor3x3x3& Tensor3x3x3::operator /= (const double& c)
{
	x/=c;y/=c;z/=c;
	return *this;
}

inline Tensor3x3x3 Tensor3x3x3::operator - () const
{
	return Tensor3x3x3(-x, -y, -z);
}

inline Tensor3x3x3 operator + (const Tensor3x3x3& m1, const Tensor3x3x3& m2)
{
	return Tensor3x3x3(m1.x+m2.x,m1.y+m2.y,m1.z+m2.z);
}

inline Tensor3x3x3 operator - (const Tensor3x3x3& m1, const Tensor3x3x3& m2)
{
	return Tensor3x3x3(m1.x-m2.x,m1.y-m2.y,m1.z-m2.z);
}

inline Tensor3x3x3 operator / (const Tensor3x3x3& m, const double& s)
{
	return Tensor3x3x3(m.x/s,m.y/s,m.z/s);
}

inline Tensor3x3x3 operator * (const double& s, const Tensor3x3x3& m)
{
	return Tensor3x3x3(s*m.x,s*m.y,s*m.z);
}

inline Tensor3x3x3 operator * (const Tensor3x3x3& m, const double& s)
{
	return Tensor3x3x3(m.x*s,m.y*s,m.z*s);
}

inline Tensor3x3 operator * (const Vector3d& v,const Tensor3x3x3& m)
{
	return Tensor3x3(v.x*m.x.x.x+v.y*m.y.x.x+v.z*m.z.x.x,
		             v.x*m.x.x.y+v.y*m.y.x.y+v.z*m.z.x.y,
					 v.x*m.x.x.z+v.y*m.y.x.z+v.z*m.z.x.z,//
					 v.x*m.y.y.x+v.y*m.y.y.x+v.z*m.z.y.x,
					 v.x*m.x.y.y+v.y*m.y.y.y+v.z*m.z.y.y,
					 v.x*m.x.y.z+v.y*m.y.y.z+v.z*m.z.y.z,//
					 v.x*m.x.z.x+v.y*m.y.z.x+v.z*m.z.z.x,
					 v.x*m.x.z.y+v.y*m.y.z.y+v.z*m.z.z.y,
					 v.x*m.x.z.z+v.y*m.y.z.z+v.z*m.z.z.z);
}

inline Tensor3x3 operator * (const Tensor3x3x3& m,const Vector3d& v)
{
	return Tensor3x3(m.x*v,m.y*v,m.z*v);
}

inline Tensor3x3x3 operator * (const Tensor3x3& m1,const Tensor3x3x3& m2)
{
	return Tensor3x3x3(m1.x.x*m2.x.x.x+m1.x.y*m2.y.x.x+m1.x.z*m2.z.x.x,
	                   m1.x.x*m2.x.x.y+m1.x.y*m2.y.x.y+m1.x.z*m2.z.x.y,
	                   m1.x.x*m2.x.x.z+m1.x.y*m2.y.x.z+m1.x.z*m2.z.x.z,//
	                   m1.x.x*m2.x.y.x+m1.x.y*m2.y.y.x+m1.x.z*m2.z.y.x,
	                   m1.x.x*m2.x.y.y+m1.x.y*m2.y.y.y+m1.x.z*m2.z.y.y,
	                   m1.x.x*m2.x.y.z+m1.x.y*m2.y.y.z+m1.x.z*m2.z.y.z,//
	                   m1.x.x*m2.x.z.x+m1.x.y*m2.y.z.x+m1.x.z*m2.z.z.x,
	                   m1.x.x*m2.x.z.y+m1.x.y*m2.y.z.y+m1.x.z*m2.z.z.y,
	                   m1.x.x*m2.x.z.z+m1.x.y*m2.y.z.z+m1.x.z*m2.z.z.z,//
	                   
	                   m1.y.x*m2.x.x.x+m1.y.y*m2.y.x.x+m1.y.z*m2.z.x.x,
	                   m1.y.x*m2.x.x.y+m1.y.y*m2.y.x.y+m1.y.z*m2.z.x.y,
	                   m1.y.x*m2.x.x.z+m1.y.y*m2.y.x.z+m1.y.z*m2.z.x.z,//
	                   m1.y.x*m2.x.y.x+m1.y.y*m2.y.y.x+m1.y.z*m2.z.y.x,
	                   m1.y.x*m2.x.y.y+m1.y.y*m2.y.y.y+m1.y.z*m2.z.y.y,
	                   m1.y.x*m2.x.y.z+m1.y.y*m2.y.y.z+m1.y.z*m2.z.y.z,//
	                   m1.y.x*m2.x.z.x+m1.y.y*m2.y.z.x+m1.y.z*m2.z.z.x,
	                   m1.y.x*m2.x.z.y+m1.y.y*m2.y.z.y+m1.y.z*m2.z.z.y,
	                   m1.y.x*m2.x.z.z+m1.y.y*m2.y.z.z+m1.y.z*m2.z.z.z,//
	                   
	                   m1.z.x*m2.x.x.x+m1.z.y*m2.y.x.x+m1.z.z*m2.z.x.x,
	                   m1.z.x*m2.x.x.y+m1.z.y*m2.y.x.y+m1.z.z*m2.z.x.y,
	                   m1.z.x*m2.x.x.z+m1.z.y*m2.y.x.z+m1.z.z*m2.z.x.z,//
	                   m1.z.x*m2.x.y.x+m1.z.y*m2.y.y.x+m1.z.z*m2.z.y.x,
	                   m1.z.x*m2.x.y.y+m1.z.y*m2.y.y.y+m1.z.z*m2.z.y.y,
	                   m1.z.x*m2.x.y.z+m1.z.y*m2.y.y.z+m1.z.z*m2.z.y.z,//
	                   m1.z.x*m2.x.z.x+m1.z.y*m2.y.z.x+m1.z.z*m2.z.z.x,
	                   m1.z.x*m2.x.z.y+m1.z.y*m2.y.z.y+m1.z.z*m2.z.z.y,
	                   m1.z.x*m2.x.z.z+m1.z.y*m2.y.z.z+m1.z.z*m2.z.z.z);//
}

inline Tensor3x3x3 operator * (const Tensor3x3x3& m1,const Tensor3x3& m2)
{
	return Tensor3x3x3(m1.x*m2,m1.y*m2,m1.z*m2);
}
inline Vector3d dot2(const Tensor3x3x3& m1, const Tensor3x3& m2){
	return Vector3d(dot2(m1.x,m2),dot2(m1.y,m2),dot2(m1.z,m2));
}

inline Tensor3x3x3 diamond(const Tensor3x3& m,const Vector3d& v)
{
	return Tensor3x3x3(m.x.x*v.x, m.x.y*v.x, m.x.z*v.x,
	                   m.x.x*v.y, m.x.y*v.y, m.x.z*v.y,
	                   m.x.x*v.z, m.x.y*v.z, m.x.z*v.z,//
	                   
					   m.y.x*v.x, m.y.y*v.x, m.y.z*v.x,
	                   m.y.x*v.y, m.y.y*v.y, m.y.z*v.y,
	                   m.y.x*v.z, m.y.y*v.z, m.y.z*v.z,//
	                   
	                   m.z.x*v.x, m.z.y*v.x, m.z.z*v.x,
	                   m.z.x*v.y, m.z.y*v.y, m.z.z*v.y,
	                   m.z.x*v.z, m.z.y*v.z, m.z.z*v.z);
}

inline Tensor3x3x3 diamond(const Vector3d& v,const Tensor3x3& m)
{
	return Tensor3x3x3(m.x.x*v.x, m.x.y*v.x, m.x.z*v.x,
	                   m.x.x*v.y, m.x.y*v.y, m.x.z*v.y,
	                   m.x.x*v.z, m.x.y*v.z, m.x.z*v.z,//
	                   
					   m.y.x*v.x, m.y.y*v.x, m.y.z*v.x,
	                   m.y.x*v.y, m.y.y*v.y, m.y.z*v.y,
	                   m.y.x*v.z, m.y.y*v.z, m.y.z*v.z,//
	                   
	                   m.z.x*v.x, m.z.y*v.x, m.z.z*v.x,
	                   m.z.x*v.y, m.z.y*v.y, m.z.z*v.y,
	                   m.z.x*v.z, m.z.y*v.z, m.z.z*v.z);
	//return left(m,v);
}

inline Tensor3x3x3 dyad(const Tensor3x3& m,const Vector3d& v)
{
	return Tensor3x3x3(m.x.x*v.x, m.x.x*v.y, m.x.x*v.z,
	                   m.x.y*v.x, m.x.y*v.y, m.x.y*v.z,
	                   m.x.z*v.x, m.x.z*v.y, m.x.z*v.z,//
	                   
	                   m.y.x*v.x, m.y.x*v.y, m.y.x*v.z,
	                   m.y.y*v.x, m.y.y*v.y, m.y.y*v.z,
	                   m.y.z*v.x, m.y.z*v.y, m.y.z*v.z,//
	                   
	                   m.z.x*v.x, m.z.x*v.y, m.z.x*v.z,
	                   m.z.y*v.x, m.z.y*v.y, m.z.y*v.z,
	                   m.z.z*v.x, m.z.z*v.y, m.z.z*v.z);
}

inline Tensor3x3x3 dyad(const Vector3d& v,const Tensor3x3& m)
{
	return Tensor3x3x3(v.x*m.x.x, v.x*m.x.y, v.x*m.x.z,
	                   v.x*m.y.x, v.x*m.y.y, v.x*m.y.z,
	                   v.x*m.z.x, v.x*m.z.y, v.x*m.z.z,//
	                   
	                   v.y*m.x.x, v.y*m.x.y, v.y*m.x.z,
	                   v.y*m.y.x, v.y*m.y.y, v.y*m.y.z,
	                   v.y*m.z.x, v.y*m.z.y, v.y*m.z.z,//
	                   
	                   v.z*m.x.x, v.z*m.x.y, v.z*m.x.z,
	                   v.z*m.y.x, v.z*m.y.y, v.z*m.y.z,
	                   v.z*m.z.x, v.z*m.z.y, v.z*m.z.z);
}

inline Tensor3x3x3 convert(const Tensor3x3& Q,const Tensor3x3x3& G)
{
	Tensor3x3 Qt=Q.transpose();
	Tensor3x3 tx=Q*G.x*Qt;
	Tensor3x3 ty=Q*G.y*Qt;
	Tensor3x3 tz=Q*G.z*Qt;
	return Q*Tensor3x3x3(tx,ty,tz);
}

inline Tensor3x3x3 operator^(const Tensor3x3x3& G,const Vector3d& r)
{
	return Tensor3x3x3(
		G.x.x.y*r.z-G.x.x.z*r.y, G.x.x.z*r.x-G.x.x.x*r.z, G.x.x.x*r.y-G.x.x.y*r.x,
		G.x.y.y*r.z-G.x.y.z*r.y, G.x.y.z*r.x-G.x.y.x*r.z, G.x.y.x*r.y-G.x.y.y*r.x,
		G.x.z.y*r.z-G.x.z.z*r.y, G.x.z.z*r.x-G.x.z.x*r.z, G.x.z.x*r.y-G.x.z.y*r.x,//
		
		G.y.x.y*r.z-G.y.x.z*r.y, G.y.x.z*r.x-G.y.x.x*r.z, G.y.x.x*r.y-G.y.x.y*r.x,
		G.y.y.y*r.z-G.y.y.z*r.y, G.y.y.z*r.x-G.y.y.x*r.z, G.y.y.x*r.y-G.y.y.y*r.x,
		G.y.z.y*r.z-G.y.z.z*r.y, G.y.z.z*r.x-G.y.z.x*r.z, G.y.z.x*r.y-G.y.z.y*r.x,//
		
		G.z.x.y*r.z-G.z.x.z*r.y, G.z.x.z*r.x-G.z.x.x*r.z, G.z.x.x*r.y-G.z.x.y*r.x,
		G.z.y.y*r.z-G.z.y.z*r.y, G.z.y.z*r.x-G.z.y.x*r.z, G.z.y.x*r.y-G.z.y.y*r.x,
		G.z.z.y*r.z-G.z.z.z*r.y, G.z.z.z*r.x-G.z.z.x*r.z, G.z.z.x*r.y-G.z.z.y*r.x
	);
	//以下たぶん間違い
	//return Tensor3x3x3(
	//	G.x.x.y*r.z-G.x.x.z*r.y, G.x.y.y*r.z-G.x.y.z*r.y, G.x.z.y*r.z-G.x.z.z*r.y,
	//	G.y.x.y*r.z-G.y.x.z*r.y, G.y.y.y*r.z-G.y.y.z*r.y, G.y.z.y*r.z-G.y.z.z*r.y,
	//	G.z.x.y*r.z-G.z.x.z*r.y, G.z.y.y*r.z-G.z.y.z*r.y, G.z.z.y*r.z-G.z.z.z*r.y,//
	//	
	//	G.x.x.z*r.x-G.x.x.x*r.z, G.x.y.z*r.x-G.x.y.x*r.z, G.x.z.z*r.x-G.x.z.x*r.z,
	//	G.y.x.z*r.x-G.y.x.x*r.z, G.y.y.z*r.x-G.y.y.x*r.z, G.y.z.z*r.x-G.y.z.x*r.z,
	//	G.z.x.z*r.x-G.z.x.x*r.z, G.z.y.z*r.x-G.z.y.x*r.z, G.z.z.z*r.x-G.z.z.x*r.z,//
	//	
	//	G.x.x.x*r.y-G.x.x.y*r.x, G.x.y.x*r.y-G.x.y.y*r.x, G.x.z.x*r.y-G.x.z.y*r.x,
	//	G.y.x.x*r.y-G.y.x.y*r.x, G.y.y.x*r.y-G.y.y.y*r.x, G.y.z.x*r.y-G.y.z.y*r.x,
	//	G.z.x.x*r.y-G.z.x.y*r.x, G.z.y.x*r.y-G.z.y.y*r.x, G.z.z.x*r.y-G.z.z.y*r.x
	//);
}

inline Tensor3x3x3 operator^(const Vector3d& r,const Tensor3x3x3& G)
{
	return Tensor3x3x3(
		r.y*G.z.x.x-r.z*G.y.x.x, r.y*G.z.x.y-r.z*G.y.x.y, r.y*G.z.x.z-r.z*G.y.x.z,
		r.y*G.z.y.x-r.z*G.y.y.x, r.y*G.z.y.y-r.z*G.y.y.y, r.y*G.z.y.z-r.z*G.y.y.z,
		r.y*G.z.z.x-r.z*G.y.z.x, r.y*G.z.z.y-r.z*G.y.z.y, r.y*G.z.z.z-r.z*G.y.z.z,

		r.z*G.x.x.x-r.x*G.z.x.x, r.z*G.x.x.y-r.x*G.z.x.y, r.z*G.x.x.z-r.x*G.z.x.z,
		r.z*G.x.y.x-r.x*G.z.y.x, r.z*G.x.y.y-r.x*G.z.y.y, r.z*G.x.y.z-r.x*G.z.y.z,
		r.z*G.x.z.x-r.x*G.z.z.x, r.z*G.x.z.y-r.x*G.z.z.y, r.z*G.x.z.z-r.x*G.z.z.z,

		r.x*G.y.x.x-r.y*G.x.x.x, r.x*G.y.x.y-r.y*G.x.x.y, r.x*G.y.x.z-r.y*G.x.x.z,
		r.x*G.y.y.x-r.y*G.x.y.x, r.x*G.y.y.y-r.y*G.x.y.y, r.x*G.y.y.z-r.y*G.x.y.z,
		r.x*G.y.z.x-r.y*G.x.z.x, r.x*G.y.z.y-r.y*G.x.z.y, r.x*G.y.z.z-r.y*G.x.z.z
	);
}

#endif //_TENSOR_3X3X3_H_

