//EulerAngle.h
//オイラー角に機能を持たせ、
//4元数中心に作ってある。

#ifndef _EULER_ANGLE_H_
#define _EULER_ANGLE_H_

#include <cmath>
#include <iostream>
using namespace std;

#include "Vector3d.h"
#include "Tensor3x3.h"

class EulerAngle{
public:
	EulerAngle(double ph=0.0,double th=0.0,double ps=0.0);
	EulerAngle(double xi,double eta,double zeta,double chi);
	EulerAngle(const EulerAngle &v);
	~EulerAngle();
	
	double xi,eta,zeta,chi;//4元数
	Vector3d x_vector()const;//回転前の(1,0,0)をメンバの4元数でベクトルを回転
	Vector3d y_vector()const;//回転前の(0,1,0)をメンバの4元数でベクトルを回転
	Vector3d z_vector()const;//回転前の(0,0,1)をメンバの4元数でベクトルを回転
	Vector3d rotate(const Vector3d& r)const;//rをメンバの角度で座標系を変換
	Vector3d rotate_(const Vector3d& r)const;//上の逆変換
	Tensor3x3 rotateMat()const;//メンバの角度の回転行列
	Tensor3x3 rotateMat_()const;//上の逆変換
	EulerAngle Q_convert(const Vector3d& angle_velocity)const;//角速度を４元数での速度に変換
	EulerAngle Q_convert(const EulerAngle& angle)const;//4元数を４元数に変換
	
	void set(double ph,double th,double ps);
	void set(double xi,double eta,double zeta,double chi);
	void set(const EulerAngle &v);

	void euler_to_four(double ph,double th,double ps);//オイラー角から4元数を計算する。
	void renormalize();//4元数の2乗を大きさ1にリスケールする
	
	EulerAngle& operator=(const EulerAngle& a);
	EulerAngle& operator+=(const EulerAngle& v);
	EulerAngle& operator-=(const EulerAngle& v);
	friend EulerAngle operator + (const EulerAngle& v1, const EulerAngle& v2);
	friend EulerAngle operator - (const EulerAngle& v1, const EulerAngle& v2);
	friend double operator * (const EulerAngle& v1, const EulerAngle& v2);
	friend EulerAngle operator * (const EulerAngle& v, const double& s);
	friend EulerAngle operator * (const double& s, const EulerAngle &v);
	friend EulerAngle operator / (const EulerAngle& v, const double& s);

	static void parenOn()	{ outputParen = true; }
	static void parenOff()	{ outputParen = false; }
	static bool parenIsOn()	{ return outputParen; }
private:
	static bool outputParen;	// default value is true
};
ostream& operator << (ostream& os, const EulerAngle& v);

inline EulerAngle& EulerAngle::operator=(const EulerAngle& a){
	xi=a.xi;eta=a.eta;zeta=a.zeta;chi=a.chi;
	return *this;
}

inline EulerAngle operator + (const EulerAngle& v1, const EulerAngle& v2)
{
	return EulerAngle(v1.xi+v2.xi,v1.eta+v2.eta,v1.zeta+v2.zeta,v1.chi+v2.chi);
}

inline EulerAngle operator - (const EulerAngle& v1, const EulerAngle& v2)
{
	return EulerAngle(v1.xi-v2.xi,v1.eta-v2.eta,v1.zeta-v2.zeta,v1.chi-v2.chi);
}

inline double operator * (const EulerAngle& v1, const EulerAngle& v2)
{
	return v1.xi*v2.xi+v1.eta*v2.eta+v1.zeta*v2.zeta+v1.chi*v2.chi;
}

inline EulerAngle& EulerAngle::operator += (const EulerAngle& v)
{
	xi+=v.xi;eta+=v.eta;zeta+=v.zeta;chi+=v.chi;
	return *this;
}

inline EulerAngle& EulerAngle::operator -= (const EulerAngle& v)
{
	xi-=v.xi;eta-=v.eta;zeta-=v.zeta;chi-=v.chi;
	return *this;
}

inline EulerAngle operator * (const EulerAngle& v, const double& s)
{
	return EulerAngle(s*v.xi, s*v.eta, s*v.zeta,s*v.chi);
}

inline EulerAngle operator * (const double& s, const EulerAngle& v)
{
	return EulerAngle(s*v.xi, s*v.eta, s*v.zeta,s*v.chi);
}

inline EulerAngle operator / (const EulerAngle& v, const double& s)
{
	return EulerAngle(v.xi/s, v.eta/s, v.zeta/s, v.chi/s);
}


#endif
