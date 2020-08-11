#include "EulerAngle.h"

bool EulerAngle::outputParen = true;

EulerAngle::EulerAngle(double ph,double th,double ps)
{
	euler_to_four(ph,th,ps);
}

EulerAngle::EulerAngle(double x,double e,double z,double c):xi(x),eta(e),zeta(z),chi(c)
{
}

EulerAngle::EulerAngle(const EulerAngle &v)
{
	xi=v.xi;eta=v.eta;zeta=v.zeta;chi=v.chi;
}

void EulerAngle::set(double ph,double th,double ps)
{
	euler_to_four(ph,th,ps);
}

void EulerAngle::set(double x,double e,double z,double c)
{
	xi=x;
	eta=e;
	zeta=z;
	chi=c;
}

void EulerAngle::set(const EulerAngle &v)
{
	xi=v.xi;
	eta=v.eta;
	zeta=v.zeta;
	chi=v.chi;
}

EulerAngle::~EulerAngle(){
}

void EulerAngle::euler_to_four(double phi,double theta,double psi){
	xi=sin(theta/2)*sin((psi-phi)/2);
	eta=sin(theta/2)*cos((psi-phi)/2);
	zeta=cos(theta/2)*sin((psi+phi)/2);
	chi=cos(theta/2)*cos((psi+phi)/2);
}

void EulerAngle::renormalize(){
	double factor;
	factor=sqrt(xi*xi+eta*eta+zeta*zeta+chi*chi);
	xi/=factor;
	eta/=factor;
	zeta/=factor;
	chi/=factor;
}

Vector3d EulerAngle::x_vector()const{
	double xx,yy,zz;
	xx=-xi*xi+eta*eta-zeta*zeta+chi*chi;
	yy=2*(-xi*eta+zeta*chi);
	zz=2*(eta*zeta+xi*chi);
	return Vector3d(xx,yy,zz);
}

Vector3d EulerAngle::y_vector()const{
	double xx,yy,zz;
	xx=-2*(zeta*chi+xi*eta);
	yy=xi*xi-eta*eta-zeta*zeta+chi*chi;
	zz=2*(-xi*zeta+eta*chi);
	return Vector3d(xx,yy,zz);
}

Vector3d EulerAngle::z_vector()const{
	double xx,yy,zz;
	xx=2*(eta*zeta-xi*chi);
	yy=-2*(eta*chi+xi*zeta);
	zz=-xi*xi-eta*eta+zeta*zeta+chi*chi;
	return Vector3d(xx,yy,zz);
}

Vector3d EulerAngle::rotate(const Vector3d& r)const{
	double xx,yy,zz;
	xx=r.x*(-xi*xi+eta*eta-zeta*zeta+chi*chi)
		+r.y*2*(zeta*chi-xi*eta)
		+r.z*2*(eta*zeta+xi*chi);
	yy=r.x*(-2*(xi*eta+zeta*chi))
		+r.y*(xi*xi-eta*eta-zeta*zeta+chi*chi)
		+r.z*2*(eta*chi-xi*zeta);
	zz=r.x*2*(eta*zeta-xi*chi)
		+r.y*(-2*(xi*zeta+eta*chi))
		+r.z*(-xi*xi-eta*eta+zeta*zeta+chi*chi);
	return Vector3d(xx,yy,zz);
}

Vector3d EulerAngle::rotate_(const Vector3d& r)const{
	double x,y,z;
	x=r.x*(-xi*xi+eta*eta-zeta*zeta+chi*chi)
		+r.y*(-2*(xi*eta+zeta*chi))
		+r.z*2*(eta*zeta-xi*chi);
	y=r.x*2*(zeta*chi-xi*eta)
		+r.y*(xi*xi-eta*eta-zeta*zeta+chi*chi)
		+r.z*(-2*(xi*zeta+eta*chi));
	z=r.x*2*(eta*zeta+xi*chi)
		+r.y*2*(eta*chi-xi*zeta)
		+r.z*(-xi*xi-eta*eta+zeta*zeta+chi*chi);
	return Vector3d(x,y,z);
}

Tensor3x3 EulerAngle::rotateMat()const{
	return Tensor3x3(
		-xi*xi+eta*eta-zeta*zeta+chi*chi,
			2*(zeta*chi-xi*eta),
				2*(eta*zeta+xi*chi),
		-2*(xi*eta+zeta*chi),
			xi*xi-eta*eta-zeta*zeta+chi*chi,
				2*(eta*chi-xi*zeta),
		2*(eta*zeta-xi*chi),
			-2*(xi*zeta+eta*chi),
				-xi*xi-eta*eta+zeta*zeta+chi*chi
	);
}

Tensor3x3 EulerAngle::rotateMat_()const{
	return Tensor3x3(
		-xi*xi+eta*eta-zeta*zeta+chi*chi,
			-2*(xi*eta+zeta*chi),
				2*(eta*zeta-xi*chi),
		2*(zeta*chi-xi*eta),
			xi*xi-eta*eta-zeta*zeta+chi*chi,
				-2*(xi*zeta+eta*chi),
		2*(eta*zeta+xi*chi),
			2*(eta*chi-xi*zeta),
				-xi*xi-eta*eta+zeta*zeta+chi*chi
	);
}

EulerAngle EulerAngle::Q_convert(const Vector3d& a_v)const{
	double xxi,eeta,zzeta,cchi;	
	xxi=0.5*(-zeta*a_v.x-chi*a_v.y+eta*a_v.z);
	eeta=0.5*(chi*a_v.x-zeta*a_v.y-xi*a_v.z);
	zzeta=0.5*(xi*a_v.x+eta*a_v.y+chi*a_v.z);
	cchi=0.5*(-eta*a_v.x+xi*a_v.y-zeta*a_v.z);
	return EulerAngle(xxi,eeta,zzeta,cchi);
}

EulerAngle EulerAngle::Q_convert(const EulerAngle& a)const{
	double xxi,eeta,zzeta,cchi;	
	xxi=0.5*(-zeta*a.xi-chi*a.eta+eta*a.zeta+xi*a.chi);
	eeta=0.5*(chi*a.xi-zeta*a.eta-xi*a.zeta+eta*a.chi);
	zzeta=0.5*(xi*a.xi+eta*a.eta+chi*a.zeta+zeta*a.chi);
	cchi=0.5*(-eta*a.xi+xi*a.eta-zeta*a.zeta+chi*a.chi);
	return EulerAngle(xxi,eeta,zzeta,cchi);
}

ostream& operator << (ostream& os, const EulerAngle& v)
{
	if( EulerAngle::parenIsOn() ) {
		os << "( " << v.xi << ", " << v.eta << ", " << v.zeta << ", " << v.chi << ") ";
	}
	else {
		os << v.xi << " " << v.eta << " " << v.zeta << " " << v.chi;
	}
	return os;
}
