#ifndef _STOKES_GREEN_FUNC_H_
#define _STOKES_GREEN_FUNC_H_

#include "../Tensor3x3x3.h"

namespace stokes_green_func{

Tensor3x3 get_free_G(const Vector3d& x,const Vector3d& y);
Tensor3x3x3 get_free_T(const Vector3d& x,const Vector3d& y);

Tensor3x3 get_free_T(const Vector3d& x,const Vector3d& y,const Vector3d& ny);

Vector3d get_free_G(const Vector3d& x,const Vector3d& y,const Vector3d& fy);
Vector3d get_free_T(const Vector3d& x,const Vector3d& y,const Vector3d& ny,const Vector3d& fy);


}

#endif // _STOKES_GREEN_FUNC_H_
