const int NGHOST_VAL = 3

#define AC_maux__mod__cparam maux__mod__cparam
#define iphF_UU F_PHIUU

#define AC_n2__mod__cparam AC_nx+NGHOST_VAL+1 
#define AC_m2__mod__cparam AC_ny+NGHOST_VAL+1 
#define AC_l2__mod__cparam AC_nz+NGHOST_VAL+1 

#define AC_NGHOST_VAL__mod__cparam NGHOST_VAL

#define AC_mx AC_mlocal.x
#define AC_my AC_mlocal.y
#define AC_mz AC_mlocal.z

#define AC_nx AC_nlocal.x
#define AC_ny AC_nlocal.y
#define AC_nz AC_nlocal.z

#define AC_nxgrid AC_ngrid.x
#define AC_nygrid AC_ngrid.y
#define AC_nzgrid AC_ngrid.z

#define AC_dsx AC_ds.x
#define AC_dsy AC_ds.y
#define AC_dsz AC_ds.z

#define n__mod__cdata (vertexIdx.z+1)
#define m__mod__cdata (vertexIdx.y+1)
#define AC_n__mod__cdata (vertexIdx.z+1)
#define AC_m__mod__cdata (vertexIdx.y+1)

#include "../stdlib/math"
#include "../stdlib/derivs.h"
#include "../stdlib/operators.h"
#define AC_NGHOST__mod__cparam nghost
#define REAL_MAX AC_REAL_MAX
//TP: nphis1 and nphis2 don't actually work. simply declared to compile the code
//
#define IN_DSL 1
#include "cparam.h"
#include "../../../cparam_c.h"
#include "../../../cparam_pencils.inc_c.h"
#include "PC_modulepardecs.h"
#include "../stdlib/integrators.h"

#include "fieldecs.h"

rk_intermediate(Field f, real df, int step_num, real dt)
{
	if(AC_itorder__mod__cdata == 1)
		return rk1_intermediate(f,df,step_num,dt)
	else if(AC_itorder__mod__cdata == 2)
		return rk2_intermediate(f,df,step_num,dt)
	else if(AC_itorder__mod__cdata == 3)
		return rk3_intermediate(f,df,step_num,dt)
	else
		return 0.0;
}
rk_intermediate(Field3 f, real3 df, int step_num, real dt)
{
	return real3(
			rk_intermediate(f.x,df.x,step_num,dt),
			rk_intermediate(f.y,df.y,step_num,dt),
			rk_intermediate(f.z,df.z,step_num,dt)
			)
}
rk_final(Field f, int step_num)
{
	if(AC_itorder__mod__cdata == 1)
		return rk1_final(f,step_num)
	else if(AC_itorder__mod__cdata == 2)
		return rk2_final(f,step_num)
	else if(AC_itorder__mod__cdata == 3)
		return rk3_final(f,step_num)
	else 
		return 0.0;
}
rk_final(Field3 f, int step_num)
{
	return real3(
			rk_final(f.x,step_num),
			rk_final(f.y,step_num),
			rk_final(f.z,step_num)
		    )
}
