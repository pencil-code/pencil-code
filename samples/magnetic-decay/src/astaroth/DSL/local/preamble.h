const int NGHOST_VAL = 3
#include "../stdlib/math"
#include "../stdlib/derivs.h"
#include "../stdlib/operators.h"
#include "fieldecs.h"
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
