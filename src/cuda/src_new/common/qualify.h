#pragma once

#define _qualified(module,name,pre,in,suf) pre##module##in##name##suf 
#define qualified(module,name,pre,in,suf) _qualified(module,name,pre,in,suf)

/*
#define density_push2c   qualified(density,push2c,MODPRE,MODIN,MODSUF)
#define hydro_push2c     qualified(hydro,push2c,MODPRE,MODIN,MODSUF)
#define viscosity_push2c qualified(viscosity,push2c,MODPRE,MODIN,MODSUF)
#define forcing_push2c   qualified(forcing,push2c,MODPRE,MODIN,MODSUF)
#define eos_push2c       qualified(equationofstate,push2c,MODPRE,MODIN,MODSUF)
*/
