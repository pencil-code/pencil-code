#pragma once

#define DX dx
#define DY dy
#define DZ dz

#define DOMAIN_SIZE_X lxyz[0]
#define DOMAIN_SIZE_Z lxyz[1]
#define DOMAIN_SIZE_Y lxyz[2]

#define XORIG xyz0[0]
#define YORIG xyz0[1]
#define ZORIG xyz0[2]

#define CDT  cdt
#define CDTV cdtv

#define T_STOP_FORCING tforce_stop
#define FORCING     force
#define KK1         4.5
#define KK2         5.5
#define KMAX        10.0
#define DKX         1.0
#define DKY         1.0
#define DKZ         1.0

#define LFORCING  lforcing
#define LSHEAR    lshear
//#define USE_HYDRO lhydro
#define LCORIOLIS 0

#define Q_SHEAR   qshear
#define OMEGA     omega
#define NU_VISC   nu
#define CS2_SOUND cs2
#define CS_SOUND  sqrt(cs2)

#define _qualified(module,name,pre,in,suf) pre##module##in##name##suf 
#define qualified(module,name,pre,in,suf) _qualified(module,name,pre,in,suf)

#define save_name        qualified(diagnostics,save_name,MODPRE,MODIN,MODSUF)

#define density_push2c   qualified(density,push2c,MODPRE,MODIN,MODSUF)
#define hydro_push2c     qualified(hydro,push2c,MODPRE,MODIN,MODSUF)
#define viscosity_push2c qualified(viscosity,push2c,MODPRE,MODIN,MODSUF)
#define forcing_push2c   qualified(forcing,push2c,MODPRE,MODIN,MODSUF)
#define eos_push2c       qualified(equationofstate,push2c,MODPRE,MODIN,MODSUF)

