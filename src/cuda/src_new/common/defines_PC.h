/*
*   Contains the defines that are extracted from Pencil Code
*   and read to the config files in config.cc
*/
#pragma once

#define PC_DSX dx
#define PC_DSY dy
#define PC_DSZ dz

#define PC_DOMAIN_SIZE_X lxyz[0]
#define PC_DOMAIN_SIZE_Z lxyz[1]
#define PC_DOMAIN_SIZE_Y lxyz[2]

#define PC_XORIG xyz0[0]
#define PC_YORIG xyz0[1]
#define PC_ZORIG xyz0[2]

#define PC_CDT  cdt
#define PC_CDTV cdtv

#define PC_T_STOP_FORCING tforce_stop
#define PC_FORCING     force
#define PC_KK1         4.5
#define PC_KK2         5.5
#define PC_KMAX        10.0
#define PC_DKX         1.0
#define PC_DKY         1.0
#define PC_DKZ         1.0

#define PC_LFORCING  lforcing
#define PC_LSHEAR    lshear
#define PC_LHYDRO    lhydro
#define PC_LCORIOLIS 0

#define PC_Q_SHEAR   qshear
#define PC_OMEGA     omega
#define PC_NU_VISC   nu
#define PC_CS2_SOUND cs2
#define PC_CS_SOUND  sqrt(cs2)

//TODO: figure out what to do with these
//#define _qualified(module,name,pre,in,suf) pre##module##in##name##suf 
//#define qualified(module,name,pre,in,suf) _qualified(module,name,pre,in,suf)

//#define save_name        qualified(diagnostics,save_name,MODPRE,MODIN,MODSUF)

//#define density_push2c   qualified(density,push2c,MODPRE,MODIN,MODSUF)
//#define hydro_push2c     qualified(hydro,push2c,MODPRE,MODIN,MODSUF)
//#define viscosity_push2c qualified(viscosity,push2c,MODPRE,MODIN,MODSUF)
//#define forcing_push2c   qualified(forcing,push2c,MODPRE,MODIN,MODSUF)
//#define eos_push2c       qualified(equationofstate,push2c,MODPRE,MODIN,MODSUF)

