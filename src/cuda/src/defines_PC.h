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

#define T_STOP_FORCING 1e99
#define FORCING     1e-5
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

