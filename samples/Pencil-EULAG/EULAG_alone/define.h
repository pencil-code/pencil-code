#define PARALLEL 2
#define PVM_IO 0
#define ANALIZE 0
#define GKS 0
#define V5D 0
#define PLOTPL 0
#define COLORPL 0
#define TURBPL 0
#define SPCTPL 0
#define VORTPL 1
#define ENERGY 0
#define ENERGY2 1
#define CRAYT3D 0
#define CRAYPVP 0
#define CRAYT3E 0
#define SGI_O2K 0
#define HP 0
#define WORKS 0
#define FUJI_VPP 0
#define IBM 0
#define LNX 0
#define CPQ 0
#define PLE 2
#define SEMILAG 0  /* 0=EULERIAN, 1=SEMI-LAGRANGIAN                   */
#define MOISTMOD 0 /* 0=DRY, 1=WARM MOIST, 2=ICE A+B (NOT READY YET)  */
#define MHD 0      /* 0=Hydro, 1=Magnetohydrodynamics       */
#define J3DIM 1    /* 0=2D MODEL, 1=3D MODEL                          */
#define SGS 0      /* 0=NO DIFFUSION, 1=OLD DISSIP, 2=DISSIP AMR      */
#define SUMR16 0   /* 0=REAL*8, 1=REAL*16 for GLOBAL SUMS             */
#define ISEND 2    /* 1=SENDRECV, 2=IRECV+SEND, 3=ISEND+IRECV    */
#define PRECF 0    /* 0=PRECON_DF or PRECON_BCZ, 1=PRECON_F           */
#define POLES   1  /* 0=NO FLOW OVER POLES, 1=FLOW OVER POLES ON      */
#define IORFLAG 0  /* 0=tape in serial mode, 1=tape in parallel mode  */
#define IOWFLAG 0  /* 0=tape in serial mode, 1=tape in parallel mode  */
#define IORTAPE 2  /* 1=default, 2=single, 3=double precision         */
#define IOWTAPE 2  /* 1=default, 2=single, 3=double precision         */
#define TIMEPLT 0  /* 1=time counts                                   */
