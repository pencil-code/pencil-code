! $Id$
!
!  Module containing global parameters (constants).
!
module Cparam
!
  implicit none
!
  integer, parameter :: ikind8=selected_int_kind(14)  ! 8-byte integer kind
  integer, parameter :: ikind4=selected_int_kind(9)   ! 4-byte integer kind
  integer, parameter :: ikind1=selected_int_kind(2)   ! 1-byte integer kind
  integer, parameter :: rkind8=selected_real_kind(12) ! 8-byte real kind
  integer, parameter :: rkind4=selected_real_kind(6)  ! 4-byte real kind
 ! integer, parameter :: rkind16 = selected_real_kind(33, 4931) ! 16-byte real kind - not accepted by all compilers
  integer, parameter :: rkind16 = rkind8
!
  include 'cparam.local'
!
  integer, parameter :: nx=nxgrid/nprocx,ny=nygrid/nprocy,nz=nzgrid/nprocz,nyz=ny*nz
  integer, parameter :: max_n = max(nx,max(ny,nz))
  integer, parameter :: nxygrid=nxgrid*nygrid,nxzgrid=nxgrid*nzgrid,nyzgrid=nygrid*nzgrid
  integer, parameter :: nprocxy=nprocx*nprocy
  integer, parameter :: nprocyz=nprocy*nprocz
  integer, parameter :: nprocxz=nprocx*nprocz
  integer, parameter :: n_forcing_cont_max=2
  integer, parameter :: ndustspec0=8
  character, dimension(3), parameter :: coornames=(/'x','y','z'/)
  character(LEN=2), dimension(12), parameter :: compnames=(/'x ','y ','z ','xx','xy','xz','yx','yy','yz','zx','zy','zz'/)
  integer, dimension(6),parameter :: compinds_6=(/1,2,3,5,6,9/)
  logical, dimension(3), parameter :: lactive_dimension = (/ nxgrid > 1, nygrid > 1, nzgrid > 1 /)
  logical, dimension(3), parameter :: linactive_dimension = (/ nxgrid == 1, nygrid == 1, nzgrid == 1 /)
  integer, parameter :: dimensionality=min(nxgrid-1,1)+min(nygrid-1,1)+min(nzgrid-1,1)
  integer, dimension(3), parameter :: grid_dims=(/nx,ny,nz/)
!
  include 'cparam.inc'
  logical, parameter :: lenergy=lentropy.or.ltemperature.or.lthermal_energy
!
  integer, parameter :: penc_name_len=16
!
  include 'cparam_pencils.inc'
!
!  Derived and fixed parameters.
!
! BEGIN CHANGE FOR DYNAMICAL ALLOCATION
  integer, parameter :: mfarray=mvar+maux+mglobal+mscratch
  integer, parameter :: mcom=mvar+maux_com
  integer, parameter :: mparray=mpvar+mpaux
  integer, parameter :: mpcom=mpvar+mpaux
  integer, parameter :: mqarray=mqvar+mqaux
! END CHANGE FOR DYNAMICAL ALLOCATION
!
  integer(KIND=ikind8), parameter :: nw=nx*ny*nz
!
!!!  integer, parameter :: mx=nx+2*nghost,l1=1+nghost,l2=mx-nghost
!!!  integer, parameter :: my=ny+2*nghost,m1=1+nghost,m2=my-nghost
!!!  integer, parameter :: mz=nz+2*nghost,n1=1+nghost,n2=mz-nghost
  integer, parameter :: mx=nx+2*nghost,l1=1+nghost
  integer, parameter :: my=ny+2*nghost,m1=1+nghost
  integer, parameter :: mz=nz+2*nghost,n1=1+nghost
  integer, parameter :: mxgrid=nxgrid+2*nghost
  integer, parameter :: mygrid=nygrid+2*nghost
  integer, parameter :: mzgrid=nzgrid+2*nghost
  integer, parameter :: mw=mx*my*mz
  integer(KIND=ikind8), parameter :: nwgrid=int(nxgrid,kind=ikind8)* &
                                            int(nygrid,kind=ikind8)* &
                                            int(nzgrid,kind=ikind8)
!
!!!  integer, parameter :: l1i=l1+nghost-1,l2i=l2-nghost+1
!!!  integer, parameter :: m1i=m1+nghost-1,m2i=m2-nghost+1
!!!  integer, parameter :: n1i=n1+nghost-1,n2i=n2-nghost+1
  integer, parameter :: l1i=l1+nghost-1
  integer, parameter :: m1i=m1+nghost-1
  integer, parameter :: n1i=n1+nghost-1


  integer, parameter :: nrcyl=nxgrid/2
  integer, parameter :: nrcylrun=max(nx/20,1)
!
!  Number of bins for Pulsar Timing Array
!
  integer, parameter :: nbin_angular=19*2
!
!  Array dimension for reduce operation (maxima and sums).
!  Use here symbol mreduce, use nreduce in call.
!
  integer, parameter :: mreduce=6
!
!  Number of slots in initlnrho etc.
!
  integer, parameter :: ninit=5
!
!  Maximum number of diagnostics in print.in
!
  integer, parameter :: max_diagnostics=100
!
!  Name:          Maximum string length of a:
!  --------------------------------------------
!  fnlen          file name
!  intlen         integer (64 bit plus sign)
!  bclen          string for boundary condition
!  labellen       label (eg. initss, initaa)
!  linelen        line to be read in
!  datelen        date-and-time string
!  max_col_width  diagnostic column
!  nscbc_len      ?
!
  integer, parameter :: fnlen=135,intlen=21,bclen=3,labellen=40,linelen=256
  integer, parameter :: datelen=30,max_col_width=30,nscbc_len=24,fmtlen=30
!
!  Significant length of random number generator state.
!  Different compilers have different lengths:
!    NAG: 1, Compaq: 2, Intel: 47, SGI: 64, NEC: 256
!
  integer, parameter :: mseed=256
!
!  Predefine maximum possible numbers.
!
  integer(KIND=ikind4), parameter :: int_sgl=0
  integer, parameter :: max_int=huge(int_sgl)
  real, parameter :: huge_real=huge(0.)
  real(KIND=rkind8), parameter :: zero_double=0., huge_double=huge(zero_double)
  real, parameter :: max_real=huge_real/10.    ! division necessary as INTEL compiler considers
                                               ! huge(0.) illegal when reading it from a namelist
!
!  Tiny and huge numbers.
!
  real, parameter :: one_real=1.0
  real, parameter :: epsi=5*epsilon(one_real),tini=5*tiny(one_real)
  real, parameter :: huge1=0.2*huge_real
  real, parameter :: min_ts=1e-99  !PAR_DOC: minimum number to be displayed in time series
!
!  A marker value that is highly unlikely ("impossible") to ever occur
!  during a meaningful run: use a very large number.
!  We use numbers ~ 2 orders of magnitude below the maximum possible
!  values, as they may still get multiplied by some moderate numbers.
!
!  This value is a marker for some variable being uninitialized, and it is
!  tempting to replace the mechanism used here by NaN.
!  This may or may not work (need to reliably create NaN [including
!  REAL_PRECISION=double], some compilers seem to trap assignment of NaN
!  values, etc.
!  Also, there is no NaN concept for integers.
!
  real, parameter :: impossible=3.9085e37
  integer, parameter :: impossible_int=-max_int/100
!
! MPI
!
  integer, parameter :: root=0
!
!  Diagnostic variable types.
!
!  Values > 0 get maxed across all processors before any
!  transformation using mpi_reduce_max;
!  values < 0 get summed over all processors before any transformation
!  using mpi_reduce_sum;
!  value 0 causes the value simply to be used from the root processor.
!
  integer, parameter :: ilabel_max=-1,ilabel_sum=1,ilabel_save=0
  integer, parameter :: ilabel_max_sqrt=-2,ilabel_sum_sqrt=2
  integer, parameter :: ilabel_sum_log10=10, ilabel_sum_masked=11
  integer, parameter :: ilabel_max_dt=-3,ilabel_max_neg=-4
  integer, parameter :: ilabel_max_reciprocal=-5
  integer, parameter :: ilabel_integrate=3,ilabel_integrate_sqrt=30, ilabel_integrate_log10=40
  integer, parameter :: ilabel_surf=4
  integer, parameter :: ilabel_sum_par=5,ilabel_sum_sqrt_par=6, ilabel_sum_log10_par=20, ilabel_sum_plain=21
  integer, parameter :: ilabel_sum_weighted=7,ilabel_sum_weighted_sqrt=8
  integer, parameter :: ilabel_sum_lim=9,ilabel_complex=100
!
  real, parameter :: lntwo=0.69314718055995d0
!
!  first zeros of Bessel functions of order 0 and 1
!  k2bessel0 is the second zero of Bessel function of order 0
!
  real, parameter :: k1bessel0=2.4048255577, k1bessel1=3.8317060
  real, parameter :: k2bessel0=5.5200781
!
!  pi and its derivatives.
!
  real, parameter :: pi=3.14159265358979323846264338327950d0
  real, parameter :: pi_1=1./pi,pi4_1=(1.0)/(pi*pi*pi*pi),pi5_1=1.0/(pi*pi*pi*pi*pi)
  real, parameter :: sqrtpi=1.77245385090551602729816748334115d0
  real, parameter :: sqrt2=1.41421356237309504880168872420970d0
  real, parameter :: sqrt21=1./sqrt2
  real, parameter :: sqrt2pi=sqrt2*sqrtpi
  real, parameter :: four_pi_over_three=4.0/3.0*pi
  real, parameter :: onethird=1./3., twothird=2./3., fourthird=4./3., onesixth=1./6.
  real, parameter :: one_over_sqrt3=0.577350269189625764509148780501958d0
  real, parameter :: twopi = 6.2831853071795864769252867665590d0
  real, parameter :: dtor = pi/180.d0
!
!  Physical constants, taken from
!  http://physics.nist.gov/cuu/Constants/index.html.
!
  real(KIND=rkind8), parameter :: hbar_cgs=1.054571596d-27  ! [erg*s]
  real(KIND=rkind8), parameter :: k_B_cgs=1.3806505d-16     ! [erg/K]
  real(KIND=rkind8), parameter :: m_u_cgs=1.66053886d-24    ! [g]
  real(KIND=rkind8), parameter :: mu0_cgs=4*pi              ! [cgs]
  ! Better express R_cgs as a derived quantity (i.e. don't define here...)
  ! (Not done yet since it breaks the interstellar test)
  !real(KIND=rkind8), parameter :: R_cgs=k_B_cgs/m_u_cgs    ! [erg/g/K]
  real(KIND=rkind8), parameter :: R_cgs=8.3144D7            ! [erg/g/K]
  ! It would be better to specify the following masses in units of m_u:
  real(KIND=rkind8), parameter :: m_p_cgs=1.67262158d-24    ! [g]
  real(KIND=rkind8), parameter :: m_e_cgs=9.10938188d-28    ! [g]
  real(KIND=rkind8), parameter :: m_H_cgs=m_e_cgs+m_p_cgs   ! [g]
  real(KIND=rkind8), parameter :: eV_cgs=1.602176462d-12    ! [erg]
  real(KIND=rkind8), parameter :: sigmaSB_cgs=5.670400d-5   ! [erg/cm^2/s/K^4]
! unclear source (probably just guessing?)
  real(KIND=rkind8), parameter :: sigmaH_cgs=4.d-17         ! [cm^2]
  real(KIND=rkind8), parameter :: kappa_es_cgs=3.4d-1       ! [cm^2/g]
  real(KIND=rkind8), parameter :: c_light_cgs=2.99792458d10 ! [cm/s]
  real(KIND=rkind8), parameter :: G_Newton_cgs=6.6742d-8    ! [cm3/g/s2]
  real(KIND=rkind8), parameter :: density_scale_cgs=1.2435d21 ![cm] 403pc Reynolds 91, etc
  real(KIND=rkind8), parameter :: N_avogadro_cgs=6.022d23 ![1/mol]
  real(KIND=rkind8), parameter :: alpha_fine=7.2973525643d-3
  real(KIND=rkind8), parameter :: sigma_Thomson_cgs=6.652458732160d-25 ![cm^2]
  real(KIND=rkind8), parameter :: e_cgs=4.8032047d-10  ![statcoulombs]
  real(KIND=rkind8), parameter :: Chypercharge=41./12. !(nondimensional)
  real(KIND=rkind8), parameter :: mass_zboson=7.48e-18 !(in Mpl, not reduced)
  real(KIND=rkind8), parameter :: mass_zboson_GeV=91.2 ![in GeV]
!
  logical, parameter :: ALWAYS_FALSE=.false.
!
!  Data structure used to gather slice information from the various modules.
!
  type slice_data
    character (LEN=labellen) :: name
    integer :: index
    logical :: ready
    real, pointer, dimension (:,:) :: xy
    real, pointer, dimension (:,:) :: xz
    real, pointer, dimension (:,:) :: xz2
    real, pointer, dimension (:,:) :: yz
    real, pointer, dimension (:,:) :: xy2
    real, pointer, dimension (:,:) :: xy3
    real, pointer, dimension (:,:) :: xy4
    real, pointer, dimension (:,:) :: r
  endtype slice_data
!
!  Data structure used to allow module specific boundary conditions.
!
  type boundary_condition
    character (len=bclen) :: bcname
    integer :: ivar
    integer :: location
    logical :: done
    real :: value1
    real :: value2
  endtype boundary_condition
!
  integer, parameter :: iBC_X_TOP=1
  integer, parameter :: iBC_X_BOT=-1
  integer, parameter :: iBC_Y_TOP=2
  integer, parameter :: iBC_Y_BOT=-2
  integer, parameter :: iBC_Z_TOP=3
  integer, parameter :: iBC_Z_BOT=-3
  integer, parameter :: BOT=1, TOP=2, BOTH=3
!
!  Indices of rho, d rho/d x, d^2 rho/d x^2, d^6 rho/d x^6, d p/d x, s, d s/d x, &
!             d^2 s/d x^2, d^6 s/d x^6, d lnrho/d z in array reference_state.
!
  integer, parameter :: iref_rho=1, iref_grho=2, iref_d2rho=3, iref_d6rho=4, &
                        iref_gp=5, iref_s=6, iref_gs=7, iref_d2s=8, iref_d6s=9
  integer, parameter :: nref_vars=9
!
!  Symbolic constants for Yin-Yang grid.
!
  integer, parameter :: BILIN=1, BIQUAD=2, BICUB=3, QUADSPLINE=4, BIQUIN=5
!
!  Symbolic constants for Cubed Sphere grid.
!  The order of the patches is the same as in MATINS.
!
  integer, parameter :: XPLUS=1, YPLUS=2, XMINUS=3, YMINUS=4, ZPLUS=5, ZMINUS=6
!
  integer, parameter :: max_threads_possible = 200
  integer, parameter :: PERF_DIAGS=1, PERF_WSNAP=2, PERF_POWERSNAP=3, PERF_WSNAP_DOWN=4
  integer, parameter :: n_helperflags=4
  integer, parameter :: n_xy_specs_max=10,nk_max=10, nz_max=10
  integer, parameter :: n_pdfs_max=10, n_cspec_max=10
  integer, parameter :: mname=100
  integer, parameter :: mname_half=20

  include 'cparam_enum.h'

endmodule Cparam
