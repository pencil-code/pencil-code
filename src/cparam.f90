! $Id$
!
!  Module containing global parameters (constants).
!
module Cparam
!
  implicit none
!
  integer, parameter :: ikind8=selected_int_kind(14)  ! 8-byte integer kind
  integer, parameter :: ikind1=selected_int_kind(2)   ! 1-byte real kind
  integer, parameter :: rkind8=selected_real_kind(12) ! 8-byte real kind
  integer, parameter :: rkind4=selected_real_kind(6)  ! 4-byte real kind
  integer, parameter :: rkind16 = selected_real_kind(33, 4931) ! 16-byte real kind
!
  include 'cparam.local'
!
  integer, parameter :: nx=nxgrid/nprocx,ny=nygrid/nprocy,nz=nzgrid/nprocz
  integer, parameter :: nxygrid=nxgrid*nygrid,nxzgrid=nxgrid*nzgrid,nyzgrid=nygrid*nzgrid
  integer, parameter :: nprocxy=nprocx*nprocy
  integer, parameter :: nprocyz=nprocy*nprocz
  integer, parameter :: nprocxz=nprocx*nprocz
  integer, parameter :: n_forcing_cont_max=2
  character, dimension(3), parameter :: coornames=(/'x','y','z'/)
  logical, dimension(3), parameter :: lactive_dimension = (/ nxgrid > 1, nygrid > 1, nzgrid > 1 /)
  integer, parameter :: dimensionality=min(nxgrid-1,1)+min(nygrid-1,1)+min(nzgrid-1,1)
  integer, dimension(3), parameter :: grid_dims=(/nx,ny,nz/)
!
  include 'cparam.inc'
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
  integer :: l2=mx-nghost
  integer, parameter :: my=ny+2*nghost,m1=1+nghost
  integer :: m2=my-nghost
  integer, parameter :: mz=nz+2*nghost,n1=1+nghost
  integer :: n2=mz-nghost
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
  integer :: l2i=mx-2*nghost+1
  integer, parameter :: m1i=m1+nghost-1
  integer :: m2i=my-2*nghost+1
  integer, parameter :: n1i=n1+nghost-1
  integer :: n2i=mz-2*nghost+1
!
  integer, parameter :: nrcyl=nxgrid/2
  integer, parameter :: nrcylrun=max(nx/20,1)
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
  integer, parameter :: max_int=huge(0)
  real, parameter :: huge_real=huge(0.)
  double precision, parameter :: huge_double=huge(0.d0)
  real, parameter :: max_real=huge_real/10.    ! division necessary as INTEL compiler considers
                                               ! huge(0.) illegal when reading it from a namelist
!
!  Tiny and huge numbers.
!
  real, parameter :: one_real=1.0
  real, parameter :: epsi=5*epsilon(one_real),tini=5*tiny(one_real)
  real, parameter :: huge1=0.2*huge_real
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
  integer, parameter :: ilabel_max=-1,ilabel_sum=1,ilabel_save=11
  integer, parameter :: ilabel_max_sqrt=-2,ilabel_sum_sqrt=2
  integer, parameter :: ilabel_sum_log10=10
  integer, parameter :: ilabel_max_dt=-3,ilabel_max_neg=-4
  integer, parameter :: ilabel_max_reciprocal=-5
  integer, parameter :: ilabel_integrate=3,ilabel_integrate_sqrt=30, ilabel_integrate_log10=40
  integer, parameter :: ilabel_surf=4
  integer, parameter :: ilabel_sum_par=5,ilabel_sum_sqrt_par=6, ilabel_sum_log10_par=20, ilabel_sum_plain=21
  integer, parameter :: ilabel_sum_weighted=7,ilabel_sum_weighted_sqrt=8
  integer, parameter :: ilabel_sum_lim=9,ilabel_complex=100
!
!  pi and its derivatives.
!
  real, parameter :: pi=3.14159265358979323846264338327950d0
  real, parameter :: pi_1=1./pi,pi4_1=pi**(-4),pi5_1=pi**(-5)
  real, parameter :: sqrtpi=1.77245385090551602729816748334115d0
  real, parameter :: sqrt2=1.41421356237309504880168872420970d0
  real, parameter :: sqrt2pi=sqrt2*sqrtpi
  real, parameter :: four_pi_over_three=4.0/3.0*pi
  real, parameter :: onethird=1./3., twothird=2./3., fourthird=4./3., onesixth=1./6.
  real, parameter :: one_over_sqrt3=0.577350269189625764509148780501958d0
  real, parameter :: twopi = 6.2831853071795864769252867665590d0
  real, parameter :: dtor = pi/180.d0
!
  real, parameter :: lntwo=0.69314718055995d0
!
!  first zeros of Bessel functions of order 0 and 1
!  k2bessel0 is the second zero of Bessel function of order 0
!
  real, parameter :: k1bessel0=2.4048255577, k1bessel1=3.8317060
  real, parameter :: k2bessel0=5.5200781
!
!  Physical constants, taken from
!  http://physics.nist.gov/cuu/Constants/index.html.
!
  double precision, parameter :: hbar_cgs=1.054571596d-27  ! [erg*s]
  double precision, parameter :: k_B_cgs=1.3806505d-16     ! [erg/K]
  double precision, parameter :: m_u_cgs=1.66053886d-24    ! [g]
  double precision, parameter :: mu0_cgs=4*pi              ! [cgs]
  ! Better express R_cgs as a derived quantity (i.e. don't define here...)
  ! (Not done yet since it breaks the interstellar test)
  !double precision, parameter :: R_cgs=k_B_cgs/m_u_cgs    ! [erg/g/K]
  double precision, parameter :: R_cgs=8.3144D7            ! [erg/g/K]
  ! It would be better to specify the following masses in units of m_u:
  double precision, parameter :: m_p_cgs=1.67262158d-24    ! [g]
  double precision, parameter :: m_e_cgs=9.10938188d-28    ! [g]
  double precision, parameter :: m_H_cgs=m_e_cgs+m_p_cgs   ! [g]
  double precision, parameter :: eV_cgs=1.602176462d-12    ! [erg]
  double precision, parameter :: sigmaSB_cgs=5.670400d-5   ! [erg/cm^2/s/K^4]
! unclear source (probably just guessing?)
  double precision, parameter :: sigmaH_cgs=4.d-17         ! [cm^2]
  double precision, parameter :: kappa_es_cgs=3.4d-1       ! [cm^2/g]
  double precision, parameter :: c_light_cgs=2.99792458d10 ! [cm/s]
  double precision, parameter :: G_Newton_cgs=6.6742d-8    ! [cm3/g/s2]
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
endmodule Cparam
