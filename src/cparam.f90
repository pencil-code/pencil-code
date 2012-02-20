! $Id$
!
!  Module containing global parameters (constants).
!
module Cparam
!
  implicit none
!
  include 'cparam.local'
!
  integer, parameter :: nx=nxgrid/nprocx,ny=nygrid/nprocy,nz=nzgrid/nprocz
  integer, parameter :: nprocxy=nprocx*nprocy
  integer, parameter :: nprocyz=nprocy*nprocz
  integer, parameter :: nprocxz=nprocx*nprocz
!
  include 'cparam.inc'
!
  include 'cparam_pencils.inc'
!
!  Derived and fixed parameters.
!
  integer, parameter :: mfarray=mvar+maux+mglobal+mscratch
  integer, parameter :: mcom=mvar+maux_com
!
  integer, parameter :: ikind8=selected_int_kind(14) ! 8-byte integer kind
  integer(KIND=ikind8), parameter :: nw=nx*ny*nz
!
  integer, parameter :: mx=nx+2*nghost,l1=1+nghost,l2=mx-nghost
  integer, parameter :: my=ny+2*nghost,m1=1+nghost,m2=my-nghost
  integer, parameter :: mz=nz+2*nghost,n1=1+nghost,n2=mz-nghost
  integer, parameter :: mxgrid=nxgrid+2*nghost
  integer, parameter :: mygrid=nygrid+2*nghost
  integer, parameter :: mzgrid=nzgrid+2*nghost
  integer, parameter :: mw=mx*my*mz
  integer(KIND=ikind8), parameter :: nwgrid=int(nxgrid,kind=ikind8)* &
                    int(nygrid,kind=ikind8)*int(nzgrid,kind=ikind8)
!
  integer, parameter :: l1i=l1+nghost-1,l2i=l2-nghost+1
  integer, parameter :: m1i=m1+nghost-1,m2i=m2-nghost+1
  integer, parameter :: n1i=n1+nghost-1,n2i=n2-nghost+1
!
  integer, parameter :: nrcyl=nx/2
  integer, parameter :: nrcylrun=max(nx/20,1)
!
!  Array dimension for reduce operation (maxima and sums).
!  Use here symbol mreduce, use nreduce in call.
!
  integer, parameter :: mreduce=6
  integer :: ip=14
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
  integer, parameter :: fnlen=135,intlen=21,bclen=3,labellen=25,linelen=256
  integer, parameter :: datelen=30,max_col_width=30,nscbc_len=24
!
!  Number of slots in initlnrho etc.
!
  integer, parameter :: ninit=5
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
  real, parameter :: max_real=huge(0.0)
!
!  A marker value that is highly unlikely ("impossible") to ever occur
!  during a meaningful run: use the highest possible number.
!  TODO: 'impossible' should be deleted. Initialization with NaNs should be
!        activated via a compiler flag. Then testing against NaN can be done
!        by using the 'is_nan' function found in the syscalls module.
!        This will require many changes in many files. Any current usage of
!        'impossible' for integers must be replaced eg. by '...=-max_int'.
!        (Bourdin.KIS)
!
  real, parameter :: impossible=3.9085e37
  integer, parameter :: impossible_int=max_int/100
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
  integer, parameter :: ilabel_max_dt=-3,ilabel_max_neg=-4
  integer, parameter :: ilabel_max_reciprocal=-5
  integer, parameter :: ilabel_integrate=3,ilabel_surf=4
  integer, parameter :: ilabel_sum_par=5,ilabel_sum_sqrt_par=6
  integer, parameter :: ilabel_sum_weighted=7,ilabel_sum_weighted_sqrt=8
  integer, parameter :: ilabel_sum_lim=9
!
!  pi and its derivatives.
!
  real, parameter :: pi=3.14159265358979323846264338327950D0
  real, parameter :: pi_1=1./pi,pi4_1=pi**(-4),pi5_1=pi**(-5)
  real, parameter :: sqrtpi=1.77245385090551602729816748334115D0
  real, parameter :: sqrt2=1.41421356237309504880168872420970D0
  real, parameter :: four_pi_over_three=4.0/3.0*pi
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
  ! Better express R_cgs as a derive quantity (i.e. don't define here...)
  ! (Not done yet since it breaks the interstellar test)
  !double precision, parameter :: R_cgs=k_B_cgs/m_u_cgs     ! [erg/K]
  double precision, parameter :: R_cgs=8.3144D7     ! [erg/K]
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
  logical, parameter :: NO_WARN=.false.
!
!  Data structure used to gather slice information from the various modules.
!
  type slice_data
    character (LEN=30) :: name
    integer :: index
    logical :: ready
    real, pointer, dimension (:,:) :: xy
    real, pointer, dimension (:,:) :: xz
    real, pointer, dimension (:,:) :: yz
    real, pointer, dimension (:,:) :: xy2
    real, pointer, dimension (:,:) :: xy3
    real, pointer, dimension (:,:) :: xy4
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
endmodule Cparam
