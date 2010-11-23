! $Id$
!
! Global particle variables
!
module Particles_cdata
!
  use Cdata
!
  implicit none
!
  public
!
  real, parameter :: npar_per_cell=npar/(1.*nwgrid)
  integer, parameter :: mspvar=mpvar+1 
!
  real :: rp_int=-impossible,rp_ext=-impossible
  real :: dsnap_par_minor=0.0, dsnap_par=0.0
  real :: rhops=1.0e10, rhops1=1.0d-10
  real :: mp_swarm=0.0, np_swarm=0.0, rhop_swarm=0.0
  real :: four_pi_rhops_over_three
  real :: np_const=0.0, rhop_const=0.0
  integer, dimension (mpar_loc) :: ipar
  integer, dimension (nspar)    :: ipar_nbody
  integer, dimension (npar_species) :: ipar_fence_species=0
!
  integer :: npvar=0, npar_loc=0, mspar=0, npar_total=0
  integer :: ixp=0, iyp=0, izp=0, ivpx=0, ivpy=0, ivpz=0, iap=0
  integer :: inpswarm=0, irhopswarm=0
  integer :: ipsx=0, ipsy=0, ipsz=0
  integer :: iupx=0, iupy=0, iupz=0
  integer :: ipviscx=0, ipviscy=0, ipviscz=0
  integer :: inp=0
  integer :: idiag_nmigmax=0, npart_radii=0
  integer, dimension(ny*nz) :: npar_imn, k1_imn, k2_imn
  logical :: linterpolate_spline=.true.
  logical :: lparticlemesh_cic=.false., lparticlemesh_tsc=.false.
  logical :: linterp_reality_check=.false., lmigration_redo=.false.
  logical :: lnocalc_np=.false., lnocalc_rhop=.false.
  logical :: lmigration_real_check=.true.
  logical :: lcheck_exact_frontier=.false.
  character (len=2*bclen+1) :: bcpx='p',bcpy='p',bcpz='p'
  character (len=2*bclen+1) :: bcspx='p',bcspy='p',bcspz='p'
!
  logical :: lshepherd_neighbour=.false.
  logical :: lrandom_particle_pencils=.false., lrandom_particle_blocks=.false.
  logical :: linsert_particles_continuously=.false.
  integer, dimension (nx) :: kshepherd
  integer, allocatable, dimension (:) :: kneighbour
!
  type quant_interp_penc
!
!  Interpolation toggles:
!
    logical :: luu,loo,lTT,lrho
!
!  Interpolation policies:
!
    integer :: pol_uu, pol_oo, pol_TT, pol_rho
  end type quant_interp_penc
!  
  type(quant_interp_penc) :: interp
!
!  Interpolated quantities: moved outside type to conform to
!  the f90 standard.
!
    real, dimension(:,:), allocatable :: interp_uu, interp_oo
    real, dimension(:), allocatable :: interp_TT, interp_rho
!
!  Interpolation policies:
!     cic := cloud in cell (linear interpolation)
!     tsc := triangular shaped cloud (spline or
!               quadratic interpolation; depends on linterpolate_spline flag)
!     ngp := nearest grid point
!
  integer, parameter :: cic=0
  integer, parameter :: tsc=1
  integer, parameter :: ngp=2
!
endmodule Particles_cdata
