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
  integer, parameter :: mspvar=mpvar+1
  integer, parameter :: lun_output=93
!
  real, parameter :: npar_per_cell=npar/(1.0*nwgrid)
  real :: rp_int=-impossible, rp_ext=-impossible
  real :: dsnap_par_minor=0.0, dsnap_par=0.0
  real :: rhopmat=1.0, rhopmat1=1.0, mpmat=0.0
  real :: mp_swarm=0.0, np_swarm=0.0, rhop_swarm=0.0
  real :: four_pi_rhopmat_over_three
  real :: np_const=0.0, rhop_const=0.0
  real :: energy_gain_shear_bcs=impossible
  real :: log_ap_min_dist=0.0, log_ap_max_dist=6.0
  real :: rsinkparticle_1=0.0
  real :: t_nextinsert=0. !The time at which new particles are going to be inserted.
!
  integer, dimension (nx) :: kshepherd
  integer, allocatable, dimension (:) :: kneighbour
  integer, dimension (mpar_loc) :: ipar
  integer, dimension (nspar) :: ipar_nbody
  integer, dimension (npar_species) :: ipar_fence_species=0
  integer, dimension(ny*nz) :: npar_imn, k1_imn, k2_imn
  integer :: npvar=0, npar_loc=0, mspar=0, npar_total=0
  integer :: ixp=0, iyp=0, izp=0, ivpx=0, ivpy=0, ivpz=0, iap=0, iaps=0
  integer :: inpswarm=0, irhopswarm=0
  integer :: ipsx=0, ipsy=0, ipsz=0
  integer :: ipss=0, ipst=0, ipxx=0, ipyy=0, ipzz=0
  integer :: iuup=0, iupx=0, iupy=0, iupz=0
  integer :: ipviscx=0, ipviscy=0, ipviscz=0
  integer :: inp=0, irhop=0, irhops=0
  integer :: idiag_nmigmax=0, npart_radii=0
  integer :: nbin_ap_dist=100
!
  logical :: linterpolate_spline=.true.
  logical :: lparticlemesh_cic=.true., lparticlemesh_tsc=.false.
  logical :: linterp_reality_check=.false., lmigration_redo=.false.
  logical :: lnocalc_np=.false., lnocalc_rhop=.false.
  logical :: lmigration_real_check=.true.
  logical :: lcheck_exact_frontier=.false.
  logical :: lshepherd_neighbour=.false.
  logical :: lrandom_particle_pencils=.false., lrandom_particle_blocks=.false.
  logical :: linsert_particles_continuously=.false.
  logical :: loutput_psize_dist=.false.
  logical :: lsinkparticle_1=.false.
  logical :: linsert_particle=.false.
  logical :: lcommunicate_rhop=.false.
  logical :: lcommunicate_np=.false.
!
  character (len=2*bclen+1) :: bcpx='p', bcpy='p', bcpz='p'
  character (len=2*bclen+1) :: bcspx='p', bcspy='p', bcspz='p'
  character (len=10), dimension(mpvar) :: pvarname
!
  type quant_interp_penc
!
!  Interpolation toggles:
!
    logical :: luu, loo, lTT, lrho, lgradTT
!
!  Interpolation policies:
!
    integer :: pol_uu, pol_oo, pol_TT, pol_rho, pol_gradTT
  end type quant_interp_penc
!  
  type(quant_interp_penc) :: interp
!
!  Interpolated quantities: moved outside type to conform to
!  the f90 standard.
!
    real, dimension(:,:), allocatable :: interp_uu, interp_oo, interp_gradTT
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
  real :: t_nextcol=0. !collision diagnostic times, set to turn-over of largest eddies
!
endmodule Particles_cdata
