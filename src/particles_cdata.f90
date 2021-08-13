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
  integer, parameter :: lun_output=93
!
  real, parameter :: npar_per_cell=npar/(1.0*nwgrid)
!
! Introduce the maximum number possible particles that can live in
! a cell.
!
  integer, parameter :: maxp=20
  integer(KIND=ikind8), parameter :: npar_maxp=npar*maxp
  integer, parameter :: max_par_per_grid=int(npar_maxp/nwgrid)+1  ! ceiling needed?
  real :: rp_int=-impossible, rp_ext=-impossible
  real :: rp_ext_width=-impossible
  real :: dsnap_par_minor=0.0, dsnap_par=0.0
  real :: rhopmat=1.0, rhopmat1=1.0, mpmat=0.0
  real :: eps_dtog = 0.01
  real :: mp_swarm=0.0, np_swarm=0.0, rhop_swarm=0.0
  real :: four_pi_rhopmat_over_three
  real :: np_const=0.0, rhop_const=0.0
  real :: particle_radius = 0.0
  real :: energy_gain_shear_bcs=impossible
  real :: log_ap_min_dist=0.0, log_ap_max_dist=6.0
  real :: rsinkparticle_1=0.0
  real :: t_nextinsert=0. !The time at which new particles are going to be inserted.
  real :: dustdensity_powerlaw=0.
  real, dimension(4) :: gab_weights
  real :: gab_width=3.0

  real :: remove_particle_at_time=-1.0, remove_particle_criteria_size=0.0
  real :: remove_particle_criteria_edtog=impossible
!
  integer, dimension(-1:1,-1:1,-1:1) :: neighbors_par = -1
  integer, dimension (nx) :: kshepherd
  integer, allocatable, dimension (:) :: kneighbour
  integer, dimension (mpar_loc) :: ipar
  integer, dimension (npar_species) :: ipar_fence_species=0
  integer, dimension(ny*nz) :: npar_imn, k1_imn, k2_imn
  integer :: npvar=0, npar_loc=0, npar_total=0, npaux=0
  integer :: ixp=0, iyp=0, izp=0, ivpx=0, ivpy=0, ivpz=0, iap=0, iaps=0, irpbeta=0
  integer :: ixp0=0, iyp0=0, izp0=0,ippersist=0
  integer :: idR11=0, idR12=0, idR13=0, idV11=0, idV12=0, idV13=0
  integer :: idR21=0, idR22=0, idR23=0, idV21=0, idV22=0, idV23=0
  integer :: idR31=0, idR32=0, idR33=0, idV31=0, idV32=0, idV33=0
  integer :: iVolp,iVelVolp
  integer :: idXpo1=0, idXpo2=0, idXpo3=0
  integer :: iuf=0,iufx=0,iufy=0,iufz=0
  integer :: isigmap11=0,isigmap12=0,isigmap13=0
  integer :: isigmap21=0,isigmap22=0,isigmap23=0
  integer :: isigmap31=0,isigmap32=0,isigmap33=0
  integer :: ilnVp=0
  integer :: iblowup=0,iPPp=0,iQQp=0,iRRp=0
  integer :: iTp=0, imp=0, iCOp=0, impinit=0, iapinit=0
  integer :: irhosurf=0
  integer :: ivpx_cart,ivpy_cart,ivpz_cart
  integer :: inpswarm=0, irhopswarm=0
  integer :: ipsx=0, ipsy=0, ipsz=0
  integer :: ipss=0, ipst=0, ipxx=0, ipyy=0, ipzz=0
  integer :: iuup=0, iupx=0, iupy=0, iupz=0
  integer :: ipviscx=0, ipviscy=0, ipviscz=0
  integer :: inp=0, irhop=0, irhops=0, ipeh=0
  integer :: iup11=0,iup12=0,iup13=0
  integer :: iup21=0,iup22=0,iup23=0
  integer :: iup31=0,iup32=0,iup33=0
  integer :: ibpx=0,ibpy=0,ibpz=0
  integer :: idiag_nmigmax=0, idiag_nmigmmax=0, npart_radii=0
  integer :: nbin_ap_dist=100
  integer :: iads=0, iads_end=0
  integer :: isurf=0,isurf_end=0
  integer :: ieffp=0
  integer :: idlncc=0
  integer :: idfg=0,idfx=0,idfy=0,idfz=0
  integer :: ivpxt=0, ivpzt=0

  integer :: npar_inserted_tot=0
  integer :: max_particles=npar
!
  logical :: linterpolate_spline=.true.
  logical :: lparticlemesh_cic=.true., lparticlemesh_tsc=.false.
  logical :: lparticlemesh_pqs_assignment=.false.
  logical :: lparticlemesh_gab = .false.
  logical :: linterp_reality_check=.false., lmigration_redo=.false.
  logical :: lnocalc_np=.false., lnocalc_rhop=.false.
  logical :: lcalc_np=.true., lcalc_rhop=.true.
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
  logical :: lparticles_radius_rpbeta=.false.
  logical :: lignore_rhop_swarm=.false.
  logical :: lnocollapse_xdir_onecell=.false.
  logical :: lnocollapse_ydir_onecell=.false.
  logical :: lnocollapse_zdir_onecell=.false.
  logical :: lswap_radius_and_number=.false.
!
  character (len=2*bclen+1) :: bcpx='p', bcpy='p', bcpz='p'
  character (len=labellen), dimension(mparray) :: pvarname
  character (len=labellen) :: particle_mesh = ''
  character (len=labellen) :: remove_particle_criteria='all'
!
  type quant_interp_penc
!
!  Interpolation toggles:
!
    logical :: luu, loo, lTT, lrho, lgradTT, lbb, lee
    logical :: lpp, lspecies, lnu
!
!  Interpolation policies:
!
    integer :: pol_uu, pol_oo, pol_TT, pol_rho, pol_gradTT, pol_bb, pol_ee
    integer :: pol_pp, pol_species, pol_nu
  end type quant_interp_penc
!
  type(quant_interp_penc) :: interp
!
!  Interpolated quantities: moved outside type to conform to
!  the f90 standard.
!
    real, dimension(:,:), allocatable :: interp_uu, interp_oo, interp_gradTT
    real, dimension(:,:), allocatable :: interp_species
    real, dimension(:), allocatable :: interp_pp, interp_TT, interp_rho
    real, dimension(:), allocatable :: interp_bb, interp_ee, interp_nu
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
! The following in necessary for particle-particle interaction. If it is used,
! the array fpwn in allocated in the module particles_potential and its
! dimensions are lpar_max and mparray respectively.
! The array fpwn contains fp of this processor with the neighbours
!  
!
! Notation: lpar_max is the maximum number of particles that
! can be accommodated in the array fpwn which includes the   
! the particles from this processors and its neighbours.
! lpar_loc is the number upto which the array fpwn is actually 
! filled, all the rest of the elements should be zero.
!
  integer :: lpar_max=0,lpar_loc=0,nslab=0
  real, allocatable, dimension(:,:) :: fpwn,fp_buffer_in,fp_buffer_out
!
! Number of grid points that are considered neighbouring grids. This is determined by
! the grid resolution and the number psigma which determines that range of interaction of
! the potential. The default is set to 0 but must be reset in post-parameter read initialization.
! We use three different numbers for three directions in case we have different resolutions
! along three directions. So far this is written only for cartesian co-ordinates.
  integer,allocatable,dimension(:,:,:,:) :: invert_ineargrid_map
  logical :: lallocated_neighbour_list=.false.
!***********************************************************************
!contains
!***********************************************************************
!    subroutine allocate_neighbour_list(neighbourx,neighboury,neighbourz,Nneighbour)
!      integer :: neighbourx,neighboury,neighbourz,Nneighbour
!
!  allocates the memory for calculation of neighbourlist
!
!      write(*,*) 'DM:neighbourx,neighboury,neighbourz',neighbourx,neighboury,neighbourz
!      allocate(nlist(1-neighbourx:nx+neighbourx,1-neighboury:nx+neighboury, &
!        1-neighbourz:nx+neighbourz,Nneighbour+1))
!      lallocated_neighbour_list=.true.
!
!    endsubroutine allocate_neighbour_list
!***********************************************************************
!
endmodule Particles_cdata
