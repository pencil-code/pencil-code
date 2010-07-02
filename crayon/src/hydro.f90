! $Id: hydro.f90 13608 2010-04-09 11:56:42Z mppiyali $
!
!  This module takes care of most of the things related to velocity.
!  Pressure, for example, is added in the energy (entropy) module,
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lhydro = .true.
!
! MVAR CONTRIBUTION 3
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED divu; oo(3); o2; ou; u2; uij(3,3); uu(3)
! PENCILS PROVIDED sij(3,3); sij2; uij5(3,3); ugu(3); ugu2; oij(3,3); qq(3)
! PENCILS PROVIDED del2u(3)
! PENCILS PROVIDED graddivu(3) 
! PENCILS PROVIDED rhougu(3)
!
!***************************************************************
module Hydro
!
  use Cdata
  use Cparam
  use Viscosity
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'hydro.h'
!
! Slice precalculation buffers
!
  real, target, dimension (nx,ny,3) :: oo_xy,oo_xy2,oo_xy3,oo_xy4
  real, target, dimension (nx,nz,3) :: oo_xz
  real, target, dimension (ny,nz,3) :: oo_yz
  real, target, dimension (nx,ny) :: divu_xy,u2_xy,o2_xy,mach_xy
  real, target, dimension (nx,ny) :: divu_xy2,u2_xy2,o2_xy2,mach_xy2
  real, target, dimension (nx,nz) :: divu_xz,u2_xz,o2_xz,mach_xz
  real, target, dimension (ny,nz) :: divu_yz,u2_yz,o2_yz,mach_yz
  real, dimension (nz,3) :: uumz,guumz=0.
  real, target, dimension (nx,ny) :: divu_xy3,divu_xy4,u2_xy3,u2_xy4,mach_xy4
  real, target, dimension (nx,ny) :: o2_xy3,o2_xy4,mach_xy3
!
!  cosine and sine function for setting test fields and analysis
!
  real, dimension (mz) :: c2z,csz,s2z,cz,sz
!
!  precession matrices
!
  real, dimension (3,3) :: mat_cori=0.,mat_cent=0.
!
! init parameters
!
  real :: widthuu=.1, radiusuu=1., urand=0., kx_uu=1., ky_uu=1., kz_uu=1.
  real :: urandi=0.
  real :: uu_left=0.,uu_right=0.,uu_lower=1.,uu_upper=1.
  real :: uy_left=0.,uy_right=0.
  real :: initpower=1.,cutoff=0.
  real, dimension (ninit) :: ampl_ux=0.0, ampl_uy=0.0, ampl_uz=0.0
  real, dimension (ninit) :: kx_ux=0.0, kx_uy=0.0, kx_uz=0.0
  real, dimension (ninit) :: ky_ux=0.0, ky_uy=0.0, ky_uz=0.0
  real, dimension (ninit) :: kz_ux=0.0, kz_uy=0.0, kz_uz=0.0
  real, dimension (ninit) :: phase_ux=0.0, phase_uy=0.0, phase_uz=0.0
  real :: omega_precession=0., alpha_precession=0.
  real, dimension (ninit) :: ampluu=0.0
  character (len=labellen), dimension(ninit) :: inituu='nothing'
  character (len=labellen) :: borderuu='nothing'
  real, dimension (3) :: uu_const=(/0.,0.,0./)
  complex, dimension (3) :: coefuu=(/0.,0.,0./)
  real :: u_out_kep=0.0, velocity_ceiling=-1.0
  real :: mu_omega=0., gap=0.
  integer :: nb_rings=0
  integer :: neddy=0
  real, dimension (5) :: om_rings=0.
  integer :: N_modes_uu=0
  logical :: ladvection_velocity=.true.
  logical :: lprecession=.false.
  logical :: lshear_rateofstrain=.false.
  logical :: luut_as_aux=.false.,loot_as_aux=.false.
  logical, target :: lpressuregradient_gas=.true.
  logical, target :: lcoriolis_force=.true.
  logical, target :: lcentrifugal_force=.false.
  logical :: lscale_tobox=.true.
  real :: incl_alpha = 0.0, rot_rr = 0.0
  real :: xsphere = 0.0, ysphere = 0.0, zsphere = 0.0
! The following is useful to debug the forcing - Dhruba
  real :: outest
  real :: ampl_Omega=0.0
  real :: omega_ini=0.0
  logical :: loutest,ldiffrot_test=.false.
  real :: r_cyl = 1.0, skin_depth = 1e-1
!
  namelist /hydro_init_pars/ &
       ampluu, ampl_ux, ampl_uy, ampl_uz, phase_ux, phase_uy, phase_uz, &
       inituu, widthuu, radiusuu, urand, urandi, lpressuregradient_gas, &
       uu_left, uu_right, uu_lower, uu_upper, kx_uu, ky_uu, kz_uu, coefuu, &
       kx_ux, ky_ux, kz_ux, kx_uy, ky_uy, kz_uy, kx_uz, ky_uz, kz_uz, &
       uy_left, uy_right,uu_const, Omega,  initpower, cutoff, &
       u_out_kep, N_modes_uu, lcoriolis_force, lcentrifugal_force, &
       ladvection_velocity, lprecession, omega_precession, alpha_precession, &
       luut_as_aux,loot_as_aux, &
       velocity_ceiling, mu_omega, nb_rings, om_rings, gap, &
       lscale_tobox, ampl_Omega,omega_ini, &
       r_cyl,skin_depth, incl_alpha, rot_rr,xsphere,ysphere,zsphere,&
       neddy
! run parameters
  real :: tdamp=0.,dampu=0.,wdamp=0.
  real :: dampuint=0.0,dampuext=0.0,rdampint=impossible,rdampext=impossible
  real :: ruxm=0.,ruym=0.,ruzm=0.
  real :: tau_damp_ruxm1=0.,tau_damp_ruym1=0.,tau_damp_ruzm1=0.
  real :: tau_damp_ruxm=0.,tau_damp_ruym=0.,tau_damp_ruzm=0.,tau_diffrot1=0.
  real :: ampl1_diffrot=0.,ampl2_diffrot=0.
  real :: Omega_int=0.,xexp_diffrot=1.,kx_diffrot=1.,kz_diffrot=0.
  real :: othresh=0.,othresh_per_orms=0.,orms=0.,othresh_scl=1.
  real :: utop=0.,ubot=0.,omega_out=0.,omega_in=0.
  real :: width_ff_uu=1.,x1_ff_uu=0.,x2_ff_uu=0.
  real :: eckmann_friction=0.0
  integer :: novec,novecmax=nx*ny*nz/4
  logical :: ldamp_fade=.false.,lOmega_int=.false.
  logical :: lfreeze_uint=.false.,lfreeze_uext=.false.
  logical :: lremove_mean_momenta=.false.
  logical :: lremove_mean_flow=.false.
  logical :: lreinitialize_uu=.false.
  logical :: lalways_use_gij_etc=.false.
  logical :: lcalc_uumean=.false.
  logical :: lcoriolis_xdep=.false.
  character (len=labellen) :: uuprof='nothing'
!
!  parameters for interior boundary conditions
!
  character (len=labellen) :: interior_bc_hydro_profile='nothing'
  logical :: lhydro_bc_interior=.false.
  logical :: lno_meridional_flow=.false.
  real :: z1_interior_bc_hydro=0.,kz_analysis=1.
!
  namelist /hydro_run_pars/ &
       Omega,theta, &
       tdamp,dampu,dampuext,dampuint,rdampext,rdampint,wdamp, &
       tau_damp_ruxm,tau_damp_ruym,tau_damp_ruzm,tau_diffrot1, &
       inituu,ampluu,kz_uu, &
       ampl1_diffrot,ampl2_diffrot,uuprof, &
       xexp_diffrot,kx_diffrot,kz_diffrot, &
       kz_analysis, &
       lreinitialize_uu,lremove_mean_momenta,lremove_mean_flow, &
       lOmega_int,Omega_int, ldamp_fade, othresh,othresh_per_orms, &
       borderuu, lfreeze_uint, lpressuregradient_gas, &
       lfreeze_uext,lcoriolis_force,lcentrifugal_force,ladvection_velocity, &
       utop,ubot,omega_out,omega_in, &
       lprecession, omega_precession, alpha_precession, lshear_rateofstrain, &
       lalways_use_gij_etc,lcalc_uumean, &
       width_ff_uu,x1_ff_uu,x2_ff_uu, &
       luut_as_aux,loot_as_aux, &
       loutest, ldiffrot_test,&
       interior_bc_hydro_profile, lhydro_bc_interior, z1_interior_bc_hydro, &
       velocity_ceiling,&
       eckmann_friction, ampl_Omega, lcoriolis_xdep,&
       lno_meridional_flow
! diagnostic variables (need to be consistent with reset list below)
  integer :: idiag_u2tm=0       ! DIAG_DOC: $\left<\uv(t)\cdot\int_0^t\uv(t')
                                ! DIAG_DOC:   dt'\right>$
  integer :: idiag_uotm=0       ! DIAG_DOC: $\left<\uv(t)\cdot\int_0^t\omv(t')
                                ! DIAG_DOC:   dt'\right>$
  integer :: idiag_outm=0       ! DIAG_DOC: $\left<\omv(t)\cdot\int_0^t\uv(t')
                                ! DIAG_DOC:   dt'\right>$
  integer :: idiag_u2m=0        ! DIAG_DOC: $\left<\uv^2\right>$
  integer :: idiag_um2=0        ! DIAG_DOC:
  integer :: idiag_uxpt=0       ! DIAG_DOC:
  integer :: idiag_uypt=0       ! DIAG_DOC:
  integer :: idiag_uzpt=0       ! DIAG_DOC:
  integer :: idiag_urms=0       ! DIAG_DOC: $\left<\uv^2\right>^{1/2}$
  integer :: idiag_umax=0       ! DIAG_DOC: $\max(|\uv|)$
  integer :: idiag_uzrms=0      ! DIAG_DOC: $\left<u_z^2\right>^{1/2}$
  integer :: idiag_uzrmaxs=0    ! DIAG_DOC:
  integer :: idiag_uxmin=0      ! DIAG_DOC: $\min(|u_x|)$
  integer :: idiag_uymin=0      ! DIAG_DOC: $\min(|u_y|)$
  integer :: idiag_uzmin=0      ! DIAG_DOC: $\min(|u_z|)$
  integer :: idiag_uxmax=0      ! DIAG_DOC: $\max(|u_x|)$
  integer :: idiag_uymax=0      ! DIAG_DOC: $\max(|u_y|)$
  integer :: idiag_uzmax=0      ! DIAG_DOC: $\max(|u_z|)$
  integer :: idiag_uxm=0        ! DIAG_DOC: $\left<u_x\right>$
  integer :: idiag_uym=0        ! DIAG_DOC: $\left<u_y\right>$
  integer :: idiag_uzm=0        ! DIAG_DOC: $\left<u_z\right>$
  integer :: idiag_ux2m=0       ! DIAG_DOC: $\left<u_x^2\right>$
  integer :: idiag_uy2m=0       ! DIAG_DOC: $\left<u_y^2\right>$
  integer :: idiag_uz2m=0       ! DIAG_DOC: $\left<u_z^2\right>$
  integer :: idiag_ux2ccm=0     ! DIAG_DOC: $\left<u_x^2\cos^2kz\right>$
  integer :: idiag_ux2ssm=0     ! DIAG_DOC: $\left<u_x^2\sin^2kz\right>$
  integer :: idiag_uy2ccm=0     ! DIAG_DOC: $\left<u_y^2\cos^2kz\right>$
  integer :: idiag_uy2ssm=0     ! DIAG_DOC: $\left<u_y^2\sin^2kz\right>$
  integer :: idiag_ux2mx=0      ! DIAG_DOC: $\left<u_x^2\right>_{yz}$
  integer :: idiag_uy2mx=0      ! DIAG_DOC: $\left<u_y^2\right>_{yz}$
  integer :: idiag_uz2mx=0      ! DIAG_DOC: $\left<u_z^2\right>_{yz}$
  integer :: idiag_ox2mx=0      ! DIAG_DOC: $\left<\omega_x^2\right>_{yz}$
  integer :: idiag_oy2mx=0      ! DIAG_DOC: $\left<\omega_y^2\right>_{yz}$
  integer :: idiag_oz2mx=0      ! DIAG_DOC: $\left<\omega_z^2\right>_{yz}$
  integer :: idiag_uxuycsm=0    ! DIAG_DOC: $\left<u_xu_y\cos kz\sin kz\right>$
  integer :: idiag_ux2my=0      ! DIAG_DOC:
  integer :: idiag_uy2my=0      ! DIAG_DOC:
  integer :: idiag_uz2my=0      ! DIAG_DOC:
  integer :: idiag_ux2mz=0      ! DIAG_DOC: $\left<u_x^2\right>_{xy}$
  integer :: idiag_uy2mz=0      ! DIAG_DOC: $\left<u_y^2\right>_{xy}$
  integer :: idiag_uz2mz=0      ! DIAG_DOC: $\left<u_z^2\right>_{xy}$
  integer :: idiag_rux2mz=0     ! DIAG_DOC: $\left<\varrho u_x^2\right>_{xy}$
  integer :: idiag_ruy2mz=0     ! DIAG_DOC: $\left<\varrho u_y^2\right>_{xy}$
  integer :: idiag_ruz2mz=0     ! DIAG_DOC: $\left<\varrho u_z^2\right>_{xy}$
  integer :: idiag_uxuym=0      ! DIAG_DOC: $\left<u_x u_y\right>$
  integer :: idiag_uxuzm=0      ! DIAG_DOC: $\left<u_x u_z\right>$
  integer :: idiag_uyuzm=0      ! DIAG_DOC: $\left<u_y u_z\right>$
  integer :: idiag_uxuymz=0     ! DIAG_DOC:
  integer :: idiag_uxuzmz=0     ! DIAG_DOC:
  integer :: idiag_uyuzmz=0     ! DIAG_DOC:
  integer :: idiag_uxuymy=0     ! DIAG_DOC:
  integer :: idiag_uxuzmy=0     ! DIAG_DOC:
  integer :: idiag_uyuzmy=0     ! DIAG_DOC:
  integer :: idiag_uxuymx=0     ! DIAG_DOC:
  integer :: idiag_uxuzmx=0     ! DIAG_DOC:
  integer :: idiag_uyuzmx=0     ! DIAG_DOC:
  integer :: idiag_uxmz=0       ! DIAG_DOC: $\left< u_x \right>_{x,y}$
                                ! DIAG_DOC:   \quad(horiz. averaged $x$
                                ! DIAG_DOC:   velocity)
  integer :: idiag_uymz=0       ! DIAG_DOC:
  integer :: idiag_uzmz=0       ! DIAG_DOC:
  integer :: idiag_oxmz=0       ! DIAG_DOC: $\left< \omega_x \right>_{xy}$
  integer :: idiag_oymz=0       ! DIAG_DOC: $\left< \omega_y \right>_{xy}$
  integer :: idiag_ozmz=0       ! DIAG_DOC: $\left< \omega_z \right>_{xy}$
  integer :: idiag_umx=0        ! DIAG_DOC: $\left< u_x \right>$
  integer :: idiag_umy=0        ! DIAG_DOC: $\left< u_y \right>$
  integer :: idiag_umz=0        ! DIAG_DOC: $\left< u_z \right>$
  integer :: idiag_uxmy=0       ! DIAG_DOC: $\left< u_x \right>_{y}$
  integer :: idiag_uymy=0       ! DIAG_DOC: $\left< u_y \right>_{y}$
  integer :: idiag_uzmy=0       ! DIAG_DOC: $\left< u_z \right>_{y}$
  integer :: idiag_u2mz=0       ! DIAG_DOC: $\left< \uv^2 \right>_{y}$
  integer :: idiag_omumz=0      ! DIAG_DOC: $\left<\left<\Wv\right>_{xy}
                                ! DIAG_DOC:   \cdot\left<\Uv\right>_{xy}
                                ! DIAG_DOC:   \right>$ \quad($xy$-averaged
                                ! DIAG_DOC:   mean cross helicity production)
  integer :: idiag_umamz=0      ! DIAG_DOC: $\left<\left<\uv\right>_{xy}\cdot\left<\Av\right>_{xy}\right>$ 
  integer :: idiag_umbmz=0      ! DIAG_DOC: $\left<\left<\Uv\right>_{xy}
                                ! DIAG_DOC:   \cdot\left<\Bv\right>_{xy}
                                ! DIAG_DOC:   \right>$ \quad($xy$-averaged
                                ! DIAG_DOC:   mean cross helicity production)
  integer :: idiag_umxbmz=0     ! DIAG_DOC: $\left<\left<\Uv\right>_{xy}
                                ! DIAG_DOC:   \times\left<\Bv\right>_{xy}
                                ! DIAG_DOC:   \right>_z$ \quad($xy$-averaged
                                ! DIAG_DOC:   mean emf)
  integer :: idiag_uxmxy=0      ! DIAG_DOC: $\left< u_x \right>_{z}$
  integer :: idiag_uymxy=0      ! DIAG_DOC: $\left< u_y \right>_{z}$
  integer :: idiag_uzmxy=0      ! DIAG_DOC: $\left< u_z \right>_{z}$
  integer :: idiag_oxmxy=0      ! DIAG_DOC: $\left< \omega_x \right>_{z}$
  integer :: idiag_oymxy=0      ! DIAG_DOC: $\left< \omega_y \right>_{z}$
  integer :: idiag_ozmxy=0      ! DIAG_DOC: $\left< \omega_z \right>_{z}$
  integer :: idiag_pvzmxy=0     ! DIAG_DOC: $\left< (\omega_z + 2\Omega)/\varrho \right>_{z}$
                                ! DIAG_DOC: \quad(z component of potential vorticity)
  integer :: idiag_ruxmxy=0     ! DIAG_DOC: $\left< \rho u_x \right>_{z}$
  integer :: idiag_ruymxy=0     ! DIAG_DOC: $\left< \rho u_y \right>_{z}$
  integer :: idiag_ruzmxy=0     ! DIAG_DOC: $\left< \rho u_z \right>_{z}$
  integer :: idiag_ux2mxy=0     ! DIAG_DOC: $\left< u_x^2 \right>_{z}$
  integer :: idiag_uy2mxy=0     ! DIAG_DOC: $\left< u_y^2 \right>_{z}$
  integer :: idiag_uz2mxy=0     ! DIAG_DOC: $\left< u_z^2 \right>_{z}$
  integer :: idiag_rux2mxy=0    ! DIAG_DOC: $\left< \rho u_x^2 \right>_{z}$
  integer :: idiag_ruy2mxy=0    ! DIAG_DOC: $\left< \rho u_y^2 \right>_{z}$
  integer :: idiag_ruz2mxy=0    ! DIAG_DOC: $\left< \rho u_z^2 \right>_{z}$
  integer :: idiag_ruxuymxy=0   ! DIAG_DOC: $\left< \rho u_x u_y \right>_{z}$
  integer :: idiag_ruxuzmxy=0   ! DIAG_DOC: $\left< \rho u_x u_z \right>_{z}$
  integer :: idiag_ruyuzmxy=0   ! DIAG_DOC: $\left< \rho u_y u_z \right>_{z}$
  integer :: idiag_rux2m=0      ! DIAG_DOC: $\left<\rho u_x^2\right>$
  integer :: idiag_ruy2m=0      ! DIAG_DOC: $\left<\rho u_y^2\right>$
  integer :: idiag_ruz2m=0      ! DIAG_DOC: $\left<\rho u_z^2\right>$
  integer :: idiag_uxmxz=0      ! DIAG_DOC: $\left< u_x \right>_{y}$
  integer :: idiag_uymxz=0      ! DIAG_DOC: $\left< u_y \right>_{y}$
  integer :: idiag_uzmxz=0      ! DIAG_DOC: $\left< u_z \right>_{y}$
  integer :: idiag_ux2mxz=0     ! DIAG_DOC: $\left< u_x^2 \right>_{y}$
  integer :: idiag_uy2mxz=0     ! DIAG_DOC: $\left< u_y^2 \right>_{y}$
  integer :: idiag_uz2mxz=0     ! DIAG_DOC: $\left< u_z^2 \right>_{y}$
  integer :: idiag_uxmx=0       ! DIAG_DOC: $\left< u_x \right>_{yz}$
  integer :: idiag_uymx=0       ! DIAG_DOC: $\left< u_y \right>_{yz}$
  integer :: idiag_uzmx=0       ! DIAG_DOC: $\left< u_z \right>_{yz}$
  integer :: idiag_divum=0      ! DIAG_DOC:$\left<{\rm div}\uv)\right>$
  integer :: idiag_divu2m=0     ! DIAG_DOC: $\left<({\rm div}\uv)^2\right>$
  integer :: idiag_u3u21m=0     ! DIAG_DOC: $\left<u_3 u_{2,1}\right>$
  integer :: idiag_u1u32m=0     ! DIAG_DOC: $\left<u_1 u_{3,2}\right>$
  integer :: idiag_u2u13m=0     ! DIAG_DOC: $\left<u_2 u_{1,3}\right>$
  integer :: idiag_u2u31m=0     ! DIAG_DOC: $\left<u_2 u_{3,1}\right>$
  integer :: idiag_u3u12m=0     ! DIAG_DOC: $\left<u_3 u_{1,2}\right>$
  integer :: idiag_u1u23m=0     ! DIAG_DOC: $\left<u_1 u_{2,3}\right>$
  integer :: idiag_u3u21mz=0    ! DIAG_DOC:
  integer :: idiag_u1u32mz=0    ! DIAG_DOC:
  integer :: idiag_u2u13mz=0    ! DIAG_DOC:
  integer :: idiag_u2u31mz=0    ! DIAG_DOC:
  integer :: idiag_u3u12mz=0    ! DIAG_DOC:
  integer :: idiag_u1u23mz=0    ! DIAG_DOC:
  integer :: idiag_u2mr=0       ! DIAG_DOC:
  integer :: idiag_urmr=0       ! DIAG_DOC:
  integer :: idiag_upmr=0       ! DIAG_DOC:
  integer :: idiag_uzmr=0       ! DIAG_DOC:
  integer :: idiag_uxfampm=0    ! DIAG_DOC:
  integer :: idiag_uyfampm=0    ! DIAG_DOC:
  integer :: idiag_uzfampm=0    ! DIAG_DOC:
  integer :: idiag_uxfampim=0   ! DIAG_DOC:
  integer :: idiag_uyfampim=0   ! DIAG_DOC:
  integer :: idiag_uzfampim=0   ! DIAG_DOC:
  integer :: idiag_ruxm=0       ! DIAG_DOC: $\left<\varrho u_x\right>$
                                ! DIAG_DOC:   \quad(mean $x$-momentum density)
  integer :: idiag_ruym=0       ! DIAG_DOC: $\left<\varrho u_y\right>$
                                ! DIAG_DOC:   \quad(mean $y$-momentum density)
  integer :: idiag_ruzm=0       ! DIAG_DOC: $\left<\varrho u_z\right>$
                                ! DIAG_DOC:   \quad(mean $z$-momentum density)
  integer :: idiag_ruxtot=0     ! DIAG_DOC: $\left<\rho |u|\right>$
                                ! DIAG_DOC:   \quad(mean absolute $x$-momentum density)
  integer :: idiag_rumax=0      ! DIAG_DOC: $\max(\varrho |\uv|)$
                                ! DIAG_DOC:   \quad(maximum modulus of momentum)
  integer :: idiag_ruxuym=0     ! DIAG_DOC: $\left<\varrho u_x u_y\right>$
                                ! DIAG_DOC:   \quad(mean Reynolds stress)
  integer :: idiag_ruxuzm=0     ! DIAG_DOC: $\left<\varrho u_x u_z\right>$
                                ! DIAG_DOC:   \quad(mean Reynolds stress)
  integer :: idiag_ruyuzm=0     ! DIAG_DOC: $\left<\varrho u_y u_z\right>$
                                ! DIAG_DOC:   \quad(mean Reynolds stress)
  integer :: idiag_ruxuymz=0    ! DIAG_DOC:
  integer :: idiag_rlxm=0       ! DIAG_DOC: $\left< \rho y u_z - z u_y \right>$
  integer :: idiag_rlym=0       ! DIAG_DOC: $\left< \rho z u_x - x u_z \right>$
  integer :: idiag_rlzm=0       ! DIAG_DOC: $\left< \rho x u_y - y u_x \right>$
  integer :: idiag_rlx2m=0      ! DIAG_DOC: $\left<(\rho y u_z-z u_y)^2\right>$
  integer :: idiag_rly2m=0      ! DIAG_DOC: $\left<(\rho z u_x-x u_z)^2\right>$
  integer :: idiag_rlz2m=0      ! DIAG_DOC: $\left<(\rho x u_y-y u_x)^2\right>$
  integer :: idiag_tot_ang_mom=0! DIAG_DOC: Total angular momentum in spherical
                                ! DIAG_DOC: coordinates about the axis.
  integer :: idiag_rufm=0       ! DIAG_DOC:
  integer :: idiag_dtu=0        ! DIAG_DOC: $\delta t/[c_{\delta t}\,\delta x
                                ! DIAG_DOC:  /\max|\mathbf{u}|]$
                                ! DIAG_DOC:  \quad(time step relative to
                                ! DIAG_DOC:   advective time step;
                                ! DIAG_DOC:   see \S~\ref{time-step})
  integer :: idiag_oum=0        ! DIAG_DOC: $\left<\boldsymbol{\omega}
                                ! DIAG_DOC:   \cdot\uv\right>$
  integer :: idiag_fum=0        ! DIAG_DOC: $\left<\fv\cdot\uv\right>$
  integer :: idiag_o2m=0        ! DIAG_DOC: $\left<\boldsymbol{\omega}^2\right>
                                ! DIAG_DOC:   \equiv \left<(\curl\uv)^2\right>$
  integer :: idiag_orms=0       ! DIAG_DOC: $\left<\boldsymbol{\omega}^2
                                ! DIAG_DOC:   \right>^{1/2}$
  integer :: idiag_omax=0       ! DIAG_DOC: $\max(|\boldsymbol{\omega}|)$
  integer :: idiag_ox2m=0       ! DIAG_DOC: $\left<\omega_x^2\right>$
  integer :: idiag_oy2m=0       ! DIAG_DOC: $\left<\omega_y^2\right>$
  integer :: idiag_oz2m=0       ! DIAG_DOC: $\left<\omega_z^2\right>$
  integer :: idiag_oxm=0        ! DIAG_DOC:
  integer :: idiag_oym=0        ! DIAG_DOC:
  integer :: idiag_ozm=0        ! DIAG_DOC:
  integer :: idiag_oxoym=0      ! DIAG_DOC: $\left<\omega_x\omega_y\right>$
  integer :: idiag_oxozm=0      ! DIAG_DOC: $\left<\omega_x\omega_z\right>$
  integer :: idiag_oyozm=0      ! DIAG_DOC: $\left<\omega_y\omega_z\right>$
  integer :: idiag_pvzm=0       ! DIAG_DOC: $\left<\omega_z + 2\Omega/\varrho\right>$
                                ! DIAG_DOC: \quad(z component of potential vorticity)
  integer :: idiag_oumx=0       ! DIAG_DOC: $\left<\boldsymbol{\omega}
                                ! DIAG_DOC: \cdot\uv\right>_{yz}$
  integer :: idiag_oumy=0       ! DIAG_DOC: $\left<\boldsymbol{\omega}
                                ! DIAG_DOC: \cdot\uv\right>_{xz}$
  integer :: idiag_oumz=0       ! DIAG_DOC: $\left<\boldsymbol{\omega}
                                ! DIAG_DOC: \cdot\uv\right>_{xy}$
  integer :: idiag_oumxy=0      ! DIAG_DOC: $\left<\boldsymbol{\omega}
                                ! DIAG_DOC: \cdot\uv\right>_{z}$
  integer :: idiag_oumxz=0      ! DIAG_DOC: $\left<\boldsymbol{\omega}
                                ! DIAG_DOC: \cdot\uv\right>_{y}$
  integer :: idiag_uguxm=0      ! DIAG_DOC:
  integer :: idiag_uguym=0      ! DIAG_DOC:
  integer :: idiag_uguzm=0      ! DIAG_DOC:
  integer :: idiag_ugu2m=0      ! DIAG_DOC:
  integer :: idiag_uguxmx=0     ! DIAG_DOC:
  integer :: idiag_uguymx=0     ! DIAG_DOC:
  integer :: idiag_uguzmx=0     ! DIAG_DOC:
  integer :: idiag_uguxmy=0     ! DIAG_DOC:
  integer :: idiag_uguymy=0     ! DIAG_DOC:
  integer :: idiag_uguzmy=0     ! DIAG_DOC:
  integer :: idiag_uguxmz=0     ! DIAG_DOC:
  integer :: idiag_uguymz=0     ! DIAG_DOC:
  integer :: idiag_uguzmz=0     ! DIAG_DOC:
  integer :: idiag_Marms=0      ! DIAG_DOC: $\left<\uv^2/\cs^2\right>$
                                ! DIAG_DOC:   \quad(rms Mach number)
  integer :: idiag_Mamax=0      ! DIAG_DOC: $\max |\uv|/\cs$
                                ! DIAG_DOC:   \quad(maximum Mach number)
  integer :: idiag_fintm=0      ! DIAG_DOC:
  integer :: idiag_fextm=0      ! DIAG_DOC:
  integer :: idiag_duxdzma=0    ! DIAG_DOC:
  integer :: idiag_duydzma=0    ! DIAG_DOC:
  integer :: idiag_ekin=0       ! DIAG_DOC: $\left<{1\over2}\varrho\uv^2\right>$
  integer :: idiag_ekintot=0    ! DIAG_DOC: $\int_V{1\over2}\varrho\uv^2\, dV$
  integer :: idiag_ekinz=0      ! DIAG_DOC:
  integer :: idiag_fmassz=0     ! DIAG_DOC:
  integer :: idiag_fkinz=0      ! DIAG_DOC: $\left<{1\over2}\varrho\uv^2 u_z\right>_{xy}$
  integer :: idiag_fkinxy=0     ! DIAG_DOC: $\left<{1\over2}\varrho\uv^2 u_x\right>_{z}$
  integer :: idiag_fxbxm=0      ! DIAG_DOC:
  integer :: idiag_fxbym=0      ! DIAG_DOC:
  integer :: idiag_fxbzm=0      ! DIAG_DOC:
  integer :: idiag_uxglnrym=0   ! DIAG_DOC: $\left<u_x\partial_y\ln\varrho\right>$
  integer :: idiag_uyglnrxm=0   ! DIAG_DOC: $\left<u_y\partial_x\ln\varrho\right>$
  integer :: idiag_uxuydivum=0  ! DIAG_DOC: $\left<u_x u_y\nabla\cdot\uv\right>$
  integer :: idiag_urmsn=0,idiag_urmss=0,idiag_urmsh=0
  integer :: idiag_ormsn=0,idiag_ormss=0,idiag_ormsh=0
  integer :: idiag_oumn=0,idiag_oums=0,idiag_oumh=0
!
  contains
!***********************************************************************
    subroutine register_hydro()
!
!  Initialise variables which should know that we solve the hydro
!  equations: iuu, etc; increase nvar accordingly.
!
!  6-nov-01/wolf: coded
!
      use FArrayManager
      use SharedVariables, only: put_shared_variable
      use Sub
!
      integer :: ierr
!
!  indices to access uu
!
      call farray_register_pde('uu',iuu,vector=3)
      iux = iuu; iuy = iuu+1; iuz = iuu+2
!
!  Share lpressuregradient_gas so Entropy module knows whether to apply
!  pressure gradient or not. But hydro wants pressure gradient only when
!  then density is computed, i.e. not even with ldensity_anelastic.
!
      if  (.not.ldensity) lpressuregradient_gas=.false.
      call put_shared_variable('lpressuregradient_gas',&
          lpressuregradient_gas,ierr)
      if (ierr/=0) call fatal_error('register_hydro',&
          'there was a problem sharing lpressuregradient_gas')
!
!  Identify version number (generated automatically by SVN).
!
      if (lroot) call svn_id( &
           "$Id: hydro.f90 13608 2010-04-09 11:56:42Z mppiyali $")
!
!  Writing files for use with IDL
!
      if (lroot) then
        if (maux == 0) then
          if (nvar < mvar) write(4,*) ',uu $'
          if (nvar == mvar) write(4,*) ',uu'
        else
          write(4,*) ',uu $'
        endif
        write(15,*) 'uu = fltarr(mx,my,mz,3)*one'
      endif
!
    endsubroutine register_hydro
!***********************************************************************
    subroutine initialize_hydro(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  24-nov-02/tony: coded
!  13-oct-03/dave: check parameters and warn (if nec.) about velocity damping
!  26-mar-10/axel: lreinitialize_uu added
!
      use FArrayManager
      use Mpicomm, only: stop_it
      use Initcond
      use SharedVariables, only: put_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mz) :: c, s
      logical :: lstarting
!
!  calculate cosz*sinz, cos^2, and sinz^2, to take moments with
!  of ux2, uxuy, etc.
!
      c=cos(kz_analysis*z)
      s=sin(kz_analysis*z)
      cz=c
      sz=s
      c2z=c**2
      s2z=s**2
      csz=c*s
!
!  Turn off advection for 0-D runs.
!
      if (nxgrid*nygrid*nzgrid==1) then
        ladvection_velocity=.false.
        print*, 'initialize_hydro: 0-D run, turned off advection of velocity'
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
      endsubroutine initialize_hydro
!***********************************************************************
    subroutine init_uu(f)
!
!  initialise velocity field ; called from start.f90
!
!  07-nov-01/wolf: coded
!  24-nov-02/tony: renamed for consistance (i.e. init_[variable name])
!
      use EquationOfState, only: cs20, gamma, beta_glnrho_scaled
      use General
      use Gravity, only: gravz_const,z1
      use Initcond
      use InitialCondition, only: initial_condition_uu
      use Mpicomm, only: stop_it
      use Sub
      use Density, only: calc_pencils_density
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      integer :: j,i
!
!  inituu corresponds to different initializations of uu (called from start).
!
      do j=1,ninit

        select case (inituu(j))

        case ('nothing'); if (lroot .and. j==1) print*,'init_uu: nothing'
        case ('zero', '0')
          if (lroot) print*,'init_uu: zero velocity'
          ! Ensure really is zero, as may have used lread_oldsnap
          f(:,:,:,iux:iuz)=0.
        case ('const_uu'); do i=1,3; f(:,:,:,iuu+i-1) = uu_const(i); enddo
        case ('mode'); call modev(ampluu(j),coefuu,f,iuu,kx_uu,ky_uu,kz_uu)
        case ('gaussian-noise'); call gaunoise(ampluu(j),f,iux,iuz)
        case ('gaussian-noise-x'); call gaunoise(ampluu(j),f,iux)
        case ('gaussian-noise-y'); call gaunoise(ampluu(j),f,iuy)
        case ('gaussian-noise-z'); call gaunoise(ampluu(j),f,iuz)
        case ('gaussian-noise-xy'); call gaunoise(ampluu(j),f,iux,iuy)
        case ('xjump')
          call jump(f,iux,uu_left,uu_right,widthuu,'x')
          call jump(f,iuy,uy_left,uy_right,widthuu,'x')
        case ('Beltrami-x'); call beltrami(ampluu(j),f,iuu,kx=kx_uu)
        case ('Beltrami-y'); call beltrami(ampluu(j),f,iuu,ky=ky_uu)
        case ('Beltrami-z'); call beltrami(ampluu(j),f,iuu,kz=kz_uu)
        case ('rolls'); call rolls(ampluu(j),f,iuu,kx_uu,kz_uu)
        case ('trilinear-x'); call trilinear(ampluu(j),f,iux)
        case ('trilinear-y'); call trilinear(ampluu(j),f,iuy)
        case ('trilinear-z'); call trilinear(ampluu(j),f,iuz)
        case ('cos-cos-sin-uz'); call cos_cos_sin(ampluu(j),f,iuz)
        case ('tor_pert'); call tor_pert(ampluu(j),f,iux)
        case ('diffrot'); call diffrot(ampluu(j),f,iuy)
        case ('olddiffrot'); call olddiffrot(ampluu(j),f,iuy)
        case ('sinwave-phase')
          call sinwave_phase(f,iux,ampl_ux(j),kx_ux(j),ky_ux(j),kz_ux(j),phase_ux(j))
          call sinwave_phase(f,iuy,ampl_uy(j),kx_uy(j),ky_uy(j),kz_uy(j),phase_uy(j))
          call sinwave_phase(f,iuz,ampl_uz(j),kx_uz(j),ky_uz(j),kz_uz(j),phase_uz(j))
        case ('coswave-phase')
          call coswave_phase(f,iux,ampl_ux(j),kx_ux(j),ky_ux(j),kz_ux(j),phase_ux(j))
          call coswave_phase(f,iuy,ampl_uy(j),kx_uy(j),ky_uy(j),kz_uy(j),phase_uy(j))
          call coswave_phase(f,iuz,ampl_uz(j),kx_uz(j),ky_uz(j),kz_uz(j),phase_uz(j))
        case ('sinwave-x'); call sinwave(ampluu(j),f,iux,kx=kx_uu)
        case ('sinwave-y'); call sinwave(ampluu(j),f,iuy,ky=ky_uu)
        case ('sinwave-z'); call sinwave(ampluu(j),f,iuz,kz=kz_uu)
        case ('sinwave-ux-kx'); call sinwave(ampluu(j),f,iux,kx=kx_uu)
        case ('sinwave-ux-ky'); call sinwave(ampluu(j),f,iux,ky=ky_uu)
        case ('sinwave-ux-kz'); call sinwave(ampluu(j),f,iux,kz=kz_uu)
        case ('sinwave-uy-kx'); call sinwave(ampluu(j),f,iuy,kx=kx_uu)
        case ('sinwave-uy-ky'); call sinwave(ampluu(j),f,iuy,ky=ky_uu)
        case ('sinwave-uy-kz'); call sinwave(ampluu(j),f,iuy,kz=kz_uu)
        case ('sinwave-uz-kx'); call sinwave(ampluu(j),f,iuz,kx=kx_uu)
        case ('sinwave-uz-ky'); call sinwave(ampluu(j),f,iuz,ky=ky_uu)
        case ('sinwave-uz-kz'); call sinwave(ampluu(j),f,iuz,kz=kz_uu)
        case ('sinwave-y-z')
          if (lroot) print*, 'init_uu: sinwave-y-z, ampluu=', ampluu(j)
          call sinwave(ampluu(j),f,iuy,kz=kz_uu)
        case ('sinwave-z-y')
          if (lroot) print*, 'init_uu: sinwave-z-y, ampluu=', ampluu(j)
          call sinwave(ampluu(j),f,iuz,ky=ky_uu)
        case ('sinwave-z-x')
          if (lroot) print*, 'init_uu: sinwave-z-x, ampluu=', ampluu(j)
          call sinwave(ampluu(j),f,iuz,kx=kx_uu)
        case ('damped_sinwave-z-x')
          if (lroot) print*, 'init_uu: damped_sinwave-z-x, ampluu=', ampluu(j)
          do m=m1,m2; do n=n1,n2
            f(:,m,n,iuz)=f(:,m,n,iuz)+ampluu(j)*sin(kx_uu*x)*exp(-10*z(n)**2)
          enddo; enddo
        case ('coswave-x'); call coswave(ampluu(j),f,iux,kx=kx_uu)
        case ('coswave-y'); call coswave(ampluu(j),f,iuy,ky=ky_uu)
        case ('coswave-z'); call coswave(ampluu(j),f,iuz,kz=kz_uu)
        case ('coswave-x-z'); call coswave(ampluu(j),f,iux,kz=kz_uu)
        case ('coswave-z-x'); call coswave(ampluu(j),f,iuz,kx=kx_uu)
        case ('x1cosycosz'); call x1_cosy_cosz(ampluu(j),f,iuy,ky=ky_uu,kz=kz_uu)
        case ('soundwave-x'); call soundwave(ampluu(j),f,iux,kx=kx_uu)
        case ('soundwave-y'); call soundwave(ampluu(j),f,iuy,ky=ky_uu)
        case ('soundwave-z'); call soundwave(ampluu(j),f,iuz,kz=kz_uu)
        case ('robertsflow'); call robertsflow(ampluu(j),f,iuu)
        case ('hawley-et-al'); call hawley_etal99a(ampluu(j),f,iuy,Lxyz)
        case ('const-ux'); f(:,:,:,iux) = ampluu(j)
        case ('const-uy'); f(:,:,:,iuy) = ampluu(j)
        case ('const-uz'); f(:,:,:,iuz) = ampluu(j)
        case default
          !
          !  Catch unknown values
          !
          if (lroot) print*, 'init_uu: No such value for inituu: ', &
            trim(inituu(j))
          call stop_it("")

        endselect
!
!  End loop over initial conditions
!
      enddo
!
!  Interface for user's own initial condition
!
      if (linitial_condition) call initial_condition_uu(f)
!
    endsubroutine init_uu
!***********************************************************************
    subroutine pencil_criteria_hydro()
!
!  All pencils that the Hydro module depends on are specified here.
!
!  20-nov-04/anders: coded
!
      use Mpicomm, only: stop_it
!
      if (ladvection_velocity) lpenc_requested(i_ugu)=.true.
      if (Omega/=0.0) lpenc_requested(i_uu)=.true.
!
!  video pencils
!
      if (dvid/=0.) then
        lpenc_video(i_oo)=.true.
        lpenc_video(i_o2)=.true.
        lpenc_video(i_u2)=.true.
      endif
!
!  diagnostic pencils
!
      lpenc_diagnos(i_uu)=.true.
      if (idiag_oumxy/=0 .or. &
          idiag_oumxz/=0) lpenc_diagnos2d(i_ou)=.true.
      if (idiag_ox2m/=0 .or. idiag_oy2m/=0 .or. idiag_oz2m/=0 .or. &
          idiag_oxm /=0 .or. idiag_oym /=0 .or. idiag_ozm /=0 .or. &
          idiag_oxoym/=0 .or. idiag_oxozm/=0 .or. idiag_oyozm/=0 .or. &
          idiag_pvzm /=0) &
          lpenc_diagnos(i_oo)=.true.
      if (idiag_orms/=0 .or. idiag_omax/=0 .or. idiag_o2m/=0 .or. &
          idiag_ormsh/=0 )  lpenc_diagnos(i_o2)=.true.
      if (idiag_oum/=0 .or. idiag_oumx/=0.or.idiag_oumy/=0.or.idiag_oumz/=0 .or. &
           idiag_oumh/=0) lpenc_diagnos(i_ou)=.true.
      if (idiag_Marms/=0 .or. idiag_Mamax/=0) lpenc_diagnos(i_Ma2)=.true.
      if (idiag_urms/=0 .or. idiag_umax/=0 .or. idiag_rumax/=0 .or. &
          idiag_u2m/=0 .or. idiag_um2/=0 .or. idiag_u2mz/=0 .or. &
          idiag_urmsh/=0) lpenc_diagnos(i_u2)=.true.
      if (idiag_duxdzma/=0 .or. idiag_duydzma/=0) lpenc_diagnos(i_uij)=.true.
      if (idiag_fmassz/=0 .or. idiag_ruxuym/=0 .or. idiag_ruxuymz/=0 .or. &
          idiag_ruxm/=0 .or. idiag_ruym/=0 .or. idiag_ruzm/=0 .or. &
          idiag_ruxuzm/=0 .or. idiag_ruyuzm/=0 .or. idiag_pvzm/=0 .or. &
          idiag_ruxtot/=0) lpenc_diagnos(i_rho)=.true.
      if (idiag_oxmxy/=0 .or. idiag_oymxy/=0 .or. idiag_ozmxy/=0 .or. &
          idiag_oxmz/=0 .or. idiag_oymz/=0 .or. idiag_ozmz/=0 .or. &
          idiag_pvzmxy/=0) &
          lpenc_diagnos2d(i_oo)=.true.
      if (idiag_pvzmxy/=0) lpenc_diagnos2d(i_rho)=.true.
       if (idiag_ekin/=0 .or. idiag_ekintot/=0 .or. idiag_fkinz/=0 .or. &
          idiag_ekinz/=0) then
        lpenc_diagnos(i_rho)=.true.
        lpenc_diagnos(i_u2)=.true.
      endif
      if (idiag_fkinxy/=0) then
        lpenc_diagnos2d(i_rho)=.true.
        lpenc_diagnos2d(i_u2)=.true.
      endif
      if (idiag_uguxm/=0 .or. idiag_uguym/=0 .or. idiag_uguzm/=0) &
          lpenc_diagnos(i_ugu)=.true.
          lpenc_diagnos(i_rhougu)=.true.
      if (idiag_ugu2m/=0) lpenc_diagnos(i_ugu2)=.true.
      if (idiag_uguxmx/=0 .or. idiag_uguymx/=0 .or. idiag_uguzmx/=0 .or. &
          idiag_uguxmy/=0 .or. idiag_uguymy/=0 .or. idiag_uguzmy/=0 .or. &
          idiag_uguxmz/=0 .or. idiag_uguymz/=0 .or. idiag_uguzmz/=0) &
          lpenc_diagnos(i_ugu)=.true.
! check whether right variables are set for half-box calculations.
      if (idiag_urmsn/=0 .or. idiag_ormsn/=0 .or. idiag_oumn/=0) then
        if ((.not.lequatory).and.(.not.lequatorz)) then
          call stop_it("You have to set either of lequatory or lequatorz to true to calculate averages over half the box")
        else
          if (lequatory) write(*,*) 'pencil-criteria_hydro: box divided along y dirn'
          if (lequatorz) write(*,*) 'pencil-criteria_hydro: box divided along z dirn'
        endif
      else
      endif
!
    endsubroutine pencil_criteria_hydro
!***********************************************************************
    subroutine pencil_interdep_hydro(lpencil_in)
!
!  Interdependency among pencils from the Hydro module is specified here.
!
!  20-nov-04/anders: coded
!
      logical, dimension (npencils) :: lpencil_in
!
      if (lpencil_in(i_u2)) lpencil_in(i_uu)=.true.
      if (lpencil_in(i_divu)) lpencil_in(i_uij)=.true.
      if (lalways_use_gij_etc) lpencil_in(i_oo)=.true.
      if (lpencil_in(i_sij)) then
        lpencil_in(i_uij)=.true.
        lpencil_in(i_divu)=.true.
      endif
      if (lalways_use_gij_etc) then
        if (lpencil_in(i_del2u))    lpencil_in(i_graddivu)=.true.
        if (lpencil_in(i_graddivu)) lpencil_in(i_oo)=.true.
      endif
      if (lpencil_in(i_oo)) lpencil_in(i_uij)=.true.
      if (lpencil_in(i_o2)) lpencil_in(i_oo)=.true.
      if (lpencil_in(i_ou)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_oo)=.true.
      endif
      if (lpencil_in(i_ugu)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_uij)=.true.
      endif
      if (lpencil_in(i_rhougu)) then
        lpencil_in(i_rho)=.true.
        lpencil_in(i_ugu)=.true.
      endif
      if (lpencil_in(i_sij2)) lpencil_in(i_sij)=.true.
!
    endsubroutine pencil_interdep_hydro
!***********************************************************************
    subroutine calc_pencils_hydro(f,p)
!
!  Calculate Hydro pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  08-nov-04/tony: coded
!  26-mar-07/axel: started using the gij_etc routine
!
      use Deriv
      use Mpicomm, only: stop_it
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
! uu
      if (lpencil(i_uu)) p%uu=f(l1:l2,m,n,iux:iuz)
! u2
      if (lpencil(i_u2)) call dot2_mn(p%uu,p%u2)
! uij
      if (lpencil(i_uij)) call gij(f,iuu,p%uij,1)
! divu
      if (lpencil(i_divu)) call div_mn(p%uij,p%divu)
! sij
      if (lpencil(i_sij)) call traceless_strain(p%uij,p%divu,p%sij)
! sij2
      if (lpencil(i_sij2)) call multm2_mn(p%sij,p%sij2)
! oo (=curlu)
      if (lpencil(i_oo)) then
        call curl_mn(p%uij,p%oo)
      endif
! o2
      if (lpencil(i_o2)) call dot2_mn(p%oo,p%o2)
! ou
      if (lpencil(i_ou)) call dot_mn(p%oo,p%uu,p%ou)
! Useful to debug forcing - Dhruba
      if (loutest.and.lpencil(i_ou))then
!      write(*,*) lpencil(i_ou)
        outest = minval(p%ou)
        if (outest.lt.(-1.0d-8))then
          write(*,*) m,n,outest,maxval(p%ou),lpencil(i_ou)
          write(*,*)'WARNING : hydro:ou has different sign than relhel'
        else
        endif
      else
      endif
! ugu
      if (lpencil(i_ugu)) then
        call u_dot_grad(iuu,p%uij,p%uu,p%ugu)
      endif

      if (lpencil(i_rhougu)) then
       p%rhougu(:,1)=p%rho*p%ugu(:,1)
       p%rhougu(:,2)=p%rho*p%ugu(:,2)
       p%rhougu(:,3)=p%rho*p%ugu(:,3)
      endif
! ugu2
      if (lpencil(i_ugu2)) call dot2_mn(p%ugu,p%ugu2)
! del2u, graddivu
      if (lalways_use_gij_etc) then
        if (lpencil(i_graddivu)) then
          if (headtt.or.ldebug) print*,'calc_pencils_hydro: call gij_etc'
          call gij_etc(f,iuu,p%oij,GRADDIV=p%graddivu)
        endif
        if (lpencil(i_del2u)) then
          call curl_mn(p%oij,p%qq)
          p%del2u=p%graddivu-p%qq
       endif
!
!   Avoid warnings from pencil consistency check...
!   WL: WHY SHOULD WE WANT TO AVOID WARNINGS FROM THE PENCIL CHECK??
!
        !if (.not. lpencil(i_uij)) p%uij=0.0
        !if (.not. lpencil(i_graddivu)) p%graddivu=0.0
      else
        if (lpencil(i_del2u)) then
          if (lpencil(i_graddivu)) then
            call del2v_etc(f,iuu,DEL2=p%del2u,GRADDIV=p%graddivu)
          else
             call del2v(f,iuu,p%del2u)
          endif
        else
          if (lpencil(i_graddivu)) call del2v_etc(f,iuu,GRADDIV=p%graddivu)
        endif
      endif
!
    endsubroutine calc_pencils_hydro
!***********************************************************************
    subroutine duu_dt(df,p)
!
!  velocity evolution
!  calculate du/dt = - u.gradu - 2Omega x u + grav + Fvisc
!  pressure gradient force added in density and entropy modules.
!
!   7-jun-02/axel: incoporated from subroutine pde
!  10-jun-02/axel+mattias: added Coriolis force
!  23-jun-02/axel: glnrho and fvisc are now calculated in here
!  17-jun-03/ulf: ux2, uy2 and uz2 added as diagnostic quantities
!  27-jun-07/dhruba: differential rotation as subroutine call
!
      use Diagnostics
      use IO
      use Mpicomm, only: stop_it
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      integer :: j
!
      intent(in) :: p
      intent(out) :: df
!
!  Identify module and boundary conditions.
!
      call timing('duu_dt','entered',mnloop=.true.)
      if (headtt.or.ldebug) print*,'duu_dt: SOLVE'
      if (headtt) then
        call identify_bcs('ux',iux)
        call identify_bcs('uy',iuy)
        call identify_bcs('uz',iuz)
      endif
!
!  Advection term.
!
      if (ladvection_velocity) then
        df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)-p%ugu
      endif
!
!  Coriolis force, -2*Omega x u (unless lprecession=T)
!  Omega=(-sin_theta, 0, cos_theta), where theta corresponds to
!  colatitude. theta=0 places the box to the north pole and theta=90
!  the equator. Cartesian coordinates (x,y,z) now correspond to
!  (theta,phi,r) i.e. (south,east,up), in spherical polar coordinates
!
      if (Omega/=0) call coriolis_cartesian(df,p%uu,iux)
!
! calculate viscous force
!
      if (lviscosity) call calc_viscous_force(df,p)
!
!  ``uu/dx'' for timestep
!
      if (lfirst.and.ldt.and.ladvection_velocity) then
          advec_uu=abs(p%uu(:,1))*dx_1(l1:l2)+ &
                   abs(p%uu(:,2))*dy_1(  m  )+ &
                   abs(p%uu(:,3))*dz_1(  n  )
      endif
!
      if (headtt.or.ldebug) print*,'duu_dt: max(advec_uu) =',maxval(advec_uu)
!
!  write slices for output in wvid in run.f90
!  This must be done outside the diagnostics loop (accessed at different times).
!  Note: ix is the index with respect to array with ghost zones.
!
      if (lvideo.and.lfirst) then
        divu_yz(m-m1+1,n-n1+1)=p%divu(ix_loc-l1+1)
        if (m.eq.iy_loc)  divu_xz(:,n-n1+1)=p%divu
        if (n.eq.iz_loc)  divu_xy(:,m-m1+1)=p%divu
        if (n.eq.iz2_loc) divu_xy2(:,m-m1+1)=p%divu
        do j=1,3
          oo_yz(m-m1+1,n-n1+1,j)=p%oo(ix_loc-l1+1,j)
          if (m==iy_loc)  oo_xz(:,n-n1+1,j)=p%oo(:,j)
          if (n==iz_loc)  oo_xy(:,m-m1+1,j)=p%oo(:,j)
          if (n==iz2_loc) oo_xy2(:,m-m1+1,j)=p%oo(:,j)
        enddo
        u2_yz(m-m1+1,n-n1+1)=p%u2(ix_loc-l1+1)
        if (m==iy_loc)  u2_xz(:,n-n1+1)=p%u2
        if (n==iz_loc)  u2_xy(:,m-m1+1)=p%u2
        if (n==iz2_loc) u2_xy2(:,m-m1+1)=p%u2
        o2_yz(m-m1+1,n-n1+1)=p%o2(ix_loc-l1+1)
        if (m==iy_loc)  o2_xz(:,n-n1+1)=p%o2
        if (n==iz_loc)  o2_xy(:,m-m1+1)=p%o2
        if (n==iz2_loc) o2_xy2(:,m-m1+1)=p%o2
        if (m.eq.iy_loc)  mach_xz(:,n-n1+1)=p%Ma2
        if (n.eq.iz_loc)  mach_xy(:,m-m1+1)=p%Ma2
        if (n.eq.iz2_loc) mach_xy2(:,m-m1+1)=p%Ma2
        mach_yz(m-m1+1,n-n1+1)=p%Ma2(ix_loc-l1+1)
        call vecout(41,trim(directory)//'/ovec',p%oo,othresh,novec)
      endif
!
!  Calculate maxima and rms values for diagnostic purposes
!
      call timing('duu_dt','just before ldiagnos',mnloop=.true.)
      if (ldiagnos) then
        if (headtt.or.ldebug) print*,'duu_dt: Calculate maxima and rms values...'
        if (idiag_dtu/=0) call max_mn_name(advec_uu/cdt,idiag_dtu,l_dt=.true.)
        if (idiag_urms/=0)   call sum_mn_name(p%u2,idiag_urms,lsqrt=.true.)
        if (idiag_urmsh/=0) then
          if (lequatory) call sum_mn_name_halfy(p%u2,idiag_urmsh)
          if (lequatorz) call sum_mn_name_halfz(p%u2,idiag_urmsh)
          fname(idiag_urmsn)=fname_half(idiag_urmsh,1)
          fname(idiag_urmss)=fname_half(idiag_urmsh,2)
          itype_name(idiag_urmsn)=ilabel_sum_sqrt
          itype_name(idiag_urmss)=ilabel_sum_sqrt
        else
        endif
        if (idiag_umax/=0)   call max_mn_name(p%u2,idiag_umax,lsqrt=.true.)
        if (idiag_uzrms/=0) &
            call sum_mn_name(p%uu(:,3)**2,idiag_uzrms,lsqrt=.true.)
        if (idiag_uzrmaxs/=0) &
            call max_mn_name(p%uu(:,3)**2,idiag_uzrmaxs,lsqrt=.true.)
        if (idiag_uxmin/=0) call max_mn_name(-p%uu(:,1),idiag_uxmin,lneg=.true.)
        if (idiag_uymin/=0) call max_mn_name(-p%uu(:,2),idiag_uymin,lneg=.true.)
        if (idiag_uzmin/=0) call max_mn_name(-p%uu(:,3),idiag_uzmin,lneg=.true.)
        if (idiag_uxmax/=0) call max_mn_name(p%uu(:,1),idiag_uxmax)
        if (idiag_uymax/=0) call max_mn_name(p%uu(:,2),idiag_uymax)
        if (idiag_uzmax/=0) call max_mn_name(p%uu(:,3),idiag_uzmax)
        if (idiag_rumax/=0) call max_mn_name(p%u2*p%rho**2,idiag_rumax,lsqrt=.true.)
        if (idiag_uguxm/=0) call sum_mn_name(p%ugu(:,1),idiag_uguxm)
        if (idiag_uguym/=0) call sum_mn_name(p%ugu(:,2),idiag_uguym)
        if (idiag_uguzm/=0) call sum_mn_name(p%ugu(:,3),idiag_uguzm)
        if (idiag_ugu2m/=0) call sum_mn_name(p%ugu2,idiag_ugu2m)
        if (idiag_u2m/=0)     call sum_mn_name(p%u2,idiag_u2m)
        if (idiag_um2/=0)     call max_mn_name(p%u2,idiag_um2)
        if (idiag_divum/=0)   call sum_mn_name(p%divu,idiag_divum)
        if (idiag_divu2m/=0)  call sum_mn_name(p%divu**2,idiag_divu2m)
        if (idiag_uxm/=0)     call sum_mn_name(p%uu(:,1),idiag_uxm)
        if (idiag_uym/=0)     call sum_mn_name(p%uu(:,2),idiag_uym)
        if (idiag_uzm/=0)     call sum_mn_name(p%uu(:,3),idiag_uzm)
        if (idiag_ux2m/=0)    call sum_mn_name(p%uu(:,1)**2,idiag_ux2m)
        if (idiag_uy2m/=0)    call sum_mn_name(p%uu(:,2)**2,idiag_uy2m)
        if (idiag_uz2m/=0)    call sum_mn_name(p%uu(:,3)**2,idiag_uz2m)
        if (idiag_ux2ccm/=0)  call sum_mn_name(c2z(n)*p%uu(:,1)**2,idiag_ux2ccm)
        if (idiag_ux2ssm/=0)  call sum_mn_name(s2z(n)*p%uu(:,1)**2,idiag_ux2ssm)
        if (idiag_uy2ccm/=0)  call sum_mn_name(c2z(n)*p%uu(:,2)**2,idiag_uy2ccm)
        if (idiag_uy2ssm/=0)  call sum_mn_name(s2z(n)*p%uu(:,2)**2,idiag_uy2ssm)
        if (idiag_uxuycsm/=0) call sum_mn_name(cz(n)*sz(n)*p%uu(:,1)*p%uu(:,2),idiag_uxuycsm)
        if (idiag_uxuym/=0)   call sum_mn_name(p%uu(:,1)*p%uu(:,2),idiag_uxuym)
        if (idiag_uxuzm/=0)   call sum_mn_name(p%uu(:,1)*p%uu(:,3),idiag_uxuzm)
        if (idiag_uyuzm/=0)   call sum_mn_name(p%uu(:,2)*p%uu(:,3),idiag_uyuzm)
        if (idiag_ruxuym/=0)  call sum_mn_name(p%rho*p%uu(:,1)*p%uu(:,2),idiag_ruxuym)
        if (idiag_ruxuzm/=0)  call sum_mn_name(p%rho*p%uu(:,1)*p%uu(:,3),idiag_ruxuzm)
        if (idiag_ruyuzm/=0)  call sum_mn_name(p%rho*p%uu(:,2)*p%uu(:,3),idiag_ruyuzm)
        if (idiag_duxdzma/=0) call sum_mn_name(abs(p%uij(:,1,3)),idiag_duxdzma)
        if (idiag_duydzma/=0) call sum_mn_name(abs(p%uij(:,2,3)),idiag_duydzma)
!
        if (idiag_rux2m/=0)    call sum_mn_name(p%rho*p%uu(:,1)**2,idiag_rux2m)
        if (idiag_ruy2m/=0)    call sum_mn_name(p%rho*p%uu(:,2)**2,idiag_ruy2m)
        if (idiag_ruz2m/=0)    call sum_mn_name(p%rho*p%uu(:,3)**2,idiag_ruz2m)
        if (idiag_ekin/=0)  call sum_mn_name(.5*p%rho*p%u2,idiag_ekin)
        if (idiag_ekintot/=0) &
            call integrate_mn_name(.5*p%rho*p%u2,idiag_ekintot)
        if (idiag_uxglnrym/=0)  call sum_mn_name(p%uu(:,1)*p%glnrho(:,2),idiag_uxglnrym)
        if (idiag_uyglnrxm/=0)  call sum_mn_name(p%uu(:,2)*p%glnrho(:,1),idiag_uyglnrxm)
        if (idiag_uxuydivum/=0) call sum_mn_name(p%uu(:,1)*p%uu(:,2)*p%divu,idiag_uxuydivum)
!
!  Kinetic field components at one point (=pt).
!
        if (lroot.and.m==mpoint.and.n==npoint) then
          if (idiag_uxpt/=0) call save_name(p%uu(lpoint-nghost,1),idiag_uxpt)
          if (idiag_uypt/=0) call save_name(p%uu(lpoint-nghost,2),idiag_uypt)
          if (idiag_uzpt/=0) call save_name(p%uu(lpoint-nghost,3),idiag_uzpt)
        endif
        if (idiag_u2mz/=0)  call xysum_mn_name_z(p%u2,idiag_u2mz)
!
!  Mean momenta.
!
        if (idiag_ruxm/=0) call sum_mn_name(p%rho*p%uu(:,1),idiag_ruxm)
        if (idiag_ruym/=0) call sum_mn_name(p%rho*p%uu(:,2),idiag_ruym)
        if (idiag_ruzm/=0) call sum_mn_name(p%rho*p%uu(:,3),idiag_ruzm)
        if (idiag_ruxtot/=0) call sum_mn_name(p%rho*abs(p%uu(:,1)),idiag_ruxtot)
!
!  Mean angular momenta.
!
        if (idiag_rlxm/=0) call sum_mn_name( &
            p%rho*(y(m)*p%uu(:,3)-z(n)*p%uu(:,2)),idiag_rlxm)
        if (idiag_rlym/=0) call sum_mn_name( &
            p%rho*(z(n)*p%uu(:,1)-x(l1:l2)*p%uu(:,3)),idiag_rlym)
        if (idiag_rlzm/=0) call sum_mn_name( &
            p%rho*(x(l1:l2)*p%uu(:,2)-y(m)*p%uu(:,1)),idiag_rlzm)
        if (idiag_rlx2m/=0) call sum_mn_name( &
            (p%rho*(y(m)*p%uu(:,3)-z(n)*p%uu(:,2)))**2,idiag_rlx2m)
        if (idiag_rly2m/=0) call sum_mn_name( &
            (p%rho*(z(n)*p%uu(:,1)-x(l1:l2)*p%uu(:,3)))**2,idiag_rly2m)
        if (idiag_rlz2m/=0) call sum_mn_name( &
            (p%rho*(x(l1:l2)*p%uu(:,2)-y(m)*p%uu(:,1)))**2,idiag_rlz2m)
!
!  Things related to vorticity.
!
        if (idiag_oum/=0) call sum_mn_name(p%ou,idiag_oum)
        if (idiag_oumh/=0) then
          if (lequatory) call sum_mn_name_halfy(p%ou,idiag_oumh)
          if (lequatorz) call sum_mn_name_halfz(p%ou,idiag_oumh)
          fname(idiag_oumn)=fname_half(idiag_oumh,1)
          fname(idiag_oums)=fname_half(idiag_oumh,2)
          itype_name(idiag_oumn)=ilabel_sum
          itype_name(idiag_oums)=ilabel_sum
        else
        endif
        if (idiag_orms/=0) call sum_mn_name(p%o2,idiag_orms,lsqrt=.true.)
        if (idiag_ormsh/=0) then
          if (lequatory) call sum_mn_name_halfy(p%o2,idiag_ormsh)
          if (lequatorz) call sum_mn_name_halfz(p%o2,idiag_ormsh)
          fname(idiag_ormsn)=fname_half(idiag_ormsh,1)
          fname(idiag_ormss)=fname_half(idiag_ormsh,2)
          itype_name(idiag_ormsn)=ilabel_sum_sqrt
          itype_name(idiag_ormss)=ilabel_sum_sqrt
        else
        endif
        if (idiag_omax/=0) call max_mn_name(p%o2,idiag_omax,lsqrt=.true.)
        if (idiag_o2m/=0)  call sum_mn_name(p%o2,idiag_o2m)
        if (idiag_ox2m/=0) call sum_mn_name(p%oo(:,1)**2,idiag_ox2m)
        if (idiag_oy2m/=0) call sum_mn_name(p%oo(:,2)**2,idiag_oy2m)
        if (idiag_oz2m/=0) call sum_mn_name(p%oo(:,3)**2,idiag_oz2m)
        if (idiag_oxm /=0) call sum_mn_name(p%oo(:,1)   ,idiag_oxm)
        if (idiag_oym /=0) call sum_mn_name(p%oo(:,2)   ,idiag_oym)
        if (idiag_ozm /=0) call sum_mn_name(p%oo(:,3)   ,idiag_ozm)
        if (idiag_oxoym/=0) call sum_mn_name(p%oo(:,1)*p%oo(:,2),idiag_oxoym)
        if (idiag_oxozm/=0) call sum_mn_name(p%oo(:,1)*p%oo(:,3),idiag_oxozm)
        if (idiag_oyozm/=0) call sum_mn_name(p%oo(:,2)*p%oo(:,3),idiag_oyozm)
        if (idiag_pvzm/=0) call sum_mn_name((p%oo(:,3) + 2.*Omega)/p%rho,idiag_pvzm)
!
!  Mach number, rms and max
!
        if (idiag_Marms/=0) call sum_mn_name(p%Ma2,idiag_Marms,lsqrt=.true.)
        if (idiag_Mamax/=0) call max_mn_name(p%Ma2,idiag_Mamax,lsqrt=.true.)
!
      endif
!
!  1d-averages. Happens at every it1d timesteps, NOT at every it1
!
      if (l1davgfirst) then
        if (idiag_fmassz/=0) call xysum_mn_name_z(p%rho*p%uu(:,3),idiag_fmassz)
        if (idiag_fkinz/=0)  call xysum_mn_name_z(.5*p%rho*p%u2*p%uu(:,3),idiag_fkinz)
        if (idiag_uxmz/=0)   call xysum_mn_name_z(p%uu(:,1),idiag_uxmz)
        if (idiag_uymz/=0)   call xysum_mn_name_z(p%uu(:,2),idiag_uymz)
        if (idiag_uzmz/=0)   call xysum_mn_name_z(p%uu(:,3),idiag_uzmz)
        if (idiag_oxmz/=0)   call xysum_mn_name_z(p%oo(:,1),idiag_oxmz)
        if (idiag_oymz/=0)   call xysum_mn_name_z(p%oo(:,2),idiag_oymz)
        if (idiag_ozmz/=0)   call xysum_mn_name_z(p%oo(:,3),idiag_ozmz)
        if (idiag_uxmy/=0)   call xzsum_mn_name_y(p%uu(:,1),idiag_uxmy)
        if (idiag_uymy/=0)   call xzsum_mn_name_y(p%uu(:,2),idiag_uymy)
        if (idiag_uzmy/=0)   call xzsum_mn_name_y(p%uu(:,3),idiag_uzmy)
        if (idiag_uxmx/=0)   call yzsum_mn_name_x(p%uu(:,1),idiag_uxmx)
        if (idiag_uymx/=0)   call yzsum_mn_name_x(p%uu(:,2),idiag_uymx)
        if (idiag_uzmx/=0)   call yzsum_mn_name_x(p%uu(:,3),idiag_uzmx)
        if (idiag_ux2mz/=0)  call xysum_mn_name_z(p%uu(:,1)**2,idiag_ux2mz)
        if (idiag_uy2mz/=0)  call xysum_mn_name_z(p%uu(:,2)**2,idiag_uy2mz)
        if (idiag_uz2mz/=0)  call xysum_mn_name_z(p%uu(:,3)**2,idiag_uz2mz)
        if (idiag_rux2mz/=0) call xysum_mn_name_z(p%rho*p%uu(:,1)**2,idiag_rux2mz)
        if (idiag_ruy2mz/=0) call xysum_mn_name_z(p%rho*p%uu(:,2)**2,idiag_ruy2mz)
        if (idiag_ruz2mz/=0) call xysum_mn_name_z(p%rho*p%uu(:,3)**2,idiag_ruz2mz)
        if (idiag_ux2my/=0)  call xzsum_mn_name_y(p%uu(:,1)**2,idiag_ux2my)
        if (idiag_uy2my/=0)  call xzsum_mn_name_y(p%uu(:,2)**2,idiag_uy2my)
        if (idiag_uz2my/=0)  call xzsum_mn_name_y(p%uu(:,3)**2,idiag_uz2my)
        if (idiag_ux2mx/=0)  call yzsum_mn_name_x(p%uu(:,1)**2,idiag_ux2mx)
        if (idiag_uy2mx/=0)  call yzsum_mn_name_x(p%uu(:,2)**2,idiag_uy2mx)
        if (idiag_uz2mx/=0)  call yzsum_mn_name_x(p%uu(:,3)**2,idiag_uz2mx)
        if (idiag_ox2mx/=0)  call yzsum_mn_name_x(p%oo(:,1)**2,idiag_ox2mx)
        if (idiag_oy2mx/=0)  call yzsum_mn_name_x(p%oo(:,2)**2,idiag_oy2mx)
        if (idiag_oz2mx/=0)  call yzsum_mn_name_x(p%oo(:,3)**2,idiag_oz2mx)
        if (idiag_uxuymz/=0) &
            call xysum_mn_name_z(p%uu(:,1)*p%uu(:,2),idiag_uxuymz)
        if (idiag_uxuzmz/=0) &
            call xysum_mn_name_z(p%uu(:,1)*p%uu(:,3),idiag_uxuzmz)
        if (idiag_uyuzmz/=0) &
            call xysum_mn_name_z(p%uu(:,2)*p%uu(:,3),idiag_uyuzmz)
        if (idiag_ruxuymz/=0) &
          call xysum_mn_name_z(p%rho*p%uu(:,1)*p%uu(:,2),idiag_ruxuymz)
        if (idiag_uxuymy/=0) &
            call xzsum_mn_name_y(p%uu(:,1)*p%uu(:,2),idiag_uxuymy)
        if (idiag_uxuzmy/=0) &
            call xzsum_mn_name_y(p%uu(:,1)*p%uu(:,3),idiag_uxuzmy)
        if (idiag_uyuzmy/=0) &
            call xzsum_mn_name_y(p%uu(:,2)*p%uu(:,3),idiag_uyuzmy)
        if (idiag_uxuymx/=0) &
            call yzsum_mn_name_x(p%uu(:,1)*p%uu(:,2),idiag_uxuymx)
        if (idiag_uxuzmx/=0) &
            call yzsum_mn_name_x(p%uu(:,1)*p%uu(:,3),idiag_uxuzmx)
        if (idiag_uyuzmx/=0) &
            call yzsum_mn_name_x(p%uu(:,2)*p%uu(:,3),idiag_uyuzmx)
        if (idiag_ekinz/=0)  call xysum_mn_name_z(.5*p%rho*p%u2,idiag_ekinz)
        if (idiag_oumx/=0)   call yzsum_mn_name_x(p%ou,idiag_oumx)
        if (idiag_oumy/=0)   call xzsum_mn_name_y(p%ou,idiag_oumy)
        if (idiag_oumz/=0)   call xysum_mn_name_z(p%ou,idiag_oumz)
        if (idiag_uguxmx/=0) call yzsum_mn_name_x(p%ugu(:,1),idiag_uguxmx)
        if (idiag_uguymx/=0) call yzsum_mn_name_x(p%ugu(:,2),idiag_uguymx)
        if (idiag_uguzmx/=0) call yzsum_mn_name_x(p%ugu(:,3),idiag_uguzmx)
        if (idiag_uguxmy/=0) call xzsum_mn_name_y(p%ugu(:,1),idiag_uguxmy)
        if (idiag_uguymy/=0) call xzsum_mn_name_y(p%ugu(:,2),idiag_uguymy)
        if (idiag_uguzmy/=0) call xzsum_mn_name_y(p%ugu(:,3),idiag_uguzmy)
        if (idiag_uguxmz/=0) call xysum_mn_name_z(p%ugu(:,1),idiag_uguxmz)
        if (idiag_uguymz/=0) call xysum_mn_name_z(p%ugu(:,2),idiag_uguymz)
        if (idiag_uguzmz/=0) call xysum_mn_name_z(p%ugu(:,3),idiag_uguzmz)
      endif
!
!  2-D averages.
!  Note that this does not necessarily happen with ldiagnos=.true.
!
      if (l2davgfirst) then
        if (idiag_uxmxz/=0) &
            call ysum_mn_name_xz(p%uu(:,1),idiag_uxmxz)
        if (idiag_uymxz/=0) &
            call ysum_mn_name_xz(p%uu(:,2),idiag_uymxz)
        if (idiag_uzmxz/=0) &
            call ysum_mn_name_xz(p%uu(:,3),idiag_uzmxz)
        if (idiag_ux2mxz/=0) &
            call ysum_mn_name_xz(p%uu(:,1)**2,idiag_ux2mxz)
        if (idiag_uy2mxz/=0) &
            call ysum_mn_name_xz(p%uu(:,2)**2,idiag_uy2mxz)
        if (idiag_uz2mxz/=0) &
            call ysum_mn_name_xz(p%uu(:,3)**2,idiag_uz2mxz)
        if (idiag_oumxz/=0) &
            call ysum_mn_name_xz(p%ou,idiag_oumxz)
        if (idiag_uxmxy/=0) call zsum_mn_name_xy(p%uu(:,1),idiag_uxmxy)
        if (idiag_uymxy/=0) call zsum_mn_name_xy(p%uu(:,2),idiag_uymxy)
        if (idiag_uzmxy/=0) call zsum_mn_name_xy(p%uu(:,3),idiag_uzmxy)
        if (idiag_oxmxy/=0) call zsum_mn_name_xy(p%oo(:,1),idiag_oxmxy)
        if (idiag_oymxy/=0) call zsum_mn_name_xy(p%oo(:,2),idiag_oymxy)
        if (idiag_ozmxy/=0) call zsum_mn_name_xy(p%oo(:,3),idiag_ozmxy)
        if (idiag_oumxy/=0) call zsum_mn_name_xy(p%ou,idiag_oumxy)
        if (idiag_pvzmxy/=0) call zsum_mn_name_xy((p%oo(:,3)+2.*Omega)/p%rho,idiag_pvzmxy)
        if (idiag_ruxmxy/=0) call zsum_mn_name_xy(p%rho*p%uu(:,1),idiag_ruxmxy)
        if (idiag_ruymxy/=0) call zsum_mn_name_xy(p%rho*p%uu(:,2),idiag_ruymxy)
        if (idiag_ruzmxy/=0) call zsum_mn_name_xy(p%rho*p%uu(:,3),idiag_ruzmxy)
        if (idiag_ux2mxy/=0) &
            call zsum_mn_name_xy(p%uu(:,1)**2,idiag_ux2mxy)
        if (idiag_uy2mxy/=0) &
            call zsum_mn_name_xy(p%uu(:,2)**2,idiag_uy2mxy)
        if (idiag_uz2mxy/=0) &
            call zsum_mn_name_xy(p%uu(:,3)**2,idiag_uz2mxy)
        if (idiag_rux2mxy/=0) &
            call zsum_mn_name_xy(p%rho*p%uu(:,1)**2,idiag_rux2mxy)
        if (idiag_ruy2mxy/=0) &
            call zsum_mn_name_xy(p%rho*p%uu(:,2)**2,idiag_ruy2mxy)
        if (idiag_ruz2mxy/=0) &
            call zsum_mn_name_xy(p%rho*p%uu(:,3)**2,idiag_ruz2mxy)
        if (idiag_ruxuymxy/=0) &
            call zsum_mn_name_xy(p%rho*p%uu(:,1)*p%uu(:,2),idiag_ruxuymxy)
        if (idiag_ruxuzmxy/=0) &
            call zsum_mn_name_xy(p%rho*p%uu(:,1)*p%uu(:,3),idiag_ruxuzmxy)
        if (idiag_ruyuzmxy/=0) &
            call zsum_mn_name_xy(p%rho*p%uu(:,2)*p%uu(:,3),idiag_ruyuzmxy)
        if (idiag_fkinxy/=0)  call zsum_mn_name_xy(.5*p%rho*p%u2*p%uu(:,1),idiag_fkinxy)
      else
!
!  idiag_uxmxy and idiag_uymxy also need to be calculated when
!  ldiagnos and idiag_umx and/or idiag_umy, so
!
!  We may need to calculate uxmxy without calculating umx. The following
!  if condition was messing up calculation of umxy_rms
!
        if (ldiagnos) then
          if (idiag_uxmxy/=0) call zsum_mn_name_xy(p%uu(:,1),idiag_uxmxy)
          if (idiag_uymxy/=0) call zsum_mn_name_xy(p%uu(:,2),idiag_uymxy)
          if (idiag_uzmxy/=0) call zsum_mn_name_xy(p%uu(:,3),idiag_uzmxy)
        endif
      endif
      call timing('duu_dt','finished',mnloop=.true.)
!
    endsubroutine duu_dt
!***********************************************************************
   subroutine coriolis_cartesian(df,uu,velind)
!
!  coriolis terms for cartesian geometry
!
!  30-oct-09/MR: outsourced, parameter velind added
!  checked to be an equivalent change by auot-test conv-slab-noequi, mdwarf
!
      real, dimension (mx,my,mz,mvar), intent(out) :: df
      real, dimension (nx,3),          intent(in)  :: uu
      integer,                         intent(in)  :: velind

! velind is start index for velocity variable to which Coriolis force corresponds
! x,y,z -components referred to by velind, velind+1, velind+2      (MR:IMMER ERF†LLT?)

      real :: c2, s2

      if (Omega==0.) return

      if (theta==0) then

        if (lcoriolis_force) then

          if (headtt) print*,'duu_dt: add Coriolis force; Omega=',Omega

          c2=2*Omega
          df(l1:l2,m,n,velind  )=df(l1:l2,m,n,velind  )+c2*uu(:,2)
          df(l1:l2,m,n,velind+1)=df(l1:l2,m,n,velind+1)-c2*uu(:,1)

        endif
!
!  add centrifugal force (doing this with periodic boundary
!  conditions in x and y would not be compatible, so it is
!  therefore usually ignored in those cases!)
!
        if (lcentrifugal_force) then

          if (headtt) print*,'duu_dt: add Centrifugal force; Omega=',Omega
          df(l1:l2,m,n,velind  )=df(l1:l2,m,n,velind  )+x(l1:l2)*Omega**2
          df(l1:l2,m,n,velind+1)=df(l1:l2,m,n,velind+1)+y(  m  )*Omega**2

        endif

      else
!
!  add Coriolis force with an angle (defined such that theta=60,
!  for example, would correspond to 30 degrees latitude).
!  Omega=(-sin_theta, 0, cos_theta).
!
        if (lcoriolis_force) then

          if (headtt) &
            print*,'duu_dt: Coriolis force; Omega, theta=', Omega, theta

          c2= 2*Omega*cos(theta*pi/180.)
          s2=-2*Omega*sin(theta*pi/180.)

          df(l1:l2,m,n,velind  )=df(l1:l2,m,n,velind  )+c2*uu(:,2)
          df(l1:l2,m,n,velind+1)=df(l1:l2,m,n,velind+1)-c2*uu(:,1)+s2*uu(:,3)
          df(l1:l2,m,n,velind+2)=df(l1:l2,m,n,velind+2)           -s2*uu(:,2)

        endif

      endif

   endsubroutine coriolis_cartesian
!***********************************************************************
    subroutine traceless_strain(uij,divu,sij)
!
!  Calculates traceless rate-of-strain tensor sij from derivative tensor uij
!  and divergence divu within each pencil;
!
!  16-oct-09/MR: carved out from calc_pencils_hydro
!
    real, dimension (nx,3,3)         :: uij, sij
    real, dimension (nx)             :: divu
!
    integer :: i,j
!
    intent(in)  :: uij, divu
    intent(out) :: sij
!
!  In-place operation is possible, i.e. uij and sij may refer to the same array.
!
    do j=1,3
      sij(:,j,j)=uij(:,j,j)
      do i=j+1,3
        sij(:,i,j)=.5*(uij(:,i,j)+uij(:,j,i))
        sij(:,j,i)=sij(:,i,j)
      enddo
      sij(:,j,j)=sij(:,j,j)-(1./3.)*divu
    enddo
!
    if (lshear) then
      if (lshear_rateofstrain) then
        sij(:,1,2)=sij(:,1,2)+Sshear
        sij(:,2,1)=sij(:,2,1)+Sshear
      endif
    endif
!
    endsubroutine traceless_strain
!***********************************************************************
    subroutine read_hydro_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=hydro_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=hydro_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_hydro_init_pars
!***********************************************************************
    subroutine write_hydro_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=hydro_init_pars)
!
    endsubroutine write_hydro_init_pars
!***********************************************************************
    subroutine read_hydro_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=hydro_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=hydro_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_hydro_run_pars
!***********************************************************************
    subroutine write_hydro_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=hydro_run_pars)
!
    endsubroutine write_hydro_run_pars
!***********************************************************************
    subroutine rprint_hydro(lreset,lwrite)
!
!  reads and registers print parameters relevant for hydro part
!
!   3-may-02/axel: coded
!  27-may-02/axel: added possibility to reset list
!
      use Diagnostics
!
      integer :: iname,inamez,inamey,inamex,ixy,ixz,iname_half
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_u2tm=0
        idiag_uotm=0
        idiag_outm=0
        idiag_u2m=0
        idiag_um2=0
        idiag_uxpt=0
        idiag_uypt=0
        idiag_uzpt=0
        idiag_urms=0
        idiag_umax=0
        idiag_uzrms=0
        idiag_uzrmaxs=0
        idiag_uxmin=0
        idiag_uymin=0
        idiag_uzmin=0
        idiag_uxmax=0
        idiag_uymax=0
        idiag_uzmax=0
        idiag_uxm=0
        idiag_uym=0
        idiag_uzm=0
        idiag_ux2m=0
        idiag_uy2m=0
        idiag_uz2m=0
        idiag_ux2ccm=0
        idiag_ux2ssm=0
        idiag_uy2ccm=0
        idiag_uy2ssm=0
        idiag_uxuycsm=0
        idiag_rux2m=0
        idiag_ruy2m=0
        idiag_ruz2m=0
        idiag_ux2mx=0
        idiag_uy2mx=0
        idiag_uz2mx=0
        idiag_ox2mx=0
        idiag_oy2mx=0
        idiag_oz2mx=0
        idiag_ux2my=0
        idiag_uy2my=0
        idiag_uz2my=0
        idiag_ux2mz=0
        idiag_uy2mz=0
        idiag_uz2mz=0
        idiag_rux2mz=0
        idiag_ruy2mz=0
        idiag_ruz2mz=0
        idiag_uxmz=0
        idiag_uymz=0
        idiag_uzmz=0
        idiag_oxmz=0
        idiag_oymz=0
        idiag_ozmz=0
        idiag_uxuym=0
        idiag_uxuzm=0
        idiag_uyuzm=0
        idiag_uxuymz=0
        idiag_uxuzmz=0
        idiag_uyuzmz=0
        idiag_uxuymz=0
        idiag_umx=0
        idiag_umy=0
        idiag_umz=0
        idiag_omumz=0
        idiag_umamz=0
        idiag_umbmz=0
        idiag_umxbmz=0
        idiag_divum=0
        idiag_divu2m=0
        idiag_uxmy=0
        idiag_uymy=0
        idiag_uzmy=0
        idiag_uxuymy=0
        idiag_uxuzmy=0
        idiag_uyuzmy=0
        idiag_u2mr=0
        idiag_urmr=0
        idiag_upmr=0
        idiag_uzmr=0
        idiag_uxfampm=0
        idiag_uyfampm=0
        idiag_uzfampm=0
        idiag_uxmxz=0
        idiag_uymxz=0
        idiag_uzmxz=0
        idiag_ux2mxz=0
        idiag_uy2mxz=0
        idiag_uz2mxz=0
        idiag_uxmxy=0
        idiag_uymxy=0
        idiag_uzmxy=0
        idiag_oxmxy=0
        idiag_oymxy=0
        idiag_ozmxy=0
        idiag_pvzmxy=0
        idiag_ruxmxy=0
        idiag_ruymxy=0
        idiag_ruzmxy=0
        idiag_ux2mxy=0
        idiag_uy2mxy=0
        idiag_uz2mxy=0
        idiag_rux2mxy=0
        idiag_ruy2mxy=0
        idiag_ruz2mxy=0
        idiag_ruxuymxy=0
        idiag_ruxuzmxy=0
        idiag_ruyuzmxy=0
        idiag_ruxm=0
        idiag_ruym=0
        idiag_ruzm=0
        idiag_ruxtot=0
        idiag_rlxm=0
        idiag_rlym=0
        idiag_rlzm=0
        idiag_rlx2m=0
        idiag_rly2m=0
        idiag_rlz2m=0
        idiag_tot_ang_mom=0
        idiag_rumax=0
        idiag_rufm=0
        idiag_dtu=0
        idiag_oum=0
        idiag_fum=0
        idiag_o2m=0
        idiag_orms=0
        idiag_omax=0
        idiag_ox2m=0
        idiag_oy2m=0
        idiag_oz2m=0
        idiag_oxm=0
        idiag_oym=0
        idiag_ozm=0
        idiag_oxoym=0
        idiag_oxozm=0
        idiag_oyozm=0
        idiag_pvzm=0
        idiag_oumx=0
        idiag_oumy=0
        idiag_oumz=0
        idiag_oumxy=0
        idiag_oumxz=0
        idiag_Marms=0
        idiag_Mamax=0
        idiag_fintm=0
        idiag_fextm=0
        idiag_duxdzma=0
        idiag_duydzma=0
        idiag_ekin=0
        idiag_ekintot=0
        idiag_ekinz=0
        idiag_fmassz=0
        idiag_fkinz=0
        idiag_fkinxy=0
        idiag_fxbxm=0
        idiag_fxbym=0
        idiag_fxbzm=0
        idiag_ruxuym=0
        idiag_ruxuzm=0
        idiag_ruyuzm=0
        idiag_ruxuymz=0
        idiag_uguxm=0
        idiag_uguym=0
        idiag_uguzm=0
        idiag_ugu2m=0
        idiag_uguxmx=0
        idiag_uguymx=0
        idiag_uguzmx=0
        idiag_uguxmy=0
        idiag_uguymy=0
        idiag_uguzmy=0
        idiag_uguxmz=0
        idiag_uguymz=0
        idiag_uguzmz=0
        idiag_uxglnrym=0
        idiag_uyglnrxm=0
        idiag_uxuydivum=0
        idiag_urmsh=0;idiag_urmsn=0;idiag_urmss=0
        idiag_ormsh=0;idiag_ormsn=0;idiag_ormss=0
        idiag_oumh=0;idiag_oumn=0;idiag_oums=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if (lroot.and.ip<14) print*,'rprint_hydro: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'ekin',idiag_ekin)
        call parse_name(iname,cname(iname),cform(iname),'ekintot',idiag_ekintot)
        call parse_name(iname,cname(iname),cform(iname),'u2tm',idiag_u2tm)
        call parse_name(iname,cname(iname),cform(iname),'uotm',idiag_uotm)
        call parse_name(iname,cname(iname),cform(iname),'outm',idiag_outm)
        call parse_name(iname,cname(iname),cform(iname),'u2m',idiag_u2m)
        call parse_name(iname,cname(iname),cform(iname),'um2',idiag_um2)
        call parse_name(iname,cname(iname),cform(iname),'o2m',idiag_o2m)
        call parse_name(iname,cname(iname),cform(iname),'oum',idiag_oum)
        call parse_name(iname,cname(iname),cform(iname),'fum',idiag_fum)
        call parse_name(iname,cname(iname),cform(iname),'oumn',idiag_oumn)
        call parse_name(iname,cname(iname),cform(iname),'oums',idiag_oums)
        call parse_name(iname,cname(iname),cform(iname),'dtu',idiag_dtu)
        call parse_name(iname,cname(iname),cform(iname),'urms',idiag_urms)
        call parse_name(iname,cname(iname),cform(iname),'urmsn',idiag_urmsn)
        call parse_name(iname,cname(iname),cform(iname),'urmss',idiag_urmss)
        call parse_name(iname,cname(iname),cform(iname),'umax',idiag_umax)
        call parse_name(iname,cname(iname),cform(iname),'uxmin',idiag_uxmin)
        call parse_name(iname,cname(iname),cform(iname),'uymin',idiag_uymin)
        call parse_name(iname,cname(iname),cform(iname),'uzmin',idiag_uzmin)
        call parse_name(iname,cname(iname),cform(iname),'uxmax',idiag_uxmax)
        call parse_name(iname,cname(iname),cform(iname),'uymax',idiag_uymax)
        call parse_name(iname,cname(iname),cform(iname),'uzmax',idiag_uzmax)
        call parse_name(iname,cname(iname),cform(iname),'uzrms',idiag_uzrms)
        call parse_name(iname,cname(iname),cform(iname),'uzrmaxs',idiag_uzrmaxs)
        call parse_name(iname,cname(iname),cform(iname),'uxm',idiag_uxm)
        call parse_name(iname,cname(iname),cform(iname),'uym',idiag_uym)
        call parse_name(iname,cname(iname),cform(iname),'uzm',idiag_uzm)
        call parse_name(iname,cname(iname),cform(iname),'ux2m',idiag_ux2m)
        call parse_name(iname,cname(iname),cform(iname),'uy2m',idiag_uy2m)
        call parse_name(iname,cname(iname),cform(iname),'uz2m',idiag_uz2m)
        call parse_name(iname,cname(iname),cform(iname),'ux2ccm',idiag_ux2ccm)
        call parse_name(iname,cname(iname),cform(iname),'ux2ssm',idiag_ux2ssm)
        call parse_name(iname,cname(iname),cform(iname),'uy2ccm',idiag_uy2ccm)
        call parse_name(iname,cname(iname),cform(iname),'uy2ssm',idiag_uy2ssm)
        call parse_name(iname,cname(iname),cform(iname),'uxuycsm',idiag_uxuycsm)
        call parse_name(iname,cname(iname),cform(iname),'rux2m',idiag_rux2m)
        call parse_name(iname,cname(iname),cform(iname),'ruy2m',idiag_ruy2m)
        call parse_name(iname,cname(iname),cform(iname),'ruz2m',idiag_ruz2m)
        call parse_name(iname,cname(iname),cform(iname),'uxuym',idiag_uxuym)
        call parse_name(iname,cname(iname),cform(iname),'uxuzm',idiag_uxuzm)
        call parse_name(iname,cname(iname),cform(iname),'uyuzm',idiag_uyuzm)
        call parse_name(iname,cname(iname),cform(iname),'ruxuym',idiag_ruxuym)
        call parse_name(iname,cname(iname),cform(iname),'ruxuzm',idiag_ruxuzm)
        call parse_name(iname,cname(iname),cform(iname),'ruyuzm',idiag_ruyuzm)
        call parse_name(iname,cname(iname),cform(iname),'ox2m',idiag_ox2m)
        call parse_name(iname,cname(iname),cform(iname),'oy2m',idiag_oy2m)
        call parse_name(iname,cname(iname),cform(iname),'oz2m',idiag_oz2m)
        call parse_name(iname,cname(iname),cform(iname),'oxm',idiag_oxm)
        call parse_name(iname,cname(iname),cform(iname),'oym',idiag_oym)
        call parse_name(iname,cname(iname),cform(iname),'ozm',idiag_ozm)
        call parse_name(iname,cname(iname),cform(iname),'oxoym',idiag_oxoym)
        call parse_name(iname,cname(iname),cform(iname),'oxozm',idiag_oxozm)
        call parse_name(iname,cname(iname),cform(iname),'oyozm',idiag_oyozm)
        call parse_name(iname,cname(iname),cform(iname),'pvzm',idiag_pvzm)
        call parse_name(iname,cname(iname),cform(iname),'orms',idiag_orms)
        call parse_name(iname,cname(iname),cform(iname),'ormsn',idiag_ormsn)
        call parse_name(iname,cname(iname),cform(iname),'ormss',idiag_ormss)
        call parse_name(iname,cname(iname),cform(iname),'omax',idiag_omax)
        call parse_name(iname,cname(iname),cform(iname),'ruxm',idiag_ruxm)
        call parse_name(iname,cname(iname),cform(iname),'ruym',idiag_ruym)
        call parse_name(iname,cname(iname),cform(iname),'ruzm',idiag_ruzm)
        call parse_name(iname,cname(iname),cform(iname),'ruxtot',idiag_ruxtot)
        call parse_name(iname,cname(iname),cform(iname),'rlxm',idiag_rlxm)
        call parse_name(iname,cname(iname),cform(iname),'rlym',idiag_rlym)
        call parse_name(iname,cname(iname),cform(iname),'rlzm',idiag_rlzm)
        call parse_name(iname,cname(iname),cform(iname),'rlx2m',idiag_rlx2m)
        call parse_name(iname,cname(iname),cform(iname),'rly2m',idiag_rly2m)
        call parse_name(iname,cname(iname),cform(iname),'rlz2m',idiag_rlz2m)
        call parse_name(iname,cname(iname),cform(iname),'tot_ang_mom',idiag_tot_ang_mom)
        call parse_name(iname,cname(iname),cform(iname),'rumax',idiag_rumax)
        call parse_name(iname,cname(iname),cform(iname),'umx',idiag_umx)
        call parse_name(iname,cname(iname),cform(iname),'umy',idiag_umy)
        call parse_name(iname,cname(iname),cform(iname),'umz',idiag_umz)
        call parse_name(iname,cname(iname),cform(iname),'omumz',idiag_omumz)
        call parse_name(iname,cname(iname),cform(iname),'umamz',idiag_umamz)
        call parse_name(iname,cname(iname),cform(iname),'umbmz',idiag_umbmz)
        call parse_name(iname,cname(iname),cform(iname),'umxbmz',idiag_umxbmz)
        call parse_name(iname,cname(iname),cform(iname),'Marms',idiag_Marms)
        call parse_name(iname,cname(iname),cform(iname),'Mamax',idiag_Mamax)
        call parse_name(iname,cname(iname),cform(iname),'divum',idiag_divum)
        call parse_name(iname,cname(iname),cform(iname),'divu2m',idiag_divu2m)
        call parse_name(iname,cname(iname),cform(iname),'uxpt',idiag_uxpt)
        call parse_name(iname,cname(iname),cform(iname),'uypt',idiag_uypt)
        call parse_name(iname,cname(iname),cform(iname),'uzpt',idiag_uzpt)
        call parse_name(iname,cname(iname),cform(iname),'fintm',idiag_fintm)
        call parse_name(iname,cname(iname),cform(iname),'fextm',idiag_fextm)
        call parse_name(iname,cname(iname),cform(iname),'duxdzma',idiag_duxdzma)
        call parse_name(iname,cname(iname),cform(iname),'duydzma',idiag_duydzma)
        call parse_name(iname,cname(iname),cform(iname),'rufm',idiag_rufm)
        call parse_name(iname,cname(iname),cform(iname),'fxbxm',idiag_fxbxm)
        call parse_name(iname,cname(iname),cform(iname),'fxbym',idiag_fxbym)
        call parse_name(iname,cname(iname),cform(iname),'fxbzm',idiag_fxbzm)
        call parse_name(iname,cname(iname),cform(iname),'uxfampm',idiag_uxfampm)
        call parse_name(iname,cname(iname),cform(iname),'uyfampm',idiag_uyfampm)
        call parse_name(iname,cname(iname),cform(iname),'uzfampm',idiag_uzfampm)
        call parse_name(iname,cname(iname),cform(iname),'uxfampim',idiag_uxfampim)
        call parse_name(iname,cname(iname),cform(iname),'uyfampim',idiag_uyfampim)
        call parse_name(iname,cname(iname),cform(iname),'uzfampim',idiag_uzfampim)
        call parse_name(iname,cname(iname),cform(iname),'uguxm',idiag_uguxm)
        call parse_name(iname,cname(iname),cform(iname),'uguym',idiag_uguym)
        call parse_name(iname,cname(iname),cform(iname),'uguzm',idiag_uguzm)
        call parse_name(iname,cname(iname),cform(iname),'ugu2m',idiag_ugu2m)
        call parse_name(iname,cname(iname),cform(iname),'uxglnrym',idiag_uxglnrym)
        call parse_name(iname,cname(iname),cform(iname),'uyglnrxm',idiag_uyglnrxm)
        call parse_name(iname,cname(iname),cform(iname),'uxuydivum',idiag_uxuydivum)
      enddo
!
! Quantities which are averaged over half (north-south) the box
!
      iname_half=name_half_max
      if ((idiag_urmsn/=0).or.(idiag_urmss/=0))then
        iname_half=iname_half+1
        idiag_urmsh=iname_half
      else
      endif
      if ((idiag_ormsn/=0).or.(idiag_ormss/=0))then
        iname_half=iname_half+1
        idiag_ormsh=iname_half
      else
      endif
      if ((idiag_oumn/=0).or.(idiag_oums/=0))then
        iname_half=iname_half+1
        idiag_oumh=iname_half
      else
      endif
      name_half_max=iname_half
!
!  Check for those quantities for which we want yz-averages.
!
      do inamex=1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'uxmx',idiag_uxmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'uymx',idiag_uymx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'uzmx',idiag_uzmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'ux2mx',idiag_ux2mx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'uy2mx',idiag_uy2mx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'uz2mx',idiag_uz2mx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'ox2mx',idiag_ox2mx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'oy2mx',idiag_oy2mx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'oz2mx',idiag_oz2mx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'uxuymx',idiag_uxuymx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'uxuzmx',idiag_uxuzmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'uyuzmx',idiag_uyuzmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'oumx',idiag_oumx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'uguxmx',idiag_uguxmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'uguymx',idiag_uguymx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'uguzmx',idiag_uguzmx)
      enddo
!
!  Check for those quantities for which we want xz-averages.
!
      do inamey=1,nnamey
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'uxmy',idiag_uxmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'uymy',idiag_uymy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'uzmy',idiag_uzmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey), &
            'ux2my',idiag_ux2my)
        call parse_name(inamey,cnamey(inamey),cformy(inamey), &
            'uy2my',idiag_uy2my)
        call parse_name(inamey,cnamey(inamey),cformy(inamey), &
            'uz2my',idiag_uz2my)
        call parse_name(inamey,cnamey(inamey),cformy(inamey), &
            'uxuymy',idiag_uxuymy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey), &
            'uxuzmy',idiag_uxuzmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey), &
            'uyuzmy',idiag_uyuzmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'oumy',idiag_oumy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey), &
            'uguxmy',idiag_uguxmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey), &
            'uguymy',idiag_uguymy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey), &
            'uguzmy',idiag_uguzmy)
      enddo
!
!  Check for those quantities for which we want xy-averages.
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uxmz',idiag_uxmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uymz',idiag_uymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uzmz',idiag_uzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'oxmz',idiag_oxmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'oymz',idiag_oymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ozmz',idiag_ozmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ux2mz',idiag_ux2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uy2mz',idiag_uy2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uz2mz',idiag_uz2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'rux2mz',idiag_rux2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ruy2mz',idiag_ruy2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ruz2mz',idiag_ruz2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'uxuymz',idiag_uxuymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'uxuzmz',idiag_uxuzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'uyuzmz',idiag_uyuzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'ruxuymz',idiag_ruxuymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'fmassz',idiag_fmassz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'fkinz',idiag_fkinz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'ekinz',idiag_ekinz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'u2mz',idiag_u2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'oumz',idiag_oumz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'uguxmz',idiag_uguxmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'uguymz',idiag_uguymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'uguzmz',idiag_uguzmz)
      enddo
!
!  Check for those quantities for which we want y-averages.
!
      do ixz=1,nnamexz
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'uxmxz',idiag_uxmxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'uymxz',idiag_uymxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'uzmxz',idiag_uzmxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'ux2mxz',idiag_ux2mxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'uy2mxz',idiag_uy2mxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'uz2mxz',idiag_uz2mxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'oumxz',idiag_oumxz)
      enddo
!
!  check for those quantities for which we want z-averages
!
      do ixy=1,nnamexy
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'uxmxy',idiag_uxmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'uymxy',idiag_uymxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'uzmxy',idiag_uzmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'oxmxy',idiag_oxmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'oymxy',idiag_oymxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'ozmxy',idiag_ozmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'oumxy',idiag_oumxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'pvzmxy',idiag_pvzmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'ruxmxy',idiag_ruxmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'ruymxy',idiag_ruymxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'ruzmxy',idiag_ruzmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'ux2mxy',idiag_ux2mxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'uy2mxy',idiag_uy2mxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'uz2mxy',idiag_uz2mxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'rux2mxy',idiag_rux2mxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'ruy2mxy',idiag_ruy2mxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'ruz2mxy',idiag_ruz2mxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'ruxuymxy',idiag_ruxuymxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'ruxuzmxy',idiag_ruxuzmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'ruyuzmxy',idiag_ruyuzmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'fkinxy',idiag_fkinxy)
      enddo
!
!  write column where which hydro variable is stored
!
      if (lwr) then
        write(3,*) 'i_ekin=',idiag_ekin
        write(3,*) 'i_ekintot=',idiag_ekintot
        write(3,*) 'i_u2tm=',idiag_u2tm
        write(3,*) 'i_uotm=',idiag_uotm
        write(3,*) 'i_outm=',idiag_outm
        write(3,*) 'i_u2m=',idiag_u2m
        write(3,*) 'i_um2=',idiag_um2
        write(3,*) 'i_o2m=',idiag_o2m
        write(3,*) 'i_oum=',idiag_oum
        write(3,*) 'i_fum=',idiag_fum
        write(3,*) 'i_oumn=',idiag_oumn
        write(3,*) 'i_oums=',idiag_oums
        write(3,*) 'i_oumh=',idiag_oumh
        write(3,*) 'i_dtu=',idiag_dtu
        write(3,*) 'i_urms=',idiag_urms
        write(3,*) 'i_urmsn=',idiag_urmsn
        write(3,*) 'i_urmss=',idiag_urmss
        write(3,*) 'i_urmsh=',idiag_urmsh
        write(3,*) 'i_umax=',idiag_umax
        write(3,*) 'i_uxmax=',idiag_uxmax
        write(3,*) 'i_uymax=',idiag_uymax
        write(3,*) 'i_uzmax=',idiag_uzmax
        write(3,*) 'i_uzrms=',idiag_uzrms
        write(3,*) 'i_uzrmaxs=',idiag_uzrmaxs
        write(3,*) 'i_ux2m=',idiag_ux2m
        write(3,*) 'i_uy2m=',idiag_uy2m
        write(3,*) 'i_uz2m=',idiag_uz2m
        write(3,*) 'i_ux2ccm=',idiag_ux2ccm
        write(3,*) 'i_ux2ssm=',idiag_ux2ssm
        write(3,*) 'i_uy2ccm=',idiag_uy2ccm
        write(3,*) 'i_uy2ssm=',idiag_uy2ssm
        write(3,*) 'i_uxuycsm=',idiag_uxuycsm
        write(3,*) 'i_rux2m=',idiag_rux2m
        write(3,*) 'i_ruy2m=',idiag_ruy2m
        write(3,*) 'i_ruz2m=',idiag_ruz2m
        write(3,*) 'i_uxuym=',idiag_uxuym
        write(3,*) 'i_uxuzm=',idiag_uxuzm
        write(3,*) 'i_uyuzm=',idiag_uyuzm
        write(3,*) 'i_ruxuym=',idiag_ruxuym
        write(3,*) 'i_ruxuzm=',idiag_ruxuzm
        write(3,*) 'i_ruyuzm=',idiag_ruyuzm
        write(3,*) 'i_ox2m=',idiag_ox2m
        write(3,*) 'i_oy2m=',idiag_oy2m
        write(3,*) 'i_oz2m=',idiag_oz2m
        write(3,*) 'i_oxm=',idiag_oxm
        write(3,*) 'i_oym=',idiag_oym
        write(3,*) 'i_ozm=',idiag_ozm
        write(3,*) 'i_oxoym=',idiag_oxoym
        write(3,*) 'i_oxozm=',idiag_oxozm
        write(3,*) 'i_oyozm=',idiag_oyozm
        write(3,*) 'i_pvzm=',idiag_pvzm
        write(3,*) 'i_orms=',idiag_orms
        write(3,*) 'i_ormsn=',idiag_ormsn
        write(3,*) 'i_ormss=',idiag_ormss
        write(3,*) 'i_ormsh=',idiag_ormsh
        write(3,*) 'i_omax=',idiag_omax
        write(3,*) 'i_ruxm=',idiag_ruxm
        write(3,*) 'i_ruym=',idiag_ruym
        write(3,*) 'i_ruzm=',idiag_ruzm
        write(3,*) 'i_ruxtot=',idiag_ruxtot
        write(3,*) 'i_rumax=',idiag_rumax
        write(3,*) 'i_umx=',idiag_umx
        write(3,*) 'i_umy=',idiag_umy
        write(3,*) 'i_umz=',idiag_umz
        write(3,*) 'i_omumz=',idiag_omumz
        write(3,*) 'i_umamz=',idiag_umamz
        write(3,*) 'i_umbmz=',idiag_umbmz
        write(3,*) 'i_umxbmz=',idiag_umxbmz
        write(3,*) 'i_Marms=',idiag_Marms
        write(3,*) 'i_Mamax=',idiag_Mamax
        write(3,*) 'i_divum=',idiag_divum
        write(3,*) 'i_divu2m=',idiag_divu2m
        write(3,*) 'i_uxfampm=',idiag_uxfampm
        write(3,*) 'i_uyfampm=',idiag_uyfampm
        write(3,*) 'i_uzfampm=',idiag_uzfampm
        write(3,*) 'i_uxfampim=',idiag_uxfampim
        write(3,*) 'i_uyfampim=',idiag_uyfampim
        write(3,*) 'i_uzfampim=',idiag_uzfampim
        write(3,*) 'i_uxpt=',idiag_uxpt
        write(3,*) 'i_uypt=',idiag_uypt
        write(3,*) 'i_uzpt=',idiag_uzpt
        write(3,*) 'i_fmassz=',idiag_fmassz
        write(3,*) 'i_fkinz=',idiag_fkinz
        write(3,*) 'i_fkinxy=',idiag_fkinxy
        write(3,*) 'i_ekinz=',idiag_ekinz
        write(3,*) 'i_uxmz=',idiag_uxmz
        write(3,*) 'i_uymz=',idiag_uymz
        write(3,*) 'i_uzmz=',idiag_uzmz
        write(3,*) 'i_oxmz=',idiag_oxmz
        write(3,*) 'i_oymz=',idiag_oymz
        write(3,*) 'i_ozmz=',idiag_ozmz
        write(3,*) 'i_ux2mz=',idiag_ux2mz
        write(3,*) 'i_uy2mz=',idiag_uy2mz
        write(3,*) 'i_uz2mz=',idiag_uz2mz
        write(3,*) 'i_rux2mz=',idiag_rux2mz
        write(3,*) 'i_ruy2mz=',idiag_ruy2mz
        write(3,*) 'i_ruz2mz=',idiag_ruz2mz
        write(3,*) 'i_uxmxy=',idiag_uxmxy
        write(3,*) 'i_uymxy=',idiag_uymxy
        write(3,*) 'i_uzmxy=',idiag_uzmxy
        write(3,*) 'i_oumxy=',idiag_oumxy
        write(3,*) 'i_oumxz=',idiag_oumxz
        write(3,*) 'i_ruxmxy=',idiag_ruxmxy
        write(3,*) 'i_ruymxy=',idiag_ruymxy
        write(3,*) 'i_ruzmxy=',idiag_ruzmxy
        write(3,*) 'i_ux2mxy=',idiag_ux2mxy
        write(3,*) 'i_uy2mxy=',idiag_uy2mxy
        write(3,*) 'i_uz2mxy=',idiag_uz2mxy
        write(3,*) 'i_rux2mxy=',idiag_rux2mxy
        write(3,*) 'i_ruy2mxy=',idiag_ruy2mxy
        write(3,*) 'i_ruz2mxy=',idiag_ruz2mxy
        write(3,*) 'i_ruxuymxy=',idiag_ruxuymxy
        write(3,*) 'i_ruxuzmxy=',idiag_ruxuzmxy
        write(3,*) 'i_ruyuzmxy=',idiag_ruyuzmxy
        write(3,*) 'i_uxmxz=',idiag_uxmxz
        write(3,*) 'i_uymxz=',idiag_uymxz
        write(3,*) 'i_uzmxz=',idiag_uzmxz
        write(3,*) 'i_ux2mxz=',idiag_ux2mxz
        write(3,*) 'i_uy2mxz=',idiag_uy2mxz
        write(3,*) 'i_uz2mxz=',idiag_uz2mxz
        write(3,*) 'i_u2mz=',idiag_u2mz
        write(3,*) 'i_fintm=',idiag_fintm
        write(3,*) 'i_fextm=',idiag_fextm
        write(3,*) 'i_duxdzma=',idiag_duxdzma
        write(3,*) 'i_duydzma=',idiag_duydzma
        write(3,*) 'i_rufm=',idiag_rufm
        write(3,*) 'i_fxbxm=',idiag_fxbxm
        write(3,*) 'i_fxbym=',idiag_fxbym
        write(3,*) 'i_fxbzm=',idiag_fxbzm
        write(3,*) 'nname=',nname
        write(3,*) 'iuu=',iuu
        write(3,*) 'iux=',iux
        write(3,*) 'iuy=',iuy
        write(3,*) 'iuz=',iuz
      endif
!
    endsubroutine rprint_hydro
!***********************************************************************
    subroutine get_slices_hydro(f,slices)
!
!  Write slices for animation of Hydro variables.
!
!  26-jul-06/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  Velocity field.
!
        case ('uu')
          if (slices%index>=3) then
            slices%ready=.false.
          else
            slices%index=slices%index+1
            slices%yz =f(slices%ix,m1:m2    ,n1:n2     ,iux-1+slices%index)
            slices%xz =f(l1:l2    ,slices%iy,n1:n2     ,iux-1+slices%index)
            slices%xy =f(l1:l2    ,m1:m2    ,slices%iz ,iux-1+slices%index)
            slices%xy2=f(l1:l2    ,m1:m2    ,slices%iz2,iux-1+slices%index)
            if (slices%index<=3) slices%ready=.true.
          endif
!
!  Divergence of velocity.
!
        case ('divu')
          slices%yz =>divu_yz
          slices%xz =>divu_xz
          slices%xy =>divu_xy
          slices%xy2=>divu_xy2
          slices%ready=.true.
!
!  Velocity squared.
!
        case ('u2')
          slices%yz =>u2_yz
          slices%xz =>u2_xz
          slices%xy =>u2_xy
          slices%xy2=>u2_xy2
          slices%ready=.true.
!
!  Vorticity.
!
        case ('oo')
          if (slices%index>=3) then
            slices%ready=.false.
          else
            slices%index=slices%index+1
            slices%yz =>oo_yz(:,:,slices%index)
            slices%xz =>oo_xz(:,:,slices%index)
            slices%xy =>oo_xy(:,:,slices%index)
            slices%xy2=>oo_xy2(:,:,slices%index)
            if (slices%index<=3) slices%ready=.true.
          endif
!
!  Vorticity squared.
!
        case ('o2')
          slices%yz =>o2_yz
          slices%xz =>o2_xz
          slices%xy =>o2_xy
          slices%xy2=>o2_xy2
          slices%ready=.true.
!
!  Mach number.
!
        case ('mach')
          slices%yz =>mach_yz
          slices%xz =>mach_xz
          slices%xy =>mach_xy
          slices%xy2=>mach_xy2
          slices%ready=.true.
!
      endselect
!
    endsubroutine get_slices_hydro
!***********************************************************************
endmodule Hydro
