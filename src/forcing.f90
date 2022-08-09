! $Id$
!
!  This module contains routines both for delta-correlated
!  and continuous forcing. The fcont pencil is only provided
!  for continuous forcing.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lforcing = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED fcont(3,n_forcing_cont_max)
!
!***************************************************************
module Forcing
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Geometrical_types
!
  implicit none
!
  include 'forcing.h'
  include 'record_types.h'
!
  interface input_persistent_forcing
     module procedure input_persist_forcing_id
     module procedure input_persist_forcing
  endinterface
!
  real :: force=0.,force2=0., force_double=0., force1_scl=1., force2_scl=1.
  real :: relhel=1., height_ff=0., r_ff=0., r_ff_hel=0., rcyl_ff=0., rel_zcomp=1.
  real :: Bconst=1., Bslope=0.
  real :: fountain=1.,width_ff=.5,nexp_ff=1.,n_hel_sin_pow=0.
  real :: crosshel=0.
  real :: radius_ff=0., k1_ff=1., kx_ff=1., ky_ff=1., kz_ff=1., z_center=0.
  real :: slope_ff=0., work_ff=0., omega_ff=0., omega_double_ff=0., n_equator_ff=1.
  real :: tforce_stop=impossible,tforce_stop2=impossible
  real :: tforce_start=0.,tforce_start2=0.
  real :: wff_ampl=0.,  xff_ampl=0.,  yff_ampl=0.,  zff_ampl=0.
  real :: wff2_ampl=0., xff2_ampl=0., yff2_ampl=0., zff2_ampl=0.
  real :: zff_hel=0.,max_force=impossible
  real :: dtforce=0., dtforce_ampl=.5, dtforce_duration=-1.0, force_strength=0.
  double precision :: tforce_ramp_down=1.1, tauforce_ramp_down=1.
  real, dimension(3) :: force_direction=(/0.,0.,0./)
  real, dimension(3,2) :: location_fixed=0.
  real, dimension(nx) :: profx_ampl=1.,profx_hel=1., profx_ampl1=0.
  real, dimension(my) :: profy_ampl=1.,profy_hel=1.
  real, dimension(mz) :: profz_ampl=1.,profz_hel=1.,qdouble_profile=1.
  integer :: kfountain=5,iff,ifx,ify,ifz,ifff,iffx,iffy,iffz,i2fff,i2ffx,i2ffy,i2ffz
  integer :: kzlarge=1
  integer :: iforcing_zsym=0, nlocation=1
  logical :: lwork_ff=.false.,lmomentum_ff=.false.
  logical :: lhydro_forcing=.true., lneutral_forcing=.false., lmagnetic_forcing=.false.
  logical :: lcrosshel_forcing=.false.,ltestfield_forcing=.false.,ltestflow_forcing=.false.
  logical :: lxxcorr_forcing=.false., lxycorr_forcing=.false.
  logical :: lhelical_test=.false., lrandom_location=.true., lrandom_time=.false.
  logical :: lwrite_psi=.false., lforce_ramp_down=.false.
  logical :: lscale_kvector_tobox=.false.,lwrite_gausspot_to_file=.false.
  logical :: lwrite_gausspot_to_file_always=.false.
  logical :: lscale_kvector_fac=.false.
  logical :: lforce_peri=.false., lforce_cuty=.false.
  logical :: lforcing2_same=.false., lforcing2_curl=.false.
  logical :: lff_as_aux = .false.
  real :: scale_kvectorx=1.,scale_kvectory=1.,scale_kvectorz=1.
  logical :: old_forcing_evector=.false., lforcing_coefs_hel_double=.false.
  character (len=labellen) :: iforce='zero', iforce2='zero'
  character (len=labellen) :: iforce_profile='nothing'
  character (len=labellen) :: iforce_tprofile='nothing'
  real :: equator=0.
  real :: kx_2df=0.,ky_2df=0.,xminf=0.,xmaxf=0.,yminf=0.,ymaxf=0.
! For helical forcing in spherical polar coordinate system
  real, allocatable, dimension(:,:,:) :: psif
  real, allocatable, dimension(:,:) :: cklist
  logical :: lfastCK=.false.,lsamesign=.true.
! allocated only if we have lfastCK=T
  real,allocatable,dimension(:,:,:) :: Zpsi_list
  real,allocatable,dimension(:,:,:) :: RYlm_list,IYlm_list
  integer :: helsign=0,nlist_ck=25
  real :: fpre = 1.0,ck_equator_gap=0.,ck_gap_step=0.
  integer :: icklist,jtest_aa0=5,jtest_uu0=1
! For random forcing
  logical :: lavoid_xymean=.false., lavoid_ymean=.false., lavoid_zmean=.false., &
             ldiscrete_phases=.false.
! For random forcing in 2d
  integer,allocatable, dimension (:,:) :: random2d_kmodes
  integer :: random2d_nmodes
  integer :: random2d_kmin=0, random2d_kmax=0
  integer :: channel_force=1
! continuous 2d forcing
  integer :: k2d
  logical :: l2dxz,l2dyz
! anisotropic forcing
  logical :: laniso_forcing_old=.true.
! Persistent stuff
  real :: tsforce=-10.
  real, dimension (3) :: location, location2
!
!  continuous forcing variables
!
  integer :: n_forcing_cont=n_forcing_cont_max
  logical :: lembed=.false.,lshearing_adjust_old=.false.
  logical, dimension(n_forcing_cont_max) :: lgentle=.false.
  character (len=labellen), dimension(n_forcing_cont_max) :: iforcing_cont='nothing'
  real, dimension(n_forcing_cont_max) :: ampl_ff=1., ampl1_ff=0., width_fcont=1., x1_fcont=0., x2_fcont=0.
  real, dimension(n_forcing_cont_max) :: kf_fcont=impossible, kf_fcont_x=impossible, kf_fcont_y=impossible, kf_fcont_z=impossible
  real, dimension(n_forcing_cont_max) :: omega_fcont=0., omegay_fcont=0., omegaz_fcont=0.
  real, dimension(n_forcing_cont_max) :: eps_fcont=0., tgentle=0., z_center_fcont=0.
  real, dimension(n_forcing_cont_max) :: ampl_bb=5.0e-2,width_bb=0.1,z_bb=0.1,eta_bb=1.0e-4
  real, dimension(n_forcing_cont_max) :: fcont_ampl=1., ABC_A=1., ABC_B=1., ABC_C=1.
  real :: ampl_diffrot=1.0,omega_exponent=1.0
  real :: omega_tidal=1.0, R0_tidal=1.0, phi_tidal=1.0, Omega_vortex=0.
  real :: cs0eff=impossible
  type(torus_rect), save :: torus
  ! GP_TC13 forcing
  real :: tcor_GP=1.,kmin_GP=1.,kmax_GP=2.,beta_GP=1.3333
  integer :: nk_GP=2
!
!  auxiliary functions for continuous forcing function
!
  real, dimension (my,n_forcing_cont_max) :: phi1_ff
  real, dimension (mx,n_forcing_cont_max) :: phi2_ff
  real, dimension (mx,n_forcing_cont_max) :: sinx,cosx,sinxt,cosxt,embedx
  real, dimension (my,n_forcing_cont_max) :: siny,cosy,sinyt,cosyt,embedy
  real, dimension (mz,n_forcing_cont_max) :: sinz,cosz,sinzt,coszt,embedz
  real, dimension (100,n_forcing_cont_max) :: xi_GP,eta_GP
!
  namelist /forcing_run_pars/ &
       tforce_start,tforce_start2,&
       iforce,force,force_double,relhel,crosshel,height_ff,r_ff,r_ff_hel, &
       rcyl_ff,width_ff,nexp_ff,lff_as_aux,Bconst,Bslope, rel_zcomp, &
       iforce2, force2, force1_scl, force2_scl, iforcing_zsym, &
       kfountain, fountain, tforce_stop, tforce_stop2, &
       radius_ff,k1_ff,kx_ff,ky_ff,kz_ff,slope_ff,work_ff,lmomentum_ff, &
       omega_ff, omega_double_ff, n_equator_ff, location_fixed, lrandom_location, nlocation, &
       lwrite_gausspot_to_file,lwrite_gausspot_to_file_always, &
       wff_ampl, xff_ampl, yff_ampl, zff_ampl, zff_hel, &
       wff2_ampl, xff2_ampl,yff2_ampl, zff2_ampl, &
       lhydro_forcing, lneutral_forcing, lmagnetic_forcing, &
       ltestfield_forcing, ltestflow_forcing, &
       lcrosshel_forcing, lxxcorr_forcing, lxycorr_forcing, &
       jtest_aa0,jtest_uu0, &
       max_force,dtforce,dtforce_duration,old_forcing_evector, &
       lforcing_coefs_hel_double, dtforce_ampl, &
       iforce_profile, iforce_tprofile, lscale_kvector_tobox, &
       force_direction, force_strength, lhelical_test, &
       lfastCK,fpre,helsign,nlist_ck,lwrite_psi,&
       ck_equator_gap,ck_gap_step,&
       ABC_A, ABC_B, ABC_C, &
       lforcing_cont,iforcing_cont, z_center_fcont, z_center, &
       lembed, k1_ff, ampl_ff, ampl1_ff, width_fcont, x1_fcont, x2_fcont, &
       kf_fcont, omega_fcont, omegay_fcont, omegaz_fcont, eps_fcont, &
       lsamesign, lshearing_adjust_old, equator, &
       lscale_kvector_fac,scale_kvectorx,scale_kvectory,scale_kvectorz, &
       lforce_peri,lforce_cuty,lforcing2_same,lforcing2_curl, &
       tgentle,random2d_kmin,random2d_kmax,l2dxz,l2dyz,k2d, &
       z_bb,width_bb,eta_bb,fcont_ampl, &
       ampl_diffrot,omega_exponent,kx_2df,ky_2df,xminf,xmaxf,yminf,ymaxf, &
       lavoid_xymean,lavoid_ymean,lavoid_zmean, ldiscrete_phases, &
       omega_tidal, R0_tidal, phi_tidal, omega_vortex, &
       lforce_ramp_down, tforce_ramp_down, tauforce_ramp_down, &
       n_hel_sin_pow, kzlarge, cs0eff, channel_force, torus, Omega_vortex, &
       lrandom_time, laniso_forcing_old, tcor_GP, kmin_GP,kmax_GP,nk_GP,beta_GP
!
! other variables (needs to be consistent with reset list below)
!
  integer :: idiag_rufm=0, idiag_rufint=0, idiag_ufm=0, idiag_ofm=0
  integer :: idiag_qfm=0, idiag_ffm=0
  integer :: idiag_ruxfxm=0, idiag_ruyfym=0, idiag_ruzfzm=0
  integer :: idiag_ruxfym=0, idiag_ruyfxm=0
  integer :: idiag_fxbxm=0, idiag_fxbym=0, idiag_fxbzm=0
!
! Auxiliaries
!
  real, dimension(:,:), pointer :: reference_state
  real, dimension(3) :: k1xyz=0.
  real, dimension(mz) :: profz_k
  logical :: lmhd_forcing
!
! Data from k.dat.
!
  integer :: nk, nk2
  real :: kav, kav2
  real, dimension(:), allocatable :: kkx,kky,kkz
  real, dimension(:), allocatable :: kkx2,kky2,kkz2
!
  real, allocatable, dimension (:,:) :: KS_k,KS_A,KS_B !or through whole field for each wavenumber?
  real, allocatable, dimension (:) :: KS_omega !or through whole field for each wavenumber?
  integer :: KS_modes = 25
!
  contains
!
!***********************************************************************
    subroutine register_forcing
!
!  add forcing in timestep
!  11-may-2002/wolf: coded
!
!  identify version number
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_forcing
!***********************************************************************
    subroutine initialize_forcing
!
!  read seed field parameters
!  nothing done from start.f90
!
!  25-sep-2014/MR: determine n_forcing_cont according to the actual selection
!  23-jan-2015/MR: reference state now fetched here and stored in module variable
!  15-feb-2015/MR: returning before entering continuous forcing section when
!                  no such forcing is requested
!  18-dec-2015/MR: minimal wavevectors k1xyz moved here from grid
!  14-Jun-2016/MR+NS: added forcing sinx*exp(-z^2)
!  11-May-2017/NS: added forcing Aycont_z
!  08-Aug-2019/MR: moved reading of k.dat, k_double.dat from individual subroutines
!                  -> nk, kav, kk[xyz] and nk2, kav2, kk2[xyz], respectively, now module 
!                  variables.
!
      use General, only: bessj,itoa
      use Mpicomm, only: stop_it
      use SharedVariables, only: get_shared_variable
      use Sub, only: step,erfunc,stepdown,register_report_aux
      use EquationOfState, only: cs0
!
      real :: zstar,rmin,rmax,a_ell,anum,adenom,jlm_ff,ylm_ff,alphar,Balpha,RYlm,IYlm
      integer :: l,m,n,i,ilread,ilm,ckno,ilist,emm,aindex,Legendrel
      logical :: lk_dot_dat_exists
!
      if (lstart) then
        if (ip<4) print*,'initialize_forcing: not needed in start'
      else
!
!  Check whether we want ordinary hydro forcing or magnetic forcing
!  Currently there is also the option of forcing just the neutrals,
!  but in future it could be combined with other forcings instead.
!
        if (lcrosshel_forcing) then
          if (.not.(lmagnetic.and.lhydro)) &
            call fatal_error('initialize_forcing','cross helicity forcing requires an MHD run')
          lmagnetic_forcing=.true.
          lhydro_forcing=.true.
          lmhd_forcing=.false.
        else
          lmhd_forcing=lmagnetic_forcing.and.lhydro_forcing .or. &
                       ltestfield_forcing.and.ltestflow_forcing
        endif 
        if (lmagnetic_forcing) then
          ifff=iaa; iffx=iax; iffy=iay; iffz=iaz
        endif 
        if (lneutral_forcing) then
          if (iuun==0) call fatal_error("lneutral_forcing","uun==0")
          ifff=iuun; iffx=iunx; iffy=iuny; iffz=iunz
        endif 
        if (lhydro_forcing) then
          if (lmagnetic_forcing) then
            i2fff=iuu; i2ffx=iux; i2ffy=iuy; i2ffz=iuz
          else
            if (iphiuu==0) then
              ifff=iuu; iffx=iux; iffy=iuy; iffz=iuz
            else
              ifff=iphiuu; iffx=iphiuu; iffy=0; iffz=0
            endif
          endif
        endif
        if (.not.(lmagnetic_forcing.or.lhydro_forcing.or.lneutral_forcing.or. &
                  ltestfield_forcing.or.ltestflow_forcing)) &
          call stop_it("initialize_forcing: No forcing function set")
        
        if (ldebug) print*,'initialize_forcing: ifff=',ifff
!
!  check whether we want constant forcing at each timestep,
!  in which case lwork_ff is set to true.
!
        if (work_ff/=0.) then
          force=1.
          lwork_ff=.true.
          if (lroot) print*,'initialize_forcing: reset force=1., because work_ff is set'
        endif
      endif
!
!  initialize location to location_fixed
!
      location=location_fixed(:,1)
      location2=location_fixed(:,2)
!
!  vertical profiles for amplitude and helicity of the forcing
!  default is constant profiles for rms velocity and helicity.
!
      if (iforce_profile=='nothing') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        profz_ampl=1.; profz_hel=1.
!
!  sine profile for amplitude of forcing 
!
      elseif (iforce_profile=='sinz') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        profz_hel=1.
        do n=1,mz
          profz_ampl(n)=sin(kzlarge*z(n))
        enddo
!
      elseif (iforce_profile=='equator') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        profz_ampl=1.
        do n=1,mz
          profz_hel(n)=sin(2*pi*z(n)*n_equator_ff/Lz)
        enddo
!
      elseif (iforce_profile=='equator_step') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        profz_ampl=1.
        do n=1,mz
          profz_hel(n)= -1.+2.*step(z(n),equator-ck_equator_gap,ck_gap_step)
        enddo
!
!  step function change in intensity of helicity at zff_hel
!
      elseif (iforce_profile=='step_ampl=z') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.; profy_hel=1.
                       profz_hel=1.
        do n=1,mz
          profz_ampl(n)=step(z(n),zff_hel,width_ff)
        enddo
!
!  sign change of helicity proportional to cosy
!
      elseif (iforce_profile=='equator_hel=cosy') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.
        do m=1,my
          profy_hel(m)=cos(y(m))
        enddo
        profz_ampl=1.; profz_hel=1.
!
! step function profile
!
      elseif (iforce_profile=='equator_hel=step') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        profz_ampl=1.
        do m=1,my
          profy_hel(m)= -1.+2.*step(y(m),yequator-ck_equator_gap,ck_gap_step)
        enddo
!
!  sign change of helicity about z=0
!
      elseif (iforce_profile=='equator_hel=z') then
        call fatal_error("initialize_forcing","use equator_hel=z/L instead")
!
!  Linear profile of helicity, normalized by ztop (or -zbot, if it is bigger)
!
      elseif (iforce_profile=='equator_hel=z/L') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        profz_ampl=1.
        do n=1,mz
          profz_hel(n)=z(n)/max(-xyz0(3),xyz1(3))
        enddo
!
!  cosine profile of helicity about z=0
!
      elseif (iforce_profile=='cos(z)') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        profz_ampl=1.
        do n=1,mz
          profz_hel(n)=cos(z(n))
        enddo
!
!  cosine profile of helicity about z=0
!
      elseif (iforce_profile=='cos(z/2)') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        profz_ampl=1.
        do n=1,mz
          profz_hel(n)=cos(.5*z(n))
        enddo
!
!  Cosine profile of helicity about z=0, with max function.
!  Should be used only if |z| < 1.5*pi.
!
      elseif (iforce_profile=='cos(z/2)_with_halo') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        profz_ampl=1.
        do n=1,mz
          profz_hel(n)=max(cos(.5*z(n)),0.)
        enddo
!
!  helicity profile proportional to z^2, but vanishing on the boundaries
!
      elseif (iforce_profile=='squared') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        profz_ampl=1.
        do n=1,mz
          profz_hel(n)=(z(n)-xyz0(3))*(xyz1(3)-z(n))
        enddo
!
!  helicity profile proportional to z^3, but vanishing on the boundaries
!
      elseif (iforce_profile=='cubic') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        profz_ampl=1.
        do n=1,mz
          profz_hel(n)=z(n)*(z(n)-xyz0(3))*(xyz1(3)-z(n))
        enddo
!
!  power law
!
      elseif (iforce_profile=='(z-z0)^n') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        profz_ampl=(abs(z-xyz0(3))/height_ff)**nexp_ff
        profz_hel=1.
        if (lroot) print*,'profz_ampl=',profz_ampl
!
!  power law with offset
!
      elseif (iforce_profile=='1+(z-z0)^n') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        if (height_ff>0) then
          zstar=xyz0(3)
        elseif (height_ff<0) then
          zstar=xyz1(3)
        else
          zstar=impossible
          call fatal_error('must have height_ff/=0','forcing')
        endif
        profz_ampl=(1.+(z-zstar)/height_ff)**nexp_ff
        profz_hel=1.
        if (lroot) print*,'profz_ampl=',profz_ampl
!
!  power law with offset
!
      elseif (iforce_profile=='1+(z-z0)^n_prof') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        if (height_ff>0) then
          zstar=xyz0(3)
        elseif (height_ff<0) then
          zstar=xyz1(3)
        else
          zstar=impossible
          call fatal_error('must have height_ff/=0','forcing')
        endif
        profz_ampl=((1.+(z-zstar)/height_ff)**nexp_ff)*.5*(1.+tanh(20.*cos(.55*z)))
        profz_hel=1.
        if (lroot) print*,'profz_ampl=',profz_ampl
!
!  exponential law
!
      elseif (iforce_profile=='exp(z/H)') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        profz_ampl=exp(z/width_ff)
        profz_hel=1.
!
      elseif (iforce_profile=='exp(-(z/H)^2)') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        profz_ampl=exp(-(z/width_ff)**2)
        profz_hel=1.
!
      elseif (iforce_profile=='xybox') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        profz_ampl=.5*(1.-erfunc(z/width_ff))
        profz_hel=1.
        profx_ampl= step(x(l1:l2),xminf,width_ff)-step(x(l1:l2),xmaxf,width_ff)
        do m=1,my
          profy_ampl(m)= step(y(m),yminf,width_ff)-step(y(m),ymaxf,width_ff)
        enddo
!
!  turn off forcing intensity above z=0
!
      elseif (iforce_profile=='surface_z') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        profz_ampl=.5*(1.-erfunc((z-r_ff)/width_ff))
        profz_hel=1.
!
!  turn off helicity of forcing above x=r_ff
!
      elseif (iforce_profile=='surface_helx') then
        profx_ampl=1.; profx_hel=.5*(1.-erfunc((x(l1:l2)-r_ff)/width_ff))
        profy_ampl=1.; profy_hel=1.
        profz_ampl=1.; profz_hel=1.
!
!  turn off helicity of forcing above x=r_ff
!
      elseif (iforce_profile=='surface_helx_cosy') then
        profx_ampl=1.; profx_hel=.5*(1.-erfunc((x(l1:l2)-r_ff)/width_ff))
       profy_ampl=1.; profy_hel=cos(y)
       profz_ampl=1.; profz_hel=1.
!
!  turn off helicity of forcing above x=r_ff
!  used in Jabbari et al. (2015)
!
      elseif (iforce_profile=='surface_helx_cosy*siny**n_hel_sin_pow') then
        profx_ampl=1.; profx_hel=.5*(1.-erfunc((x(l1:l2)-r_ff)/width_ff))
        profy_ampl=1.; profy_hel=cos(y)*sin(y)**n_hel_sin_pow
        profz_ampl=1.; profz_hel=1.
!
!  turn off helicity of forcing above x=r_ff
!  but with step function in radius, as in Warnecke et al. (2011)
!
      elseif (iforce_profile=='surface_stepx_cosy*siny**n_hel_sin_pow') then
        profx_ampl=.5*(1.-erfunc((x(l1:l2)-r_ff)/width_ff))
        profx_hel=1.
        profy_ampl=1.; profy_hel=cos(y)*sin(y)**n_hel_sin_pow
        profz_ampl=1.; profz_hel=1.
!
!  turn off helicity of forcing above z=r_ff
!
      elseif (iforce_profile=='surface_helz') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        profz_ampl=1.
        profz_hel=.5*(1.-erfunc((z-r_ff)/width_ff))
!
!  turn off helicity of forcing above z=0
!
      elseif (iforce_profile=='surface_amplhelz') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        profz_ampl=.5*(1.-erfunc((z-r_ff)/width_ff))
        profz_hel=.5*(1.-erfunc((z-r_ff_hel)/width_ff))
!
!  turn off forcing intensity above z=z0, and
!  stepx profile of helicity
!
      elseif (iforce_profile=='surface_z_stepx') then
        profx_ampl=1.
        profy_ampl=1.; profy_hel=1.
        do l=l1,l2
          profx_hel(l-nghost)= -1.+2.*step(x(l),0.,width_ff)
        enddo
        profz_ampl=.5*(1.-erfunc((z-r_ff)/width_ff))
        profz_hel=1.
!
!  turn off forcing intensity above z=z0, and
!  stepy profile of helicity
!
      elseif (iforce_profile=='surface_z_stepy') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.
        do m=1,my
          profy_hel(m)= -1.+2.*step(y(m),0.,width_ff)
        enddo
        profz_ampl=.5*(1.-erfunc((z-r_ff)/width_ff))
        profz_hel=1.
!
! turn on forcing in a certain layer (generalization of 'forced_convection')
!
      elseif (iforce_profile=='xyzbox') then
        profx_ampl=.5*(1.+erfunc((x(l1:l2)-xff_ampl)/wff_ampl)) &
                  -.5*(1.+erfunc((x(l1:l2)-xff2_ampl)/wff2_ampl))
        profy_ampl=.5*(1.+erfunc((y-yff_ampl)/wff_ampl)) &
                  -.5*(1.+erfunc((y-yff2_ampl)/wff2_ampl))
        profz_ampl=.5*(1.+erfunc((z-zff_ampl)/wff_ampl)) &
                  -.5*(1.+erfunc((z-zff2_ampl)/wff2_ampl))
        profx_hel=1.
        profy_hel=1.
        profz_hel=1.
!
! turn on forcing in a certain layer (generalization of 'forced_convection')
!
      elseif (iforce_profile=='zlayer') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        profz_ampl=.5*(1.+erfunc((z-zff_ampl)/wff_ampl)) &
                  -.5*(1.+erfunc((z-zff2_ampl)/wff2_ampl))
!        profz_hel=1.
        do n=1,mz
          profz_hel(n)=z(n)/max(-xyz0(3),xyz1(3))
        enddo
!
! turn on forcing in the bulk of the convection zone
!
      elseif (iforce_profile=='forced_convection') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        profz_ampl=.5*(1.+ erfunc((z)/width_ff)) - 0.5*(1.+ erfunc((z-1.)/(2.*width_ff)))
        profz_hel=1.
!
!  cosx modulation
!
      elseif (iforce_profile=='cosx2') then
        print*,'--------------------'
        print*,'using iforce profile ', iforce_profile
        print*,'--------------------'
        profx_ampl=cos(kx_ff*x(l1:l2))**2
        profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        profz_ampl=1.; profz_hel=1.
!
!  cosx modulation
!
      elseif (iforce_profile=='cosx^nexp') then
        profx_ampl=fountain*cos(kx_ff*x(l1:l2))**nexp_ff
        profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        profz_ampl=1.; profz_hel=1.
!
!  turn off forcing intensity above x=x0, and
!  cosy profile of helicity
!
      elseif (iforce_profile=='surface_x_cosy') then
        profx_ampl=.5*(1.-erfunc((x(l1:l2)-r_ff)/width_ff))
        profx_hel=1.
        profy_ampl=1.
        do m=1,my
          profy_hel(m)=cos(y(m))
        enddo
        profz_ampl=1.; profz_hel=1.
!
!  turn off forcing intensity above x=x0, and
!  stepy profile of helicity, used in Warnecke et al. (2011)
!
      elseif (iforce_profile=='surface_x_stepy') then
        profx_ampl=.5*(1.-erfunc((x(l1:l2)-r_ff)/width_ff))
        profx_hel=1.
        profy_ampl=1.
        do m=1,my
          profy_hel(m)= -1.+2.*step(y(m),pi/2.,width_ff)
        enddo
        profz_ampl=1.; profz_hel=1.
!
!  turn off forcing intensity above x=x0
!
      elseif (iforce_profile=='surface_x') then
        profx_ampl=.5*(1.-erfunc((x(l1:l2)-r_ff)/width_ff))
        profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        profz_ampl=1.; profz_hel=1.
!
!  turn off forcing intensity above x=r_ff and y=r_ff
!
      elseif (iforce_profile=='surface_xy') then
        profx_ampl=.5*(1.-erfunc((x(l1:l2)-r_ff)/width_ff)); profx_hel=1.
        profy_ampl=.5*(1.-erfunc((y-r_ff)/width_ff)); profy_hel=1.
        profz_ampl=1.; profz_hel=1.
!
!  turn off forcing intensity above x=r_ff and y=r_ff
!
      elseif (iforce_profile=='surface_xz') then
        profx_ampl=.5*(1.-erfunc((x(l1:l2)-r_ff)/width_ff)); profx_hel=1.
        profz_ampl=.5*(1.-erfunc((z-r_ff)/width_ff)); profz_hel=1.
        profy_ampl=1.; profy_hel=1.
!
!  turn on forcing intensity above x=r_ff
!
      elseif (iforce_profile=='above_x0') then
        profx_ampl=.5*(1.+erfunc((x(l1:l2)-r_ff)/width_ff))
        profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        profz_ampl=1.; profz_hel=1.
!
!  turn on forcing intensity above x=x0 with cosy profile
!
      elseif (iforce_profile=='above_x0_cosy') then
        profx_ampl=.5*(1.+erfunc((x(l1:l2)-r_ff)/width_ff))
        profy_ampl=1.
        profx_hel=1.
        do m=1,my
          profy_hel(m)=cos(y(m));
        enddo
        profz_ampl=1.; profz_hel=1.
!
!  just a change in intensity in the z direction
!
      elseif (iforce_profile=='intensity') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        profz_hel=1.
        do n=1,mz
          profz_ampl(n)=.5+.5*cos(z(n))
        enddo
!
!  Galactic profile both for intensity and helicity
!
      elseif (iforce_profile=='galactic-old') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        do n=1,mz
          if (abs(z(n))<zff_ampl) profz_ampl(n)=.5*(1.-cos(z(n)))
          if (abs(z(n))<zff_hel ) profz_hel (n)=.5*(1.+cos(z(n)/2.))
        enddo
!
!  Galactic profile both for intensity and helicity
!  Normally one would put zff_ampl=zff_hel=pi
!
      elseif (iforce_profile=='galactic') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        do n=1,mz
          profz_ampl(n)=.0
          profz_hel (n)=.0
          if (abs(z(n))<zff_ampl) profz_ampl(n)=.5*(1.+cos(z(n)))
          if (abs(z(n))<zff_hel ) profz_hel (n)=sin(z(n))
        enddo
!
! Galactic helicity profile for helicity
      elseif (iforce_profile=='gal-wind') then
         profx_ampl=1.; profx_hel=1.
         profy_ampl=1.; profy_hel=1.
         do n=1,mz
            profz_hel(n)=sin(pi*z(n)/Lz)
         enddo
!
      elseif (iforce_profile=='gal-wind2') then
         profx_ampl=1.; profx_hel=1.
         profy_ampl=1.; profy_hel=1.
         do n=1,mz
            profz_hel(n)=sin(2.*pi*z(n)/Lz)
         enddo
!
!  corona
!
      elseif (iforce_profile=='diffrot_corona') then
        profz_ampl=1.; profz_hel=1.
        profy_ampl=1.; profy_hel=1.
        profx_ampl=.5*(1.-tanh((x(l1:l2)-xff_ampl)/wff_ampl))
        profx_hel=1.
!
      else
        call fatal_error('initialize_forcing','iforce_profile value does not exist')
      endif
!
!  Turn on forcing intensity for force_double above z=0.
!  Adding a z-profile is the default. To turn it off, put width_ff=0.
!  If it is turned off, we put qdouble_profile=.5, to give both 50%,
!  if force_double=force.
!
      if (lforcing_coefs_hel_double) then
        if (width_ff==0.) then
          qdouble_profile=.5
        else
          qdouble_profile=.5*(1.+erfunc((z-r_ff)/width_ff))
        endif
!
!  Debug output
!
        if (ip<17.and.lroot) then
          print*,'r_ff,width_ff=',r_ff,width_ff
          print*,'qdouble_profile=',qdouble_profile
        endif
      endif
!
!  at the first step, the sin and cos functions are calculated for all
!  x,y,z points and are then saved and used for all subsequent steps
!  and pencils
!
      if (ip<=6) print*,'forcing_cont:','lforcing_cont=',lforcing_cont,iforcing_cont

      if (lstart) return

      if (iforce=='2'           .or.iforce=='helical'     .or.iforce=='helical_kprof'.or. &
          iforce=='helical_both'.or.iforce=='irrotational'.or.iforce=='hel_smooth'   .or. &
          iforce=='noshear') then

        inquire(FILE="k.dat", EXIST=lk_dot_dat_exists)

        if (lk_dot_dat_exists) then
          if (lroot.and.ip<14) print*,'initialize_forcing: opening k.dat'
          open(9,file='k.dat',status='old')
          read(9,*) nk,kav
          if (lroot.and.ip<14) print*,'initialize_forcing: average k=',kav
          if (allocated(kkx)) deallocate(kkx,kky,kkz)
          allocate(kkx(nk),kky(nk),kkz(nk))
          read(9,*) kkx
          read(9,*) kky
          read(9,*) kkz
          close(9)
        else
          call inevitably_fatal_error('initialize_forcing', 'you must give an input k.dat file')
        endif

        if (lforcing_coefs_hel_double) then
!
!  Read k_double.dat once at first call.
!
          inquire(FILE="k_double.dat", EXIST=lk_dot_dat_exists)
          if (lk_dot_dat_exists) then
            if (lroot.and.ip<14) print*,'initialize_forcing: opening k_double.dat'
            open(9,file='k_double.dat',status='old')
            read(9,*) nk2,kav2
            if (lroot.and.ip<14) print*,'initialize_forcing: average k(2)=',kav2
            if (allocated(kkx2)) deallocate(kkx2,kky2,kkz2)
            allocate(kkx2(nk2),kky2(nk2),kkz2(nk2))
            read(9,*) kkx2
            read(9,*) kky2
            read(9,*) kkz2
            close(9)
          else
            call inevitably_fatal_error ('initialize_forcing', &
                                         'you must give an input k_double.dat file')
          endif
        endif
!
      elseif (iforce=='chandra_kendall'.or.iforce=='cktest') then
!
        if (.not. lspherical_coords) call fatal_error('initialize_forcing', &
                        'Chandrasekhar-Kendall forcing works only in spherical coordinates!')
!
        if (abs(helsign)/=1) &
          call fatal_error("initialize_forcing", "CK forcing: helsign must be +1 or -1")
!
        if (iforce=='cktest') then
          lhelical_test=.true.
          icklist=0
        endif

        if (lroot) print*,'Chandrasekhar-Kendall forcing in spherical polar coordinates'
        if (.not.allocated(psif)) allocate(psif(mx,my,mz),cklist(nlist_ck,5))
!
! Read the list of values for emm, ell and alpha from file "alpha_in.dat". This code is designed for 25 such.
!
        open(unit=76,file="alpha_in.dat",status="old")
        read(76,*) ckno,rmin,rmax
        if (ckno/=nlist_ck) &
          call fatal_error("initialize_forcing", &
               "Number of entries in alpha_in.dat "//trim(itoa(ckno))//" unequal nlist_ck="//trim(itoa(nlist_ck))//".")

        do ilread=1,nlist_ck
          read(76,*) (cklist(ilread,ilm),ilm=1,5)
        enddo
        close(76)

        if (ck_equator_gap/=0) &
          profy_ampl=1.-step(y,pi/2.-ck_equator_gap,ck_gap_step)+step(y,pi/2.+ck_equator_gap,ck_gap_step)

        if (lfastCK) then

          if (.not.allocated(Zpsi_list)) then
            allocate(Zpsi_list(mx,nlist_ck,3))
            allocate(RYlm_list(my,mz,nlist_ck),IYlm_list(my,mz,nlist_ck))
          endif

          do ilist=1,nlist_ck
            emm = cklist(ilist,1)
            Legendrel = cklist(ilist,2)
            do n=1,mz
              do m=1,my
                call sp_harm_real(RYlm,Legendrel,emm,y(m),z(n))
                call sp_harm_imag(IYlm,Legendrel,emm,y(m),z(n))
                RYlm_list(m,n,ilist)=RYlm
                IYlm_list(m,n,ilist)=IYlm
              enddo
            enddo

            do aindex=1,3
              Balpha = cklist(ilist,2+aindex)
              call sp_bessely_l(anum,Legendrel,Balpha*x(l1))
              call sp_besselj_l(adenom,Legendrel,Balpha*x(l1))
              a_ell = -anum/adenom
              do l=1,mx
                alphar=Balpha*x(l)
                call sp_besselj_l(jlm_ff,Legendrel,alphar)
                call sp_bessely_l(ylm_ff,Legendrel,alphar)
                Zpsi_list(l,ilist,aindex) = (a_ell*jlm_ff+ylm_ff)
              enddo
            enddo
          enddo

        endif 
      endif

      if (lff_as_aux) call register_report_aux('ff', iff, ifx, ify, ifz)
!
!  Get reference_state. Requires that density is initialized before forcing.
!
      if (lreference_state) &
        call get_shared_variable('reference_state',reference_state,caller='initialize_forcing')
!
!  Make sure that testfields/testflows are registered.  
!    
      ltestfield_forcing = ltestfield_forcing.and.iaatest>0
      ltestflow_forcing = ltestflow_forcing.and.iuutest>0
!
!  At the moment, cs0 is used for normalization.
!  It is saved in param2.nml.
!
      if (cs0eff==impossible) then
        if (cs0==impossible) then
          cs0eff=1.
          if (headt) print*,'initialize_forcing: for normalization, use cs0eff=',cs0eff
        else
          cs0eff=cs0
        endif
      endif

      if (.not.lforcing_cont) return
      
      do i=1,n_forcing_cont_max
        if ( iforcing_cont(i)=='nothing' ) then
          n_forcing_cont=i-1
          exit
        else
          if (kf_fcont(i)  ==impossible) kf_fcont(i)  =k1_ff
          if (kf_fcont_x(i)==impossible) kf_fcont_x(i)=kx_ff
          if (kf_fcont_y(i)==impossible) kf_fcont_y(i)=ky_ff
          if (kf_fcont_z(i)==impossible) kf_fcont_z(i)=kz_ff
          if (z_center_fcont(i)==0.) z_center_fcont(i)=z_center
        endif
!
!  all 3 sin and cos functions needed
!
        if (iforcing_cont(i)=='ABC' .or. &
            iforcing_cont(i)=='Straining') then
          if (lroot) print*,'forcing_cont: ABC--calc sinx, cosx, etc'
          sinx(:,i)=sin(kf_fcont(i)*x); cosx(:,i)=cos(kf_fcont(i)*x)
          siny(:,i)=sin(kf_fcont(i)*y); cosy(:,i)=cos(kf_fcont(i)*y)
          sinz(:,i)=sin(kf_fcont(i)*z); cosz(:,i)=cos(kf_fcont(i)*z)
        elseif (iforcing_cont(i)=='RobertsFlow' .or. &
                iforcing_cont(i)=='RobertsFlowII' .or. &
                iforcing_cont(i)=='RobertsFlowMask' .or. &
                iforcing_cont(i)=='RobertsFlow_exact' ) then
          if (lroot) print*,'forcing_cont: Roberts Flow'
          sinx(:,i)=sin(kf_fcont(i)*x); cosx(:,i)=cos(kf_fcont(i)*x)
          siny(:,i)=sin(kf_fcont(i)*y); cosy(:,i)=cos(kf_fcont(i)*y)
          if (iforcing_cont(i)=='RobertsFlowMask') then
            profx_ampl=exp(-.5*x(l1:l2)**2/radius_ff**2)
            profy_ampl=exp(-.5*y**2/radius_ff**2)
          endif
        elseif (iforcing_cont(i)=='RobertsFlow-zdep'.or.iforcing_cont(i)=='Roberts-for-SSD') then
          if (lroot) print*,'forcing_cont: z-dependent Roberts Flow'
          sinx(:,i)=sin(kf_fcont(i)*x); cosx(:,i)=cos(kf_fcont(i)*x)
          siny(:,i)=sin(kf_fcont(i)*y); cosy(:,i)=cos(kf_fcont(i)*y)
        elseif (iforcing_cont(i)=='nocos') then
          if (lroot) print*,'forcing_cont: nocos flow'
          sinx(:,i)=sin(kf_fcont(i)*x)
          siny(:,i)=sin(kf_fcont(i)*y)
          sinz(:,i)=sin(kf_fcont(i)*z)
        elseif (iforcing_cont(i)=='KolmogorovFlow-x') then
          if (lroot) print*,'forcing_cont: Kolmogorov flow'
          cosx(:,i)=cos(kf_fcont(i)*x)
        elseif (iforcing_cont(i)=='KolmogorovFlow-z') then
          if (lroot) print*,'forcing_cont: Kolmogorov flow z'
          cosz(:,i)=cos(kf_fcont(i)*z)
        elseif (iforcing_cont(i)=='TG') then
          if (lroot) print*,'forcing_cont: TG'
          sinx(:,i)=sin(kf_fcont(i)*x); cosx(:,i)=cos(kf_fcont(i)*x)
          siny(:,i)=sin(kf_fcont(i)*y); cosy(:,i)=cos(kf_fcont(i)*y)
          cosz(:,i)=cos(kf_fcont(i)*z)
          if (lembed) then
            embedx(:,i)=.5+.5*tanh(x/width_fcont(i))*tanh((pi-x)/width_fcont(i))
            embedy(:,i)=.5+.5*tanh(y/width_fcont(i))*tanh((pi-y)/width_fcont(i))
            embedz(:,i)=.5+.5*tanh(z/width_fcont(i))*tanh((pi-z)/width_fcont(i))
            sinx(:,i)=embedx(:,i)*sinx(:,i); cosx(:,i)=embedx(:,i)*cosx(:,i)
            siny(:,i)=embedy(:,i)*siny(:,i); cosy(:,i)=embedy(:,i)*cosy(:,i)
            cosz(:,i)=embedz(:,i)*cosz(:,i)
          endif
        elseif (iforcing_cont(i)=='cosx*cosy*cosz') then
          if (lroot) print*,'forcing_cont: cosx(:,i)*cosy(:,i)*cosz(:,i)'
          sinx(:,i)=sin(kf_fcont(i)*x); cosx(:,i)=cos(kf_fcont(i)*x)
          siny(:,i)=sin(kf_fcont(i)*y); cosy(:,i)=cos(kf_fcont(i)*y)
          sinz(:,i)=sin(kf_fcont(i)*z); cosz(:,i)=cos(kf_fcont(i)*z)
        elseif (iforcing_cont(i)=='sinx') then
          sinx(:,i)=sin(kf_fcont(i)*x)
          if (tgentle(i) > 0.) then
            lgentle(i)=.true.
            if (lroot) print *, 'initialize_forcing: gentle forcing till t = ', tgentle
          endif
        elseif (iforcing_cont(i)=='(0,sinx,0)') then
          sinx(:,i)=sin(kf_fcont(i)*x)
        elseif (iforcing_cont(i)=='(0,0,cosx)') then
          cosx(:,i)=cos(kf_fcont(i)*x)
        elseif (iforcing_cont(i)=='(0,0,cosxcosy)') then
          cosx(:,i)=cos(kf_fcont(i)*x)
          cosy(:,i)=cos(kf_fcont(i)*y)
        elseif (iforcing_cont(i)=='B=(0,0,cosxcosy)') then
          sinx(:,i)=sin(kf_fcont(i)*x)
          siny(:,i)=sin(kf_fcont(i)*y)
          cosx(:,i)=cos(kf_fcont(i)*x)
          cosy(:,i)=cos(kf_fcont(i)*y)
        elseif (iforcing_cont(i)=='(sinz,cosz,0)') then
          sinz(:,i)=sin(kf_fcont(i)*z)
          cosz(:,i)=cos(kf_fcont(i)*z)
        elseif (iforcing_cont(i)=='(0,cosx*cosz,0)' &
           .or. iforcing_cont(i)=='(0,cosx*cosz,0)_Lor') then
          cosx(:,i)=cos(kf_fcont_x(i)*x)
          cosz(:,i)=cos(kf_fcont_z(i)*z)
          sinx(:,i)=sin(kf_fcont_x(i)*x)
          sinz(:,i)=sin(kf_fcont_z(i)*z)
          profy_ampl=exp(-(kf_fcont_y(i)*y)**2)
        elseif (iforcing_cont(i)=='(0,sinx*exp(-z^2),0)') then
          sinx(:,i)=sin(kf_fcont_x(i)*x)
          cosx(:,i)=cos(kf_fcont_x(i)*x)
        elseif (iforcing_cont(i)=='J0_k1x') then
          do l=l1,l2
            profx_ampl (l-l1+1)=ampl_ff(i) *bessj(0,k1bessel0*x(l))
            profx_ampl1(l-l1+1)=ampl1_ff(i)*bessj(1,k1bessel0*x(l))
          enddo
        elseif (iforcing_cont(i)=='gaussian-z') then
          profx_ampl=exp(-.5*x(l1:l2)**2/radius_ff**2)*ampl_ff(i)
          profy_ampl=exp(-.5*y**2/radius_ff**2)
          profz_ampl=exp(-.5*z**2/radius_ff**2)
        elseif (iforcing_cont(i)=='fluxring_cylindrical') then
          if (lroot) print*,'forcing_cont: fluxring cylindrical'
        elseif (iforcing_cont(i)=='vortex') then
          call torus_init(torus)
        elseif (iforcing_cont(i)=='KS') then
        !call periodic_KS_setup(-5./3.) !Kolmogorov spec. periodic KS
        !call random_isotropic_KS_setup(-5./3.,1.,(nxgrid)/2.) !old form
        call random_isotropic_KS_setup_test !Test KS model code with 3 specific
        !modes.
        endif
      enddo
      if (lroot .and. n_forcing_cont==0) &
        call warning('forcing','no valid continuous iforcing_cont specified')
!
!  minimal wavenumbers
!
      where( Lxyz/=0.) k1xyz=2.*pi/Lxyz
!
      if (r_ff /=0. .or. rcyl_ff/=0.) profz_k = tanh(z/width_ff)
!
    endsubroutine initialize_forcing
!***********************************************************************
    subroutine addforce(f)
!
!  Add forcing at the end of each time step.
!  Since forcing is constant during one time step,
!  this can be added as an Euler 1st order step
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      logical, save :: lfirstforce=.true., lfirstforce2=.true.
      logical, save :: llastforce=.true., llastforce2=.true.
!
!  Turn off forcing if t<tforce_start or t>tforce_stop.
!  This can be useful for producing good initial conditions
!  for turbulent decay experiments.
!
      call timing('addforce','entered')
!
      if ( (t>tforce_stop) .or. (t<tforce_start) ) then
        if ( (t>tforce_stop) .and. llastforce) then
          if (iforce/='zero' .and. lroot) print*, 'addforce: t>tforce_stop; no forcing'
          llastforce=.false.
        endif
      else
        if ( iforce/='zero' .and. lfirstforce .and. lroot ) &
!--         print*, 'addforce: addforce started'
        lfirstforce=.false.
!
!  calculate and add forcing function
!

        select case (iforce)
        case ('2drandom_xy');     call forcing_2drandom_xy(f)
        case ('2drxy_simple');    call forcing_2drandom_xy_simple(f)
        case ('ABC');             call forcing_ABC(f)
        case ('MHD_mode');        call forcing_mhd_mode(f)
        case ('blobs');           call forcing_blobs(f)
        case ('blobHS_random');   call forcing_blobHS_random(f)
        case ('chandra_kendall','cktest'); call forcing_chandra_kendall(f)
        case ('diffrot');         call forcing_diffrot(f,force)
        case ('fountain', '3');   call forcing_fountain(f)
        case ('hillrain');        call forcing_hillrain(f,force)
        case ('gaussianpot');     call forcing_gaussianpot(f,force)
        case ('white_noise');     call forcing_white_noise(f)
        case ('GP');              call forcing_GP(f)
        case ('Galloway-Proctor-92'); call forcing_GP92(f)
        case ('irrotational');    call forcing_irro(f,force)
        case ('helical', '2');    call forcing_hel(f)
        case ('helical_both');    call forcing_hel_both(f)
        case ('helical_kprof');   call forcing_hel_kprof(f)
        case ('hel_smooth');      call forcing_hel_smooth(f)
        case ('horiz-shear');     call forcing_hshear(f)
        case ('nocos');           call forcing_nocos(f)
        case ('TG');              call forcing_TG(f)
        case ('twist');           call forcing_twist(f)
        case ('tidal');           call forcing_tidal(f,force)
        case ('zero'); if (headt.and.ip<10) print*,'addforce: No forcing'
        case default
          if (lroot) print*,'addforce: No such forcing iforce=',trim(iforce)
        endselect
      endif
!
!  add *additional* forcing function
!
      if ( (t>tforce_stop2) .or. (t<tforce_start2) ) then
        if ( (t>tforce_stop2) .and. llastforce2 .and. lroot) &
            print*,'addforce: t>tforce_stop2; no forcing'
        if (t>tforce_stop2) llastforce2=.false.
      else
        if ( (iforce2/='zero') .and. lfirstforce2 .and. lroot) &
            print*, 'addforce: addforce2 started'
        lfirstforce2=.false.
!
        select case (iforce2)
        case ('diffrot');      call forcing_diffrot(f,force2)
        case ('fountain');     call forcing_fountain(f)
        case ('helical');      call forcing_hel(f)
        case ('helical_both'); call forcing_hel_both(f)
        case ('horiz-shear');  call forcing_hshear(f)
        case ('irrotational'); call forcing_irro(f,force2)
        case ('tidal');        call forcing_tidal(f,force2)
        case ('zero');
          if (headtt .and. lroot) print*,'addforce: No additional forcing'
        case default;
          if (lroot) print*,'addforce: No such forcing iforce2=',trim(iforce2)
        endselect
!
        if (headtt.or.ldebug) print*,'addforce: done addforce'
      endif
      call timing('addforce','finished')
!
    endsubroutine addforce
!***********************************************************************
    subroutine forcing_2drandom_xy(f)
!
!  Random force in two dimensions (x and y) limited to bands of wave-vector
!  in space.
!
!  14-feb-2011/ dhruba : coded
!
      use EquationOfState, only: cs0
      use General, only: random_number_wrapper
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical, save :: lfirst_call=.true.
      integer :: ikmodes,iran1,iran2,kx1,ky1,kx2,ky2
      real, dimension(nx,3) :: forcing_rhs
      real, dimension(nx) :: xkx1,xkx2
      real,dimension(4) :: fran
      real :: phase1,phase2,pi_over_Lx,force_norm
!
      if (lfirst_call) then
! If this is the first time this routine is being called select a set of
! k-vectors in two dimensions according to inputparameters:
        if (lroot) print*,'forcing_2drandom_xy: selecting k vectors'
        call get_2dmodes (.true.)
        allocate(random2d_kmodes (2,random2d_nmodes))
        call get_2dmodes (.false.)
        if (lroot) then
! The root processors also write out the forced modes.
          open(unit=10,file=trim(datadir)//'/2drandomk.out',status='unknown')
          do ikmodes = 1, random2d_nmodes
            write(10,*) random2d_kmodes(1,ikmodes),random2d_kmodes(2,ikmodes)
          enddo
        endif
        lfirst_call=.false.
      endif
!
! force = xhat [ cos (k_1 x + \phi_1 ) + cos (k_2 y + \phi_1) ] +
!         yhat [ cos (k_3 x + \phi_2 ) + cos (k_4 y + \phi_2) ]
! where k_1 -- k_2 and k_3 -- k_4 are two randomly chose pairs in the list
! random2d_kmodes and  \phi_1 and \phi_2 are two random phases.
!
      call random_number_wrapper(fran,CHANNEL=channel_force)
      phase1=pi*(2*fran(1)-1.)
      phase2=pi*(2*fran(2)-1.)
      iran1=random2d_nmodes*(.9999*fran(3))+1
      iran2=random2d_nmodes*(.9999*fran(4))+1
!
!  normally we want to use the wavevectors as they are,
!  but in some cases, e.g. when the box is bigger than 2pi,
!  we want to rescale k so that k=1 now corresponds to a smaller value.
!
      if (lscale_kvector_fac) then
        kx1=random2d_kmodes(1,iran1)*scale_kvectorx
        ky1=random2d_kmodes(2,iran1)*scale_kvectory
        kx2=random2d_kmodes(1,iran2)*scale_kvectorx
        ky2=random2d_kmodes(2,iran2)*scale_kvectory
        pi_over_Lx=0.5
      elseif (lscale_kvector_tobox) then
        kx1=random2d_kmodes(1,iran1)*(2.*pi/Lxyz(1))
        ky1=random2d_kmodes(2,iran1)*(2.*pi/Lxyz(2))
        kx2=random2d_kmodes(1,iran2)*(2.*pi/Lxyz(1))
        ky2=random2d_kmodes(2,iran2)*(2.*pi/Lxyz(2))
        pi_over_Lx=pi/Lxyz(1)
      else
        kx1=random2d_kmodes(1,iran1)
        ky1=random2d_kmodes(2,iran1)
        kx2=random2d_kmodes(1,iran2)
        ky2=random2d_kmodes(2,iran2)
        pi_over_Lx=0.5
      endif
!
! Now add the forcing
!
      force_norm = force*cs0*cs0*sqrt(dt)
      do n=n1,n2
        do m=m1,m2
          xkx1 = x(l1:l2)*kx1+phase1
          xkx2 = x(l1:l2)*kx2+phase2
          forcing_rhs(:,1) = force_norm*( cos(xkx1) + cos(y(m)*ky1+phase1) )
          forcing_rhs(:,2) = force_norm*( cos(xkx2) + cos(y(m)*ky2+phase2) )
          forcing_rhs(:,3) = 0.
          if (lhelical_test) then
            f(l1:l2,m,n,iuu:iuu+2)=forcing_rhs(:,1:3)
          else
            f(l1:l2,m,n,iuu:iuu+2)=f(l1:l2,m,n,iuu:iuu+2)+forcing_rhs(:,1:3)
          endif
        enddo
      enddo
!
      if (ip<=9) print*,'forcing_2drandom_xy: forcing OK'
!
    endsubroutine forcing_2drandom_xy
!***********************************************************************
    subroutine get_2dmodes (lonly_total)
!
      integer :: ik1,ik2,modk
      integer :: imode=1
      logical :: lonly_total
!
!  15-feb-2011/dhruba: coded
!
      imode=1
      do ik1=0,random2d_kmax
        do ik2 = 0, random2d_kmax
          modk = nint(sqrt ( float(ik1*ik1) + float (ik2*ik2) ))
          if ( (modk<=random2d_kmax).and.(modk>=random2d_kmin)) then
            if (.not.lonly_total) then
              random2d_kmodes(1,imode) = ik1
              random2d_kmodes(2,imode) = ik2
            endif
            imode=imode+1
          endif
        enddo
      enddo
      if (lonly_total) random2d_nmodes = imode
!
    endsubroutine get_2dmodes
!***********************************************************************
    subroutine forcing_2drandom_xy_simple(f)
!
!  Random force in two dimensions (x and y) limited to bands of wave-vector
!  in space.
!
!  14-feb-2011/ dhruba : coded
!
      use EquationOfState, only: cs0
      use General, only: random_number_wrapper, itoa
      use Sub, only: step
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension(nx,3) :: forcing_rhs
      real, dimension(nx) :: xkx
      real,dimension(4) :: fran
      real :: phase1,phase2,force_norm
!
! force = xhat [ cos (k_2 y + \phi_1) ] +
!         yhat [ cos (k_3 x + \phi_2 ) ]
! where k_1 -- k_2 and k_3 -- k_4 are two randomly chose pairs in the list
! random2d_kmodes and  \phi_1 and \phi_2 are two random phases.
!
      call random_number_wrapper(fran,CHANNEL=channel_force)
      phase1=pi*(2*fran(1)-1.)
      phase2=pi*(2*fran(2)-1.)
!
!  normally we want to use the wavevectors as they are,
!  but in some cases, e.g. when the box is bigger than 2pi,
!  we want to rescale k so that k=1 now corresponds to a smaller value.
!

!
! Now add the forcing
!
      force_norm = force*cs0*cs0*sqrt(dt)
      do n=n1,n2
        do m=m1,m2
          xkx = x(l1:l2)*kx_2df+phase1
          forcing_rhs(:,1) = force_norm*(cos(y(m)*ky_2df+phase2) )*profx_ampl*profy_ampl(m)
          forcing_rhs(:,2) = force_norm*(sin(xkx) )*profx_ampl*profy_ampl(m)
          forcing_rhs(:,3) = 0.
          if (lhelical_test) then
            f(l1:l2,m,n,iuu:iuu+2)=forcing_rhs(:,1:3)
          else
            f(l1:l2,m,n,iuu:iuu+2)=f(l1:l2,m,n,iuu:iuu+2)+forcing_rhs(:,1:3)
          endif
        enddo
      enddo
!
      if (ip<=9) print*,'forcing_2drandom_xy_simple: forcing OK'
!
    endsubroutine forcing_2drandom_xy_simple
!***********************************************************************
    subroutine forcing_irro(f,force_ampl)
!
!  add acoustic forcing function, using a set of precomputed wavevectors
!  This forcing drives pressure waves
!
!  10-sep-01/axel: coded
!   6-feb-13/MR: discard wavevectors [0,0,kz] if lavoid_yxmean
!  06-dec-13/nishant: made kkx etc allocatable
!  20-aug-14/MR: discard wavevectors [kx,0,kz] if lavoid_ymean, [kx,ky,0] if lavoid_zmean
!  21-jan-15/MR: changes for use of reference state.
!
      use General, only: random_number_wrapper
      use Mpicomm, only: mpireduce_sum,mpibcast_real
      use Sub, only: del2v_etc,dot
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: kx0,kx,ky,kz,force_ampl,pi_over_Lx
      real :: phase,ffnorm,iqfm,fsum,fsum_tmp
      real, dimension (2) :: fran
      real, dimension (nx) :: rho1,qf
      real, dimension (nx,3) :: forcing_rhs,curlo
      complex, dimension (mx) :: fx
      complex, dimension (my) :: fy
      complex, dimension (mz) :: fz
      complex, dimension (3) :: ikk
      integer :: ik,j,jf
!
      call keep_compiler_quiet(force_ampl)
!
      do
        call random_number_wrapper(fran,CHANNEL=channel_force)
        phase=pi*(2*fran(1)-1.)
        ik=nk*(.9999*fran(2))+1
!
!  if lavoid_xymean=T and wavevector is close enough to [0,0,kz] discard it
!  and look for a new one
!
        if ( lavoid_xymean ) then
          if ( abs(kkx(ik))>.9*k1xyz(1) .or. abs(kky(ik))>.9*k1xyz(2) ) exit
        elseif ( lavoid_ymean ) then
          if ( abs(kky(ik))>.9*k1xyz(2) ) exit
        elseif ( lavoid_zmean ) then
          if ( abs(kkz(ik))>.9*k1xyz(3) ) exit
        else
          exit
        endif
      enddo
!
      if (ip<=6) print*,'forcing_irro: ik,phase,kk=',ik,phase,kkx(ik),kky(ik),kkz(ik),dt
!
!  normally we want to use the wavevectors as they are,
!  but in some cases, e.g. when the box is bigger than 2pi,
!  we want to rescale k so that k=1 now corresponds to a smaller value.
!
      if (lscale_kvector_fac) then
        kx0=kkx(ik)*scale_kvectorx
        ky=kky(ik)*scale_kvectory
        kz=kkz(ik)*scale_kvectorz
        pi_over_Lx=0.5
      elseif (lscale_kvector_tobox) then
        kx0=kkx(ik)*(2.*pi/Lxyz(1))
        ky=kky(ik)*(2.*pi/Lxyz(2))
        kz=kkz(ik)*(2.*pi/Lxyz(3))
        pi_over_Lx=pi/Lxyz(1)
      else
        kx0=kkx(ik)
        ky=kky(ik)
        kz=kkz(ik)
        pi_over_Lx=0.5
      endif
!
!  in the shearing sheet approximation, kx = kx0 - St*k_y.
!  Here, St=-deltay/Lx. However, to stay near kx0, we ignore
!  integer shifts.
!
      if (Sshear==0.) then
        kx=kx0
      else
        if (lshearing_adjust_old) then
          kx=kx0+ky*deltay/Lx
        else
          kx=kx0+mod(ky*deltay/Lx-pi_over_Lx,2.*pi_over_Lx)+pi_over_Lx
        endif
      endif
!
!  Need to multiply by dt (for Euler step), but it also needs to be
!  divided by sqrt(dt), because square of forcing is proportional
!  to a delta function of the time difference.
!  Divide also by kav, so it would be ffnorm=force*sqrt(kav/dt)*dt/kav,
!  but this can be simplified to give:
!
      ffnorm=force*sqrt(dt/kav)
!
!  pre-calculate for the contributions e^(ikx*x), e^(iky*y), e^(ikz*z),
!  as well as the phase factor.
!
      
      fx=exp(cmplx(0.,kx*x+phase))*ffnorm
      fy=exp(cmplx(0.,ky*y))
      fz=exp(cmplx(0.,kz*z))
!
!  write i*k as vector
!
      ikk(1)=cmplx(0.,kx)
      ikk(2)=cmplx(0.,ky)
      ikk(3)=cmplx(0.,kz)
!
!  Loop over all directions, but skip over directions with no extent.
!
      iqfm=0.
      do n=n1,n2
      do m=m1,m2
!
!  Can force either velocity (default) or momentum (perhaps more physical).
!
        if (lmomentum_ff) then
          if (ldensity_nolog) then
            if (lreference_state) then
              rho1=1./(f(l1:l2,m,n,irho)+reference_state(:,iref_rho))
            else
              rho1=1./f(l1:l2,m,n,irho)
            endif
          else
            rho1=exp(-f(l1:l2,m,n,ilnrho))
          endif
        else
          rho1=1.
        endif
        do j=1,3
          if (lactive_dimension(j)) then
            jf=j+ifff-1
            forcing_rhs(:,j)=profx_ampl*profy_ampl(m)*profz_ampl(n) &
                *rho1*real(ikk(j)*fx(l1:l2)*fy(m)*fz(n))
            f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+forcing_rhs(:,j)
          endif
        enddo
        if (lout) then
          if (idiag_qfm/=0) then
            call del2v_etc(f,iuu,curlcurl=curlo)
            call dot(curlo,forcing_rhs,qf)
            iqfm=iqfm+sum(qf)
          endif
        endif
      enddo
      enddo
!
!  For printouts
!  On different processors, irufm needs to be communicated
!  to other processors.
!
      if (lout) then
        if (idiag_qfm/=0) then
          fsum_tmp=iqfm/nwgrid
          call mpireduce_sum(fsum_tmp,fsum)
          iqfm=fsum
          call mpibcast_real(iqfm)
          fname(idiag_qfm)=iqfm
          itype_name(idiag_qfm)=ilabel_sum
        endif
!
      endif
!
    endsubroutine forcing_irro
!***********************************************************************
    subroutine forcing_coefs_hel(coef1,coef2,coef3,fda,fx,fy,fz)
!
!  Calculates position-independent and 1D coefficients for helical forcing.
!
!  4-oct-17/MR: outsourced from forcing_hel.
!               Spotted bug: for old_forcing_evector=T, kk and ee remain undefined
!                            - needs to be fixed
      use Sub
!
      real,    dimension (3), intent(out) :: coef1,coef2,coef3,fda
      complex, dimension (mx),intent(out) :: fx
      complex, dimension (my),intent(out) :: fy
      complex, dimension (mz),intent(out) :: fz

      real :: phase, fact 
      real, dimension(3) :: kk
!
!  The routine fconst_coefs_hel generates various coefficients, including the phase.
!
      call fconst_coefs_hel(force,kkx,kky,kkz,nk,kav,coef1,coef2,coef3,kk,phase,fact,fda)
      call fxyz_coefs_hel(coef1,coef2,coef3,kk,phase,fact,fx,fy,fz)

    endsubroutine forcing_coefs_hel
!***********************************************************************
    subroutine forcing_pars_hel(coef1,coef2,coef3,fda,kk,phase,fact)
!
!  Calculates position-independent and 1D coefficients for helical forcing.
!  But this routine doesn't seem to be called from anywhere anymore.
!
!  4-oct-17/MR: outsourced from forcing_hel.
!               Spotted bug: for old_forcing_evector=T, kk and ee remain undefined
!                            - needs to be fixed
      use Sub
!
      real, dimension (3), intent(out) :: coef1,coef2,coef3,fda,kk
      real, intent(out) :: phase, fact 
!
      call fconst_coefs_hel(force,kkx,kky,kkz,nk,kav,coef1,coef2,coef3,kk,phase,fact,fda)

    endsubroutine forcing_pars_hel
!***********************************************************************
    subroutine forcing_coefs_hel2(coef1,coef2,coef3,fda,fx,fy,fz)
!
!  Modified copy of forcing_coefs_hel
!  Calculates position-independent and 1D coefficients for helical forcing.
!
!  10-nov-18/axel: adapted from forcing_coefs_hel to read also k_double.dat.
!
      use Sub
!
      real,    dimension (3), intent(out) :: coef1,coef2,coef3,fda
      complex, dimension (mx),intent(out) :: fx
      complex, dimension (my),intent(out) :: fy
      complex, dimension (mz),intent(out) :: fz
!
      real, dimension (3) :: kk
      real :: phase, fact 
!
!  This one also calculates coefficients, which are new ones, independent
!  of those in subroutine forcing_coefs_hel
!
      call fconst_coefs_hel(force_double,kkx2,kky2,kkz2,nk2,kav2,coef1,coef2,coef3,kk,phase,fact,fda)
      call fxyz_coefs_hel(coef1,coef2,coef3,kk,phase,fact,fx,fy,fz)
!
    endsubroutine forcing_coefs_hel2
!***********************************************************************
    subroutine fconst_coefs_hel(force_fact,kkx,kky,kkz,nk,kav,coef1,coef2,coef3,kk,phase,fact,fda)
!
!  This routine is can be called with any values of kkx,kky,kkz
!  to produce coef1,coef2,coef3,kk,phase,fact and fda.
!
!  08-aug-19/MR: modified to provide kk,phase,fact instead of fx,fy,fz
! 
      use General, only: random_number_wrapper,random_seed_wrapper
      use Sub
      use Mpicomm, only: stop_it
!      use Cdata, only: seed, nseed
!
      real,                   intent(in ) :: force_fact, kav
      integer,                intent(in ) :: nk
      real,    dimension (nk),intent(in ) :: kkx,kky,kkz
      real,    dimension (3), intent(out) :: coef1,coef2,coef3,kk,fda
      real,                   intent(out) :: phase,fact
!
      real :: ffnorm
      real, dimension (2) :: fran
      integer :: ik
      real :: kx0,kx,ky,kz,k2,k,pi_over_Lx
      real :: ex,ey,ez,kde,kex,key,kez,kkex,kkey,kkez
      real, dimension(3) :: e1,e2,ee
      real :: norm,phi
      real :: fd,fd2
!
!  generate random coefficients -1 < fran < 1
!  ff=force*Re(exp(i(kx+phase)))
!  |k_i| < akmax
!
      do
        call random_number_wrapper(fran,CHANNEL=channel_force)
!        call random_seed_wrapper(GET=seed)
!if (lroot) write(20,*) 'forcing: seed=', seed(1:nseed)
        phase=pi*(2*fran(1)-1.)
!AB: to add time-dependence XX
        ik=nk*(.9999*fran(2))+1
!
!  if lavoid_xymean=T and wavevector is close enough to [0,0,kz] discard it
!  and look for a new one
!
        if ( lavoid_xymean ) then
          if ( abs(kkx(ik))>.9*k1xyz(1) .or. abs(kky(ik))>.9*k1xyz(2) ) exit
        elseif ( lavoid_ymean ) then
          if ( abs(kky(ik))>.9*k1xyz(2) ) exit
        elseif ( lavoid_zmean ) then
          if ( abs(kkz(ik))>.9*k1xyz(3) ) exit
        else
          exit
        endif
      enddo
!
      if (ip<=6) then
        print*,'fconst_coefs_hel: ik,phase=',ik,phase
        print*,'fconst_coefs_hel: kx,ky,kz=',kkx(ik),kky(ik),kkz(ik)
      endif
!
!  normally we want to use the wavevectors as they are,
!  but in some cases, e.g. when the box is bigger than 2pi,
!  we want to rescale k so that k=1 now corresponds to a smaller value.
!
      if (lscale_kvector_fac) then
        kx0=kkx(ik)*scale_kvectorx
        ky=kky(ik)*scale_kvectory
        kz=kkz(ik)*scale_kvectorz
        pi_over_Lx=0.5
      elseif (lscale_kvector_tobox) then
        kx0=kkx(ik)*(2.*pi/Lxyz(1))
        ky=kky(ik)*(2.*pi/Lxyz(2))
        kz=kkz(ik)*(2.*pi/Lxyz(3))
        pi_over_Lx=pi/Lxyz(1)
      else
        kx0=kkx(ik)
        ky=kky(ik)
        kz=kkz(ik)
        pi_over_Lx=0.5
      endif
!
!  in the shearing sheet approximation, kx = kx0 - St*k_y.
!  Here, St=-deltay/Lx. However, to stay near kx0, we ignore
!  integer shifts.
!
      if (Sshear==0.) then
        kx=kx0
      else
        if (lshearing_adjust_old) then
          kx=kx0+ky*deltay/Lx
        else
          kx=kx0+mod(ky*deltay/Lx-pi_over_Lx,2.*pi_over_Lx)+pi_over_Lx
        endif
      endif
!
!  compute k^2 and output wavenumbers
!
      k2=kx**2+ky**2+kz**2
      k=sqrt(k2)
      if (ip<4) then
        open(89,file='forcing_hel_output.dat',position='append')
        write(89,'(6f10.5)') k,kx0,kx,ky,kz,deltay
        close(89)
      endif
!
!  Find e-vector:
!  Start with old method (not isotropic) for now.
!  Pick e1 if kk not parallel to ee1. ee2 else.
!
      if ((ky==0).and.(kz==0)) then
        ex=0; ey=1; ez=0
      else
        ex=1; ey=0; ez=0
      endif
!
      if (old_forcing_evector) then
!!! kk, ee not defined!
      else
!
!  Isotropize ee in the plane perp. to kk by
!  (1) constructing two basis vectors for the plane perpendicular
!      to kk, and
!  (2) choosing a random direction in that plane (angle phi)
!  Need to do this in order for the forcing to be isotropic.
!
        kk = (/kx, ky, kz/)
        ee = (/ex, ey, ez/)
        call cross(kk,ee,e1)
        call dot2(e1,norm); e1=e1/sqrt(norm) ! e1: unit vector perp. to kk
        call cross(kk,e1,e2)
        call dot2(e2,norm); e2=e2/sqrt(norm) ! e2: unit vector perp. to kk, e1
        call random_number_wrapper(phi,CHANNEL=channel_force); phi = phi*2*pi
        ee = cos(phi)*e1 + sin(phi)*e2
        ex=ee(1); ey=ee(2); ez=ee(3)
      endif
!
!  k.e
!
      call dot(kk,ee,kde)
!
!  k x e
!
      kex=ky*ez-kz*ey
      key=kz*ex-kx*ez
      kez=kx*ey-ky*ex
!
!  k x (k x e)
!
      kkex=ky*kez-kz*key
      kkey=kz*kex-kx*kez
      kkez=kx*key-ky*kex
!
!  ik x (k x e) + i*phase
!
!  Normalize ff; since we don't know dt yet, we finalize this
!  within timestep where dt is determined and broadcast.
!
!  This does already include the new sqrt(2) factor (missing in B01).
!  So, in order to reproduce the 0.1 factor mentioned in B01
!  we have to set force=0.07.
!
!  Furthermore, for |relhel| < 1, sqrt(2) should be replaced by
!  sqrt(1.+relhel**2). This is done now (9-nov-02).
!  This means that the previous value of force=0.07 (for relhel=0)
!  should now be replaced by 0.05.
!
!  Note: kav is not to be scaled with k1_ff (forcing should remain
!  unaffected when changing k1_ff).
!
      ffnorm = sqrt(1.+relhel**2) &
              *k*sqrt(k2-kde**2)/sqrt(kav*cs0eff**3)*(k/kav)**slope_ff
      if (ip<=9) then
        print*,'fconst_coefs_hel: k,kde,kav=',k,kde,kav
        print*,'fconst_coefs_hel: cs0eff,ffnorm=',cs0eff,ffnorm
        print*,'fconst_coefs_hel: k*sqrt(k2-kde**2)=',k*sqrt(k2-kde**2)
      endif
!
!  need to multiply by dt (for Euler step), but it also needs to be
!  divided by sqrt(dt), because square of forcing is proportional
!  to a delta function of the time difference
!
      fact=force_fact/ffnorm*sqrt(dt)
!
!  prefactor; treat real and imaginary parts separately (coef1 and coef2),
!  so they can be multiplied by different profiles below.
!
      coef1=k*(/kex,key,kez/)
      coef2=relhel*(/kkex,kkey,kkez/)
!
!  possibly multiply forcing by sgn(z) and radial profile
!
      if (rcyl_ff==0.) coef3=crosshel*k*(/kkex,kkey,kkez/)

      !if (ip<=5) print*,'fconst_coefs_hel: coef=',coef1,coef2
!
! An attempt to implement anisotropic forcing using direction
! dependent forcing amplitude. Activated only if force_strength,
! describing the anisotropic part of the forcing, is
! nonzero. force_direction, which is a vector, defines the preferred
! direction of forcing. The expression for the forcing amplitude used
! at the moment is:
!
!  f(i)=f0*[1+epsilon(delta_ij*(k(i)*fd(j))/(|k||fd|))^2*fd(i)/|fd|]
!
! here f0 and fd are shorthand for force and forcing_direction,
! respectively, and epsilon=force_strength/force.
!
! 09-feb-22/hongzhe: possibility of a new implementation. If
!                    laniso_forcing_old=F then fda is the same in all
!                    three directions, and fda=sqrt(1+epsilon*mu^2)
!                    where mu = cos of the angle between kk and fd
!
      if (force_strength/=0.) then
        call dot(force_direction,force_direction,fd2)
        fd=sqrt(fd2)
        if (laniso_forcing_old) then
          fda = 1. + (force_strength/force) &
                *(kk*force_direction/(k*fd))**2 * force_direction/fd
        else
          fda = sqrt(1. + (force_strength/force) &
                *(dot_product(kk,force_direction)/(k*fd))**2)
        endif
      else
        fda = 1.
      endif
!
    endsubroutine fconst_coefs_hel
!**************************************************************************
    subroutine fxyz_coefs_hel(coef1,coef2,coef3,kk,phase,fact,fx,fy,fz)
!
!  08-aug-19/MR: carved out to produce fx,fy,fz from the other parameters.
! 
      use Mpicomm, only: stop_it

      real,    dimension (3), intent(in ) :: coef1,coef2,coef3,kk
      real,                   intent(in ) :: phase, fact
      complex, dimension (mx),intent(out) :: fx
      complex, dimension (my),intent(out) :: fy
      complex, dimension (mz),intent(out) :: fz

      real :: kx, kz, phase_x

      if (ldiscrete_phases) then
        kx=kk(1)*k1_ff
        phase_x=nint(phase/(kx*dx))*dx
        fx=cmplx(cos(kx*(x+phase_x)),sin(kx*(x+phase_x)))*fact
      else
        fx=exp(cmplx(0.,kk(1)*k1_ff*x+phase))*fact
      endif

      fy=exp(cmplx(0.,kk(2)*k1_ff*y))
!
!  symmetry of forcing function about z direction
!
      kz = kk(3)
      select case (iforcing_zsym)
        case(0); fz=exp(cmplx(0.,kz*k1_ff*z))
        case(1); fz=cos(kz*k1_ff*z)
        case(-1); fz=sin(kz*k1_ff*z)
        case default; call stop_it('fxyz_coefs_hel: incorrect iforcing_zsym')
      endselect
!
!  possibly multiply forcing by z-profile
!  (This stuff is now supposed to be done in initialize; keep for now)
!
!-    if (height_ff/=0.) then
!-      if (lroot) print*,'forcing_hel: include z-profile'
!-      tmpz=(z/height_ff)**2
!-      fz=fz*exp(-tmpz**5/max(1.-tmpz,1e-5))
!-    endif
!
! need to discuss with axel
!
      if (rcyl_ff/=0.) then
!       if (lroot) &
!         print*,'fxyz_coefs_hel: applying sgn(z)*xi(r) profile'
        !
        ! only z-dependent part can be done here; radial stuff needs to go
!       ! into the loop
        !
        fz = fz*profz_k
      endif

      if (ip<=5) then
        print*,'fxyz_coefs_hel: fx=',fx
        print*,'fxyz_coefs_hel: fy=',fy
        print*,'fxyz_coefs_hel: fz=',fz
      endif
!
    endsubroutine fxyz_coefs_hel
!***********************************************************************
    subroutine forcing_hel(f)
!
!  Add helical forcing function, using a set of precomputed wavevectors.
!  The relative helicity of the forcing function is determined by the factor
!  sigma, called here also relhel. If it is +1 or -1, the forcing is a fully
!  helical Beltrami wave of positive or negative helicity. For |relhel| < 1
!  the helicity less than maximum. For relhel=0 the forcing is nonhelical.
!  The forcing function is now normalized to unity (also for |relhel| < 1).
!
!  10-apr-00/axel: coded
!   3-sep-02/axel: introduced k1_ff, to rescale forcing function if k1/=1.
!  25-sep-02/axel: preset force_ampl to unity (in case slope is not controlled)
!   9-nov-02/axel: corrected normalization factor for the case |relhel| < 1.
!  23-feb-10/axel: added helicity profile with finite second derivative.
!  13-jun-13/axel: option of symmetry of forcing function about z direction
!  06-dec-13/nishant: made kkx etc allocatable
!  20-aug-14/MR: discard wavevectors analogously to forcing_irrot
!  21-jan-15/MR: changes for use of reference state.
!  26-nov-15/AB: The cs0eff factor used for normalization can now be set.
!   4-oct-17/MR: outsourced forcing_coefs_hel for calculation of coefficients
!                which can be obtained outside an mn-loop. 
!                Moved density calculation into mn-loop where it belongs.
!                Spotted possible bug: for lforcing2_curl=T, calculation of forcing_rhs 
!                seems not to be meaningful.
!  12-jan-18/axel: added periodic forcing for omega_ff /= 0.
!   3-aug-22/axel: added omega_double_ff for second forcing function
!
      use Mpicomm
      use Sub
      use EquationOfState, only: rho0
      use DensityMethods, only: getrho1
!
      real, dimension (mx,my,mz,mfarray), intent(INOUT) :: f

      real :: irufm,iruxfxm,iruxfym,iruyfxm,iruyfym,iruzfzm,fsum_tmp,fsum
      real, dimension (nx) :: rho1,ruf,rho,force_ampl
      real, dimension (nx,3) :: variable_rhs,forcing_rhs,forcing_rhs2
      real, dimension (nx,3) :: forcing_rhs_old,forcing_rhs2_old
      real, dimension (nx,3) :: force_all
      real, dimension (3), save :: fda, fda2, fda_old, fda2_old
      complex, dimension (mx), save :: fx=0., fx2=0., fx_old=0., fx2_old=0.
      complex, dimension (my), save :: fy=0., fy2=0., fy_old=0., fy2_old=0.
      complex, dimension (mz), save :: fz=0., fz2=0., fz_old=0., fz2_old=0.
      real, dimension (3), save :: coef1,coef2,coef3,coef1b,coef2b,coef3b
      integer :: j,jf,j2f
      complex, dimension (nx) :: fxyz, fxyz2, fxyz_old, fxyz2_old
      real :: profyz, pforce, qforce, aforce, tmp
      real, dimension(3) :: profyz_hel_coef2, profyz_hel_coef2b
!
!  Compute forcing coefficients.
!
      if (t>=tsforce) then
        fx_old=fx
        fy_old=fy
        fz_old=fz
        fda_old=fda
        call forcing_coefs_hel(coef1,coef2,coef3,fda,fx,fy,fz)
!
!  Possibility of reading in data for second forcing function and
!  computing the relevant coefficients (fx2,fy2,fz2,fda2) here.
!  Unlike coef[1-3], where we add the letter b, we add a 2 for the
!  other coefficients for better readibility.
!
        if (lforcing_coefs_hel_double) then
          fx2_old=fx2
          fy2_old=fy2
          fz2_old=fz2
          fda2_old=fda2
          call forcing_coefs_hel2(coef1b,coef2b,coef3b,fda2,fx2,fy2,fz2)
        elseif (lmhd_forcing) then
          fx2_old=fx2
          fy2_old=fy2
          fz2_old=fz2
          fda2_old=fda2
          call forcing_coefs_hel(coef1b,coef2b,coef3b,fda2,fx2,fy2,fz2)
        endif
!
!  By default, dtforce=0, so new forcing is applied at every time step.
!  Alternatively, it can be set to any other time, so the forcing is
!  updated every dtforce.
!
        tsforce=t+dtforce
      endif
!
!  weight factor, pforce is the factor for the new term,
!  and qforce is that for the outgoing term.
!  If no outphasing is involved, pforce=1. and qforce=0.
!
      if (dtforce==0.) then
        pforce=1.
        qforce=0.
      else
        aforce=sin(   pi*(tsforce-t)/dtforce)**2*dtforce_ampl+(1.-dtforce_ampl)
        pforce=cos(.5*pi*(tsforce-t)/dtforce)**2*aforce
        qforce=sin(.5*pi*(tsforce-t)/dtforce)**2*aforce
      endif
!
!  loop the two cases separately, so we don't check for r_ff during
!  each loop cycle which could inhibit (pseudo-)vectorisation
!  calculate energy input from forcing; must use lout (not ldiagnos)
!
      irufm=0; iruxfxm=0; iruxfym=0; iruyfxm=0; iruyfym=0; iruzfzm=0
!
!  Here standard case. The case rcyl_ff==0 is to be removed sometime.
!
      if (rcyl_ff == 0) then
        do n=n1,n2
          do m=m1,m2
!
!  Compute useful shorthands for primary forcing function.
!
            profyz=profy_ampl(m)*profz_ampl(n)
            profyz_hel_coef2=profy_hel(m)*profz_hel(n)*coef2
!
!  Compute the combined complex forcing function fxyz.
!  If lforcing_coefs_hel_double is set, we also add a
!  contribution from fxyz2=fx2(l1:l2)*fy2(m)*fz2(n),
!  which is weighted with factor qdouble_profile(n).
!  qdouble_profile turns on fxyz2 in the upper parts.
!
            fxyz=fx(l1:l2)*fy(m)*fz(n)
            fxyz_old=fx_old(l1:l2)*fy_old(m)*fz_old(n)
            force_ampl=profx_ampl*profyz
!
!  Do the same for secondary forcing function.
!
            if (lforcing_coefs_hel_double.or.lmhd_forcing) then
              profyz_hel_coef2b=profy_hel(m)*profz_hel(n)*coef2b
              fxyz2=fx2(l1:l2)*fy2(m)*fz2(n)
              fxyz2_old=fx2_old(l1:l2)*fy2_old(m)*fz2_old(n)
            endif
!
!  Possibility of compute work done by forcing.
!
            if (lwork_ff) &
              force_ampl=force_ampl*calc_force_ampl(f,fx,fy,fz,profyz*cmplx(coef1,profyz_hel_coef2))
!
!  In the past we always forced the du/dt, but in some cases
!  it may be better to force rho*du/dt (if lmomentum_ff=.true.)
!  For compatibility with earlier results, lmomentum_ff=.false. by default.
!
            if (ldensity) then
              if (lmomentum_ff.or.lout) then
                call getrho1(f(:,m,n,ilnrho),rho1)
                if (lmomentum_ff) force_ampl=force_ampl*rho1
              endif
            endif

            do j=1,3
              if (lactive_dimension(j)) then
!
!  Primary forcing function: assemble here forcing_rhs(:,j).
!  Add here possibility of periodic forcing proportional to cos(om*t).
!  By default, omega_ff=0.
!
                forcing_rhs(:,j) = force_ampl*fda(j)*cos(omega_ff*t) &
                                  *real(cmplx(coef1(j),profx_hel*profyz_hel_coef2(j))*fxyz)

                if (qforce/=0.) &
                  forcing_rhs_old(:,j) = force_ampl*fda_old(j)*cos(omega_ff*t) &
                                        *real(cmplx(coef1(j),profx_hel*profyz_hel_coef2(j))*fxyz_old)

                if (lmhd_forcing) then
                  forcing_rhs2(:,j) = force_ampl*fda2(j)*cos(omega_ff*t) & 
                                     *real(cmplx(coef1b(j),profx_hel*profyz_hel_coef2b(j))*fxyz2)
                  if (qforce/=0.) &
                    forcing_rhs2_old(:,j) = force_ampl*fda2_old(j)*cos(omega_ff*t) &
                                           *real(cmplx(coef1b(j),profx_hel*profyz_hel_coef2b(j))*fxyz2_old)
                endif
!
!  Possibility of adding second forcing function.
! 
                if (lforcing_coefs_hel_double) then
                  forcing_rhs(:,j) = (1.-qdouble_profile(n))*forcing_rhs(:,j) &
                                    +qdouble_profile(n)*force_ampl*fda2(j)*cos(omega_double_ff*t) &
                                    *real(cmplx(coef1b(j),profx_hel*profyz_hel_coef2b(j))*fxyz2)
!
                  if (qforce/=0.) &
                    forcing_rhs_old(:,j) = (1.-qdouble_profile(n))*forcing_rhs_old(:,j) &
                                          +qdouble_profile(n)*force_ampl*fda2_old(j)*cos(omega_double_ff*t) &
                                          *real(cmplx(coef1b(j),profx_hel*profyz_hel_coef2b(j))*fxyz2_old)
                endif
!
!  assemble new combination
!
                if (qforce/=0.) &
                  forcing_rhs(:,j) = pforce*forcing_rhs(:,j) + qforce*forcing_rhs_old(:,j)
!
! put force into auxiliary variable, if requested
!
                if (lff_as_aux) f(l1:l2,m,n,iff+j-1) = f(l1:l2,m,n,iff+j-1)+forcing_rhs(:,j)
!
!  Compute additional forcing function (used for velocity if crosshel=1).
!  It can optionally be the same. Alternatively, one has to set crosshel=1.
!
                if (lforcing2_same) then
                  forcing_rhs2(:,j)=forcing_rhs(:,j)
                elseif (lforcing2_curl) then
!
!  Something seems to be missing here as real(cmplx(0.,coef3(j))) is zero.
!
                  forcing_rhs2(:,j) = force_ampl*real(cmplx(0.,coef3(j)))*fda(j) 
                else
                  forcing_rhs2(:,j) = force_ampl*real(cmplx(0.,coef3(j))*fxyz)*fda(j)
                endif
!
!  Choice of different possibilities.
!  Here is the decisive step where forcing is applied.
!  Added possibility of linearly ramping down the forcing in time.
!
                if (ifff/=0) then

                  jf=j+ifff-1
                  j2f=j+i2fff-1
!
                  if (lhelical_test) then
                    f(l1:l2,m,n,jf)=forcing_rhs(:,j)
                  else
                    if (lforce_ramp_down) then
                      tmp=max(0d0,1d0+min(0d0,(tforce_ramp_down-t)/tauforce_ramp_down))
                      f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+forcing_rhs(:,j)*force1_scl*tmp
                    else
                      f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+forcing_rhs(:,j)*force1_scl
                    endif
!
!  Allow here for forcing both in u and in b=curla. In that case one sets
!  lhydro_forcing=T and lmagnetic_forcing=T.
!
                    if (lmhd_forcing) then 
                      f(l1:l2,m,n,j2f)=f(l1:l2,m,n,j2f)+forcing_rhs2(:,j)*force2_scl
                    elseif (lcrosshel_forcing) then
                      f(l1:l2,m,n,j2f)=f(l1:l2,m,n,j2f)+forcing_rhs(:,j)*force2_scl
                    endif
!
!  Forcing with enhanced xx correlation.
!
                    if (j==1.and.lxxcorr_forcing) &
                      f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+forcing_rhs(:,1)*force2_scl
                  endif
                endif
!
!  If one of the testfield methods is used, we need to add a forcing term
!  in one of the auxiliary equations. Their location is denoted by jtest_aa0
!  and jtest_uu0 for the testfield and testflow equations, respectively.
!  In the testflow module, jtest_uu0=1 is used, while in the testfield_nonlinear_z
!  module jtest_uu0=5 is used, so for now we give them by hand.
!
                if (ltestfield_forcing) then
                  jf=j+iaatest-1+3*(jtest_aa0-1)
                  f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+forcing_rhs(:,j)*force1_scl
                endif
                if (ltestflow_forcing) then
                  jf=j+iuutest-1+3*(jtest_uu0-1)
                  if (ltestfield_forcing) then
                    f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+forcing_rhs2(:,j)*force2_scl
                  else
                    f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+forcing_rhs(:,j)*force2_scl
                  endif
                endif
              endif
            enddo
!
!  Forcing with enhanced xy correlation.
!  Can only come outside the previous j loop.
!  This is apparently not yet applied to testfield or testflow forcing.
!
            if (ifff/=0.and..not.lhelical_test.and.lxycorr_forcing) then
              do j=1,3
                if (lactive_dimension(j)) then
                  jf=j+ifff-1
                  if (j==1) f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf) &
                           +forcing_rhs(:,2)*force2_scl
                  if (j==2) f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf) &
                           +forcing_rhs(:,1)*force2_scl
                endif
              enddo
            endif
!
!  Sum up.
!
            if (lout) then
              if (idiag_rufm/=0 .or. &
                  idiag_ruxfxm/=0 .or. idiag_ruxfym/=0 .or. &
                  idiag_ruyfxm/=0 .or. idiag_ruyfym/=0 .or. &
                  idiag_ruzfzm/=0) then
!
!  Compute density.
!
                if (ldensity) then
                  rho=1./rho1
                else
                  rho=rho0
                endif
!
!  Evaluate the various choices.
!
                if (idiag_rufm/=0) then
!
!  Compute rhs.
!
                  variable_rhs=f(l1:l2,m,n,iffx:iffz)
                  call multsv_mn(rho/dt,forcing_rhs,force_all)
                  call dot_mn(variable_rhs,force_all,ruf)
                  irufm=irufm+sum(ruf)
                endif
                if (idiag_ruxfxm/=0) iruxfxm=iruxfxm+sum( &
                    rho*f(l1:l2,m,n,iux)*forcing_rhs(:,1))
                if (idiag_ruxfym/=0) iruxfym=iruxfym+sum( &
                    rho*f(l1:l2,m,n,iux)*forcing_rhs(:,2))
                if (idiag_ruyfxm/=0) iruyfxm=iruyfxm+sum( &
                    rho*f(l1:l2,m,n,iuy)*forcing_rhs(:,1))
                if (idiag_ruyfym/=0) iruyfym=iruyfym+sum( &
                    rho*f(l1:l2,m,n,iuy)*forcing_rhs(:,2))
                if (idiag_ruzfzm/=0) iruzfzm=iruzfzm+sum( &
                    rho*f(l1:l2,m,n,iuz)*forcing_rhs(:,3))
              endif
            endif
!
!  End of mn loop.
!
          enddo
        enddo
      else
!
!  Radial profile, but this is old fashioned and probably no longer used.
!
call fatal_error('forcing_hel','check that radial profile with rcyl_ff/=0. works ok')
        do j=1,3
          if (lactive_dimension(j)) then
            jf=j+ifff-1
            do n=n1,n2
              do m=m1,m2
                if (ldensity) then
                  call getrho1(f(:,m,n,ilnrho),rho1)
                else
                  rho1=1./rho0
                endif
                forcing_rhs(:,j)=rho1*profx_ampl*profy_ampl(m)*profz_ampl(n) &
                  *real(cmplx(coef1(j),profy_hel(m)*profz_k(n)*coef2(j)) &
                  *fx(l1:l2)*fy(m)*fz(n))
                  if (lhelical_test) then
                    f(l1:l2,m,n,jf)=forcing_rhs(:,j)
                  else
                    f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+forcing_rhs(:,j)
                  endif
              enddo
            enddo
          endif
        enddo
      endif
!
!  For printouts
!  On different processors, irufm needs to be communicated
!  to other processors.
!
      if (lout) then
        if (idiag_rufm/=0) then
          fsum_tmp=irufm
          call mpireduce_sum(fsum_tmp,fsum)
          irufm=fsum
          call mpibcast_real(irufm)
          fname(idiag_rufm)=irufm
          itype_name(idiag_rufm)=ilabel_sum
        endif
        if (idiag_ruxfxm/=0) then
          fsum_tmp=iruxfxm
          call mpireduce_sum(fsum_tmp,fsum)
          iruxfxm=fsum
          call mpibcast_real(iruxfxm)
          fname(idiag_ruxfxm)=iruxfxm
          itype_name(idiag_ruxfxm)=ilabel_sum
        endif
        if (idiag_ruxfym/=0) then
          fsum_tmp=iruxfym
          call mpireduce_sum(fsum_tmp,fsum)
          iruxfym=fsum
          call mpibcast_real(iruxfym)
          fname(idiag_ruxfym)=iruxfym
          itype_name(idiag_ruxfym)=ilabel_sum
        endif
        if (idiag_ruyfxm/=0) then
          fsum_tmp=iruyfxm
          call mpireduce_sum(fsum_tmp,fsum)
          iruyfxm=fsum
          call mpibcast_real(iruyfxm)
          fname(idiag_ruyfxm)=iruyfxm
          itype_name(idiag_ruyfxm)=ilabel_sum
        endif
        if (idiag_ruyfym/=0) then
          fsum_tmp=iruyfym
          call mpireduce_sum(fsum_tmp,fsum)
          iruyfym=fsum
          call mpibcast_real(iruyfym)
          fname(idiag_ruyfym)=iruyfym
          itype_name(idiag_ruyfym)=ilabel_sum
        endif
        if (idiag_ruzfzm/=0) then
          fsum_tmp=iruzfzm
          call mpireduce_sum(fsum_tmp,fsum)
          iruzfzm=fsum
          call mpibcast_real(iruzfzm)
          fname(idiag_ruzfzm)=iruzfzm
          itype_name(idiag_ruzfzm)=ilabel_sum
        endif
      endif
!
      if (ip<=9) print*,'forcing_hel: forcing OK'
!
    endsubroutine forcing_hel
!***********************************************************************
    subroutine forcing_hel_kprof(f)
!
!  Add helical forcing function, using a set of precomputed wavevectors.
!  The relative helicity of the forcing function is determined by the factor
!  sigma, called here also relhel. If it is +1 or -1, the forcing is a fully
!  helical Beltrami wave of positive or negative helicity. For |relhel| < 1
!  the helicity less than maximum. For relhel=0 the forcing is nonhelical.
!  The forcing function is now normalized to unity (also for |relhel| < 1).
!
!  30-jan-11/axel: adapted from forcing_hel and added z-dependent scaling of k
!  06-dec-13/nishant: made kkx etc allocatable
!
      use Diagnostics, only: sum_mn_name
      use EquationOfState, only: cs0,rho0
      use General, only: random_number_wrapper
      use Sub, only: del2v_etc,curl,cross,dot,dot2
      use DensityMethods, only: getrho1, getrho
!
      real :: phase,ffnorm
      real, dimension (2) :: fran
      real, dimension (nx) :: rho1,ff,uf,of,qf,rho
      real, dimension (mz) :: kfscl
      real, dimension (nx,3) :: variable_rhs,forcing_rhs,forcing_rhs2
      real, dimension (nx,3) :: fda,uu,oo,bb,fxb,curlo
      real, dimension (mx,my,mz,mfarray) :: f
      complex, dimension (mx) :: fx
      complex, dimension (my) :: fy
      complex, dimension (mz) :: fz
      real, dimension (3) :: coef1,coef2,coef3
      integer :: ik,j,jf,j2f
      real :: kx0,kx,ky,kz,k2,k,force_ampl,pi_over_Lx
      real :: ex,ey,ez,kde,sig,fact,kex,key,kez,kkex,kkey,kkez
      real, dimension(3) :: e1,e2,ee,kk
      real :: norm,phi
      real :: fd,fd2
!
!  generate random coefficients -1 < fran < 1
!  ff=force*Re(exp(i(kx+phase)))
!  |k_i| < akmax
!
      call random_number_wrapper(fran,CHANNEL=channel_force)
      phase=pi*(2*fran(1)-1.)
      ik=nk*(.9999*fran(2))+1
      if (ip<=6) then
        print*,'forcing_hel_kprof: ik,phase=',ik,phase
        print*,'forcing_hel_kprof: kx,ky,kz=',kkx(ik),kky(ik),kkz(ik)
        print*,'forcing_hel_kprof: dt=',dt
      endif
!
!  compute z-dependent scaling factor, regard kav as kmax and write
!    kinv=1./kav+(1.-1./kav)*(ztop-z)/(ztop-zbot)
!    kfscl=exp((z(n)-xyz0(3))/Lxyz(3))/kav
!  i.e. divide by kav, so that k=1 when z=zbot.
!
      kfscl=1./(1.+(kav-1.)*(xyz0(3)+Lxyz(3)-z)/Lxyz(3))
      !kfscl=exp((z-xyz0(3))/Lxyz(3))/kav
!
!  Loop over all k, but must have the call to random_number_wrapper
!  outside this loop.
!
      call random_number_wrapper(phi,CHANNEL=channel_force); phi = phi*2*pi
      do n=n1,n2
!
!  normally we want to use the wavevectors as they are,
!  but in some cases, e.g. when the box is bigger than 2pi,
!  we want to rescale k so that k=1 now corresponds to a smaller value.
!
        if (lscale_kvector_fac) then
          kx0=kkx(ik)*scale_kvectorx
          ky=kky(ik)*scale_kvectory
          kz=kkz(ik)*scale_kvectorz
          pi_over_Lx=0.5
        elseif (lscale_kvector_tobox) then
          kx0=kkx(ik)*(2.*pi/Lxyz(1))
          ky=kky(ik)*(2.*pi/Lxyz(2))
          kz=kkz(ik)*(2.*pi/Lxyz(3))
          pi_over_Lx=pi/Lxyz(1)
        else
          kx0=kkx(ik)
          ky=kky(ik)
          kz=kkz(ik)
          pi_over_Lx=0.5
        endif
!
!  scale all components of k
!
        kx0=kx0*kfscl(n)
        ky=ky*kfscl(n)
        kz=kz*kfscl(n)
!
!  in the shearing sheet approximation, kx = kx0 - St*k_y.
!  Here, St=-deltay/Lx. However, to stay near kx0, we ignore
!  integer shifts.
!
        if (Sshear==0.) then
          kx=kx0
        else
          if (lshearing_adjust_old) then
            kx=kx0+ky*deltay/Lx
          else
            kx=kx0+mod(ky*deltay/Lx-pi_over_Lx,2.*pi_over_Lx)+pi_over_Lx
          endif
        endif
!
!  compute k^2 and output wavenumbers
!
        k2=kx**2+ky**2+kz**2
        k=sqrt(k2)
        if (ip<4) then
          open(89,file='forcing_hel_kprof_output.dat',position='append')
          write(89,'(6f10.5)') k,kx0,kx,ky,kz,deltay
          close(89)
        endif
!
!  Find e-vector:
!  Start with old method (not isotropic) for now.
!  Pick e1 if kk not parallel to ee1. ee2 else.
!
        if ((ky==0).and.(kz==0)) then
          ex=0; ey=1; ez=0
        else
          ex=1; ey=0; ez=0
        endif
        if (.not. old_forcing_evector) then
   !
   !  Isotropize ee in the plane perp. to kk by
   !  (1) constructing two basis vectors for the plane perpendicular
   !      to kk, and
   !  (2) choosing a random direction in that plane (angle phi)
   !  Need to do this in order for the forcing to be isotropic.
   !
          kk = (/kx, ky, kz/)
          ee = (/ex, ey, ez/)
          call cross(kk,ee,e1)
          call dot2(e1,norm); e1=e1/sqrt(norm) ! e1: unit vector perp. to kk
          call cross(kk,e1,e2)
          call dot2(e2,norm); e2=e2/sqrt(norm) ! e2: unit vector perp. to kk, e1
          ee = cos(phi)*e1 + sin(phi)*e2
          ex=ee(1); ey=ee(2); ez=ee(3)
        endif
!
!  k.e
!
        call dot(kk,ee,kde)
!
!  k x e
!
        kex=ky*ez-kz*ey
        key=kz*ex-kx*ez
        kez=kx*ey-ky*ex
!
!  k x (k x e)
!
        kkex=ky*kez-kz*key
        kkey=kz*kex-kx*kez
        kkez=kx*key-ky*kex
!
!  ik x (k x e) + i*phase
!
!  Normalize ff; since we don't know dt yet, we finalize this
!  within timestep where dt is determined and broadcast.
!
!  This does already include the new sqrt(2) factor (missing in B01).
!  So, in order to reproduce the 0.1 factor mentioned in B01
!  we have to set force=0.07.
!
!  Furthermore, for |relhel| < 1, sqrt(2) should be replaced by
!  sqrt(1.+relhel**2). This is done now (9-nov-02).
!  This means that the previous value of force=0.07 (for relhel=0)
!  should now be replaced by 0.05.
!
!  Note: kav is not to be scaled with k1_ff (forcing should remain
!  unaffected when changing k1_ff).
!
        ffnorm=sqrt(1.+relhel**2) &
          *k*sqrt(k2-kde**2)/sqrt(kav*cs0**3)*(k/kav)**slope_ff
        if (ip<=9) then
          print*,'forcing_hel_kprof: k,kde,ffnorm,kav=',k,kde,ffnorm,kav
          print*,'forcing_hel_kprof: k*sqrt(k2-kde**2)=',k*sqrt(k2-kde**2)
        endif
!
!  need to multiply by dt (for Euler step), but it also needs to be
!  divided by sqrt(dt), because square of forcing is proportional
!  to a delta function of the time difference
!
        fact=force/ffnorm*sqrt(dt)
!
!  The wavevector is for the case where Lx=Ly=Lz=2pi. If that is not the
!  case one needs to scale by 2pi/Lx, etc.
!
        fx=exp(cmplx(0.,kx*k1_ff*x+phase))*fact
        fy=exp(cmplx(0.,ky*k1_ff*y))
        fz=exp(cmplx(0.,kz*k1_ff*z))
!
!  possibly multiply forcing by sgn(z) and radial profile
!
        if (rcyl_ff/=0.) then
          if (lroot) print*,'forcing_hel_kprof: applying sgn(z)*xi(r) profile'
          !
          ! only z-dependent part can be done here; radial stuff needs to go
          ! into the loop
          !
          fz = fz*profz_k
        endif
!
        if (ip<=5) then
          print*,'forcing_hel_kprof: fx=',fx
          print*,'forcing_hel_kprof: fy=',fy
          print*,'forcing_hel_kprof: fz=',fz
        endif
!
!  prefactor; treat real and imaginary parts separately (coef1 and coef2),
!  so they can be multiplied by different profiles below.
!
        coef1(1)=k*kex; coef2(1)=relhel*kkex; coef3(1)=crosshel*k*kkex
        coef1(2)=k*key; coef2(2)=relhel*kkey; coef3(2)=crosshel*k*kkey
        coef1(3)=k*kez; coef2(3)=relhel*kkez; coef3(3)=crosshel*k*kkez
        if (ip<=5) print*,'forcing_hel_kprof: coef=',coef1,coef2
!
! An attempt to implement anisotropic forcing using direction
! dependent forcing amplitude. Activated only if force_strength,
! describing the anisotropic part of the forcing, is
! nonzero. force_direction, which is a vector, defines the preferred
! direction of forcing. The expression for the forcing amplitude used
! at the moment is:
!
!  f(i)=f0*[1+epsilon(delta_ij*(k(i)*fd(j))/(|k||fd|))^2*fd(i)/|fd|]
!
! here f0 and fd are shorthand for force and forcing_direction,
! respectively, and epsilon=force_strength/force.
!
        if (force_strength/=0.) then
          call dot(force_direction,force_direction,fd2)
          fd=sqrt(fd2)
          do j=1,3
             fda(:,j) = 1. + (force_strength/force) &
                  *(kk(j)*force_direction(j)/(k*fd))**2 &
                  *force_direction(j)/fd
          enddo
        else
          fda = 1.
        endif
!
!  loop the two cases separately, so we don't check for r_ff during
!  each loop cycle which could inhibit (pseudo-)vectorisation
!  calculate energy input from forcing; must use lout (not ldiagnos)
!
        force_ampl=1.0
        if (rcyl_ff == 0) then       ! no radial profile
          do m=m1,m2
!
!  In the past we always forced the du/dt, but in some cases
!  it may be better to force rho*du/dt (if lmomentum_ff=.true.)
!  For compatibility with earlier results, lmomentum_ff=.false. by default.
!
            if (ldensity.and.lmomentum_ff) then
              call getrho1(f(:,m,n,ilnrho),rho1)
            else
              rho1=1./rho0
            endif
            if (lwork_ff) force_ampl=calc_force_ampl(f,fx,fy,fz,profy_ampl(m)*profz_ampl(n) &
                                                     *cmplx(coef1,profy_hel(m)*profz_hel(n)*coef2))
            variable_rhs=f(l1:l2,m,n,iffx:iffz)
            do j=1,3
              if (lactive_dimension(j)) then
!
                forcing_rhs(:,j)=rho1*profx_ampl*profy_ampl(m)*profz_ampl(n)*force_ampl &
                  *real(cmplx(coef1(j),profx_hel*profy_hel(m)*profz_hel(n)*coef2(j)) &
                  *fx(l1:l2)*fy(m)*fz(n))*fda(:,j)
!
                forcing_rhs2(:,j)=rho1*profx_ampl*profy_ampl(m)*profz_ampl(n)*force_ampl &
                  *real(cmplx(0.,coef3(j)) &
                  *fx(l1:l2)*fy(m)*fz(n))*fda(:,j)
!
                if (ifff/=0) then
!
                  jf=j+ifff-1
                  j2f=j+i2fff-1
!
                  if (lhelical_test) then
                    f(l1:l2,m,n,jf)=forcing_rhs(:,j)
                  else
                    f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+forcing_rhs(:,j)*force1_scl
!
!  allow here for forcing both in u and in b=curla. In that case one sets
!  lhydro_forcing=T and lmagnetic_forcing=T.
!
                    if (lmhd_forcing) then
                      f(l1:l2,m,n,j2f)=f(l1:l2,m,n,j2f)+forcing_rhs2(:,j)*force2_scl
                    elseif (lcrosshel_forcing) then
                      f(l1:l2,m,n,j2f)=f(l1:l2,m,n,j2f)+forcing_rhs(:,j)*force2_scl
                    endif
                  endif
!
                endif
!
!  If one of the testfield methods is used, we need to add a forcing term
!  in one of the auxiliary equations. Their location is denoted by jtest_aa0
!  and jtest_uu0 for the testfield and testflow equations, respectively.
!  In the testflow module, jtest_uu0=1 is used, while in the testfield_nonlinear
!  module jtest_uu0=5 is used, so for now we give them by hand.
!
                if (ltestfield_forcing) then
                  jf=j+iaatest-1+3*(jtest_aa0-1)
                  f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+forcing_rhs(:,j)*force1_scl
                endif
                if (ltestflow_forcing) then
                  jf=j+iuutest-1+3*(jtest_uu0-1)
                  if (ltestfield_forcing) then
                    f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+forcing_rhs2(:,j)*force2_scl
                  else
                    f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+forcing_rhs(:,j)*force2_scl
                  endif
                endif
              endif
            enddo
!
!  For printouts:
!
            if (lout) then
              if (ldensity.and.idiag_rufm/=0) then

                if (lmomentum_ff) then
                  rho=1./rho1
                else
                  call getrho(f(:,m,n,ilnrho),rho)
                endif

                uu=f(l1:l2,m,n,iux:iuz)
                call dot(uu,forcing_rhs,uf)
                call sum_mn_name(rho*uf,idiag_rufm)

              endif
              if (idiag_ufm/=0) then
                uu=f(l1:l2,m,n,iux:iuz)
                call dot(uu,forcing_rhs,uf)
                call sum_mn_name(uf,idiag_ufm)
              endif
              if (idiag_ofm/=0) then
                call curl(f,iuu,oo)
                call dot(oo,forcing_rhs,of)
                call sum_mn_name(of,idiag_ofm)
              endif
              if (idiag_qfm/=0) then
                call del2v_etc(f,iuu,curlcurl=curlo)
                call dot(curlo,forcing_rhs,qf)
                call sum_mn_name(of,idiag_qfm)
              endif
              if (idiag_ffm/=0) then
                call dot2(forcing_rhs,ff)
                call sum_mn_name(ff,idiag_ffm)
              endif
              if (lmagnetic) then
                if (idiag_fxbxm/=0.or.idiag_fxbym/=0.or.idiag_fxbzm/=0) then
                  call curl(f,iaa,bb)
                  call cross(forcing_rhs,bb,fxb)
                  call sum_mn_name(fxb(:,1),idiag_fxbxm)
                  call sum_mn_name(fxb(:,2),idiag_fxbym)
                  call sum_mn_name(fxb(:,3),idiag_fxbzm)
                endif
              endif
            endif
          enddo

        else
!
!  Radial profile, but this is old fashioned and probably no longer used.
!
          call fatal_error('forcing_hel_kprof','check that radial profile with rcyl_ff works ok')
          do j=1,3
            if (lactive_dimension(j)) then
              jf=j+ifff-1
              sig = relhel*profz_k(n)
              coef1(1)=cmplx(k*kex,sig*kkex)
              coef1(2)=cmplx(k*key,sig*kkey)
              coef1(3)=cmplx(k*kez,sig*kkez)
              do m=m1,m2
!
!  In the past we always forced the du/dt, but in some cases
!  it may be better to force rho*du/dt (if lmomentum_ff=.true.)
!  For compatibility with earlier results, lmomentum_ff=.false. by default.
!
                if (ldensity.and.lmomentum_ff) then
                  call getrho1(f(:,m,n,ilnrho),rho1)
                else
                  rho1=1./rho0
                endif

                forcing_rhs(:,j)=rho1 &
                  *profx_ampl*profy_ampl(m)*profz_ampl(n) &
                  *real(cmplx(coef1(j),profy_hel(m)*coef2(j)) &
                  *fx(l1:l2)*fy(m)*fz(n))
                if (lhelical_test) then
                  f(l1:l2,m,n,jf)=forcing_rhs(:,j)
                else
                  f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+forcing_rhs(:,j)
                endif
              enddo
            endif
          enddo
        endif
!
      enddo
!
      if (ip<=9) print*,'forcing_hel_kprof: forcing OK'
!
    endsubroutine forcing_hel_kprof
!***********************************************************************
    subroutine forcing_hel_both(f)
!
!  Add helical forcing function, using a set of precomputed wavevectors.
!  The relative helicity of the forcing function is determined by the factor
!  sigma, called here also relhel. If it is +1 or -1, the forcing is a fully
!  helical Beltrami wave of positive or negative helicity. For |relhel| < 1
!  the helicity less than maximum. For relhel=0 the forcing is nonhelical.
!  The forcing function is now normalized to unity (also for |relhel| < 1).
!  This adds positive helical forcing to the "northern hemisphere" (y above the
!  midplane and negative helical forcing to the "southern hemisphere". The
!  two forcing are merged at the "equator" (midplane) where both are smoothly
!  set to zero.
!
!  22-sep-08/dhruba: adapted from forcing_hel
!   6-oct-09/MR: according to Axel, this routine is now superseded by forcing_hel and should be deleted
!  06-dec-13/nishant: made kkx etc allocatable
!
      use EquationOfState, only: cs0
      use General, only: random_number_wrapper
      use Sub
!
      real :: phase,ffnorm
      real, dimension (2) :: fran
      real, dimension (nx) :: rho1
      real, dimension (nx,3) :: variable_rhs,forcing_rhs
      real, dimension (mx,my,mz,mfarray) :: f
      complex, dimension (mx) :: fx
      complex, dimension (my) :: fy
      complex, dimension (mz) :: fz
      real, dimension (3) :: coef1,coef2
      integer :: ik,j,jf
      real :: kx0,kx,ky,kz,k2,k
      real :: ex,ey,ez,kde,fact,kex,key,kez,kkex,kkey,kkez
      real, dimension(3) :: e1,e2,ee,kk
      real :: norm,phi,pi_over_Lx
!
!  additional stuff for test fields
!
      integer :: jtest
!
!  generate random coefficients -1 < fran < 1
!  ff=force*Re(exp(i(kx+phase)))
!  |k_i| < akmax
!
      call random_number_wrapper(fran,CHANNEL=channel_force)
      phase=pi*(2*fran(1)-1.)
      ik=nk*(.9999*fran(2))+1
      if (ip<=6) then 
        print*,'forcing_hel_both: ik,phase=',ik,phase
        print*,'forcing_hel_both: kx,ky,kz=',kkx(ik),kky(ik),kkz(ik)
        print*,'forcing_hel_both: dt=',dt
      endif
!
!  normally we want to use the wavevectors as they are,
!  but in some cases, e.g. when the box is bigger than 2pi,
!  we want to rescale k so that k=1 now corresponds to a smaller value.
!
      if (lscale_kvector_fac) then
        kx0=kkx(ik)*scale_kvectorx
        ky=kky(ik)*scale_kvectory
        kz=kkz(ik)*scale_kvectorz
        pi_over_Lx=0.5
      elseif (lscale_kvector_tobox) then
        kx0=kkx(ik)*(2.*pi/Lxyz(1))
        ky=kky(ik)*(2.*pi/Lxyz(2))
        kz=kkz(ik)*(2.*pi/Lxyz(3))
        pi_over_Lx=pi/Lxyz(1)
      else
        kx0=kkx(ik)
        ky=kky(ik)
        kz=kkz(ik)
        pi_over_Lx=0.5
      endif
!
!  in the shearing sheet approximation, kx = kx0 - St*k_y.
!  Here, St=-deltay/Lx
!
      if (Sshear==0.) then
        kx=kx0
      else
        kx=kx0+ky*deltay/Lx
      endif
!
      if (headt.or.ip<5) print*, 'forcing_hel_both: kx0,kx,ky,kz=',kx0,kx,ky,kz
      k2=kx**2+ky**2+kz**2
      k=sqrt(k2)
!
! Find e-vector
!
      !
      ! Start with old method (not isotropic) for now.
      ! Pick e1 if kk not parallel to ee1. ee2 else.
      !
      if ((ky==0).and.(kz==0)) then
        ex=0; ey=1; ez=0
      else
        ex=1; ey=0; ez=0
      endif
      if (.not. old_forcing_evector) then
        !
        !  Isotropize ee in the plane perp. to kk by
        !  (1) constructing two basis vectors for the plane perpendicular
        !      to kk, and
        !  (2) choosing a random direction in that plane (angle phi)
        !  Need to do this in order for the forcing to be isotropic.
        !
        kk = (/kx, ky, kz/)
        ee = (/ex, ey, ez/)
        call cross(kk,ee,e1)
        call dot2(e1,norm); e1=e1/sqrt(norm) ! e1: unit vector perp. to kk
        call cross(kk,e1,e2)
        call dot2(e2,norm); e2=e2/sqrt(norm) ! e2: unit vector perp. to kk, e1
        call random_number_wrapper(phi,CHANNEL=channel_force); phi = phi*2*pi
        ee = cos(phi)*e1 + sin(phi)*e2
        ex=ee(1); ey=ee(2); ez=ee(3)
      endif
!
!  k.e
!
      call dot(kk,ee,kde)
!
!  k x e
!
      kex=ky*ez-kz*ey
      key=kz*ex-kx*ez
      kez=kx*ey-ky*ex
!
!  k x (k x e)
!
      kkex=ky*kez-kz*key
      kkey=kz*kex-kx*kez
      kkez=kx*key-ky*kex
!
!  ik x (k x e) + i*phase
!
!  Normalize ff; since we don't know dt yet, we finalize this
!  within timestep where dt is determined and broadcast.
!  For further details on normalization see forcing_hel subroutine.
!
      ffnorm=sqrt(1.+relhel**2) &
        *k*sqrt(k2-kde**2)/sqrt(kav*cs0**3)*(k/kav)**slope_ff
      if (ip<=9) then
        print*,'forcing_hel_both: k,kde,ffnorm,kav=',k,kde,ffnorm,kav
        print*,'forcing_hel_both: k*sqrt(k2-kde**2)=',k*sqrt(k2-kde**2)
      endif
!
!  need to multiply by dt (for Euler step), but it also needs to be
!  divided by sqrt(dt), because square of forcing is proportional
!  to a delta function of the time difference
!
      fact=force/ffnorm*sqrt(dt)
!
!  The wavevector is for the case where Lx=Ly=Lz=2pi. If that is not the
!  case one needs to scale by 2pi/Lx, etc.
!
      fx=exp(cmplx(0.,kx*k1_ff*x+phase))*fact
! only sines, to make sure that the force goes to zero at the equator
      fy=cmplx(0.,sin(ky*k1_ff*y))
      fz=exp(cmplx(0.,kz*k1_ff*z))
!
      if (ip<=5) then
        print*,'forcing_hel_both: fx=',fx
        print*,'forcing_hel_both: fy=',fy
        print*,'forcing_hel_both: fz=',fz
      endif
!
!  prefactor; treat real and imaginary parts separately (coef1 and coef2),
!  so they can be multiplied by different profiles below.
!
      coef1(1)=k*kex; coef2(1)=relhel*kkex
      coef1(2)=k*key; coef2(2)=relhel*kkey
      coef1(3)=k*kez; coef2(3)=relhel*kkez
      if (ip<=5) print*,'forcing_hel_both: coef=',coef1,coef2
!
!  loop the two cases separately, so we don't check for r_ff during
!  each loop cycle which could inhibit (pseudo-)vectorisation
!  calculate energy input from forcing; must use lout (not ldiagnos)
!
     rho1=1.
      do n=n1,n2
        do m=m1,m2
          variable_rhs=f(l1:l2,m,n,iffx:iffz)
          do j=1,3
            if (lactive_dimension(j)) then
              jf=j+ifff-1
              if (y(m)>equator)then
                forcing_rhs(:,j)=rho1*real(cmplx(coef1(j),coef2(j)) &
                                    *fx(l1:l2)*fy(m)*fz(n))&
                        *(1.-step(y(m),equator-ck_equator_gap,ck_gap_step)+&
                            step(y(m),equator+ck_equator_gap,ck_gap_step))
              else
                forcing_rhs(:,j)=rho1*real(cmplx(coef1(j),-coef2(j)) &
                  *fx(l1:l2)*fy(m)*fz(n))&
                  *(1.-step(y(m),equator-ck_equator_gap,ck_gap_step)+&
                      step(y(m),equator+ck_equator_gap,ck_gap_step))
              endif
              if (lhelical_test) then
                f(l1:l2,m,n,jf)=forcing_rhs(:,j)
              else
                f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+forcing_rhs(:,j)
              endif
              if (ltestfield_forcing) then
                jf=j+iaatest-1+3*(jtest_aa0-1)
                f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+forcing_rhs(:,j)
              endif
              if (ltestflow_forcing) then
                jf=j+iuutest-1+3*(jtest_uu0-1)
                f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+forcing_rhs(:,j)
              endif
            endif
          enddo
        enddo
      enddo
!
      if (ip<=9) print*,'forcing_hel_both: forcing OK'
!
    endsubroutine forcing_hel_both
!***********************************************************************
    subroutine forcing_chandra_kendall(f)
!
!  Add helical forcing function in spherical polar coordinate system.
!  25-jul-07/dhruba: adapted from forcing_hel
!
      use EquationOfState, only: cs0
      use General, only: random_number_wrapper
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      real, dimension(3) :: ee
      real, dimension(nx,3) :: capitalT,capitalS,capitalH,psi
      real, dimension(nx,3,3) :: psi_ij,Tij
      integer :: emm,l,j,jf,Legendrel,lmindex,aindex
      real :: a_ell,anum,adenom,jlm_ff,ylm_ff,rphase1,fnorm,alphar,Balpha,&
              psilm,RYlm,IYlm
      real :: rz,rindex,ralpha,rphase2
      real, dimension(mx) :: Z_psi
!
! This is designed for 5 emm values and for each one 5 ell values. Total 25 values.
!
      if (lhelical_test) then
        if (icklist==nlist_ck) &
          call fatal_error("forcing_chandra_kendall","CK testing: no more values in list")
        icklist=icklist+1
        lmindex=icklist
      else
        call random_number_wrapper(rindex,CHANNEL=channel_force)
        lmindex=nint(rindex*(nlist_ck-1))+1
      endif
!
      emm = cklist(lmindex,1)
      Legendrel = cklist(lmindex,2)
!
      call random_number_wrapper(ralpha,CHANNEL=channel_force)
      aindex=nint(ralpha*2)
      Balpha = cklist(lmindex,3+aindex)
!
! Now calculate the "potential" for the helical forcing. The expression
! is taken from Chandrasekhar and Kendall.
! Now construct Z_psi(r)
!
      call random_number_wrapper(rphase1,CHANNEL=channel_force)
      rphase1=rphase1*2*pi

      if (lfastCK) then
        do n=1,mz
          do m=1,my
            psilm = RYlm_list(m,n,lmindex)*cos(rphase1)- &
                    IYlm_list(m,n,lmindex)*sin(rphase1)
            psif(:,m,n) = psilm*Zpsi_list(:,lmindex,aindex+1)
            if (ck_equator_gap/=0) psif(:,m,n)=psif(:,m,n)*profy_ampl(m)
          enddo
        enddo

      else
        call sp_bessely_l(anum,Legendrel,Balpha*x(l1))
        call sp_besselj_l(adenom,Legendrel,Balpha*x(l1))
        a_ell = -anum/adenom
!        write(*,*) 'dhruba:',anum,adenom,Legendrel,Bessel_alpha,x(l1)
        do l=1,mx
          alphar=Balpha*x(l)
          call sp_besselj_l(jlm_ff,Legendrel,alphar)
          call sp_bessely_l(ylm_ff,Legendrel,alphar)
          Z_psi(l) = (a_ell*jlm_ff+ylm_ff)
        enddo
!
        do n=1,mz
          do m=1,my
            call sp_harm_real(RYlm,Legendrel,emm,y(m),z(n))
            call sp_harm_imag(IYlm,Legendrel,emm,y(m),z(n))
            psilm = RYlm*cos(rphase1)-IYlm*sin(rphase1)
            psif(:,m,n) = Z_psi*psilm
            if (ck_equator_gap/=0) psif(:,m,n)=psif(:,m,n)*profy_ampl(m)
          enddo
        enddo
      endif
!
! ----- Now calculate the force from the potential and add this to velocity
! get a random unit vector with three components ee_r, ee_theta, ee_phi
! psi at present is just Z_{ell}^m. We next do a sum over random coefficients
! get random psi.
! ----------now generate and add the force ------------
!
      call random_number_wrapper(rz,CHANNEL=channel_force)
      ee(3) = rz
      call random_number_wrapper(rphase2,CHANNEL=channel_force)
      rphase2 = pi*rphase2
      ee(1) = sqrt(1.-rz*rz)*cos(rphase2)
      ee(2) = sqrt(1.-rz*rz)*sin(rphase2)
      fnorm = fpre*cs0*cs0*sqrt(1./(cs0*Balpha))*sqrt(dt)
!     write(*,*) 'dhruba:',fnorm*sqrt(dt),dt,ee(1),ee(2),ee(3)

      do n=n1,n2
        do m=m1,m2
          psi(:,1) = psif(l1:l2,m,n)*ee(1)
          psi(:,2) = psif(l1:l2,m,n)*ee(2)
          psi(:,3) = psif(l1:l2,m,n)*ee(3)
          call gij_psi(psif,ee,psi_ij)
          call curl_mn(psi_ij,capitalT,psi)
          call gij_psi_etc(psif,ee,psi,psi_ij,Tij)
          call curl_mn(Tij,capitalS,capitalT)
          capitalS = float(helsign)*(1./Balpha)*capitalS
          if ((y(m)<pi/2.).or.(lsamesign)) then
            capitalH = capitalT + capitalS
          else
            capitalH = capitalT - capitalS
          endif
          do j=1,3
            jf = iuu+j-1
            if (r_ff /= 0.) capitalH(:,j)=profx_ampl*capitalH(:,j)
            if (lhelical_test) then
              if (lwrite_psi) then
                f(l1:l2,m,n,jf) = psif(l1:l2,m,n)
              else
                f(l1:l2,m,n,jf) = fnorm*capitalH(:,j)
              endif
            else
!
! stochastic euler scheme of integration [sqrt(dt) is already included in fnorm]
!
              f(l1:l2,m,n,jf) = f(l1:l2,m,n,jf) + fnorm*capitalH(:,j)
            endif
          enddo
        enddo
      enddo
!
    endsubroutine forcing_chandra_kendall
!***********************************************************************
    subroutine forcing_GP(f)
!
!  Add Galloway-Proctor forcing function.
!
!  24-jul-06/axel: coded
!
      use Mpicomm
      use Sub
      use DensityMethods, only: getrho
!
      real :: irufm
      real, dimension (nx) :: ruf,rho
      real, dimension (nx,3) :: variable_rhs,forcing_rhs,force_all
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx) :: cosx,sinx
      real :: cost,sint,cosym,sinym,fsum_tmp,fsum
      integer :: j,jf
      real :: fact
!
!  The default for omega_ff is now (12-jan-2018) changed from 1 to 0.
!
      if (ip<=6) print*,'forcing_GP: t=',t
      cost=cos(omega_ff*t)
      sint=sin(omega_ff*t)
!
!  Normalize ff; since we don't know dt yet, we finalize this
!  within timestep where dt is determined and broadcast.
!
!  need to multiply by dt (for Euler step), but it also needs to be
!  divided by sqrt(dt), because square of forcing is proportional
!  to a delta function of the time difference
!
      fact=sqrt(1.5)*force*sqrt(dt)
!
!  loop the two cases separately, so we don't check for r_ff during
!  each loop cycle which could inhibit (pseudo-)vectorisation
!  calculate energy input from forcing; must use lout (not ldiagnos)
!
      irufm=0
      do m=m1,m2
        cosx=cos(k1_ff*x+cost)
        sinx=sin(k1_ff*x+cost)
        cosym=cos(k1_ff*y(m)+sint)
        sinym=sin(k1_ff*y(m)+sint)
        forcing_rhs(:,1)=-fact*sinym
        forcing_rhs(:,2)=-fact*cosx(l1:l2)
        forcing_rhs(:,3)=+fact*(sinx(l1:l2)+cosym)
        do n=n1,n2
          variable_rhs=f(l1:l2,m,n,iffx:iffz)
          do j=1,3
            jf=j+ifff-1
            f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+forcing_rhs(:,j)
          enddo
          if (lout) then
            if (idiag_rufm/=0) then
              call getrho(f(:,m,n,ilnrho),rho)
              call multsv_mn(rho/dt,forcing_rhs,force_all)
              call dot_mn(variable_rhs,force_all,ruf)
              irufm=irufm+sum(ruf)
            endif
          endif
        enddo
      enddo
      !
      ! For printouts
      !
      if (lout) then
        if (idiag_rufm/=0) then
          irufm=irufm/(nwgrid)
          !
          !  on different processors, irufm needs to be communicated
          !  to other processors
          !
          fsum_tmp=irufm
          call mpireduce_sum(fsum_tmp,fsum)
          irufm=fsum
          call mpibcast_real(irufm)
          !
          fname(idiag_rufm)=irufm
          itype_name(idiag_rufm)=ilabel_sum
        endif
      endif
!
      if (ip<=9) print*,'forcing_GP: forcing OK'
!
    endsubroutine forcing_GP
!***********************************************************************
    subroutine forcing_GP92(f)
!
!  Add Galloway-Proctor (1992) forcing function.
!
!  23-mar-20/axel: adapted from forcing_GP92
!
      use Mpicomm
      use Sub
      use DensityMethods, only: getrho
!
      real :: irufm
      real, dimension (nx) :: ruf,rho
      real, dimension (nx,3) :: variable_rhs,forcing_rhs,force_all
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx) :: cosx,sinx
      real :: cost,sint,cosym,sinym,fsum_tmp,fsum
      integer :: j,jf
      real :: fact
!
!  The default for omega_ff is now (12-jan-2018) changed from 1 to 0.
!
      if (ip<=6) print*,'forcing_GP: t=',t
      cost=cos(omega_ff*t)
      sint=sin(omega_ff*t)
!
!  Normalize ff; since we don't know dt yet, we finalize this
!  within timestep where dt is determined and broadcast.
!
!  need to multiply by dt (for Euler step), but it also needs to be
!  divided by sqrt(dt), because square of forcing is proportional
!  to a delta function of the time difference
!
      fact=sqrt(1.5)*force*sqrt(dt)
!
!  loop the two cases separately, so we don't check for r_ff during
!  each loop cycle which could inhibit (pseudo-)vectorisation
!  calculate energy input from forcing; must use lout (not ldiagnos)
!
      irufm=0
      do m=m1,m2
        cosx=cos(k1_ff*x+cost)
        sinx=sin(k1_ff*x+cost)
        cosym=cos(k1_ff*y(m)+sint)
        sinym=sin(k1_ff*y(m)+sint)
        forcing_rhs(:,1)=-fact*sinym
        forcing_rhs(:,2)=-fact*cosx(l1:l2)
        forcing_rhs(:,3)=+fact*(sinx(l1:l2)+cosym)
        do n=n1,n2
          variable_rhs=f(l1:l2,m,n,iffx:iffz)
          do j=1,3
            jf=j+ifff-1
            f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+forcing_rhs(:,j)
          enddo
          if (lout) then
            if (idiag_rufm/=0) then
              call getrho(f(:,m,n,ilnrho),rho)
              call multsv_mn(rho/dt,forcing_rhs,force_all)
              call dot_mn(variable_rhs,force_all,ruf)
              irufm=irufm+sum(ruf)
            endif
          endif
        enddo
      enddo
      !
      ! For printouts
      !
      if (lout) then
        if (idiag_rufm/=0) then
          irufm=irufm/(nwgrid)
          !
          !  on different processors, irufm needs to be communicated
          !  to other processors
          !
          fsum_tmp=irufm
          call mpireduce_sum(fsum_tmp,fsum)
          irufm=fsum
          call mpibcast_real(irufm)
          !
          fname(idiag_rufm)=irufm
          itype_name(idiag_rufm)=ilabel_sum
        endif
      endif
!
      if (ip<=9) print*,'forcing_GP: forcing OK'
!
    endsubroutine forcing_GP92
!***********************************************************************
    subroutine forcing_TG(f)
!
!  Add Taylor-Green forcing function.
!
!   9-oct-04/axel: coded
!
      use Mpicomm
      use Sub
      use DensityMethods, only: getrho
!
      real :: irufm,fsum_tmp,fsum
      real, dimension (nx) :: ruf,rho
      real, dimension (nx,3) :: variable_rhs,forcing_rhs,force_all
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx), save :: sinx,cosx
      real, dimension (my), save :: siny,cosy
      real, dimension (mz), save :: cosz
      logical, save :: lfirst_call=.true.
      integer :: j,jf
      real :: fact
!
      if (lfirst_call) then
        if (lroot) print*,'forcing_TG: calculate sinx,cosx,siny,cosy,cosz'
        sinx=sin(k1_ff*x)
        cosx=cos(k1_ff*x)
        siny=sin(k1_ff*y)
        cosy=cos(k1_ff*y)
        cosz=cos(k1_ff*z)
        lfirst_call=.false.
      endif
!
      if (ip<=6) print*,'forcing_TG: dt, lfirst_call=',dt,lfirst_call
!
!  Normalize ff; since we don't know dt yet, we finalize this
!  within timestep where dt is determined and broadcast.
!
!  need to multiply by dt (for Euler step), but it also needs to be
!  divided by sqrt(dt), because square of forcing is proportional
!  to a delta function of the time difference
!
      fact=2*force*sqrt(dt)
!
!  loop the two cases separately, so we don't check for r_ff during
!  each loop cycle which could inhibit (pseudo-)vectorisation
!  calculate energy input from forcing; must use lout (not ldiagnos)
!
      irufm=0
      do n=n1,n2
        do m=m1,m2
          variable_rhs=f(l1:l2,m,n,iffx:iffz)
          forcing_rhs(:,1)=+fact*sinx(l1:l2)*cosy(m)*cosz(n)
          forcing_rhs(:,2)=-fact*cosx(l1:l2)*siny(m)*cosz(n)
          forcing_rhs(:,3)=0.
          do j=1,3
            if (lactive_dimension(j)) then
              jf=j+ifff-1
              f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+forcing_rhs(:,j)
            endif
          enddo
          if (lout) then
            if (idiag_rufm/=0) then
              call getrho(f(:,m,n,ilnrho),rho)
              call multsv_mn(rho/dt,forcing_rhs,force_all)
              call dot_mn(variable_rhs,force_all,ruf)
              irufm=irufm+sum(ruf)
            endif
          endif
        enddo
      enddo
      !
      ! For printouts
      !
      if (lout) then
        if (idiag_rufm/=0) then
          irufm=irufm/(nwgrid)
          !
          !  on different processors, irufm needs to be communicated
          !  to other processors
          !
          fsum_tmp=irufm
          call mpireduce_sum(fsum_tmp,fsum)
          irufm=fsum
          call mpibcast_real(irufm)
          !
          fname(idiag_rufm)=irufm
          itype_name(idiag_rufm)=ilabel_sum
        endif
      endif
!
      if (ip<=9) print*,'forcing_TG: forcing OK'
!
    endsubroutine forcing_TG
!***********************************************************************
    subroutine forcing_ABC(f)
!
!  Added ABC forcing function
!
!  17-jul-06/axel: coded
!
      use Diagnostics
      use DensityMethods, only: getrho
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      real :: irufm,fsum_tmp,fsum
      real, dimension (nx) :: ruf,rho
      real, dimension (nx,3) :: variable_rhs,forcing_rhs,force_all,bb,fxb
      real, dimension (mx), save :: sinx,cosx
      real, dimension (my), save :: siny,cosy
      real, dimension (mz), save :: sinz,cosz
      logical, save :: lfirst_call=.true.
      integer :: j,jf
      real :: fact
!
!  at the first step, the sin and cos functions are calculated for all
!  x,y,z points and are then saved and used for all subsequent steps
!  and pencils
!
      if (ip<=6) print*,'forcing_ABC: lfirst_call=',lfirst_call
      if (lfirst_call) then
        if (lroot) print*,'forcing_ABC: calculate sinx,cosx,siny,cosy,sinz,cosz'
        sinx=sin(k1_ff*x); cosx=cos(k1_ff*x)
        siny=sin(k1_ff*y); cosy=cos(k1_ff*y)
        sinz=sin(k1_ff*z); cosz=cos(k1_ff*z)
        lfirst_call=.false.
      endif
      if (ip<=6) print*,'forcing_ABC: dt, lfirst_call=',dt,lfirst_call
!
!  Normalize ff; since we don't know dt yet, we finalize this
!  within timestep where dt is determined and broadcast.
!
!  need to multiply by dt (for Euler step), but it also needs to be
!  divided by sqrt(dt), because square of forcing is proportional
!  to a delta function of the time difference
!
      fact=2*force*sqrt(dt)
!
!  loop the two cases separately, so we don't check for r_ff during
!  each loop cycle which could inhibit (pseudo-)vectorisation
!  calculate energy input from forcing; must use lout (not ldiagnos)
!
      irufm=0
      do n=n1,n2
        do m=m1,m2
          variable_rhs=f(l1:l2,m,n,iffx:iffz)
          forcing_rhs(:,1)=fact*(sinz(n    )+cosy(m)    )
          forcing_rhs(:,2)=fact*(sinx(l1:l2)+cosz(n)    )
          forcing_rhs(:,3)=fact*(siny(m    )+cosx(l1:l2))
          do j=1,3
            jf=j+ifff-1
            f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+forcing_rhs(:,j)
          enddo
          if (lout) then
            if (idiag_rufm/=0) then
              call getrho(f(:,m,n,ilnrho),rho)
              call multsv_mn(rho/dt,forcing_rhs,force_all)
              call dot_mn(variable_rhs,force_all,ruf)
              irufm=irufm+sum(ruf)
            endif
          endif
        enddo
      enddo
      !
      ! For printouts
      !
      if (lout) then
        if (idiag_rufm/=0) then
          irufm=irufm/(nwgrid)
          !
          !  on different processors, irufm needs to be communicated
          !  to other processors
          !
          fsum_tmp=irufm
          call mpireduce_sum(fsum_tmp,fsum)
          irufm=fsum
          call mpibcast_real(irufm)
          !
          fname(idiag_rufm)=irufm
          itype_name(idiag_rufm)=ilabel_sum
        endif
        if (lmagnetic) then
          if (idiag_fxbxm/=0.or.idiag_fxbym/=0.or.idiag_fxbzm/=0) then
            call curl(f,iaa,bb)
            call cross(forcing_rhs,bb,fxb)
            call sum_mn_name(fxb(:,1),idiag_fxbxm)
            call sum_mn_name(fxb(:,2),idiag_fxbym)
            call sum_mn_name(fxb(:,3),idiag_fxbzm)
          endif
        endif
      endif
!
      if (ip<=9) print*,'forcing_ABC: forcing OK'
!
    endsubroutine forcing_ABC
!***********************************************************************
    subroutine forcing_mhd_mode(f)
!
!  Added ABC forcing function
!
!  17-jul-06/axel: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      real :: fact
      integer :: m,n
!
      fact=force*sqrt(dt)
      do n=n1,n2
        do m=m1,m2
          f(:,m,n,iuy)=f(:,m,n,iuy)+fact*sin(k1_ff*x)
          f(:,m,n,iAy)=f(:,m,n,iAy)+fact*sin(k1_ff*x)
        enddo
      enddo
      !
      if (ip<=9) print*,'forcing_mhd: forcing OK'
!
    endsubroutine forcing_mhd_mode
!***********************************************************************
    subroutine forcing_nocos(f)
!
!  Add no-cosine forcing function.
!
!  27-oct-04/axel: coded
!
      use DensityMethods, only: getrho
      use Mpicomm
      use Sub
!
      real :: irufm,fsum_tmp,fsum
      real, dimension (nx) :: ruf,rho
      real, dimension (nx,3) :: variable_rhs,forcing_rhs,force_all
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx), save :: sinx
      real, dimension (my), save :: siny
      real, dimension (mz), save :: sinz
      logical, save :: lfirst_call=.true.
      integer :: j,jf
      real :: fact
!
      if (lfirst_call) then
        if (lroot) print*,'forcing_nocos: calculate sinx,siny,sinz'
        sinx=sin(k1_ff*x)
        siny=sin(k1_ff*y)
        sinz=sin(k1_ff*z)
        lfirst_call=.false.
      endif
!
      if (ip<=6) print*,'forcing_hel_noshear: dt, lfirst_call=',dt,lfirst_call
!
!  Normalize ff; since we don't know dt yet, we finalize this
!  within timestep where dt is determined and broadcast.
!
!  need to multiply by dt (for Euler step), but it also needs to be
!  divided by sqrt(dt), because square of forcing is proportional
!  to a delta function of the time difference
!
      fact=force*sqrt(dt)
!
!  loop the two cases separately, so we don't check for r_ff during
!  each loop cycle which could inhibit (pseudo-)vectorisation
!  calculate energy input from forcing; must use lout (not ldiagnos)
!
      irufm=0
      do n=n1,n2
        do m=m1,m2
          variable_rhs=f(l1:l2,m,n,iffx:iffz)
          forcing_rhs(:,1)=fact*sinz(n)
          forcing_rhs(:,2)=fact*sinx(l1:l2)
          forcing_rhs(:,3)=fact*siny(m)
          do j=1,3
            if (lactive_dimension(j)) then
              jf=j+ifff-1
              f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+forcing_rhs(:,j)
            endif
          enddo
          if (lout) then
            if (idiag_rufm/=0) then
              call getrho(f(:,m,n,ilnrho),rho)
              call multsv_mn(rho/dt,forcing_rhs,force_all)
              call dot_mn(variable_rhs,force_all,ruf)
              irufm=irufm+sum(ruf)
            endif
          endif
        enddo
      enddo
      !
      ! For printouts
      !
      if (lout) then
        if (idiag_rufm/=0) then
          irufm=irufm/(nwgrid)
          !
          !  on different processors, irufm needs to be communicated
          !  to other processors
          !
          fsum_tmp=irufm
          call mpireduce_sum(fsum_tmp,fsum)
          irufm=fsum
          call mpibcast_real(irufm)
          !
          fname(idiag_rufm)=irufm
          itype_name(idiag_rufm)=ilabel_sum
        endif
      endif
!
      if (ip<=9) print*,'forcing_nocos: forcing OK'
!
    endsubroutine forcing_nocos
!***********************************************************************
    subroutine forcing_gaussianpot(f,force_ampl)
!
!  gradient of gaussians as forcing function
!
!  19-dec-05/tony: coded, adapted from forcing_nocos
!  14-jul-10/axel: in less then 3-D, project forcing to computational domain
!
      use DensityMethods, only: getrho
      use EquationOfState, only: cs0
      use General, only: random_number_wrapper
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: force_ampl, force_tmp
!
      real, dimension (3) :: fran
      real, dimension (nx) :: radius2, gaussian, gaussian_fact, ruf, rho
      real, dimension (nx,3) :: variable_rhs,force_all,delta
      integer :: j, jf, ilocation
      real :: irufm,fact,width_ff21,fsum_tmp,fsum
!
!  check length of time step
!
      if (ip<=6) print*,'forcing_gaussianpot: dt=',dt
!
!  generate random numbers
!
      if (t>tsforce) then
        if (lrandom_location) then
          call random_number_wrapper(fran,CHANNEL=channel_force)
          location=fran*Lxyz+xyz0
        else
          location=location_fixed(:,1)
        endif
!
!  It turns out that in the presence of shear, and even for weak shear,
!  vorticity is being produced. In order to check whether the shearing
!  periodic boundaries are to blame, we can cut the y extent of forcing
!  locations by half.
!
        if (lforce_cuty) location(2)=location(2)*.5
!
!  reset location(i) to x, y, or z, so forcing does still occur,
!  but it is restricted to the plane or line that is being solved.
!
        if (.not.lactive_dimension(1)) location(1)=x(1)
        if (.not.lactive_dimension(2)) location(2)=y(1)
        if (.not.lactive_dimension(3)) location(3)=z(1)
!
!  write location to file; this is off by default.
!
        if (lroot .and. lwrite_gausspot_to_file) then
          open(1,file=trim(datadir)//'/gaussian_pot_forcing.dat', &
              status='unknown',position='append')
          write(1,'(4f14.7)') t, location
          close (1)
        endif
!
!  set next forcing time
!
        tsforce=t+dtforce
        if (ip<=6) print*,'forcing_gaussianpot: location=',location
      endif
!
!  Set forcing amplitude (same value for each location by default).
!  For iforce_tprofile=='sin^2', we have a temporal sine profile.
!  By default, iforce_tprofile='nothing'.
!
      if (iforce_tprofile=='nothing') then
        force_tmp=force_ampl
      elseif (iforce_tprofile=='sin^2') then
        force_tmp=force_ampl*sin(pi*(tsforce-t)/dtforce)**2
      else
        call  fatal_error('forcing_gaussianpot','iforce_tprofile not good')
      endif
!
!  Possibility of outputting data at each time step (for testing).
!
      if (lroot .and. lwrite_gausspot_to_file_always) then
        open(1,file=trim(datadir)//'/gaussian_pot_forcing_always.dat', &
            status='unknown',position='append')
        write(1,'(5f14.7)') t, location, force_tmp
        close (1)
      endif
!
!  Let explosion last dtforce_duration or, by default, until next explosion.
!  By default, dtforce_duration=1 (this is a special value, different from zero),
!  but when t>dtforce, forcing is off again, regardless of dtforce_duration.
!
      if ( (dtforce_duration<0.0) .or. &
           (t-(tsforce-dtforce))<=dtforce_duration ) then
!
!  Normalize ff; since we don't know dt yet, we finalize this
!  within timestep where dt is determined and broadcast.
!
!  We multiply the forcing term by dt and add to the right-hand side
!  of the momentum equation for an Euler step, but it also needs to be
!  divided by sqrt(dt), because square of forcing is proportional
!  to a delta function of the time difference.
!  When dtforce is finite, take dtforce+.5*dt.
!  The 1/2 factor takes care of round-off errors.
!  Also define width_ff21 = 1/width^2.
!  The factor 2 in fact takes care of the 2 in 2*xi/R^2.
!
        width_ff21=1./width_ff**2
        fact=2.*width_ff21*force_tmp*dt*sqrt(cs0*width_ff/max(dtforce+.5*dt,dt))
!
!  Loop the two cases separately, so we don't check for r_ff during
!  each loop cycle which could inhibit (pseudo-)vectorization.
!  Calculate energy input from forcing; must use lout (not ldiagnos).
!
        irufm=0
!
!  loop over all pencils
!
        do n=n1,n2
          do m=m1,m2
!
!  loop over multiple locations
!
            do ilocation=1,nlocation
!
!  Obtain distance (=delta) to center of blob
!
              if (ilocation==1) then
                delta(:,1)=x(l1:l2)-location(1)
                delta(:,2)=y(m)-location(2)
                delta(:,3)=z(n)-location(3)
              else
                delta(:,1)=x(l1:l2)-location2(1)
                delta(:,2)=y(m)-location2(2)
                delta(:,3)=z(n)-location2(3)
              endif
              do j=1,3
                if (lperi(j)) then
                  if (lforce_peri) then
                    if (j==2) then
                      delta(:,2)=2*atan(tan(.5*(delta(:,2) &
                          +2.*deltay*atan(1000.*tan(.25* &
                          (pi+x(l1:l2)-location(1)))))))
                    else
                      delta(:,j)=2*atan(tan(.5*delta(:,j)))
                    endif
                  else
                    where (delta(:,j) >  Lxyz(j)/2.) delta(:,j)=delta(:,j)-Lxyz(j)
                    where (delta(:,j) < -Lxyz(j)/2.) delta(:,j)=delta(:,j)+Lxyz(j)
                  endif
              endif
              if (.not.lactive_dimension(j)) delta(:,j)=0.
              enddo
!
!  compute gaussian blob and set to forcing function
!
              radius2=delta(:,1)**2+delta(:,2)**2+delta(:,3)**2
              gaussian=exp(-radius2*width_ff21)
              gaussian_fact=gaussian*fact
              variable_rhs=f(l1:l2,m,n,iffx:iffz)
              if (iphiuu==0) then
                do j=1,3
                  if (lactive_dimension(j)) then
                    jf=j+ifff-1
                    if (iforce_profile=='nothing') then
                      f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+gaussian_fact*delta(:,j)
                    else
                      f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+gaussian_fact*delta(:,j) &
                        *profx_ampl*profy_ampl(m)*profz_ampl(n)
                    endif
                  endif
                enddo
              else
!
!  add possibility of modulation (but only if iphiuu/=0)
!
                if (iforce_profile=='nothing') then
                  f(l1:l2,m,n,ifff)=f(l1:l2,m,n,ifff)+gaussian
                else
                  f(l1:l2,m,n,ifff)=f(l1:l2,m,n,ifff)+gaussian*profx_ampl*profy_ampl(m)*profz_ampl(n)
                endif
              endif
!
!  test
!
!--           if (icc/=0) f(l1:l2,m,n,icc)=f(l1:l2,m,n,icc)+gaussian
!
!  diagnostics (currently without this possibility of a modulation applied above)
!
              if (lout) then
                if (idiag_rufm/=0) then
                  call getrho(f(:,m,n,ilnrho),rho)
                  call multsv_mn(rho/dt,spread(gaussian,2,3)*delta,force_all)
                  call dot_mn(variable_rhs,force_all,ruf)
                  irufm=irufm+sum(ruf)
                endif
              endif
            enddo
          enddo
        enddo
      endif
!
!  For printouts
!
      if (lout) then
        if (idiag_rufm/=0) then
          irufm=irufm/(nwgrid)
!
!  on different processors, irufm needs to be communicated
!  to other processors
!
          fsum_tmp=irufm
          call mpireduce_sum(fsum_tmp,fsum)
          irufm=fsum
          call mpibcast_real(irufm)
!
          fname(idiag_rufm)=irufm
          itype_name(idiag_rufm)=ilabel_sum
        endif
      endif
!
      if (ip<=9) print*,'forcing_gaussianpot: forcing OK'
!
    endsubroutine forcing_gaussianpot
!***********************************************************************
    subroutine forcing_hillrain(f,force_ampl)
!
!  gradient of gaussians as forcing function
!
!  29-sep-15/axel: adapted from forcing_gaussianpot
!
      use DensityMethods, only: getrho
      use EquationOfState, only: cs0
      use General, only: random_number_wrapper
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: force_ampl, force_tmp
!
      real, dimension (3) :: fran
      real, dimension (nx) :: r, r2, r3, r5, pom2, ruf, rho
      real, dimension (nx,3) :: variable_rhs, force_all, delta
      integer :: j,jf
      real :: irufm, fact, fsum_tmp, fsum
      real :: a_hill, a2_hill, a3_hill
!
!  check length of time step
!
      if (ip<=6) print*,'forcing_hillrain: dt=',dt
!
!  generate random numbers
!
      if (t>tsforce) then
        if (lrandom_location) then
          call random_number_wrapper(fran,CHANNEL=channel_force)
          location=fran*Lxyz+xyz0
        else
          location=location_fixed(:,1)
        endif
!
!  It turns out that in the presence of shear, and even for weak shear,
!  vorticitity is being produced. In order to check whether the shearing
!  periodic boundaries are to blame, we can cut the y extent of forcing
!  locations by half.
!
        if (lforce_cuty) location(2)=location(2)*.5
!
!  reset location(i) to x or y, and keep location(3)=0 fixed
!
        if (.not.lactive_dimension(1)) location(1)=x(1)
        if (.not.lactive_dimension(2)) location(2)=y(1)
        location(3)=0.
!
!  write location to file
!
        if (lroot .and. lwrite_gausspot_to_file) then
          open(1,file=trim(datadir)//'/hillrain_forcing.dat', &
              status='unknown',position='append')
          write(1,'(4f14.7)') t, location
          close (1)
        endif
!
!  set next forcing time
!
        tsforce=t+dtforce
        if (ip<=6) print*,'forcing_hillrain: location=',location
      endif
!
!  Set forcing amplitude (same value for each location by default)
!  We multiply the forcing term by dt and add to the right-hand side
!  of the momentum equation for an Euler step, but it also needs to be
!  divided by sqrt(dt), because square of forcing is proportional
!  to a delta function of the time difference.
!  When dtforce is finite, take dtforce+.5*dt.
!  The 1/2 factor takes care of round-off errors.
!
!
      if (iforce_tprofile=='nothing') then
        force_tmp=force_ampl
        fact=force_tmp*dt*sqrt(cs0*radius_ff/max(dtforce+.5*dt,dt))
      elseif (iforce_tprofile=='sin^2') then
        force_tmp=force_ampl*sin(pi*(tsforce-t)/dtforce)**2
        fact=force_tmp*dt*sqrt(cs0*radius_ff/max(dtforce+.5*dt,dt))
      elseif (iforce_tprofile=='delta') then
        if (tsforce-t <= dt) then
          force_tmp=force_ampl
        else
          force_tmp=0.
        endif
        fact=force_tmp
      else
        call  fatal_error('forcing_hillrain','iforce_tprofile not good')
      endif
!
!  Let explosion last dtforce_duration or, by default, until next explosion.
!
      if ( force_tmp /= 0. .and. ((dtforce_duration<0.0) .or. &
           (t-(tsforce-dtforce))<=dtforce_duration) ) then
!
!  loop the two cases separately, so we don't check for r_ff during
!  each loop cycle which could inhibit (pseudo-)vectorisation
!  calculate energy input from forcing; must use lout (not ldiagnos)
!
        irufm=0
!
!  loop over all pencils
!
        do n=n1,n2; do m=m1,m2
!
!  Obtain distance to center of blob
!
          delta(:,1)=x(l1:l2)-location(1)
          delta(:,2)=y(m)-location(2)
          delta(:,3)=z(n)-location(3)
          do j=1,3
            if (lperi(j)) then
              if (lforce_peri) then
                if (j==2) then
                  delta(:,2)=2*atan(tan(.5*(delta(:,2) &
                      +2.*deltay*atan(1000.*tan(.25* &
                      (pi+x(l1:l2)-location(1)))))))
                else
                  delta(:,j)=2*atan(tan(.5*delta(:,j)))
                endif
              else
                where (delta(:,j) >  Lxyz(j)/2.) delta(:,j)=delta(:,j)-Lxyz(j)
                where (delta(:,j) < -Lxyz(j)/2.) delta(:,j)=delta(:,j)+Lxyz(j)
              endif
            endif
            if (.not.lactive_dimension(j)) delta(:,j)=0.
          enddo
!
!  compute factors for Hill vortex
!
          r2=delta(:,1)**2+delta(:,2)**2+delta(:,3)**2
          pom2=delta(:,1)**2+delta(:,2)**2
          r=sqrt(r2)
          r3=r2*r
          r5=r2*r3
          a_hill=radius_ff
          a2_hill=a_hill**2
          a3_hill=a_hill*a2_hill
!
!  compute Hill vortex
!
          where (r<=a_hill)
            variable_rhs(:,1)=-1.5*delta(:,1)*delta(:,3)/a2_hill
            variable_rhs(:,2)=-1.5*delta(:,2)*delta(:,3)/a2_hill
            variable_rhs(:,3)=-2.5+1.5*(pom2+r2)/a2_hill
          elsewhere
            variable_rhs(:,1)=-1.5*delta(:,1)*delta(:,3)*a3_hill/r5
            variable_rhs(:,2)=-1.5*delta(:,2)*delta(:,3)*a3_hill/r5
            variable_rhs(:,3)=-a3_hill/r3+1.5*pom2*a3_hill/r5
          endwhere
!
!  add to the relevant variable
!
          do j=1,3
            if (lactive_dimension(j)) then
              jf=j+ifff-1
              f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+fact*variable_rhs(:,j)
            endif
          enddo
!
!  passive scalar as a possible test
!
!--         if (icc/=0) f(l1:l2,m,n,icc)=f(l1:l2,m,n,icc)+gaussian
!
!  diagnostics
!
          if (lout) then
            if (idiag_rufm/=0) then
              call getrho(f(:,m,n,ilnrho),rho)
              call multsv_mn(rho/dt,f(l1:l2,m,n,iux:iuz),force_all)
              call dot_mn(variable_rhs,force_all,ruf)
              irufm=irufm+sum(ruf)
            endif
          endif
!
!  For printouts
!
          if (lout) then
            if (idiag_rufm/=0) then
              irufm=irufm/(nwgrid)
!
!  on different processors, irufm needs to be communicated
!  to other processors
!
              fsum_tmp=irufm
              call mpireduce_sum(fsum_tmp,fsum)
              irufm=fsum
              call mpibcast_real(irufm)
!
              fname(idiag_rufm)=irufm
              itype_name(idiag_rufm)=ilabel_sum
            endif
          endif
        enddo; enddo
      endif
!
      if (ip<=9) print*,'forcing_hillrain: forcing OK'
!
    endsubroutine forcing_hillrain
!***********************************************************************
    subroutine forcing_white_noise(f)
!
!  gradient of gaussians as forcing function
!
!  19-dec-13/axel: added
!
      use DensityMethods, only: getrho
      use EquationOfState, only: cs0
      use General, only: random_number_wrapper
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,fsum_tmp,fsum
!
      real, dimension (nx) :: r,p,tmp,rho,ruf
      real, dimension (nx,3) :: force_all,variable_rhs,forcing_rhs
      integer :: j,jf
      real :: irufm
!
!  check length of time step
!
      if (ip<=6) print*,'forcing_white_noise: dt=',dt
!
!  extent
!
      if (.not.lactive_dimension(1)) location(1)=x(1)
      if (.not.lactive_dimension(2)) location(2)=y(1)
      if (.not.lactive_dimension(3)) location(3)=z(1)
      if (ip<=6) print*,'forcing_white_noise: location=',location
!
!  Normalize ff; since we don't know dt yet, we finalize this
!  within timestep where dt is determined and broadcast.
!
!  We multiply the forcing term by dt and add to the right-hand side
!  of the momentum equation for an Euler step, but it also needs to be
!  divided by sqrt(dt), because square of forcing is proportional
!  to a delta function of the time difference.
!  When dtforce is finite, take dtforce+.5*dt.
!  The 1/2 factor takes care of round-off errors.
!  Also define width_ff21 = 1/width^2
!
      ampl=force*sqrt(dt*cs0)*cs0
!
!  loop the two cases separately, so we don't check for r_ff during
!  each loop cycle which could inhibit (pseudo-)vectorisation
!  calculate energy input from forcing; must use lout (not ldiagnos)
!
      irufm=0
!
!  loop over all pencils
!
      do n=n1,n2
        do m=m1,m2
          do j=1,3
            if (lactive_dimension(j)) then
              jf=j+ifff-1
              if (modulo(j-1,2)==0) then
                call random_number_wrapper(r,CHANNEL=channel_force)
                call random_number_wrapper(p,CHANNEL=channel_force)
                tmp=sqrt(-2*log(r))*sin(2*pi*p)
              else
                tmp=sqrt(-2*log(r))*cos(2*pi*p)
              endif
              f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+ampl*tmp
              forcing_rhs(:,j)=ampl*tmp
            endif
          enddo
!
!  diagnostics
!
          if (lout) then
            if (idiag_rufm/=0) then
              variable_rhs=f(l1:l2,m,n,iffx:iffz)
              call getrho(f(:,m,n,ilnrho),rho)
              call multsv_mn(rho/dt,forcing_rhs,force_all)
              call dot_mn(variable_rhs,force_all,ruf)
              irufm=irufm+sum(ruf)
            endif
          endif
        enddo
      enddo
!
!  For printouts
!
      if (lout) then
        if (idiag_rufm/=0) then
          irufm=irufm/(nwgrid)
!
!  on different processors, irufm needs to be communicated
!  to other processors
!
          fsum_tmp=irufm
          call mpireduce_sum(fsum_tmp,fsum)
          irufm=fsum
          call mpibcast_real(irufm)
!
          fname(idiag_rufm)=irufm
          itype_name(idiag_rufm)=ilabel_sum
        endif
      endif
!
      if (ip<=9) print*,'forcing_white_noise: forcing OK'
!
    endsubroutine forcing_white_noise
!***********************************************************************
    function calc_force_ampl(f,fx,fy,fz,coef) result(force_ampl)
!
!  calculates the coefficient for a forcing that satisfies
!  <rho*u*f> = constant.
!
!   7-sep-02/axel: coded
!
      use DensityMethods, only: getrho
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3) :: uu
      real, dimension (nx) :: rho,udotf
      complex, dimension (mx) :: fx
      complex, dimension (my) :: fy
      complex, dimension (mz) :: fz
      complex, dimension (3) :: coef
      real :: rho_uu_ff,force_ampl,fsum_tmp,fsum
      integer :: j,m,n
!
      rho_uu_ff=0.
      do n=n1,n2
        do m=m1,m2
          call getrho(f(:,m,n,ilnrho),rho)
          uu=f(l1:l2,m,n,iffx:iffz)
          udotf=0.
          do j=1,3
            udotf=udotf+uu(:,j)*real(coef(j)*fx(l1:l2)*fy(m)*fz(n))
          enddo
          rho_uu_ff=rho_uu_ff+sum(rho*udotf)
        enddo
      enddo
!
!  on different processors, this result needs to be communicated
!  to other processors
!
      fsum_tmp=rho_uu_ff
      call mpireduce_sum(fsum_tmp,fsum)
      if (lroot) rho_uu_ff=fsum/nwgrid
!      if (lroot) rho_uu_ff=rho_uu_ff/nwgrid
      call mpibcast_real(rho_uu_ff)
!
!  scale forcing function
!  but do this only when rho_uu_ff>0.; never allow it to change sign
!
        if (headt) print*,'calc_force_ampl: divide forcing function by rho_uu_ff=',rho_uu_ff
        !      force_ampl=work_ff/(.1+max(0.,rho_uu_ff))
        force_ampl=work_ff/rho_uu_ff
        if (force_ampl > max_force) force_ampl=max_force
        if (force_ampl < -max_force) force_ampl=-max_force
!
    endfunction calc_force_ampl
!***********************************************************************
    subroutine forcing_hel_noshear(f)
!
!  add helical forcing function, using a set of precomputed wavevectors
!
!  10-apr-00/axel: coded
!  06-dec-13/nishant: made kkx etc allocatable
!
      use EquationOfState, only: cs0
      use General, only: random_number_wrapper
      use Mpicomm
      use Sub
!
      real :: phase,ffnorm
      real, dimension (2) :: fran
      real, dimension (nx) :: radius,tmpx
!      real, dimension (mz) :: tmpz
      real, dimension (mx,my,mz,mfarray) :: f
      complex, dimension (mx) :: fx
      complex, dimension (my) :: fy
      complex, dimension (mz) :: fz
      complex, dimension (3) :: coef
      integer :: ik,j,jf,kx,ky,kz,kex,key,kez,kkex,kkey,kkez
      real :: k2,k,ex,ey,ez,kde,sig,fact
      real, dimension(3) :: e1,e2,ee,kk
      real :: norm,phi
!
!  generate random coefficients -1 < fran < 1
!  ff=force*Re(exp(i(kx+phase)))
!  |k_i| < akmax
!
      call random_number_wrapper(fran,CHANNEL=channel_force)
      phase=pi*(2*fran(1)-1.)
      ik=nk*.9999*fran(2)+1
      if (ip<=6) print*,'force_hel_noshear: ik,phase,kk=',ik,phase,kkx(ik),kky(ik),kkz(ik),dt
!
      kx=kkx(ik)
      ky=kky(ik)
      kz=kkz(ik)
      if (ip<=4) print*, 'force_hel_noshear: kx,ky,kz=',kx,ky,kz
!
      k2=float(kx**2+ky**2+kz**2)
      k=sqrt(k2)
!
! Find e-vector
!
      !
      ! Start with old method (not isotropic) for now.
      ! Pick e1 if kk not parallel to ee1. ee2 else.
      !
      if ((ky==0).and.(kz==0)) then
        ex=0; ey=1; ez=0
      else
        ex=1; ey=0; ez=0
      endif
      if (.not. old_forcing_evector) then
        !
        !  Isotropize ee in the plane perp. to kk by
        !  (1) constructing two basis vectors for the plane perpendicular
        !      to kk, and
        !  (2) choosing a random direction in that plane (angle phi)
        !  Need to do this in order for the forcing to be isotropic.
        !
        kk = (/kx, ky, kz/)
        ee = (/ex, ey, ez/)
        call cross(kk,ee,e1)
        call dot2(e1,norm); e1=e1/sqrt(norm) ! e1: unit vector perp. to kk
        call cross(kk,e1,e2)
        call dot2(e2,norm); e2=e2/sqrt(norm) ! e2: unit vector perp. to kk, e1
        call random_number_wrapper(phi,CHANNEL=channel_force); phi = phi*2*pi
        ee = cos(phi)*e1 + sin(phi)*e2
        ex=ee(1); ey=ee(2); ez=ee(3)
      endif
!
!  k.e
!
      call dot(kk,ee,kde)
!
!  k x e
!
      kex=ky*ez-kz*ey
      key=kz*ex-kx*ez
      kez=kx*ey-ky*ex
!
!  k x (k x e)
!
      kkex=ky*kez-kz*key
      kkey=kz*kex-kx*kez
      kkez=kx*key-ky*kex
!
!  ik x (k x e) + i*phase
!
!  Normalise ff; since we don't know dt yet, we finalize this
!  within timestep where dt is determined and broadcast.
!
!  This does already include the new sqrt(2) factor (missing in B01).
!  So, in order to reproduce the 0.1 factor mentioned in B01
!  we have to set force=0.07.
!
      ffnorm=sqrt(2.)*k*sqrt(k2-kde**2)/sqrt(kav*cs0**3)
      if (ip<=12) then
        print*,'force_hel_noshear: k,kde,ffnorm,kav,dt,cs0=',k,kde,ffnorm,kav,dt,cs0
        print*,'force_hel_noshear: k*sqrt(k2-kde**2)=',k*sqrt(k2-kde**2)
      endif
      !!(debug...) write(21,'(f10.4,3i3,f7.3)') t,kx,ky,kz,phase
!
!  need to multiply by dt (for Euler step), but it also needs to be
!  divided by sqrt(dt), because square of forcing is proportional
!  to a delta function of the time difference
!
      fact=force/ffnorm*sqrt(dt)
!
!  The wavevector is for the case where Lx=Ly=Lz=2pi. If that is not the
!  case one needs to scale by 2pi/Lx, etc.
!
      fx=exp(cmplx(0.,2*pi/Lx*kx*x+phase))*fact
      fy=exp(cmplx(0.,2*pi/Ly*ky*y))
      fz=exp(cmplx(0.,2*pi/Lz*kz*z))
!
!  possibly multiply forcing by z-profile
!
!-    if (height_ff/=0.) then
!-      if (lroot) print*,'forcing_hel_noshear: include z-profile'
!-      tmpz=(z/height_ff)**2
!-      fz=fz*exp(-tmpz**5/max(1.-tmpz,1e-5))
!-    endif
!
!  need to discuss with axel
!
!  possibly multiply forcing by sgn(z) and radial profile
!
!      if (r_ff/=0.) then
!        if (lroot) &
!             print*,'forcing_hel_noshear: applying sgn(z)*xi(r) profile'
!        !
!        ! only z-dependent part can be done here; radial stuff needs to go
!        ! into the loop
!        !
!        tmpz = tanh(z/width_ff)
!        fz = fz*tmpz
!      endif
!
      if (ip<=5) then
        print*,'force_hel_noshear: fx=',fx
        print*,'force_hel_noshear: fy=',fy
        print*,'force_hel_noshear: fz=',fz
      endif
!
!  prefactor
!
      sig=relhel
      coef(1)=cmplx(k*float(kex),sig*float(kkex))
      coef(2)=cmplx(k*float(key),sig*float(kkey))
      coef(3)=cmplx(k*float(kez),sig*float(kkez))
      if (ip<=5) print*,'force_hel_noshear: coef=',coef
!
! loop the two cases separately, so we don't check for r_ff during
! each loop cycle which could inhibit (pseudo-)vectorisation
!
      if (r_ff == 0) then       ! no radial profile
        do j=1,3
          jf=j+ifff-1
          do n=n1,n2
            do m=m1,m2
              f(l1:l2,m,n,jf) = &
                   f(l1:l2,m,n,jf)+real(coef(j)*fx(l1:l2)*fy(m)*fz(n))
            enddo
          enddo
        enddo
      else                      ! with radial profile
        do j=1,3
          jf=j+ifff-1
          do n=n1,n2
!---        sig = relhel*tmpz(n)
!AB: removed tmpz factor
              sig = relhel
call fatal_error('forcing_hel_noshear','radial profile should be quenched')
            coef(1)=cmplx(k*float(kex),sig*float(kkex))
            coef(2)=cmplx(k*float(key),sig*float(kkey))
            coef(3)=cmplx(k*float(kez),sig*float(kkez))
            do m=m1,m2
              radius = sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
              tmpx = 0.5*(1.-tanh((radius-r_ff)/width_ff))
              f(l1:l2,m,n,jf) = &
                   f(l1:l2,m,n,jf) + real(coef(j)*tmpx*fx(l1:l2)*fy(m)*fz(n))
            enddo
          enddo
        enddo
      endif
!
      if (ip<=12) print*,'force_hel_noshear: forcing OK'
!
    endsubroutine forcing_hel_noshear
!***********************************************************************
    subroutine forcing_roberts(f)
!
!  add some artificial fountain flow
!  (to check for example small scale magnetic helicity loss)
!
!  30-may-02/axel: coded
!
      use Mpicomm
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: sxx,cxx
      real, dimension (mx) :: sx,cx
      real, dimension (my) :: sy,cy
      real, dimension (mz) :: sz,cz,tmpz,gz,gg,ss,gz1
      real :: kx,ky,kz,ffnorm,fac
!
!  identify ourselves
!
      if (headtt.or.ldebug) print*,'forcing_roberts: ENTER'
!
!  need to multiply by dt (for Euler step), but it also needs to be
!  divided by sqrt(dt), because square of forcing is proportional
!  to a delta function of the time difference
!
      kx=kfountain
      ky=kfountain
      kz=1.
!
      sx=sin(kx*x); cx=cos(kx*x)
      sy=sin(ky*y); cy=cos(ky*y)
      sz=sin(kz*z); cz=cos(kz*z)
!
!  abbreviation
!
      sxx=sx(l1:l2)
      cxx=cx(l1:l2)
!
!  g(z) and g'(z)
!  use z-profile to cut off
!
      if (height_ff/=0.) then
        tmpz=(z/height_ff)**2
        gz=sz*exp(-tmpz**5/max(1.-tmpz,1e-5))
!
        fac=1./(60.*dz)
        gg(1:3)=0.; gg(mz-2:mz)=0. !!(border pts to zero)
        gg(4:mz-3)=fac*(45.*(gz(5:mz-2)-gz(3:mz-4)) &
                        -9.*(gz(6:mz-1)-gz(2:mz-5)) &
                           +(gz(7:mz)  -gz(1:mz-6)))
      else
        gz=0
        gg=0
      endif
!
!  make sign antisymmetric
!
      where(z<0)
        ss=-1.
      elsewhere
        ss=1.
      endwhere
      gz1=-ss*gz !!(negative for z>0)
!
!AB: removed nu dependence here. This whole routine is probably not
!AB: needed at the moment, because it is superseded by continuous forcing
!AB: in hydro.f90
!
      !ffnorm=fountain*nu*dt
      ffnorm=fountain*dt
!
!  set forcing function
!
      do n=n1,n2
      do m=m1,m2
        f(l1:l2,m,n,iffx)=f(l1:l2,m,n,iffx)+ffnorm*(+sxx*cy(m)*gz1(n)+cxx*sy(m)*gg(n))
        f(l1:l2,m,n,iffy)=f(l1:l2,m,n,iffy)+ffnorm*(-cxx*sy(m)*gz1(n)+sxx*cy(m)*gg(n))
        f(l1:l2,m,n,iffz)=f(l1:l2,m,n,iffz)+ffnorm*sxx*sy(m)*gz(n)*2.
      enddo
      enddo
!
    endsubroutine forcing_roberts
!***********************************************************************
    subroutine forcing_fountain(f)
!
!  add some artificial fountain flow
!  (to check for example small scale magnetic helicity loss)
!
!  30-may-02/axel: coded
!
      use Mpicomm
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: sxx,cxx
      real, dimension (mx) :: sx,cx
      real, dimension (my) :: sy,cy
      real, dimension (mz) :: sz,cz,tmpz,gz,gg,ss,gz1
      real :: kx,ky,kz,ffnorm,fac
!
!  identify ourselves
!
      if (headtt.or.ldebug) print*,'forcing_fountain: ENTER'
!
!  need to multiply by dt (for Euler step), but it also needs to be
!  divided by sqrt(dt), because square of forcing is proportional
!  to a delta function of the time difference
!
      kx=kfountain
      ky=kfountain
      kz=1.
!
      sx=sin(kx*x); cx=cos(kx*x)
      sy=sin(ky*y); cy=cos(ky*y)
      sz=sin(kz*z); cz=cos(kz*z)
!
!  abbreviation
!
      sxx=sx(l1:l2)
      cxx=cx(l1:l2)
!
!  g(z) and g'(z)
!  use z-profile to cut off
!
      if (height_ff/=0.) then
        tmpz=(z/height_ff)**2
        gz=sz*exp(-tmpz**5/max(1.-tmpz,1e-5))
!
        fac=1./(60.*dz)
        gg(1:3)=0.; gg(mz-2:mz)=0. !!(border pts to zero)
        gg(4:mz-3)=fac*(45.*(gz(5:mz-2)-gz(3:mz-4)) &
                        -9.*(gz(6:mz-1)-gz(2:mz-5)) &
                           +(gz(7:mz)  -gz(1:mz-6)))
      else
        gz=0
        gg=0
      endif
!
!  make sign antisymmetric
!
      where(z<0)
        ss=-1.
      elsewhere
        ss=1.
      endwhere
      gz1=-ss*gz !!(negative for z>0)
!
!AB: removed nu dependence here. This whole routine is probably not
!AB: needed at the moment, because it should be superseded by continuous
!AB: forcing in hydro.f90
!
      !ffnorm=fountain*nu*kfountain**2*dt
      ffnorm=fountain*kfountain**2*dt
!
!  set forcing function
!
      do n=n1,n2
      do m=m1,m2
        f(l1:l2,m,n,iffx)=f(l1:l2,m,n,iffx)+ffnorm*(cxx*sy(m)*gg(n))
        f(l1:l2,m,n,iffy)=f(l1:l2,m,n,iffy)+ffnorm*(sxx*cy(m)*gg(n))
        f(l1:l2,m,n,iffz)=f(l1:l2,m,n,iffz)+ffnorm*sxx*sy(m)*gz(n)*2.
      enddo
      enddo
!
    endsubroutine forcing_fountain
!***********************************************************************
    subroutine forcing_hshear(f)
!
!  add horizontal shear
!
!  19-jun-02/axel+bertil: coded
!
      use Mpicomm
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: fx
      real, dimension (mz) :: fz
      real :: kx,ffnorm
!
!  need to multiply by dt (for Euler step).
!  Define fz with ghost zones, so fz(n) is the correct position
!  Comment: Brummell et al have a polynomial instead!
!
      kx=2*pi/Lx
      fx=cos(kx*x(l1:l2))
      fz=1./cosh(z/width_ff)**2
      ffnorm=force*dt  !(dt for the timestep)
!
!  add to velocity (here only y-component)
!
      do n=n1,n2
      do m=m1,m2
        f(l1:l2,m,n,iuy)=f(l1:l2,m,n,iuy)+ffnorm*fx*fz(n)
      enddo
      enddo
!
    endsubroutine forcing_hshear
!***********************************************************************
    subroutine forcing_twist(f)
!
!  add circular twisting motion, (ux, 0, uz)
!
!  19-jul-02/axel: coded
!
      use Mpicomm
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,nz) :: xx,zz,r2,tmp,fx,fz
      real :: ffnorm,ry2,fy,ytwist1,ytwist2
!
!  identifier
!
      if (headt) print*,'forcing_twist: r_ff,width_ff=',r_ff,width_ff
!
!  need to multiply by dt (for Euler step).
!
      ffnorm=force*dt  !(dt for the timestep)
!
!  add to velocity
!  calculate r2=(x^2+z^2)/r^2
!
      xx=spread(x(l1:l2),2,nz)
      zz=spread(z(n1:n2),1,nx)
      if (r_ff==0.) then
        if (lroot) print*,'forcing_twist: division by r_ff=0!!'
      endif
      r2=(xx**2+zz**2)/r_ff**2
      tmp=exp(-r2/max(1.-r2,1e-5))*ffnorm
      fx=-zz*tmp
      fz=+xx*tmp
!
!  have opposite twists at
!
      y0=xyz0(2)
      ytwist1=y0+0.25*Ly
      ytwist2=y0+0.75*Ly
!
      do m=m1,m2
        !
        ! first twister
        !
        ry2=((y(m)-ytwist1)/width_ff)**2
        fy=exp(-ry2/max(1.-ry2,1e-5))
        f(l1:l2,m,n1:n2,iffx)=f(l1:l2,m,n1:n2,iffx)+fy*fx
        f(l1:l2,m,n1:n2,iffz)=f(l1:l2,m,n1:n2,iffz)+fy*fz
        !
        ! second twister
        !
        ry2=((y(m)-ytwist2)/width_ff)**2
        fy=exp(-ry2/max(1.-ry2,1e-5))
        f(l1:l2,m,n1:n2,iffx)=f(l1:l2,m,n1:n2,iffx)-fy*fx
        f(l1:l2,m,n1:n2,iffz)=f(l1:l2,m,n1:n2,iffz)-fy*fz
      enddo
!
    endsubroutine forcing_twist
!***********************************************************************
    subroutine forcing_diffrot(f,force_ampl)
!
!  Add differential rotation. This routine does not employ continuous
!  forcing, which would be better. It is therefore inferior to the
!  differential rotation procedure implemented directly in hydro.
!
!  26-jul-02/axel: coded
!
      use Mpicomm
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,nz) :: fx,fz,tmp
      real :: force_ampl,ffnorm,ffnorm2
!
!  identifier
!
      if (headt) print*,'forcing_diffrot: ENTER'
!
!  need to multiply by dt (for Euler step).
!
      ffnorm=force_ampl*dt  !(dt for the timestep)
!
!  prepare velocity, Uy=cosx*cosz
!
      fx=spread(cos(x(l1:l2)),2,nz)
      fz=spread(cos(z(n1:n2)),1,nx)
!
!  this forcing term is balanced by diffusion operator;
!  need to multiply by nu*k^2, but k=sqrt(1+1) for the forcing
!
!AB: removed nu dependence here. This whole routine is probably not
!AB: needed at the moment, because it should be superseded by continuous
!AB: forcing in hydro.f90
!
      !ffnorm2=ffnorm*nu*2
      ffnorm2=ffnorm*2
      tmp=ffnorm2*fx*fz
!
!  add
!
      do m=m1,m2
        f(l1:l2,m,n1:n2,iuy)=f(l1:l2,m,n1:n2,iuy)+tmp
      enddo
!
    endsubroutine forcing_diffrot
!***********************************************************************
    subroutine forcing_blobs(f)
!
!  add blobs in entropy every dtforce time units
!
!  28-jul-02/axel: coded
!
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, save :: tforce=0.
      logical, save :: lfirst_call=.true.
      integer, save :: nforce=0
      logical :: lforce
      character (len=intlen) :: ch
      character (len=fnlen) :: file
!
!  identifier
!
      if (headt) print*,'forcing_blobs: ENTER'
!
!  the last forcing time is recorded in tforce.dat
!
      file=trim(datadir)//'/tforce.dat'
      if (lfirst_call) then
        call read_snaptime(trim(file),tforce,nforce,dtforce,t)
        lfirst_call=.false.
      endif
!
!  Check whether we want to do forcing at this time.
!
      call update_snaptime(file,tforce,nforce,dtforce,t,lforce,ch)
      if (lforce) then
        call blob(force,f,iss,radius_ff,location(1),location(2),location(3))
      endif
!
    endsubroutine forcing_blobs
!***********************************************************************
    subroutine forcing_blobHS_random(f)
!
!  add blobs in HS in entropy and density
!  location & time can be random
!
!  08-nov-20/boris: coded
!
      use Sub
      use General, only: random_number_wrapper
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, save :: t_next_blob=1.
      logical, save :: lfirst_call=.true.
      integer, save :: t_interval_blobs=50.
      logical :: lforce
      character (len=intlen) :: ch
      character (len=fnlen) :: file
      real, dimension (3) :: fran
      real :: scaled_interval
!
!  identifier
!
      if (headt) print*,'forcing_blobHS_random: ENTER'
!
      if (lfirst_call) then
        t_next_blob=t+1.
        lfirst_call=.false.
      endif
!
!  Check whether we want to do forcing at this time.
!
      if (t>=t_next_blob) then
        if (lrandom_location) then
          call random_number_wrapper(fran,CHANNEL=channel_force)
          location=fran*Lxyz+xyz0
        else
          location=location_fixed(:,1)
        endif
        open (111,file=trim(datadir)//'/tblobs.dat',position="append")
        write (111,'(f12.3,3f8.4)') t_next_blob, location
        close (111)
!
!  Add a blob in HS equilibrium (entropy & density at the same time)
!
        call blob(force,f,iss,radius_ff,location(1),location(2),location(3))
        call blob(-force,f,ilnrho,radius_ff,location(1),location(2),location(3))
!
        if (lrandom_time) then
          call random_number_wrapper(fran)
          scaled_interval=-log(fran(1))*t_interval_blobs
          t_next_blob=t+scaled_interval
        else
          t_next_blob=t+dtforce
        endif
        print*,'t_next_blob=', t_next_blob
      endif
!
    endsubroutine forcing_blobHS_random
!***********************************************************************
    subroutine forcing_hel_smooth(f)
!
      use DensityMethods, only: getrho
      use General, only: random_number_wrapper
      use Mpicomm
      use Sub
!
!  06-dec-13/nishant: made kkx etc allocatable
!  23-dec-18/axel: forcing_helicity has now similar capabilities
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,3) :: force1,force2,force_vec
      real, dimension (nx) :: ruf,rho
      real, dimension (nx,3) :: variable_rhs,forcing_rhs,force_all
      real :: phase1,phase2,p_weight
      real :: kx01,ky1,kz1,kx02,ky2,kz2
      real :: mulforce_vec,irufm,fsum_tmp,fsum
      integer :: ik1,ik2,ik
!
!  Re-calculate forcing wave numbers if necessary
!
      !tsforce is set to -10 in cdata.f90. It should also be saved in a file
      !so that it will be read again on restarts.
      if (t > tsforce) then
        if (tsforce < 0) then
          call random_number_wrapper(fran1,CHANNEL=channel_force)
        else
          fran1=fran2
        endif
        call random_number_wrapper(fran2,CHANNEL=channel_force)
        tsforce=t+dtforce
      endif
      phase1=pi*(2*fran1(1)-1.)
      ik1=nk*.9999*fran1(2)+1
      kx01=kkx(ik1)
      ky1=kky(ik1)
      kz1=kkz(ik1)
      phase2=pi*(2*fran2(1)-1.)
      ik2=nk*.9999*fran2(2)+1
      kx02=kkx(ik2)
      ky2=kky(ik2)
      kz2=kkz(ik2)
!
!  Calculate forcing function
!
      call hel_vec(f,kx01,ky1,kz1,phase1,kav,force1)
      call hel_vec(f,kx02,ky2,kz2,phase2,kav,force2)
!
!  Determine weight parameter
!
      p_weight=(tsforce-t)/dtforce
      force_vec=p_weight*force1+(1-p_weight)*force2
!
! Find energy input
!
      if (lout .or. lwork_ff) then
        if (idiag_rufm/=0 .or. lwork_ff) then
          irufm=0
          do n=n1,n2
            do m=m1,m2
              forcing_rhs=force_vec(l1:l2,m,n,:)
              variable_rhs=f(l1:l2,m,n,iffx:iffz)!-force_vec(l1:l2,m,n,:)
              call getrho(f(:,m,n,ilnrho),rho)
              call multsv_mn(rho/dt,forcing_rhs,force_all)
              call dot_mn(variable_rhs,force_all,ruf)
              irufm=irufm+sum(ruf)
              !call sum_mn_name(ruf/nwgrid,idiag_rufm)
            enddo
          enddo
        endif
      endif
      irufm=irufm/nwgrid
!
! If we want to make energy input constant
!
      if (lwork_ff) then
!
!
!  on different processors, irufm needs to be communicated
!  to other processors
!
        fsum_tmp=irufm
        call mpireduce_sum(fsum_tmp,fsum)
        irufm=fsum
        call mpibcast_real(irufm)
!
! What should be added to force_vec in order to make the energy
! input equal to work_ff?
!
        mulforce_vec=work_ff/irufm
        if (mulforce_vec > max_force)  mulforce_vec=max_force
      else
        mulforce_vec = 1.0
      endif
!
!  Add forcing
!
      f(l1:l2,m1:m2,n1:n2,1:3)= &
           f(l1:l2,m1:m2,n1:n2,1:3)+force_vec(l1:l2,m1:m2,n1:n2,:)*mulforce_vec
!
! Save for printouts
!
      if (lout) then
        if (idiag_rufm/=0) then
          fname(idiag_rufm)=irufm*mulforce_vec
          itype_name(idiag_rufm)=ilabel_sum
        endif
      endif
!
    endsubroutine forcing_hel_smooth
!***********************************************************************
    subroutine forcing_tidal(f,force_ampl)
!
!  Added tidal forcing function
!  NB: This is still experimental. Use with care!
!
!  20-Nov-12/simon: coded
!  21-jan-15/MR: changes for use of reference state.
!
      use Diagnostics
      use DensityMethods, only: getrho, getrho1
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      real :: force_ampl
      real :: irufm
      real, dimension (nx) :: ruf,rho,rho1
      real, dimension (nx,3) :: variable_rhs,forcing_rhs,force_all
!      real, dimension (nx,3) :: bb,fxb
      integer :: j,jf,l
      real :: fact, dist3
!
!  Normalize ff; since we don't know dt yet, we finalize this
!  within timestep where dt is determined and broadcast.
!
!  need to multiply by dt (for Euler step), but it also needs to be
!  divided by sqrt(dt), because square of forcing is proportional
!  to a delta function of the time difference
!
      fact=2*force_ampl*sqrt(dt)
!
!  loop the two cases separately, so we don't check for r_ff during
!  each loop cycle which could inhibit (pseudo-)vectorisation
!  calculate energy input from forcing; must use lout (not ldiagnos)
!
      irufm=0
      do n=n1,n2
        do m=m1,m2
          variable_rhs=f(l1:l2,m,n,iffx:iffz)
          if (lmomentum_ff) then
            if (ldensity_nolog) then
              if (lreference_state) then
                rho1=1./(f(l1:l2,m,n,irho)+reference_state(:,iref_rho))
              else
                rho1=1./f(l1:l2,m,n,irho)
              endif
            else
              rho1=exp(-f(l1:l2,m,n,ilnrho))
              call getrho1(f(:,m,n,ilnrho),rho1)
            endif
          else
            rho1=1.
          endif
          do l=l1,l2
            dist3 = sqrt( &
                    (R0_tidal*cos(omega_tidal*t)*cos(phi_tidal)-x(l))**2 + &
                    (R0_tidal*sin(omega_tidal*t)               -y(m))**2 + &
                    (R0_tidal*cos(omega_tidal*t)*sin(phi_tidal)-z(n))**2)**3
            forcing_rhs(l-l1+1,1) = &
                rho1(l)*fact*(R0_tidal*cos(omega_tidal*t)*cos(phi_tidal)-x(l))/dist3
            forcing_rhs(l-l1+1,2) = &
                rho1(l)*fact*(R0_tidal*sin(omega_tidal*t)-y(m))/dist3
            forcing_rhs(l-l1+1,3) = &
                rho1(l)*fact*(R0_tidal*cos(omega_tidal*t)*sin(phi_tidal)-z(n))/dist3
          enddo
          do j=1,3
            jf=j+ifff-1
            f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+forcing_rhs(:,j)
          enddo
          if (lout) then
            if (idiag_rufm/=0) then
              call getrho(f(:,m,n,ilnrho),rho)
              call multsv_mn(rho/dt,forcing_rhs,force_all)
              call dot_mn(variable_rhs,force_all,ruf)
              irufm=irufm+sum(ruf)
            endif
          endif
        enddo
      enddo
      !
      ! For printouts
      !
!       if (lout) then
!         if (idiag_rufm/=0) then
!           irufm=irufm/(nwgrid)
!           !
!           !  on different processors, irufm needs to be communicated
!           !  to other processors
!           !
!           fsum_tmp=irufm
!           call mpireduce_sum(fsum_tmp,fsum)
!           irufm=fsum
!           call mpibcast_real(irufm)
!           !
!           fname(idiag_rufm)=irufm
!           itype_name(idiag_rufm)=ilabel_sum
!         endif
!         if (lmagnetic) then
!           if (idiag_fxbxm/=0.or.idiag_fxbym/=0.or.idiag_fxbzm/=0) then
!             call curl(f,iaa,bb)
!             call cross(forcing_rhs,bb,fxb)
!             call sum_mn_name(fxb(:,1),idiag_fxbxm)
!             call sum_mn_name(fxb(:,2),idiag_fxbym)
!             call sum_mn_name(fxb(:,3),idiag_fxbzm)
!           endif
!         endif
!       endif
!
      if (ip<=9) print*,'forcing_tidal: forcing OK'
!
    endsubroutine forcing_tidal
!***********************************************************************
    subroutine hel_vec(f,kx0,ky,kz,phase,kav,force1)
!
!  Add helical forcing function, using a set of precomputed wavevectors.
!  The relative helicity of the forcing function is determined by the factor
!  sigma, called here also relhel. If it is +1 or -1, the forcing is a fully
!  helical Beltrami wave of positive or negative helicity. For |relhel| < 1
!  the helicity less than maximum. For relhel=0 the forcing is nonhelical.
!  The forcing function is now normalized to unity (also for |relhel| < 1).
!
!  10-apr-00/axel: coded
!   3-sep-02/axel: introduced k1_ff, to rescale forcing function if k1/=1.
!  25-sep-02/axel: preset force_ampl to unity (in case slope is not controlled)
!   9-nov-02/axel: corrected normalization factor for the case |relhel| < 1.
!  17-jan-03/nils: adapted from forcing_hel
!
      use EquationOfState, only: cs0
      use General, only: random_number_wrapper
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: kx0,ky,kz
      real :: phase
      real :: kav
      real, dimension (mx,my,mz,3) :: force1
!
      real, dimension (nx) :: radius,tmpx
      complex, dimension (mx) :: fx
      complex, dimension (my) :: fy
      complex, dimension (mz) :: fz
      complex, dimension (3) :: coef
      integer :: j,jf
      real :: kx,k2,k,force_ampl,ffnorm
      real :: ex,ey,ez,kde,sig,fact,kex,key,kez,kkex,kkey,kkez
      real, dimension(3) :: e1,e2,ee,kk
      real :: norm,phi
!
!  in the shearing sheet approximation, kx = kx0 - St*k_y.
!  Here, St=-deltay/Lx
!
      if (Sshear==0.) then
        kx=kx0
      else
        kx=kx0+ky*deltay/Lx
      endif
!
      if (headt.or.ip<5) print*, 'hel_vec: kx0,kx,ky,kz=',kx0,kx,ky,kz
      k2=kx**2+ky**2+kz**2
      k=sqrt(k2)
!
! Find e-vector
!
      !
      ! Start with old method (not isotropic) for now.
      ! Pick e1 if kk not parallel to ee1. ee2 else.
      !
      if ((ky==0).and.(kz==0)) then
        ex=0; ey=1; ez=0
      else
        ex=1; ey=0; ez=0
      endif
      if (.not. old_forcing_evector) then
        !
        !  Isotropize ee in the plane perp. to kk by
        !  (1) constructing two basis vectors for the plane perpendicular
        !      to kk, and
        !  (2) choosing a random direction in that plane (angle phi)
        !  Need to do this in order for the forcing to be isotropic.
        !
        kk = (/kx, ky, kz/)
        ee = (/ex, ey, ez/)
        call cross(kk,ee,e1)
        call dot2(e1,norm); e1=e1/sqrt(norm) ! e1: unit vector perp. to kk
        call cross(kk,e1,e2)
        call dot2(e2,norm); e2=e2/sqrt(norm) ! e2: unit vector perp. to kk, e1
        call random_number_wrapper(phi,CHANNEL=channel_force); phi = phi*2*pi
        ee = cos(phi)*e1 + sin(phi)*e2
        ex=ee(1); ey=ee(2); ez=ee(3)
      endif
!
!  k.e
!
      call dot(kk,ee,kde)
!
!  k x e
!
      kex=ky*ez-kz*ey
      key=kz*ex-kx*ez
      kez=kx*ey-ky*ex
!
!  k x (k x e)
!
      kkex=ky*kez-kz*key
      kkey=kz*kex-kx*kez
      kkez=kx*key-ky*kex
!
!  ik x (k x e) + i*phase
!
!  Normalize ff; since we don't know dt yet, we finalize this
!  within timestep where dt is determined and broadcast.
!
!  This does already include the new sqrt(2) factor (missing in B01).
!  So, in order to reproduce the 0.1 factor mentioned in B01
!  we have to set force=0.07.
!
!  Furthermore, for |relhel| < 1, sqrt(2) should be replaced by
!  sqrt(1.+relhel**2). This is done now (9-nov-02).
!  This means that the previous value of force=0.07 (for relhel=0)
!  should now be replaced by 0.05.
!
!  Note: kav is not to be scaled with k1_ff (forcing should remain
!  unaffected when changing k1_ff).
!
      ffnorm=sqrt(1.+relhel**2) &
        *k*sqrt(k2-kde**2)/sqrt(kav*cs0**3)*(k/kav)**slope_ff
      if (ip<=9) then
        print*,'hel_vec: k,kde,ffnorm,kav,dt,cs0=',k,kde,ffnorm,kav,dt,cs0
        print*,'hel_vec: k*sqrt(k2-kde**2)=',k*sqrt(k2-kde**2)
      endif
!
!  need to multiply by dt (for Euler step), but it also needs to be
!  divided by sqrt(dt), because square of forcing is proportional
!  to a delta function of the time difference
!
      fact=force/ffnorm*sqrt(dt)
!
!  The wavevector is for the case where Lx=Ly=Lz=2pi. If that is not the
!  case one needs to scale by 2pi/Lx, etc.
!
      fx=exp(cmplx(0.,kx*k1_ff*x+phase))*fact
      fy=exp(cmplx(0.,ky*k1_ff*y))
      fz=exp(cmplx(0.,kz*k1_ff*z))
!
!  possibly multiply forcing by z-profile
!
!-    if (height_ff/=0.) then
!-      if (lroot) print*,'hel_vec: include z-profile'
!-      tmpz=(z/height_ff)**2
!-      fz=fz*exp(-tmpz**5/max(1.-tmpz,1e-5))
!-    endif
!
!  possibly multiply forcing by sgn(z) and radial profile
!
      if (r_ff/=0.) then
        if (lroot) &
             print*,'hel_vec: applying sgn(z)*xi(r) profile'
        !
        ! only z-dependent part can be done here; radial stuff needs to go
        ! into the loop
        !
        fz = fz*profz_k
      endif
!
      if (ip<=5) then
        print*,'hel_vec: fx=',fx
        print*,'hel_vec: fy=',fy
        print*,'hel_vec: fz=',fz
      endif
!
!  prefactor
!
      coef(1)=cmplx(k*kex,relhel*kkex)
      coef(2)=cmplx(k*key,relhel*kkey)
      coef(3)=cmplx(k*kez,relhel*kkez)
      if (ip<=5) print*,'hel_vec: coef=',coef
!
! loop the two cases separately, so we don't check for r_ff during
! each loop cycle which could inhibit (pseudo-)vectorisation
!
      if (r_ff == 0) then       ! no radial profile
        if (lwork_ff) then
          force_ampl=calc_force_ampl(f,fx,fy,fz,coef)
        else
          force_ampl=1.
        endif
        do j=1,3
          jf=j+ifff-1
          do n=n1,n2
            do m=m1,m2
               force1(l1:l2,m,n,jf) = &
                +force_ampl*real(coef(j)*fx(l1:l2)*fy(m)*fz(n))
            enddo
          enddo
        enddo
      else                      ! with radial profile
        do j=1,3
          jf=j+ifff-1
          do n=n1,n2
!--         sig = relhel*tmpz(n)
!AB: removed tmpz factor
              sig = relhel
call fatal_error('hel_vec','radial profile should be quenched')
            coef(1)=cmplx(k*kex,sig*kkex)
            coef(2)=cmplx(k*key,sig*kkey)
            coef(3)=cmplx(k*kez,sig*kkez)
            do m=m1,m2
              radius = sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
              tmpx = 0.5*(1.-tanh((radius-r_ff)/width_ff))
              force1(l1:l2,m,n,jf) =  real(coef(j)*tmpx*fx(l1:l2)*fy(m)*fz(n))
            enddo
          enddo
        enddo
      endif
!
      if (ip<=9) print*,'hel_vec: forcing OK'
!
    endsubroutine hel_vec
!***********************************************************************
    subroutine calc_fluxring_cylindrical(force,i)
!
!   4-aug-11/dhruba+axel: adapted from fluxring_cylindrical
!
      real, dimension (nx,3), intent(out):: force
      integer,                intent(in) :: i
!
      real, dimension (nx) :: argum,s
      real ::  b0=1., s0=2., width=.2
!
!  density for the magnetic flux flux ring
!
      s=x(l1:l2)
      argum=(s-s0)/width
      force(:,1)=2.*s*(b0/(width*s0))**2*exp(-2*argum**2)*(width**2-s**2+s*s0)
      force(:,2)=0.
      force(:,3)=0.
      force = fcont_ampl(i)*force
!
    endsubroutine calc_fluxring_cylindrical
!***********************************************************************
    subroutine calc_counter_centrifugal(force,i)
!
!   4-aug-11/dhruba+axel: adapted from fluxring_cylindrical
! Calculates the force required to counter the centrifugal force coming
! from an existing differential rotation.
!
      real, dimension (nx,3), intent(out):: force
      integer,                intent(in) :: i
!
      real, dimension (nx) :: vphi,omega_diffrot
!
      omega_diffrot = ampl_diffrot*x(l1:l2)**(omega_exponent)
      vphi =  x(l1:l2)*omega_diffrot
      force(:,1)=-vphi**2/x(l1:l2)
      force(:,2)=0.
      force(:,3)=0.
      force = fcont_ampl(i)*force
!
    endsubroutine calc_counter_centrifugal
!***********************************************************************
    subroutine calc_GP_TC13(i,sp)
!
!   4-apr-22/hongzhe: The Galloway-Proctor forcing used in Tobias &
! Cattaneo (2013). See also Appendix A of Pongkitiwanichakul+2016.
!
      use General, only: random_number_wrapper
!
      integer,  intent(in) :: i
      character (len=*), intent(in) :: sp
!
      integer :: ii
      real :: k,knorm,omega,tcor,R_GP,xi,eta
!
      sinxt(:,i)=0.
      cosxt(:,i)=0.
      sinyt(:,i)=0.
      cosyt(:,i)=0.
      sinzt(:,i)=0.
      coszt(:,i)=0.
!
      if (nk_GP>100) call fatal_error('calc_GP_TC13','large nk_GP not implemented')
      do ii=1,nk_GP
        if (nk_GP==1) then
          k = kmin_GP
        else
          k = kmin_GP + (ii-1) * (kmax_GP-kmin_GP)/(nk_GP-1)
        endif
        knorm = k/kmin_GP
        omega = omega_fcont(i)*( knorm**(2-beta_GP) )
        tcor = tcor_GP*( knorm**(beta_GP-2) )
        !
        !  Refresh random parameters
        !
        if (lfirst) then
          call random_number_wrapper(R_GP,CHANNEL=channel_force)
          if ( it==1 .or. R_GP<=(dt/tcor) ) then
            call random_number_wrapper(xi_GP(ii,i), CHANNEL=channel_force)
            call random_number_wrapper(eta_GP(ii,i),CHANNEL=channel_force)
          endif
        endif
        xi = 2*pi*( xi_GP(ii,i)-0.5)
        eta= 2*pi*( eta_GP(ii,i)-0.5)
        !
        !  Compute needed functions
        !
        select case (sp)
        case('xyz')
          sinxt(:,i) = sinxt(:,i) + k*( knorm**(-beta_GP) ) * &
              sin( k * ( x-xi  + cos(omega*t) ) )
          cosxt(:,i) = cosxt(:,i) + k*( knorm**(-beta_GP) ) * &
              cos( k * ( x-xi  + cos(omega*t) ) )
          sinyt(:,i) = sinyt(:,i) + k*( knorm**(-beta_GP) ) * &
              sin( k * ( y-eta + sin(omega*t) ) )
          cosyt(:,i) = sinyt(:,i) + k*( knorm**(-beta_GP) ) * &
              cos( k * ( y-eta + sin(omega*t) ) )
        case('yzx')
          sinzt(:,i) = sinzt(:,i) + k*( knorm**(-beta_GP) ) * &
              sin( k * ( z-xi  + cos(omega*t) ) )
          coszt(:,i) = coszt(:,i) + k*( knorm**(-beta_GP) ) * &
              cos( k * ( z-xi  + cos(omega*t) ) )
          sinxt(:,i) = sinxt(:,i) + k*( knorm**(-beta_GP) ) * &
              sin( k * ( x-eta + sin(omega*t) ) )
          cosxt(:,i) = sinxt(:,i) + k*( knorm**(-beta_GP) ) * &
              cos( k * ( x-eta + sin(omega*t) ) )
        case default
          call fatal_error('calc_GP_TC13','no valid iforcing_cont specified')
        endselect
      enddo
!
    endsubroutine calc_GP_TC13
!***********************************************************************
    subroutine forcing_after_boundary(f)
!
      real, dimension (mx,my,mz,mfarray),intent(OUT) :: f
!
      if (lforcing_cont) call forcing_cont_after_boundary
!
      if (lff_as_aux) f(l1:l2,m1:m2,n1:n2,ifx:ifz)=0.

    endsubroutine forcing_after_boundary
!***********************************************************************
    subroutine forcing_cont_after_boundary
!
!  precalculate parameters that are new at each timestep,
!  but the same for all pencils
!
      integer :: i
      real :: ecost,esint,ecoxt,ecoyt,ecozt
!
!  for the AKA effect, calculate auxiliary functions phi1_ff and phi2_ff
!
      do i=1,n_forcing_cont
        if (iforcing_cont(i)=='AKA') then
          phi1_ff(:,i)=cos(kf_fcont(i)*y+omega_fcont(i)*t)
          phi2_ff(:,i)=cos(kf_fcont(i)*x-omega_fcont(i)*t)
        elseif (iforcing_cont(i)=='GP') then
          ecost=eps_fcont(i)*cos(omega_fcont(i)*t)
          esint=eps_fcont(i)*sin(omega_fcont(i)*t)
          sinxt(:,i)=sin(kf_fcont(i)*x+ecost)
          cosxt(:,i)=cos(kf_fcont(i)*x+ecost)
          sinyt(:,i)=sin(kf_fcont(i)*y+esint)
          cosyt(:,i)=cos(kf_fcont(i)*y+esint)
        elseif (iforcing_cont(i)=='GP_TC13') then
          call calc_GP_TC13(i,'xyz')
        elseif (iforcing_cont(i)=='GP_TC13_yzx') then
          call calc_GP_TC13(i,'yzx')
        elseif (iforcing_cont(i)=='ABCtdep') then
          ecoxt=eps_fcont(i)*cos( omega_fcont(i)*t)
          ecoyt=eps_fcont(i)*cos(omegay_fcont(i)*t)
          ecozt=eps_fcont(i)*cos(omegaz_fcont(i)*t)
          sinxt(:,i)=sin(kf_fcont(i)*x+ecoxt); cosxt(:,i)=cos(kf_fcont(i)*x+ecoxt)
          sinyt(:,i)=sin(kf_fcont(i)*y+ecoyt); cosyt(:,i)=cos(kf_fcont(i)*y+ecoyt)
          sinzt(:,i)=sin(kf_fcont(i)*z+ecozt); coszt(:,i)=cos(kf_fcont(i)*z+ecozt)
        endif
      enddo
!
    endsubroutine forcing_cont_after_boundary
!***********************************************************************
    subroutine pencil_criteria_forcing
!
!  All pencils that the Forcing module depends on are specified here.
!
!  24-mar-08/axel: adapted from density.f90
!
      if (lforcing_cont) lpenc_requested(i_fcont)=.true.
      if (iforcing_cont(1)=='(0,cosx*cosz,0)_Lor') &
         lpenc_requested(i_rho1)=.true.
      if (lmomentum_ff) lpenc_requested(i_rho1)=.true.
      if (idiag_qfm/=0) lpenc_requested(i_curlo)=.true.
!
    endsubroutine pencil_criteria_forcing
!***********************************************************************
    subroutine pencil_interdep_forcing(lpencil_in)
!
!  Interdependency among pencils from the Forcing module is specified here.
!
!  24-mar-08/axel: adapted from density.f90
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_forcing
!***********************************************************************
    subroutine calc_pencils_forcing(f,p)
!
!  Calculate forcing pencils.
!
!  24-mar-08/axel: adapted from density.f90
!   6-feb-09/axel: added epsilon factor in ABC flow (eps_fcont=1. -> nocos)
!
      use Sub, only: multsv_mn
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(inout) :: f,p
!
      integer :: i
!
!  calculate forcing
!
      if (lpencil(i_fcont)) then
        if (headtt .and. lroot) print*,'forcing: add continuous forcing'

        do i=1,n_forcing_cont
          call forcing_cont(i,p%fcont(:,:,i),rho1=p%rho1)
          ! put force into auxiliary variable, if requested
          if (lff_as_aux) then
            if (i == 1) then
              f(l1:l2,m,n,ifx:ifz) = p%fcont(:,:,i)
            else
              f(l1:l2,m,n,ifx:ifz) = f(l1:l2,m,n,ifx:ifz) + p%fcont(:,:,i)
            endif
          endif
!
!  divide by rho if lmomentum_ff=T
!  MR: better to place it in hydro
!
          if (i==1 .and. lmomentum_ff) call multsv_mn(p%rho1,p%fcont(:,:,1),p%fcont(:,:,1))
        enddo
      endif
!
    endsubroutine calc_pencils_forcing
!***********************************************************************
    subroutine random_isotropic_KS_setup(initpower,kmin,kmax)
!
!  Produces random, isotropic field from energy spectrum following the
!  KS method (Malik and Vassilicos, 1999.)
!
!  More to do; unsatisfactory so far - at least for a steep power-law
!  energy spectrum.
!
!  27-may-05/tony: modified from snod's KS hydro initial
!  03-feb-06/weezy: Attempted rewrite to guarantee periodicity of
!                    KS modes.
!
      use Sub, only: cross, dot2
      use General, only: random_number_wrapper
!
      integer :: modeN
!
      real, dimension (3) :: k_unit
      real, dimension (3) :: ee,e1,e2
      real, dimension (6) :: r
      real :: initpower,kmin,kmax
      real, dimension(KS_modes) :: k,dk,energy,ps
      real :: theta,phi,alpha,beta
      real :: ex,ey,ez,norm,a
!
      if (.not.allocated(KS_k)) then
        allocate(KS_k(3,KS_modes))
        allocate(KS_A(3,KS_modes))
        allocate(KS_B(3,KS_modes))
        allocate(KS_omega(KS_modes))
      endif
!
      kmin=2.*pi      !/(1.0*Lxyz(1))
      kmax=128.*pi    !nx*pi
      a=(kmax/kmin)**(1./(KS_modes-1.))
!
!  Loop over all modes.
!
      do modeN=1,KS_modes
!
!  Pick wavenumber.
!
        k=kmin*(a**(modeN-1.))
!
!  Pick 4 random angles for each mode.
!
        call random_number_wrapper(r);
        theta=pi*(r(1) - 0.)
        phi=pi*(2*r(2) - 0.)
        alpha=pi*(2*r(3) - 0.)
        beta=pi*(2*r(4) - 0.)
!
!  Make a random unit vector by rotating fixed vector to random position
!  (alternatively make a random transformation matrix for each k).
!
        k_unit(1)=sin(theta)*cos(phi)
        k_unit(2)=sin(theta)*sin(phi)
        k_unit(3)=cos(theta)
!
        energy=(((k/kmin)**2. +1.)**(-11./6.))*(k**2.) &
            *exp(-0.5*(k/kmax)**2.)
!
!  Make a vector KS_k of length k from the unit vector for each mode.
!
        KS_k(:,modeN)=k*k_unit(:)
        KS_omega(:)=sqrt(energy(:)*(k(:)**3.))
!
!  Construct basis for plane having rr normal to it
!  (bit of code from forcing to construct x', y').
!
      if ((k_unit(2)==0).and.(k_unit(3)==0)) then
          ex=0.; ey=1.; ez=0.
        else
          ex=1.; ey=0.; ez=0.
        endif
        ee = (/ex, ey, ez/)
!
        call cross(k_unit(:),ee,e1)
!  e1: unit vector perp. to KS_k
        call dot2(e1,norm); e1=e1/sqrt(norm)
        call cross(k_unit(:),e1,e2)
!  e2: unit vector perp. to KS_k, e1
        call dot2(e2,norm); e2=e2/sqrt(norm)
!
!  Make two random unit vectors KS_B and KS_A in the constructed plane.
!
        KS_A(:,modeN) = cos(alpha)*e1 + sin(alpha)*e2
        KS_B(:,modeN) = cos(beta)*e1  + sin(beta)*e2
!
!  Make sure dk is set.
!
        call error('random_isotropic_KS_setup', 'Using uninitialized dk')
        dk=0.                     ! to make compiler happy
!
        ps=sqrt(2.*energy*dk)   !/3.0)
!
!  Give KS_A and KS_B length ps.
!
        KS_A(:,modeN)=ps*KS_A(:,modeN)
        KS_B(:,modeN)=ps*KS_B(:,modeN)
!
      enddo
!
!  Form RA = RA x k_unit and RB = RB x k_unit.
!  Note: cannot reuse same vector for input and output.
!
      do modeN=1,KS_modes
        call cross(KS_A(:,modeN),k_unit(:),KS_A(:,modeN))
        call cross(KS_B(:,modeN),k_unit(:),KS_B(:,modeN))
      enddo
!
      call keep_compiler_quiet(initpower)
!
    endsubroutine random_isotropic_KS_setup
!***********************************************************************
    subroutine random_isotropic_KS_setup_test
!
!  Produces random, isotropic field from energy spectrum following the
!  KS method (Malik and Vassilicos, 1999.)
!  This test case only uses 3 very specific modes (useful for comparison
!  with Louise's kinematic dynamo code.
!
!  03-feb-06/weezy: modified from random_isotropic_KS_setup
!
      use Sub, only: cross
      use General, only: random_number_wrapper
!
      integer :: modeN
!
      real, dimension (3,KS_modes) :: k_unit
      real, dimension(KS_modes) :: k,dk,energy,ps
      real :: initpower,kmin,kmax
!
      if (.not.allocated(KS_k)) then
        allocate(KS_k(3,KS_modes))
        allocate(KS_A(3,KS_modes))
        allocate(KS_B(3,KS_modes))
        allocate(KS_omega(KS_modes))
      endif
!
      initpower=-5./3.
      kmin=10.88279619
      kmax=23.50952672
!
      KS_k(1,1)=2.00*pi
      KS_k(2,1)=-2.00*pi
      KS_k(3,1)=2.00*pi
!
      KS_k(1,2)=-4.00*pi
      KS_k(2,2)=0.00*pi
      KS_k(3,2)=2.00*pi
!
      KS_k(1,3)=4.00*pi
      KS_k(2,3)=2.00*pi
      KS_k(3,3)=-6.00*pi
!
      KS_k(1,1)=+1; KS_k(2,1)=-1; KS_k(3,1)=1
      KS_k(1,2)=+0; KS_k(2,2)=-2; KS_k(3,2)=1
      KS_k(1,3)=+0; KS_k(2,3)=-0; KS_k(3,3)=1
!
      k(1)=kmin
      k(2)=14.04962946
      k(3)=kmax
!
      do modeN=1,KS_modes
        k_unit(:,modeN)=KS_k(:,modeN)/k(modeN)
      enddo
!
      kmax=k(KS_modes)
      kmin=k(1)
!
      do modeN=1,KS_modes
        if (modeN==1) dk(modeN)=(k(modeN+1)-k(modeN))/2.
        if (modeN>1.and.modeN<KS_modes) &
            dk(modeN)=(k(modeN+1)-k(modeN-1))/2.
        if (modeN==KS_modes) dk(modeN)=(k(modeN)-k(modeN-1))/2.
      enddo
!
      do modeN=1,KS_modes
         energy(modeN)=((k(modeN)**2 +1.)**(-11./6.))*(k(modeN)**2) &
             *exp(-0.5*(k(modeN)/kmax)**2)
      enddo
!
      ps=sqrt(2.*energy*dk)
!
      KS_A(1,1)=1.00/sqrt(2.00)
      KS_A(2,1)=-1.00/sqrt(2.00)
      KS_A(3,1)=0.00
!
      KS_A(1,2)=1.00/sqrt(3.00)
      KS_A(2,2)=1.00/sqrt(3.00)
      KS_A(3,2)=-1.00/sqrt(3.00)
!
      KS_A(1,3)=-1.00/2.00
      KS_A(2,3)=-1.00/2.00
      KS_A(3,3)=1.00/sqrt(2.00)
!
      KS_B(1,3)=1.00/sqrt(2.00)
      KS_B(2,3)=-1.00/sqrt(2.00)
      KS_B(3,3)=0.00
!
      KS_B(1,1)=1.00/sqrt(3.00)
      KS_B(2,1)=1.00/sqrt(3.00)
      KS_B(3,1)=-1.00/sqrt(3.00)
!
      KS_B(1,2)=-1.00/2.00
      KS_B(2,2)=-1.00/2.00
      KS_B(3,2)=1.00/sqrt(2.00)
!
      do modeN=1,KS_modes
        KS_A(:,modeN)=ps(modeN)*KS_A(:,modeN)
        KS_B(:,modeN)=ps(modeN)*KS_B(:,modeN)
      enddo
!
!  Form RA = RA x k_unit and RB = RB x k_unit.
!
       do modeN=1,KS_modes
         call cross(KS_A(:,modeN),k_unit(:,modeN),KS_A(:,modeN))
         call cross(KS_B(:,modeN),k_unit(:,modeN),KS_B(:,modeN))
       enddo
!
    endsubroutine random_isotropic_KS_setup_test
!***********************************************************************
    subroutine forcing_cont(i,force,rho1)
!
!   9-apr-10/MR: added RobertsFlow_exact forcing, compensates \nu\nabla^2 u
!                and u.grad u for Roberts geometry
!   4-nov-11/MR: now also compensates Coriolis force
!
!  Note: It is not enough to set lforcing_cont = T in input parameters of
!  forcing one must also set  lforcing_cont_uu = T in hydro for the
!  continuous-in-time forcing to be added to velocity.
!  Alternatively, one can set lforcing_cont_bb = T in magnetic.
!
      use Gravity, only: gravz
      use Mpicomm, only: stop_it
      use Sub, only: quintic_step, quintic_der_step, step, vortex
      use Viscosity, only: getnu
!
      integer,                intent(in) :: i
      real, dimension (nx,3), intent(out):: force
      real, dimension (nx), optional, intent(in) :: rho1
!
      real, dimension (nx) :: tmp
      real :: fact, fact1, fact2, fpara, dfpara, sqrt21k1
      real :: kf, kx, ky, kz, nu, arg, ecost, esint
      integer :: i2d1=1,i2d2=2,i2d3=3,modeN
      real, dimension(nx) :: kdotxwt, cos_kdotxwt, sin_kdotxwt
!
        select case (iforcing_cont(i))
        case('Fy=const')
          force(:,1)=0.
          force(:,2)=ampl_ff(i)
          force(:,3)=0.
        case('Fz=const')
          force(:,1)=0.
          force(:,2)=0.
          force(:,3)=ampl_ff(i)
        case('ABC')
          fact2=relhel
          fact1=.5*(1.+relhel**2)
          fact=ampl_ff(i)/sqrt(fact1*(ABC_A(i)**2+ABC_B(i)**2+ABC_C(i)**2))
          force(:,1)=fact*(ABC_C(i)*sinz(n    ,i)+fact2*ABC_B(i)*cosy(m    ,i))
          force(:,2)=fact*(ABC_A(i)*sinx(l1:l2,i)+fact2*ABC_C(i)*cosz(n    ,i))
          force(:,3)=fact*(ABC_B(i)*siny(m    ,i)+fact2*ABC_A(i)*cosx(l1:l2,i))
        case('ABCtdep')
          fact2=relhel
          fact=ampl_ff(i)/sqrt(ABC_A(i)**2+ABC_B(i)**2+ABC_C(i)**2)
          force(:,1)=fact*(ABC_C(i)*sinzt(n    ,i)+fact2*ABC_B(i)*cosyt(m    ,i))
          force(:,2)=fact*(ABC_A(i)*sinxt(l1:l2,i)+fact2*ABC_C(i)*coszt(n    ,i))
          force(:,3)=fact*(ABC_B(i)*sinyt(m    ,i)+fact2*ABC_A(i)*cosxt(l1:l2,i))
        case ('AKA')
          fact=sqrt(2.)*ampl_ff(i)
          force(:,1)=fact*phi1_ff(m    ,i)
          force(:,2)=fact*phi2_ff(l1:l2,i)
          force(:,3)=fact*(phi1_ff(m,i)+phi2_ff(l1:l2,i))
        case ('grav_z')
          force(:,1)=0
          force(:,2)=0
          force(:,3)=gravz*ampl_ff(i)*cos(omega_ff*t)
        case ('uniform_vorticity')
          if (lcylindrical_coords) then
            force(:,1)=-x(l1:l2)*y(m)*ampl_ff(i)*cos(omega_ff*t)
            force(:,2)=x(l1:l2)*ampl_ff(i)*cos(omega_ff*t)
            force(:,3)=0
          else
            force(:,1)=z(n)*ampl_ff(i)*cos(omega_ff*t)
            force(:,2)=0
            force(:,3)=-x(l1:l2)*ampl_ff(i)*cos(omega_ff*t)
          endif
        case('KolmogorovFlow-x')
          fact=ampl_ff(i)
          force(:,1)=0
          force(:,2)=fact*cosx(l1:l2,i)
          force(:,3)=0
        case('KolmogorovFlow-z')
          fact=ampl_ff(i)
          force(:,1)=fact*cosz(n,i)
          force(:,2)=0.
          force(:,3)=0.
        case ('nocos')
          fact=ampl_ff(i)
          force(:,1)=fact*sinz(n    ,i)
          force(:,2)=fact*sinx(l1:l2,i)
          force(:,3)=fact*siny(m    ,i)
!
!  Amplitude is given by ampl_ff/nu*k^2, and k^2=2 in a 2-D box
!  of size (2pi)^2, for example.
!
        case('Straining')
          fact=ampl_ff(i)
          fact2=-(dimensionality-1)*ampl_ff(i)
          force(:,1)=fact *sinx(l1:l2,i)*cosy(m,i)*cosz(n,i)
          force(:,2)=fact *cosx(l1:l2,i)*siny(m,i)*cosz(n,i)
          force(:,3)=fact2*cosx(l1:l2,i)*cosy(m,i)*sinz(n,i)
!
!  Amplitude is given by ampl_ff/nu*k^2, and k^2=2 in a 2-D box
!  of size (2pi)^2, for example.
!
        case('StrainingExcact')
          call fatal_error("forcing_cont","not checked yet")
          call getnu(nu_input=nu)
          kx=k1_ff; kz=k1_ff
          kf=sqrt(kx**2+kz**2)
          fact=ampl_ff(i)*kf**2*nu
          fact1=ampl_ff(i)*ampl_ff(i)*kx*kz
          fact2=-(dimensionality-1)*ampl_ff(i)*kf**2*nu
          force(:,1)=fact *sinx(l1:l2,i)*cosy(m,i)*cosz(n,i)-fact1*kz*cosx(l1:l2,i)*cosy(m,i)*sinz(n,i)
          force(:,2)=fact *cosx(l1:l2,i)*siny(m,i)*cosz(n,i)
          force(:,3)=fact2*cosx(l1:l2,i)*cosy(m,i)*sinz(n,i)-fact1*kx*sinx(l1:l2,i)*cosy(m,i)*cosz(n,i)
        case('RobertsFlow')
          fact=ampl_ff(i)
          force(:,1)=-fact*cosx(l1:l2,i)*siny(m,i)
          force(:,2)=+fact*sinx(l1:l2,i)*cosy(m,i)
          force(:,3)=+relhel*fact*cosx(l1:l2,i)*cosy(m,i)*sqrt2
        case('RobertsFlowII')
          fact=ampl_ff(i)
          force(:,1)=-fact*cosx(l1:l2,i)*siny(m,i)
          force(:,2)=+fact*sinx(l1:l2,i)*cosy(m,i)
          force(:,3)=+relhel*fact*sinx(l1:l2,i)*siny(m,i)*sqrt2
        case('RobertsFlowMask')
          tmp=ampl_ff(i)*profx_ampl*profy_ampl(m)
          force(:,1)=-tmp*cosx(l1:l2,i)*siny(m,i)
          force(:,2)=+tmp*sinx(l1:l2,i)*cosy(m,i)
          force(:,3)=+relhel*tmp*cosx(l1:l2,i)*cosy(m,i)*sqrt2
        case('RobertsFlow2d')
          fact=ampl_ff(i)
          i2d1=1;i2d2=2;i2d3=3
          if (l2dxz) then
            i2d1=2;i2d2=1;i2d3=2
          endif
          if (l2dyz) then
            i2d1=3;i2d2=2;i2d3=1
          endif
          force(:,i2d1)=-fact*cos(k2d*x(l1:l2))*sin(k2d*y(m))
          force(:,i2d2)=+fact*sin(k2d*x(l1:l2))*cos(k2d*y(m))
          force(:,i2d3)= 0.
        case ('RobertsFlow_exact')
          kx=kf_fcont(i); ky=kf_fcont(i)
          kf=sqrt(kx*kx+ky*ky)
          call getnu(nu_input=nu)
          fact=ampl_ff(i)*kf*kf*nu
          fact2=ampl_ff(i)*ampl_ff(i)*kx*ky
!
!!print*, 'forcing: kx, ky, kf, fact, fact2=', kx, ky, kf, fact, fact2
!
          force(:,1)=-fact*ky*cosx(l1:l2,i)*siny(m,i) - fact2*ky*sinx(l1:l2,i)*cosx(l1:l2,i)
          force(:,2)=+fact*kx*sinx(l1:l2,i)*cosy(m,i) - fact2*kx*siny(m,i)*cosy(m,i)
          force(:,3)=+fact*relhel*kf*cosx(l1:l2,i)*cosy(m,i)
!
          if ( Omega/=0. .and. theta==0. ) then              ! Obs, only implemented for rotation axis in z direction.
            fact = 2.*ampl_ff(i)*Omega
            force(:,1)= force(:,1)-fact*kx*sinx(l1:l2,i)*cosy(m,i)
            force(:,2)= force(:,2)-fact*ky*cosx(l1:l2,i)*siny(m,i)
          endif
!
        case ('RobertsFlow-zdep')
          if (headtt) print*,'z-dependent Roberts flow; eps_fcont=',eps_fcont
          fpara=quintic_step(z(n),-1.+eps_fcont(i),eps_fcont(i)) &
               -quintic_step(z(n),+1.-eps_fcont(i),eps_fcont(i))
          dfpara=quintic_der_step(z(n),-1.+eps_fcont(i),eps_fcont(i))&
                -quintic_der_step(z(n),+1.-eps_fcont(i),eps_fcont(i))
!
!  abbreviations
!
          sqrt21k1=1./(sqrt2*k1_ff)
!
!  amplitude factor missing in upper lines
!
          force(:,1)=-ampl_ff(i)*cosx(l1:l2,i)*siny(m,i) &
              -dfpara*ampl_ff(i)*sinx(l1:l2,i)*cosy(m,i)*sqrt21k1
          force(:,2)=+ampl_ff(i)*sinx(l1:l2,i)*cosy(m,i) &
              -dfpara*ampl_ff(i)*cosx(l1:l2,i)*siny(m,i)*sqrt21k1
          force(:,3)=+fpara*ampl_ff(i)*cosx(l1:l2,i)*cosy(m,i)*sqrt2
!
!  f=(sinx,0,0)
!
        case('Roberts-for-SSD')
          fact=ampl_ff(i)
          force(:,1)=-fact*cosx(l1:l2,i)*siny(m,i)
          force(:,2)=+fact*sinx(l1:l2,i)*cosy(m,i)
          force(:,3)= fact*rel_zcomp*(cosx(l1:l2,i)**2-0.5)*(cosy(m,i)**2-0.5)
!
        case ('sinx')
          if (lgentle(i).and.t<tgentle(i)) then
            fact=.5*ampl_ff(i)*(1.-cos(pi*t/tgentle(i)))
          else
            fact=ampl_ff(i)
          endif
          force(:,1)=fact*sinx(l1:l2,i)
          force(:,2)=0.
          force(:,3)=0.
!
!  f=(0,0,cosx)
!
        case ('(0,0,cosx)')
          force(:,1)=0.
          force(:,2)=0.
          force(:,3)=ampl_ff(i)*cosx(l1:l2,i)
!
!  f=(0,0,cosx)
!
        case ('(0,0,cosxcosy)')
          force(:,1)=0.
          force(:,2)=0.
          force(:,3)=ampl_ff(i)*cosx(l1:l2,i)*cosy(m,i)
!
!  f=B=(0,0,cosx*cosy)
!
        case ('B=(0,0,cosxcosy)')
          force(:,1)=-0.5*ampl_ff(i)*cosx(l1:l2,i)*siny(m,i)
          force(:,2)= 0.5*ampl_ff(i)*sinx(l1:l2,i)*cosy(m,i)
          force(:,3)= 0.
!
!  f=(0,x,0)
!
        case ('(0,x,0)')
          force(:,1)=0.
          force(:,2)=ampl_ff(i)*x(l1:l2)
          force(:,3)=0.
!
!  f=(0,sinx,0)
!
        case ('(0,sinx,0)')
          force(:,1)=0.
          force(:,2)=ampl_ff(i)*sinx(l1:l2,i)
          force(:,3)=0.
!
!  f=(0,cosx*cosz,0), modulated by profy_ampl
!
        case ('(0,cosx*cosz,0)')
          force(:,1)=0.
          force(:,2)=ampl_ff(i)*cosx(l1:l2,i)*cosz(n,i)*profy_ampl(m)
          force(:,3)=0.
!
!  f=(0,sinx*exp(-z^2),0)
!
        case ('(0,sinx*exp(-z^2),0)')
          force(:,1)=0.
          force(:,2)=ampl_ff(i)*sinx(l1:l2,i)*exp(-((z(n)-r_ff)/width_ff)**2) &
                     *(step(x(l1:l2),xminf,2.)-step(x(l1:l2),xmaxf,2.))
          force(:,3)=0.
!
!  f=(0,Aycont_z,0)
!  This ensures vanishing Ay at both boundaries if Bslope=-2Bconst/Lz
!  (making it suitable for perfect conductor boundary condition)
!
        case ('(0,Aycont_z,0)')
          force(:,1)=0.
          force(:,2)=-ampl_ff(i)*((Bconst*(z(n)-xyz0(3))) &
                      +(0.5*Bslope*((z(n)-xyz0(3))**2.)))
          force(:,3)=0.
!
!  f=(sinz,cosz,0)
!
        case ('(sinz,cosz,0)')
          force(:,1)=ampl_ff(i)*sinz(n,i)
          force(:,2)=ampl_ff(i)*cosz(n,i)
          force(:,3)=0.
!
!  Taylor-Green forcing
!
        case ('TG')
          fact=2.*ampl_ff(i)
          force(:,1)=+fact*sinx(l1:l2,i)*cosy(m,i)*cosz(n,i)
          force(:,2)=-fact*cosx(l1:l2,i)*siny(m,i)*cosz(n,i)
          force(:,3)=0.
!
!  Compressive u=-grad(phi) with phi=cos(x+y+z) forcing
!
        case ('cosx*cosy*cosz')
          fact=-ampl_ff(i)
          force(:,1)=fact*sinx(l1:l2,i)*cosy(m,i)*cosz(n,i)
          force(:,2)=fact*cosx(l1:l2,i)*siny(m,i)*cosz(n,i)
          force(:,3)=fact*cosx(l1:l2,i)*cosy(m,i)*sinz(n,i)
!
!  Galloway-Proctor (GP) forcing
!
        case ('GP')
          fact=ampl_ff(i)/sqrt(2.)
          force(:,1)=-fact*                sinyt(m,i)
          force(:,2)=+fact* sinxt(l1:l2,i)
          force(:,3)=-fact*(cosxt(l1:l2,i)+cosyt(m,i))
!
!  Galloway-Proctor (GP) forcing
!
        case ('Galloway-Proctor-92')
          fact=ampl_ff(i)
          ecost=eps_fcont(i)*cos(omega_fcont(i)*t)
          esint=eps_fcont(i)*sin(omega_fcont(i)*t)
          force(:,1)=fact*(sin(kf_fcont(i)*z(n)+esint)+cos(kf_fcont(i)*y(m)+ecost))
          force(:,2)=fact* cos(kf_fcont(i)*z(n)+esint)
          force(:,3)=fact* sin(kf_fcont(i)*y(m)+ecost)
!
!  Galloway-Proctor (GP) forcing in Tobias & Cattaneo (2013)
!
        case ('GP_TC13')
          fact=ampl_ff(i)
          force(:,1) = -fact*sinyt(m,i)
          force(:,2) = -fact*cosxt(l1:l2,i)
          force(:,3) = +fact*( sinxt(l1:l2,i)+cosyt(m,i) )
!
!  Same as GP_TC13 with a rotation of axes y->z->x->y.
!  i.e., the new z axis is the old y axis, etc.
!  This is for the convenience of testfield_z.
!
        case ('GP_TC13_yzx')
          fact=ampl_ff(i)
          force(:,1) = -fact*coszt(n,i)
          force(:,2) = +fact*( sinzt(n,i)+cosxt(l1:l2,i) )
          force(:,3) = -fact*sinxt(l1:l2,i)
!
! Continuous emf required in Induction equation for the Mag Buoy Inst
!
        case('mbi_emf')
          fact=2.*ampl_bb(i)*eta_bb(i)/width_bb(i)**2
          arg=(z(n)-z_bb(i))/width_bb(i)
          force(:,1)=fact*tanh(arg)/(cosh(arg))**2
          force(:,2)=0.0
          force(:,3)=0.0
!
!  Bessel function forcing
!
        case('J0_k1x')
          force(:,1)=0.0
          force(:,2)=profx_ampl1
          force(:,3)=profx_ampl
!
!  Gaussian blob in the z direction
!
        case('gaussian-z')
          force(:,1)=0.0
          force(:,2)=0.0
          force(:,3)=profx_ampl*profy_ampl(m)*profz_ampl(n)
!
!  fluxring_cylindrical
!
        case('fluxring_cylindrical')
          call calc_fluxring_cylindrical(force,i)
        case('counter_centrifugal')
          call calc_counter_centrifugal(force,i)
        case('vortex')
          call vortex(torus,Omega_vortex,force)
!
!  blob-like disturbance (gradient of gaussian)
!
        case('blob')
          fact=ampl_ff(i)/radius_ff**2
          tmp=fact*exp(-.5*((x(l1:l2)-location_fixed(1,1))**2 &
                           +(y(m)    -location_fixed(2,1))**2 &
                           +(z(n)    -location_fixed(3,1))**2)/radius_ff**2)
          force(:,1)=tmp*(x(l1:l2)-location_fixed(1,1))
          force(:,2)=tmp*(y(m)    -location_fixed(2,1))
          force(:,3)=tmp*(z(n)    -location_fixed(3,1))
!
!  blob-like disturbance (gradient of gaussian)
!
        case('zblob')
          tmp=ampl_ff(i)*exp(-.5*(x(l1:l2)**2+y(m)**2+z(n)**2)/radius_ff**2)
          force(:,1)=0.
          force(:,2)=0.
          force(:,3)=tmp
!
!  blob-like vertical field patch
!
        case('vert_field_blob')
          tmp=x(l1:l2)**2+y(m)**2+z(n)**2
          force(:,1)=0.
          where (tmp<=1.)
            force(:,2)=ampl_ff(i)*x(l1:l2)
          elsewhere
            force(:,2)=0.
          endwhere
          force(:,3)=0.
!
!  KS-flow
!
      case ('KS')
        force=0.
        do modeN=1,KS_modes  ! sum over KS_modes modes
          kdotxwt=KS_k(1,modeN)*x(l1:l2)+(KS_k(2,modeN)*y(m)+KS_k(3,modeN)*z(n))+KS_omega(modeN)*t
          cos_kdotxwt=cos(kdotxwt) ;  sin_kdotxwt=sin(kdotxwt)
          force(:,1) = force(:,1) + cos_kdotxwt*KS_A(1,modeN) + &
                                    sin_kdotxwt*KS_B(1,modeN)
          force(:,2) = force(:,2) + cos_kdotxwt*KS_A(2,modeN) + &
                                    sin_kdotxwt*KS_B(2,modeN)
          force(:,3) = force(:,3) + cos_kdotxwt*KS_A(3,modeN) + &
                                    sin_kdotxwt*KS_B(3,modeN)
        enddo
        force=ampl_ff(i)*force
!
!  possibility of putting zero, e.g., for purely magnetic forcings
!
        case('zero')
          force=0.
!
!  nothing (but why not?)
!
        case ('nothing')
          if (i==1) call stop_it('forcing: no valid continuous iforcing_cont specified')
          return
!
        case default
          call stop_it('forcing: no valid continuous iforcing_cont specified')
        endselect
!
    endsubroutine forcing_cont
!***********************************************************************
    subroutine calc_diagnostics_forcing(p)
!
!  add a continuous forcing term (used to be in hydro.f90)
!
!  17-sep-06/axel: coded
!
      use Diagnostics
      use Sub
!
      type (pencil_case), intent(IN) :: p

      real, dimension (nx,3) :: fxb
      real, dimension (nx) :: tmp
!
!  diagnostics
!
      if (ldiagnos) then
        if (idiag_rufm/=0) then
          call dot_mn(p%uu,p%fcont(:,:,1),tmp)
          call sum_mn_name(p%rho*tmp,idiag_rufm)
        endif
!
        if (idiag_rufint/=0) then
          call dot_mn(p%uu,p%fcont(:,:,1),tmp)
          call integrate_mn_name(p%rho*tmp,idiag_rufint)
        endif
!
        if (idiag_ufm/=0) then
          call dot_mn(p%uu,p%fcont(:,:,1),tmp)
          call sum_mn_name(tmp,idiag_ufm)
        endif
!
        if (idiag_ofm/=0) then
          call dot_mn(p%oo,p%fcont(:,:,1),tmp)
          call sum_mn_name(tmp,idiag_ofm)
        endif
!
        if (idiag_qfm/=0) then
          call dot_mn(p%curlo,p%fcont(:,:,1),tmp)
          call sum_mn_name(tmp,idiag_qfm)
        endif
!
        if (lmagnetic) then
          if (idiag_fxbxm/=0.or.idiag_fxbym/=0.or.idiag_fxbzm/=0) then
            call cross(p%fcont(:,:,2),p%bb,fxb)
            call sum_mn_name(fxb(:,1),idiag_fxbxm)
            call sum_mn_name(fxb(:,2),idiag_fxbym)
            call sum_mn_name(fxb(:,3),idiag_fxbzm)
          endif
        endif
      endif
!
    endsubroutine calc_diagnostics_forcing
!***********************************************************************
    subroutine read_forcing_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=forcing_run_pars, IOSTAT=iostat)
!
    endsubroutine read_forcing_run_pars
!***********************************************************************
    subroutine write_forcing_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=forcing_run_pars)
!
    endsubroutine write_forcing_run_pars
!***********************************************************************
    subroutine input_persist_forcing_id(id,done)
!
!  Read in the stored time of the next SNI
!
!  21-dec-05/tony: coded
!  13-Dec-2011/Bourdin.KIS: reworked
!
      use IO, only: read_persist
!
      integer :: id
      logical :: done
!
      select case (id)
        case (id_record_FORCING_LOCATION)
          if (read_persist ('FORCING_LOCATION', location)) return
          done = .true.
        case (id_record_FORCING_TSFORCE)
          if (read_persist ('FORCING_TSFORCE', tsforce)) return
          done = .true.
        case (id_record_FORCING_TORUS)
          if (read_persist ('FORCING_TORUS', torus)) return
          done = .true.
      endselect
!
      if (lroot) print *, 'input_persist_forcing: ', location, tsforce
!
    endsubroutine input_persist_forcing_id
!***********************************************************************
    subroutine input_persist_forcing
!
!  Read in the persistent forcing variables.
!
!  11-Oct-2019/PABourdin: coded
!
      use IO, only: read_persist
!
      logical :: error
!
      error = read_persist ('FORCING_LOCATION', location)
      if (lroot .and. .not. error) print *, 'input_persist_forcing: location: ', location
!
      error = read_persist ('FORCING_TSFORCE', tsforce)
      if (lroot .and. .not. error) print *, 'input_persist_forcing: tsforce: ', tsforce
!
      !error = read_persist ('FORCING_TORUS', torus)
      !if (lroot .and. .not. error) print *, 'input_persist_forcing: torus: ', torus
!
    endsubroutine input_persist_forcing
!***********************************************************************
    logical function output_persistent_forcing()
!
!  This is used, for example, for forcing functions with temporal
!  memory, such as in the paper by Mee & Brandenburg (2006, MNRAS)
!
!  21-dec-05/tony: coded
!  13-Dec-2011/Bourdin.KIS: reworked
!
      use IO, only: write_persist
!
      if (ip<=6.and.lroot .and. (tsforce>=0.)) &
          print *, 'output_persistent_forcing: ', location, tsforce
!
!  write details
!
      output_persistent_forcing = .true.
!
      if (write_persist ('FORCING_LOCATION', id_record_FORCING_LOCATION, location)) return
      if (write_persist ('FORCING_TSFORCE', id_record_FORCING_TSFORCE, tsforce)) return
      if (torus%thick>0.) then
        if (write_persist ('FORCING_TORUS', id_record_FORCING_TORUS, torus)) return
      endif
!
      output_persistent_forcing = .false.
!
    endfunction output_persistent_forcing
!***********************************************************************
    subroutine rprint_forcing(lreset,lwrite)
!
!  reads and registers print parameters relevant for hydro part
!
!  26-jan-04/axel: coded
!
      use Diagnostics
!
      integer :: iname
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
        idiag_rufm=0; idiag_rufint=0; idiag_ufm=0; idiag_ofm=0; idiag_qfm=0; idiag_ffm=0
        idiag_ruxfxm=0; idiag_ruyfym=0; idiag_ruzfzm=0
        idiag_ruxfym=0; idiag_ruyfxm=0
        idiag_fxbxm=0; idiag_fxbym=0; idiag_fxbzm=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if (lroot.and.ip<14) print*,'rprint_forcing: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'rufm',idiag_rufm)
        call parse_name(iname,cname(iname),cform(iname),'rufint',idiag_rufint)
        call parse_name(iname,cname(iname),cform(iname),'ruxfxm',idiag_ruxfxm)
        call parse_name(iname,cname(iname),cform(iname),'ruxfym',idiag_ruxfym)
        call parse_name(iname,cname(iname),cform(iname),'ruyfxm',idiag_ruyfxm)
        call parse_name(iname,cname(iname),cform(iname),'ruyfym',idiag_ruyfym)
        call parse_name(iname,cname(iname),cform(iname),'ruzfzm',idiag_ruzfzm)
        call parse_name(iname,cname(iname),cform(iname),'ufm',idiag_ufm)
        call parse_name(iname,cname(iname),cform(iname),'ofm',idiag_ofm)
        call parse_name(iname,cname(iname),cform(iname),'qfm',idiag_qfm)
        call parse_name(iname,cname(iname),cform(iname),'ffm',idiag_ffm)
        call parse_name(iname,cname(iname),cform(iname),'fxbxm',idiag_fxbxm)
        call parse_name(iname,cname(iname),cform(iname),'fxbym',idiag_fxbym)
        call parse_name(iname,cname(iname),cform(iname),'fxbzm',idiag_fxbzm)
      enddo
!
!  write column where which forcing variable is stored
!
      if (lwr) then
!
      endif
!
    endsubroutine rprint_forcing
!***********************************************************************
    subroutine forcing_clean_up
!
!  deallocate the variables needed for forcing
!
!   12-aug-09/dhruba: coded
!
      if (iforce=='chandra-kendall' .or. iforce=='cktest') then
        deallocate(psif,cklist)
        if (lfastCK) deallocate(Zpsi_list,RYlm_list,IYlm_list)
      endif
!
    endsubroutine forcing_clean_up
!***********************************************************************
    subroutine pushdiags2c(p_diag)

    integer, parameter :: n_diags=0
    integer(KIND=ikind8), dimension(:) :: p_diag

    call keep_compiler_quiet(p_diag)

    endsubroutine pushdiags2c
!***********************************************************************
    subroutine pushpars2c(p_par)

    use Syscalls, only: copy_addr

    integer, parameter :: n_pars=9
    integer(KIND=ikind8), dimension(n_pars) :: p_par

    call copy_addr(k1_ff,p_par(1))
    call copy_addr(tforce_stop,p_par(2))
    call copy_addr(iforcing_zsym,p_par(3))
    call copy_addr(profx_ampl,p_par(4))  ! (nx)
    call copy_addr(profy_ampl,p_par(5))  ! (my)
    call copy_addr(profz_ampl,p_par(6))  ! (mz)
    call copy_addr(profx_hel,p_par(7))   ! (nx)
    call copy_addr(profy_hel,p_par(8))   ! (my)
    call copy_addr(profz_hel,p_par(9))   ! (mz)

    endsubroutine pushpars2c
!*******************************************************************
endmodule Forcing
