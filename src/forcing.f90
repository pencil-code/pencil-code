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
! PENCILS PROVIDED fcont(3)
!
!***************************************************************
!
module Forcing
!
  use Cdata
  use General
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'record_types.h'
  include 'forcing.h'
!
  real :: force=0.,force2=0., force1_scl=1., force2_scl=1.
  real :: relhel=1.,height_ff=0.,r_ff=0.,rcyl_ff=0.
  real :: fountain=1.,width_ff=.5,nexp_ff=1.
  real :: crosshel=0.
  real :: dforce=0.,radius_ff=0.,k1_ff=1.,slope_ff=0.,work_ff=0.
  real :: omega_ff=1.
  real :: tforce_stop=impossible,tforce_stop2=impossible
  real :: tforce_start=0.,tforce_start2=0.
  real :: wff_ampl=0.,xff_ampl=0.,zff_ampl=0.,zff_hel=0.,max_force=impossible
  real :: dtforce=0., dtforce_duration=-1.0, force_strength=0.
  real, dimension(3) :: force_direction=(/0.,0.,0./)
  real, dimension(3) :: location_fixed=(/0.,0.,0./)
  real, dimension(nx) :: profx_ampl=1.,profx_hel=1., profx_ampl1=0.
  real, dimension(my) :: profy_ampl=1.,profy_hel=1.
  real, dimension(mz) :: profz_ampl=1.,profz_hel=1.
  integer :: kfountain=5,ifff,iffx,iffy,iffz,i2fff,i2ffx,i2ffy,i2ffz
  integer :: itestflow_forcing_offset=0,itestfield_forcing_offset=0
  logical :: lwork_ff=.false.,lmomentum_ff=.false.
  logical :: lhydro_forcing=.true.,lmagnetic_forcing=.false.
  logical :: lcrosshel_forcing=.false.,ltestfield_forcing=.false.,ltestflow_forcing=.false.
  logical :: lxxcorr_forcing=.false., lxycorr_forcing=.false.
  logical :: lhelical_test=.false.,lrandom_location=.true.
  logical :: lwrite_psi=.false.
  logical :: lscale_kvector_tobox=.false.,lwrite_gausspot_to_file=.false.
  logical :: lwrite_gausspot_to_file_always=.false.
  logical :: lscale_kvector_fac=.false.
  logical :: lforce_peri=.false., lforce_cuty=.false.
  logical :: lforcing2_same=.false., lforcing2_curl=.false.
  real :: scale_kvectorx=1.,scale_kvectory=1.,scale_kvectorz=1.
  logical :: old_forcing_evector=.false.
  character (len=labellen) :: iforce='zero', iforce2='zero'
  character (len=labellen) :: iforce_profile='nothing'
  character (len=labellen) :: iforce_tprofile='nothing'
  real :: equator=0.
! For helical forcing in spherical polar coordinate system
  real,allocatable,dimension(:,:,:) :: psif
  real,allocatable,dimension(:,:) :: cklist
  logical :: lfastCK=.false.,lsamesign=.true.
! allocated only if we have lfastCK=T
  real,allocatable,dimension(:,:,:) :: Zpsi_list
  real,allocatable,dimension(:,:,:) :: RYlm_list,IYlm_list
  integer :: helsign=0,nlist_ck=25
  real :: fpre = 1.0,ck_equator_gap=0.,ck_gap_step=0.
  integer :: icklist,jtest_aa0=5,jtest_uu0=1
! For random forcing in 2d
  integer,allocatable, dimension (:,:) :: random2d_kmodes
  integer :: random2d_nmodes
  integer :: random2d_kmin,random2d_kmax
! continious 2d forcing
  integer :: k2d
  logical :: l2dxz,l2dyz
! Persistent stuff
  real :: tsforce=-10.
  real, dimension (3) :: location
!
!  continuous forcing variables
!
  logical :: lembed=.false.,lshearing_adjust_old=.false.
  logical :: lgentle=.false.
  character (len=labellen) :: iforcing_cont='ABC'
  real :: ampl_ff=1., ampl1_ff=0., width_fcont=1., x1_fcont=0., x2_fcont=0.
  real :: kf_fcont=0.,omega_fcont=0.,eps_fcont=0.
  real :: tgentle=0.
  real :: ampl_bb=5.0e-2,width_bb=0.1,z_bb=0.1,eta_bb=1.0e-4
  real :: fcont_ampl
!
!  auxiliary functions for continuous forcing function
!
  real, dimension (my) :: phi1_ff
  real, dimension (mx) :: phi2_ff
  real, dimension (mx) :: sinx,cosx,embedx
  real, dimension (my) :: siny,cosy,embedy
  real, dimension (mz) :: sinz,cosz,embedz
!
  integer :: dummy              ! We cannot define empty namelists
  namelist /forcing_init_pars/ dummy
!
  namelist /forcing_run_pars/ &
       tforce_start,tforce_start2,&
       iforce,force,relhel,crosshel,height_ff,r_ff,rcyl_ff,width_ff,nexp_ff, &
       iforce2, force2, force1_scl, force2_scl, &
       kfountain,fountain,tforce_stop,tforce_stop2, &
       dforce,radius_ff,k1_ff,slope_ff,work_ff,lmomentum_ff, &
       omega_ff,location_fixed,lrandom_location, &
       lwrite_gausspot_to_file,lwrite_gausspot_to_file_always, &
       wff_ampl,xff_ampl,zff_ampl,zff_hel, &
       lhydro_forcing,lmagnetic_forcing,lcrosshel_forcing,ltestfield_forcing, &
       lxxcorr_forcing, lxycorr_forcing, &
       ltestflow_forcing,jtest_aa0,jtest_uu0, &
       max_force,dtforce,dtforce_duration,old_forcing_evector, &
       iforce_profile, iforce_tprofile, lscale_kvector_tobox, &
       force_direction, force_strength, lhelical_test, &
       lfastCK,fpre,helsign,nlist_ck,lwrite_psi,&
       ck_equator_gap,ck_gap_step,&
       lforcing_cont,iforcing_cont, &
       lembed, k1_ff, ampl_ff, ampl1_ff, width_fcont, x1_fcont, x2_fcont, &
       kf_fcont,omega_fcont,eps_fcont,lsamesign,&
       lshearing_adjust_old,equator,&
       lscale_kvector_fac,scale_kvectorx,scale_kvectory,scale_kvectorz, &
       lforce_peri,lforce_cuty,lforcing2_same,lforcing2_curl, &
       tgentle,random2d_kmin,random2d_kmax,l2dxz,l2dyz,k2d,&
       z_bb,width_bb,eta_bb,fcont_ampl
! other variables (needs to be consistent with reset list below)
  integer :: idiag_rufm=0, idiag_ufm=0, idiag_ofm=0, idiag_ffm=0
  integer :: idiag_ruxfxm=0, idiag_ruyfym=0, idiag_ruzfzm=0
  integer :: idiag_ruxfym=0, idiag_ruyfxm=0
  integer :: idiag_fxbxm=0, idiag_fxbym=0, idiag_fxbzm=0
!
  contains
!
!***********************************************************************
    subroutine register_forcing()
!
!  add forcing in timestep()
!  11-may-2002/wolf: coded
!
      use Mpicomm
      use Sub
!
!  identify version number
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_forcing
!***********************************************************************
    subroutine initialize_forcing(lstarting)
!
!  read seed field parameters
!  nothing done from start.f90 (lstarting=.true.)
!
      use General, only: bessj
      use Mpicomm, only: stop_it
      use Sub, only: inpui,step_scalar,erfunc
      real :: zstar
      integer :: l
!
      logical :: lstarting
!
      if (lstarting) then
        if (ip<4) print*,'initialize_forcing: not needed in start'
      else
!
!  check whether we want ordinary hydro forcing or magnetic forcing
!
        if (lmagnetic_forcing) then
          ifff=iaa; iffx=iax; iffy=iay; iffz=iaz
        elseif (lhydro_forcing) then
          ifff=iuu; iffx=iux; iffy=iuy; iffz=iuz
        elseif (lcrosshel_forcing) then
          ifff=iaa; iffx=iax; iffy=iay; iffz=iaz
          i2fff=iuu; i2ffx=iux; i2ffy=iuy; i2ffz=iuz
        else
          call stop_it("initialize_forcing: No forcing function set")
        endif
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
!  vertical profiles for amplitude and helicity of the forcing
!  default is constant profiles for rms velocity and helicity.
!
      if (iforce_profile=='nothing') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        profz_ampl=1.; profz_hel=1.
      elseif (iforce_profile=='equator') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        profz_ampl=1.
        do n=1,mz
          profz_hel(n)=sin(2*pi*z(n)/Lz)
        enddo
      elseif (iforce_profile=='equator_step') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        profz_ampl=1.
        do n=1,mz
          profz_hel(n)= -1.+2.*step_scalar(z(n),equator-ck_equator_gap,ck_gap_step)
        enddo
!
!  step function change in intensity of helicity at zff_hel
!
      elseif (iforce_profile=='step_ampl=z') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.; profy_hel=1.
                       profz_hel=1.
        do n=1,mz
          profz_ampl(n)=step_scalar(z(n),zff_hel,width_ff)
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
          profy_hel(m)= -1.+2.*step_scalar(y(m),yequator-ck_equator_gap,ck_gap_step)
        enddo
!
!  sign change of helicity about z=0
!
      elseif (iforce_profile=='equator_hel=z') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        profz_ampl=1.
        do n=1,mz
          profz_hel(n)=z(n)
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
!  turn off forcing intensity above z=0
!
      elseif (iforce_profile=='surface_z') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        profz_ampl=.5*(1.-erfunc(z/width_ff))
        profz_hel=1.
!
! turn on forcing in the bulk of the convection zone
     elseif (iforce_profile=='forced_convection') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        profz_ampl=.5*(1.+ erfunc((z)/width_ff)) - 0.5*(1.+ erfunc((z-1.)/(2.*width_ff)))
        profz_hel=1.
!
!  turn off forcing intensity above x=x0, and
!  cosy profile of helicity
!
      elseif (iforce_profile=='surface_x_cosy') then
        profx_ampl=.5*(1.-erfunc((x-r_ff)/width_ff))
        profx_hel=1.
        profy_ampl=1.
        do m=1,my
          profy_hel(m)=cos(y(m))
        enddo
        profz_ampl=1.; profz_hel=1.
!
!  turn off forcing intensity above x=x0, and
!  stepy profile of helicity
!
      elseif (iforce_profile=='surface_x_stepy') then
        profx_ampl=.5*(1.-erfunc((x-r_ff)/width_ff))
        profx_hel=1.
        profy_ampl=1.
        do m=1,my
          profy_hel(m)= -1.+2.*step_scalar(y(m),pi/2.,width_ff)
        enddo
        profz_ampl=1.; profz_hel=1.
!
!  turn off forcing intensity above x=x0
!
      elseif (iforce_profile=='surface_x') then
        profx_ampl=.5*(1.-erfunc((x-r_ff)/width_ff))
        profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        profz_ampl=1.; profz_hel=1.
!
!  turn off forcing intensity above x=r_ff and y=r_ff
!
      elseif (iforce_profile=='surface_xy') then
        profx_ampl=.5*(1.-erfunc((x-r_ff)/width_ff)); profx_hel=1.
        profy_ampl=.5*(1.-erfunc((y-r_ff)/width_ff)); profy_hel=1.
        profz_ampl=1.; profz_hel=1.
!
!  turn off forcing intensity above x=r_ff and y=r_ff
!
      elseif (iforce_profile=='surface_xz') then
        profx_ampl=.5*(1.-erfunc((x-r_ff)/width_ff)); profx_hel=1.
        profz_ampl=.5*(1.-erfunc((z-r_ff)/width_ff)); profz_hel=1.
        profy_ampl=1.; profy_hel=1.
!
!  turn on forcing intensity above x=r_ff
!
      elseif (iforce_profile=='above_x0') then
        profx_ampl=.5*(1.+erfunc((x-r_ff)/width_ff))
        profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        profz_ampl=1.; profz_hel=1.
!
!  turn on forcing intensity above x=x0 with cosy profile
!
      elseif (iforce_profile=='above_x0_cosy') then
        profx_ampl=.5*(1.+erfunc((x-r_ff)/width_ff))
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
      elseif (iforce_profile=='galactic') then
        profx_ampl=1.; profx_hel=1.
        profy_ampl=1.; profy_hel=1.
        do n=1,mz
          if (abs(z(n))<zff_ampl) profz_ampl(n)=.5*(1.-cos(z(n)))
          if (abs(z(n))<zff_hel ) profz_hel (n)=.5*(1.+cos(z(n)/2.))
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
      else
        call fatal_error('initialize_forcing','iforce_profile value does not exist')
      endif
!
!  at the first step, the sin and cos functions are calculated for all
!  x,y,z points and are then saved and used for all subsequent steps
!  and pencils
!
      if (ip<=6) print*,'forcing_cont:','lforcing_cont=',lforcing_cont,iforcing_cont
      if (iforcing_cont=='ABC') then
        if (lroot) print*,'forcing_cont: ABC--calc sinx, cosx, etc'
        sinx=sin(k1_ff*x); cosx=cos(k1_ff*x)
        siny=sin(k1_ff*y); cosy=cos(k1_ff*y)
        sinz=sin(k1_ff*z); cosz=cos(k1_ff*z)
      elseif (iforcing_cont=='RobertsFlow' .or. iforcing_cont=='RobertsFlow_exact' ) then
        if (lroot) print*,'forcing_cont: Roberts Flow'
        sinx=sin(k1_ff*x); cosx=cos(k1_ff*x)
        siny=sin(k1_ff*y); cosy=cos(k1_ff*y)
      elseif (iforcing_cont=='RobertsFlow-zdep') then
        if (lroot) print*,'forcing_cont: z-dependent Roberts Flow'
        sinx=sin(k1_ff*x); cosx=cos(k1_ff*x)
        siny=sin(k1_ff*y); cosy=cos(k1_ff*y)
      elseif (iforcing_cont=='nocos') then
        if (lroot) print*,'forcing_cont: nocos flow'
        sinx=sin(k1_ff*x)
        siny=sin(k1_ff*y)
        sinz=sin(k1_ff*z)
     elseif (iforcing_cont=='KolmogorovFlow-x') then
        if (lroot) print*,'forcing_cont: Kolmogorov flow'
        cosx=cos(k1_ff*x)
     elseif (iforcing_cont=='KolmogorovFlow-z') then
        if (lroot) print*,'forcing_cont: Kolmogorov flow z'
        cosz=cos(k1_ff*z)
      elseif (iforcing_cont=='TG') then
        if (lroot) print*,'forcing_cont: TG'
        sinx=sin(k1_ff*x); cosx=cos(k1_ff*x)
        siny=sin(k1_ff*y); cosy=cos(k1_ff*y)
        cosz=cos(k1_ff*z)
        if (lembed) then
          embedx=.5+.5*tanh(x/width_fcont)*tanh((pi-x)/width_fcont)
          embedy=.5+.5*tanh(y/width_fcont)*tanh((pi-y)/width_fcont)
          embedz=.5+.5*tanh(z/width_fcont)*tanh((pi-z)/width_fcont)
          sinx=embedx*sinx; cosx=embedx*cosx
          siny=embedy*siny; cosy=embedy*cosy
          cosz=embedz*cosz
        endif
      elseif (iforcing_cont=='sinx') then
        sinx=sin(k1_ff*x)
        if (tgentle > 0.) then
          lgentle=.true.
          if (lroot) print *, 'initialize_forcing: gentle forcing till t = ', tgentle
        endif
      elseif (iforcing_cont=='(0,0,cosx)') then
        cosx=cos(k1_ff*x)
      elseif (iforcing_cont=='J0_k1x') then
        do l=l1,l2
          profx_ampl(l-l1+1)=ampl_ff*bessj(0,k1bessel0*x(l))
          profx_ampl1(l-l1+1)=ampl1_ff*bessj(1,k1bessel0*x(l))
        enddo
      elseif (iforcing_cont=='fluxring_cylindrical') then
        if (lroot) print*,'forcing_cont: fluxring cylindrical'
        !nothing...
      endif
!
    endsubroutine initialize_forcing
!***********************************************************************
    subroutine addforce(f)
!
!  add forcing at the end of each time step
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
      if ( (t>tforce_stop) .or. (t<tforce_start) ) then
        if ( (t>tforce_stop) .and. llastforce .and. lroot) &
            print*, 'addforce: t>tforce_stop; no forcing'
        if (t>tforce_stop) llastforce=.false.
      else
        if ( iforce/='zero' .and. lfirstforce .and. lroot ) &
!--         print*, 'addforce: addforce started'
        lfirstforce=.false.
!
!  calculate and add forcing function
!
        select case (iforce)
        case ('2drandom_xy');     call forcing_2drandom_xy(f)
        case ('ABC');             call forcing_ABC(f)
        case ('blobs');           call forcing_blobs(f)
        case ('chandra_kendall'); call forcing_chandra_kendall(f)
        case ('cktest');          call forcing_cktest(f)
        case ('diffrot');         call forcing_diffrot(f,force)
        case ('fountain', '3');   call forcing_fountain(f)
        case ('gaussianpot');     call forcing_gaussianpot(f,force)
        case ('GP');              call forcing_GP(f)
        case ('irrotational');    call forcing_irro(f,force)
        case ('helical', '2');    call forcing_hel(f)
        case ('helical_both');    call forcing_hel_both(f)
        case ('helical_kprof');   call forcing_hel_kprof(f)
        case ('hel_smooth');      call forcing_hel_smooth(f)
        case ('horiz-shear');     call forcing_hshear(f)
        case ('nocos');           call forcing_nocos(f)
        case ('TG');              call forcing_TG(f)
        case ('twist');           call forcing_twist(f)
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
      integer, save :: ifirst=0
      integer :: ikmodes,iran1,iran2,kx1,ky1,kx2,ky2
      real, dimension(nx,3) :: forcing_rhs
      real, dimension(nx) :: xkx1,xkx2
      real,dimension(4) :: fran
      real :: phase1,phase2,pi_over_Lx,force_norm
!
      if (ifirst==0) then
! If this is the first time this routine is being called select a set of
! k-vectors in two dimensions according to input parameters:
        if (lroot) print*,'forcing_2drandom_xy: selecting k vectors'
        call get_2dmodes (.true.)
        allocate(random2d_kmodes (2,random2d_nmodes))
        call get_2dmodes (.false.)
        if (lroot) then
! The root processors also write out the forced modes.
          open(unit=10,file=trim(datadir)//'2drandomk.out',status='unknown')
          do ikmodes = 1, random2d_nmodes
            write(10,*) random2d_kmodes(1,ikmodes),random2d_kmodes(2,ikmodes)
          enddo
        endif
      endif
      ifirst=ifirst+1
!
! force = xhat [ cos (k_1 x + \phi_1 ) + cos (k_2 y + \phi_1) ] +
!         yhat [ cos (k_3 x + \phi_2 ) + cos (k_4 y + \phi_2) ]
! where k_1 -- k_2 and k_3 -- k_4 are two randomly chose pairs in the list
! random2d_kmodes and  \phi_1 and \phi_2 are two random phases.
!
      call random_number_wrapper(fran)
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
!   15-feb-2011 : dhruba/coded
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
    subroutine forcing_irro(f,force_ampl)
!
!  add acoustic forcing function, using a set of precomputed wavevectors
!  This forcing drives pressure waves
!
!  10-sep-01/axel: coded
!
      use Mpicomm, only: mpifinalize
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: force_ampl
!
      real :: phase,ffnorm
      real, save :: kav
      real, dimension (2) :: fran
      real, dimension (nx) :: rho1
      complex, dimension (mx) :: fx
      complex, dimension (my) :: fy
      complex, dimension (mz) :: fz
      complex, dimension (3) :: ikk
      logical, dimension (3), save :: extent
      integer, parameter :: mk=3000
      real, dimension(mk), save :: kkx,kky,kkz
      integer, save :: ifirst=0,nk
      integer :: ik,j,jf
      logical :: lk_dot_dat_exists
!
      if (ifirst==0) then
        if (lroot) print*,'forcing_irro: opening k.dat'
        inquire(FILE="k.dat", EXIST=lk_dot_dat_exists)
        if (lk_dot_dat_exists) then
          open(9,file='k.dat',status='old')
          read(9,*) nk,kav
          if (lroot) print*,'forcing_irro: average k=',kav
          if (nk>mk) then
            if (lroot) print*,'forcing_irro: dimension mk in forcing_irro is insufficient'
            print*,'nk=',nk,'mk=',mk
            call mpifinalize
          endif
          read(9,*) (kkx(ik),ik=1,nk)
          read(9,*) (kky(ik),ik=1,nk)
          read(9,*) (kkz(ik),ik=1,nk)
          close(9)
        else
          call inevitably_fatal_error ('forcing_irro:', &
              'you must give an input k.dat file')
        endif
        extent(1)=nx/=1
        extent(2)=ny/=1
        extent(3)=nz/=1
      endif
      ifirst=ifirst+1
!
      call random_number_wrapper(fran)
      phase=pi*(2*fran(1)-1.)
      ik=nk*(.9999*fran(2))+1
      if (ip<=6) print*,'forcing_irro: ik,phase,kk=',ik,phase,kkx(ik),kky(ik),kkz(ik),dt,ifirst
!
!  need to multiply by dt (for Euler step), but it also needs to be
!  divided by sqrt(dt), because square of forcing is proportional
!  to a delta function of the time difference
!
      ffnorm=force*sqrt(kav/dt)*dt
      fx=exp(cmplx(0.,kkx(ik)*x+phase))*ffnorm
      fy=exp(cmplx(0.,kky(ik)*y))
      fz=exp(cmplx(0.,kkz(ik)*z))
!
      ikk(1)=cmplx(0.,kkx(ik))
      ikk(2)=cmplx(0.,kky(ik))
      ikk(3)=cmplx(0.,kkz(ik))
!
      do j=1,3
        if (extent(j)) then
          jf=j+ifff-1
          do n=n1,n2
          do m=m1,m2
!
!  Can force either velocity (default) or momentum (slightly more physical).
!
            if (lmomentum_ff) then
              if (ldensity_nolog) then
                rho1=1/f(l1:l2,m,n,ilnrho)
              else
                rho1=exp(-f(l1:l2,m,n,ilnrho))
              endif
            else
              rho1=1.
            endif
            f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+rho1*real(ikk(j)*fx(l1:l2)*fy(m)*fz(n))
          enddo
          enddo
        endif
      enddo
!
      call keep_compiler_quiet(force_ampl)
!
    endsubroutine forcing_irro
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
!
      use Diagnostics
      use EquationOfState, only: cs0
      use General
      use Mpicomm
      use Sub
!
      real :: phase,ffnorm,irufm,iruxfxm,iruxfym,iruyfxm,iruyfym,iruzfzm
      real, save :: kav
      real, dimension (1) :: fsum_tmp,fsum
      real, dimension (2) :: fran
      real, dimension (nx) :: rho1,ruf,rho
      real, dimension (mz) :: tmpz
      real, dimension (nx,3) :: variable_rhs,forcing_rhs,forcing_rhs2
      real, dimension (nx,3) :: force_all
      real, dimension (nx,3) :: fda
      real, dimension (mx,my,mz,mfarray) :: f
      complex, dimension (mx) :: fx
      complex, dimension (my) :: fy
      complex, dimension (mz) :: fz
      real, dimension (3) :: coef1,coef2,coef3
      logical, dimension (3), save :: extent
!      integer, parameter :: mk=3000
      integer, parameter :: mk=6000
      real, dimension(mk), save :: kkx,kky,kkz
      integer, save :: ifirst=0,nk
      integer :: ik,j,jf,j2f
      real :: kx0,kx,ky,kz,k2,k,force_ampl,pi_over_Lx
      real :: ex,ey,ez,kde,sig,fact,kex,key,kez,kkex,kkey,kkez
      real, dimension(3) :: e1,e2,ee,kk
      real :: norm,phi
      real :: fd,fd2
      logical :: lk_dot_dat_exists
!
!  additional stuff for test fields
!
      if (ifirst==0) then
        if (lroot.and.ip<14) print*,'forcing_hel: opening k.dat'
        inquire(FILE="k.dat", EXIST=lk_dot_dat_exists)
        if (lk_dot_dat_exists) then
          open(9,file='k.dat',status='old')
          read(9,*) nk,kav
          if (lroot.and.ip<14) print*,'forcing_hel: average k=',kav
          if (nk>mk) then
            if (lroot) print*,'forcing_hel: mk in forcing_hel is set too small'
            print*,'nk=',nk,'mk=',mk
            call mpifinalize
          endif
          read(9,*) (kkx(ik),ik=1,nk)
          read(9,*) (kky(ik),ik=1,nk)
          read(9,*) (kkz(ik),ik=1,nk)
          close(9)
        else
          call inevitably_fatal_error ('forcing_hel:', &
              'you must give an input k.dat file')
        endif
        extent(1)=nx/=1
        extent(2)=ny/=1
        extent(3)=nz/=1
      endif
      ifirst=ifirst+1
!
!  generate random coefficients -1 < fran < 1
!  ff=force*Re(exp(i(kx+phase)))
!  |k_i| < akmax
!
      call random_number_wrapper(fran)
      phase=pi*(2*fran(1)-1.)
      ik=nk*(.9999*fran(2))+1
      if (ip<=6) print*,'forcing_hel: ik,phase=',ik,phase
      if (ip<=6) print*,'forcing_hel: kx,ky,kz=',kkx(ik),kky(ik),kkz(ik)
      if (ip<=6) print*,'forcing_hel: dt, ifirst=',dt,ifirst
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
      if (ip<4) write(88,'(6f10.5)') k,kx0,kx,ky,kz,deltay
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
        call random_number_wrapper(phi); phi = phi*2*pi
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
      if (ip<=9) print*,'forcing_hel: k,kde,ffnorm,kav=',k,kde,ffnorm,kav
      if (ip<=9) print*,'forcing_hel: k*sqrt(k2-kde**2)=',k*sqrt(k2-kde**2)
!
!  need to multiply by dt (for Euler step), but it also needs to be
!  divided by sqrt(dt), because square of forcing is proportional
!  to a delta function of the time difference
!
      fact=force/ffnorm*sqrt(dt)
!
!  The wavevector is for the case where Lx=Ly=Lz=2pi. If that is not the
!  case one needs to scale by 2pi/Lx, etc.
!DM: the previous comments seems to be wrong because we have already used
! scale_kvector_to_box option above. Axel, what do you think ?
!
      fx=exp(cmplx(0.,kx*k1_ff*x+phase))*fact
      fy=exp(cmplx(0.,ky*k1_ff*y))
      fz=exp(cmplx(0.,kz*k1_ff*z))
!
!  possibly multiply forcing by z-profile
!  (This stuff is now supposed to be done in initialize; keep for now)
!
!-    if (height_ff/=0.) then
!-      if (lroot .and. ifirst==1) print*,'forcing_hel: include z-profile'
!-      tmpz=(z/height_ff)**2
!-      fz=fz*exp(-tmpz**5/max(1.-tmpz,1e-5))
!-    endif
!
! need to discuss with axel
!
!  possibly multiply forcing by sgn(z) and radial profile
!
       if (rcyl_ff/=0.) then
         if (lroot .and. ifirst==1) &
              print*,'forcing_hel: applying sgn(z)*xi(r) profile'
         !
         ! only z-dependent part can be done here; radial stuff needs to go
!        ! into the loop
         !
         tmpz = tanh(z/width_ff)
         fz = fz*tmpz
       endif
!
      if (ip<=5) print*,'forcing_hel: fx=',fx
      if (ip<=5) print*,'forcing_hel: fy=',fy
      if (ip<=5) print*,'forcing_hel: fz=',fz
!
!  prefactor; treat real and imaginary parts separately (coef1 and coef2),
!  so they can be multiplied by different profiles below.
!
      coef1(1)=k*kex; coef2(1)=relhel*kkex; coef3(1)=crosshel*k*kkex
      coef1(2)=k*key; coef2(2)=relhel*kkey; coef3(2)=crosshel*k*kkey
      coef1(3)=k*kez; coef2(3)=relhel*kkez; coef3(3)=crosshel*k*kkez
      if (ip<=5) print*,'forcing_hel: coef=',coef1,coef2
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
!  In the past we always forced the du/dt, but in some cases
!  it may be better to force rho*du/dt (if lmomentum_ff=.true.)
!  For compatibility with earlier results, lmomentum_ff=.false. by default.
!
      if (ldensity) then
        if (lmomentum_ff) then
          rho1=exp(-f(l1:l2,m,n,ilnrho))
          rho=1./rho1
        else
          rho1=1.
          rho=exp(f(l1:l2,m,n,ilnrho))
        endif
      else
        rho1=1.
        rho=1.
      endif
!
!  loop the two cases separately, so we don't check for r_ff during
!  each loop cycle which could inhibit (pseudo-)vectorisation
!  calculate energy input from forcing; must use lout (not ldiagnos)
!
      force_ampl=1.0
      irufm=0; iruxfxm=0; iruxfym=0; iruyfxm=0; iruyfym=0; iruzfzm=0
      if (rcyl_ff == 0) then       ! no radial profile
        do n=n1,n2
          do m=m1,m2
            if (lwork_ff) call calc_force_ampl(f,fx,fy,fz,profy_ampl(m)*profz_ampl(n) &
              *cmplx(coef1,profy_hel(m)*profz_hel(n)*coef2),force_ampl)
            variable_rhs=f(l1:l2,m,n,iffx:iffz)
            do j=1,3
              if (extent(j)) then
!
!  Primary forcing function.
!
                forcing_rhs(:,j)=rho1*profx_ampl*profy_ampl(m)*profz_ampl(n)*force_ampl &
                  *real(cmplx(coef1(j),profx_hel*profy_hel(m)*profz_hel(n)*coef2(j)) &
                  *fx(l1:l2)*fy(m)*fz(n))*fda(:,j)
!
!  Compute additional forcing function (used for velocity if crosshel=1).
!  It can optionally be the same. Alterantively, one has to set crosshel=1.
!
                if (lforcing2_same) then
                  forcing_rhs2(:,j)=forcing_rhs(:,j)
                elseif (lforcing2_curl) then
                  forcing_rhs2(:,j)=rho1*profx_ampl*profy_ampl(m)*profz_ampl(n)*force_ampl &
                    *real(cmplx(coef3(j),-profx_hel*profy_hel(m)*profz_hel(n)*coef2(j)) &
                    *fx(l1:l2)*fy(m)*fz(n))*fda(:,j)
                else
                  forcing_rhs2(:,j)=rho1*profx_ampl*profy_ampl(m)*profz_ampl(n)*force_ampl &
                    *real(cmplx(0.,coef3(j)) &
                    *fx(l1:l2)*fy(m)*fz(n))*fda(:,j)
                endif
!
!  Choice of different possibilities.
!
                if (ifff/=0) then
                  jf=j+ifff-1
                  j2f=j+i2fff-1
!
                  if (lhelical_test) then
                    f(l1:l2,m,n,jf)=forcing_rhs(:,j)
                  else
                    f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+forcing_rhs(:,j)*force1_scl
!
!  Allow here for forcing both in u and in b=curla. In that case one sets
!  lhydro_forcing=F, lmagnetic_forcing=F, lcrosshel_forcing=T
!
                    if (lcrosshel_forcing) then
                      f(l1:l2,m,n,j2f)=f(l1:l2,m,n,j2f)+forcing_rhs2(:,j)*force2_scl
                    endif
!
!  Forcing with enhanced xx correlation.
!
                    if (lxxcorr_forcing) then
                      if (j==1) f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf) &
                        +forcing_rhs(:,1)*force2_scl
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
                  f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+forcing_rhs(:,j)
                endif
                if (ltestflow_forcing) then
                  jf=j+iuutest-1+3*(jtest_uu0-1)
                  f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+forcing_rhs(:,j)
                endif
              endif
            enddo
!
!  Forcing with enhanced xy correlation.
!  Can only come outside the previous j loop.
!
            do j=1,3
              if (extent(j)) then
                if (ifff/=0) then
                  jf=j+ifff-1
                  j2f=j+i2fff-1
                  if (.not.lhelical_test) then
                    if (lxycorr_forcing) then
                      if (j==1) f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf) &
                        +forcing_rhs(:,2)*force2_scl
                      if (j==2) f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf) &
                        +forcing_rhs(:,1)*force2_scl
                    endif
                  endif
                endif
              endif
            enddo
!
!  Sum up.
!
            if (lout) then
              if (idiag_rufm/=0 .or. &
                  idiag_ruxfxm/=0 .or. idiag_ruxfym/=0 .or. &
                  idiag_ruyfxm/=0 .or. idiag_ruyfym/=0 .or. &
                  idiag_ruzfzm/=0) then
!
!  Compute rhs and density.
!
                variable_rhs=f(l1:l2,m,n,iffx:iffz)
                if (ldensity_nolog) then
                  rho=f(l1:l2,m,n,irho)
                else
                  rho=exp(f(l1:l2,m,n,ilnrho))
                endif
!
!  Evaluate the various choices.
!
                if (idiag_rufm/=0) then
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
        do j=1,3
          if (extent(j)) then
            jf=j+ifff-1
            do n=n1,n2
              sig = relhel*tmpz(n)
call fatal_error('forcing_hel','check that radial profile with rcyl_ff works ok')
              coef1(1)=cmplx(k*kex,sig*kkex)
              coef1(2)=cmplx(k*key,sig*kkey)
              coef1(3)=cmplx(k*kez,sig*kkez)
              do m=m1,m2
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
          fsum_tmp(1)=irufm
          call mpireduce_sum(fsum_tmp,fsum,1)
          irufm=fsum(1)
          call mpibcast_real(irufm,1)
          fname(idiag_rufm)=irufm
          itype_name(idiag_rufm)=ilabel_sum
        endif
        if (idiag_ruxfxm/=0) then
          fsum_tmp(1)=iruxfxm
          call mpireduce_sum(fsum_tmp,fsum,1)
          iruxfxm=fsum(1)
          call mpibcast_real(iruxfxm,1)
          fname(idiag_ruxfxm)=iruxfxm
          itype_name(idiag_ruxfxm)=ilabel_sum
        endif
        if (idiag_ruxfym/=0) then
          fsum_tmp(1)=iruxfym
          call mpireduce_sum(fsum_tmp,fsum,1)
          iruxfym=fsum(1)
          call mpibcast_real(iruxfym,1)
          fname(idiag_ruxfym)=iruxfym
          itype_name(idiag_ruxfym)=ilabel_sum
        endif
        if (idiag_ruyfxm/=0) then
          fsum_tmp(1)=iruyfxm
          call mpireduce_sum(fsum_tmp,fsum,1)
          iruyfxm=fsum(1)
          call mpibcast_real(iruyfxm,1)
          fname(idiag_ruyfxm)=iruyfxm
          itype_name(idiag_ruyfxm)=ilabel_sum
        endif
        if (idiag_ruyfym/=0) then
          fsum_tmp(1)=iruyfym
          call mpireduce_sum(fsum_tmp,fsum,1)
          iruyfym=fsum(1)
          call mpibcast_real(iruyfym,1)
          fname(idiag_ruyfym)=iruyfym
          itype_name(idiag_ruyfym)=ilabel_sum
        endif
        if (idiag_ruzfzm/=0) then
          fsum_tmp(1)=iruzfzm
          call mpireduce_sum(fsum_tmp,fsum,1)
          iruzfzm=fsum(1)
          call mpibcast_real(iruzfzm,1)
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
!
      use Diagnostics
      use EquationOfState, only: cs0
      use General
      use Mpicomm
      use Sub
!
      real :: phase,ffnorm,irufm
      real, save :: kav
      real, dimension (2) :: fran
      real, dimension (nx) :: rho1,ff,uf,of,rho
      real, dimension (mz) :: tmpz,kfscl
      real, dimension (nx,3) :: variable_rhs,forcing_rhs,forcing_rhs2
      real, dimension (nx,3) :: fda,uu,oo,bb,fxb
      real, dimension (mx,my,mz,mfarray) :: f
      complex, dimension (mx) :: fx
      complex, dimension (my) :: fy
      complex, dimension (mz) :: fz
      real, dimension (3) :: coef1,coef2,coef3
      logical, dimension (3), save :: extent
      integer, parameter :: mk=6000
      real, dimension(mk), save :: kkx,kky,kkz
      integer, save :: ifirst=0,nk
      integer :: ik,j,jf,j2f
      real :: kx0,kx,ky,kz,k2,k,force_ampl,pi_over_Lx
      real :: ex,ey,ez,kde,sig,fact,kex,key,kez,kkex,kkey,kkez
      real, dimension(3) :: e1,e2,ee,kk
      real :: norm,phi
      real :: fd,fd2
      logical :: lk_dot_dat_exists
!
!  additional stuff for test fields
!
      if (ifirst==0) then
        if (lroot) print*,'forcing_hel_kprof: opening k.dat'
        inquire(FILE="k.dat", EXIST=lk_dot_dat_exists)
        if (lk_dot_dat_exists) then
          if (lroot.and.ip<14) print*,'forcing_hel_kprof: opening k.dat'
          open(9,file='k.dat',status='old')
          read(9,*) nk,kav
          if (lroot.and.ip<14) print*,'forcing_hel_kprof: average k=',kav
          if (nk>mk) then
            if (lroot) print*,'forcing_hel_kprof: mk in forcing_hel_kprof is set too small'
            print*,'nk=',nk,'mk=',mk
            call mpifinalize
          endif
          read(9,*) (kkx(ik),ik=1,nk)
          read(9,*) (kky(ik),ik=1,nk)
          read(9,*) (kkz(ik),ik=1,nk)
          close(9)
        else
          call inevitably_fatal_error ('forcing_hel_kprof:', &
              'you must give an input k.dat file')
        endif
        extent(1)=nx/=1
        extent(2)=ny/=1
        extent(3)=nz/=1
      endif
      ifirst=ifirst+1
!
!  generate random coefficients -1 < fran < 1
!  ff=force*Re(exp(i(kx+phase)))
!  |k_i| < akmax
!
      call random_number_wrapper(fran)
      phase=pi*(2*fran(1)-1.)
      ik=nk*(.9999*fran(2))+1
      if (ip<=6) print*,'forcing_hel_kprof: ik,phase=',ik,phase
      if (ip<=6) print*,'forcing_hel_kprof: kx,ky,kz=',kkx(ik),kky(ik),kkz(ik)
      if (ip<=6) print*,'forcing_hel_kprof: dt, ifirst=',dt,ifirst
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
      call random_number_wrapper(phi); phi = phi*2*pi
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
      if (ip<4) write(88,'(6f10.5)') k,kx0,kx,ky,kz,deltay
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
      if (ip<=9) print*,'forcing_hel_kprof: k,kde,ffnorm,kav=',k,kde,ffnorm,kav
      if (ip<=9) print*,'forcing_hel_kprof: k*sqrt(k2-kde**2)=',k*sqrt(k2-kde**2)
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
         if (lroot .and. ifirst==1) &
              print*,'forcing_hel_kprof: applying sgn(z)*xi(r) profile'
         !
         ! only z-dependent part can be done here; radial stuff needs to go
!        ! into the loop
         !
         tmpz = tanh(z/width_ff)
         fz = fz*tmpz
       endif
!
      if (ip<=5) print*,'forcing_hel_kprof: fx=',fx
      if (ip<=5) print*,'forcing_hel_kprof: fy=',fy
      if (ip<=5) print*,'forcing_hel_kprof: fz=',fz
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
!  In the past we always forced the du/dt, but in some cases
!  it may be better to force rho*du/dt (if lmomentum_ff=.true.)
!  For compatibility with earlier results, lmomentum_ff=.false. by default.
!
      if (ldensity) then
        if (lmomentum_ff) then
          rho1=exp(-f(l1:l2,m,n,ilnrho))
          rho=1./rho1
        else
          rho1=1.
          rho=exp(f(l1:l2,m,n,ilnrho))
        endif
      else
        rho1=1.
        rho=1.
      endif
!
!  loop the two cases separately, so we don't check for r_ff during
!  each loop cycle which could inhibit (pseudo-)vectorisation
!  calculate energy input from forcing; must use lout (not ldiagnos)
!
      force_ampl=1.0
      irufm=0
      if (rcyl_ff == 0) then       ! no radial profile
          do m=m1,m2
            if (lwork_ff) call calc_force_ampl(f,fx,fy,fz,profy_ampl(m)*profz_ampl(n) &
              *cmplx(coef1,profy_hel(m)*profz_hel(n)*coef2),force_ampl)
            variable_rhs=f(l1:l2,m,n,iffx:iffz)
            do j=1,3
              if (extent(j)) then
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
                    f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+forcing_rhs(:,j)
!
!  allow here for forcing both in u and in b=curla. In that case one sets
!  lhydro_forcing=F, lmagnetic_forcing=F, lcrosshel_forcing=T
!
                    if (lcrosshel_forcing) then
                      f(l1:l2,m,n,j2f)=f(l1:l2,m,n,j2f)+forcing_rhs2(:,j)
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
                  f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+forcing_rhs(:,j)
                endif
                if (ltestflow_forcing) then
                  jf=j+iuutest-1+3*(jtest_uu0-1)
                  f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+forcing_rhs(:,j)
                endif
              endif
            enddo
          enddo
      else
!
!  Radial profile, but this is old fashioned and probably no longer used.
!
        do j=1,3
          if (extent(j)) then
            jf=j+ifff-1
              sig = relhel*tmpz(n)
call fatal_error('forcing_hel_kprof','check that radial profile with rcyl_ff works ok')
              coef1(1)=cmplx(k*kex,sig*kkex)
              coef1(2)=cmplx(k*key,sig*kkey)
              coef1(3)=cmplx(k*kez,sig*kkez)
              do m=m1,m2
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
!  enddo
!
      enddo
!
!  For printouts:
!
      if (lout) then
        if (idiag_rufm/=0) then
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
!
      use Mpicomm
      use General
      use Sub
      use EquationOfState, only: cs0
!
      real :: phase,ffnorm
      real, save :: kav
      real, dimension (2) :: fran
      real, dimension (nx) :: rho1
      real, dimension (nx,3) :: variable_rhs,forcing_rhs
      real, dimension (mx,my,mz,mfarray) :: f
      complex, dimension (mx) :: fx
      complex, dimension (my) :: fy
      complex, dimension (mz) :: fz
      real, dimension (3) :: coef1,coef2
      logical, dimension (3), save :: extent
      integer, parameter :: mk=3000
      real, dimension(mk), save :: kkx,kky,kkz
      integer, save :: ifirst=0,nk
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
      if (ifirst==0) then
        if (lroot) print*,'forcing_hel_both: opening k.dat'
        if (lroot) print*,'Equator',equator
        open(9,file='k.dat',status='old')
        read(9,*) nk,kav
        if (lroot) print*,'forcing_hel_both: average k=',kav
        if (nk>mk) then
          if (lroot) print*,'forcing_hel_both: mk in forcing_hel is set too small'
          print*,'nk=',nk,'mk=',mk
          call mpifinalize
        endif
        read(9,*) (kkx(ik),ik=1,nk)
        read(9,*) (kky(ik),ik=1,nk)
        read(9,*) (kkz(ik),ik=1,nk)
        close(9)
        extent(1)=nx/=1
        extent(2)=ny/=1
        extent(3)=nz/=1
      endif
      ifirst=ifirst+1
!
!  generate random coefficients -1 < fran < 1
!  ff=force*Re(exp(i(kx+phase)))
!  |k_i| < akmax
!
      call random_number_wrapper(fran)
      phase=pi*(2*fran(1)-1.)
      ik=nk*(.9999*fran(2))+1
      if (ip<=6) print*,'forcing_hel_both: ik,phase=',ik,phase
      if (ip<=6) print*,'forcing_hel_both: kx,ky,kz=',kkx(ik),kky(ik),kkz(ik)
      if (ip<=6) print*,'forcing_hel_both: dt, ifirst=',dt,ifirst
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
        call random_number_wrapper(phi); phi = phi*2*pi
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
      if (ip<=9) print*,'forcing_hel_both: k,kde,ffnorm,kav=',k,kde,ffnorm,kav
      if (ip<=9) print*,'forcing_hel_both: k*sqrt(k2-kde**2)=',k*sqrt(k2-kde**2)
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
      if (ip<=5) print*,'forcing_hel: fx=',fx
      if (ip<=5) print*,'forcing_hel: fy=',fy
      if (ip<=5) print*,'forcing_hel: fz=',fz
!
!  prefactor; treat real and imaginary parts separately (coef1 and coef2),
!  so they can be multiplied by different profiles below.
!
      coef1(1)=k*kex; coef2(1)=relhel*kkex
      coef1(2)=k*key; coef2(2)=relhel*kkey
      coef1(3)=k*kez; coef2(3)=relhel*kkez
      if (ip<=5) print*,'forcing_hel: coef=',coef1,coef2
!
!  loop the two cases separately, so we don't check for r_ff during
!  each loop cycle which could inhibit (pseudo-)vectorisation
!  calculate energy input from forcing; must use lout (not ldiagnos)
!
     rho1=1
      do n=n1,n2
        do m=m1,m2
          variable_rhs=f(l1:l2,m,n,iffx:iffz)
          do j=1,3
            if (extent(j)) then
              jf=j+ifff-1
              if (y(m)>equator)then
                forcing_rhs(:,j)=rho1*real(cmplx(coef1(j),coef2(j)) &
                                    *fx(l1:l2)*fy(m)*fz(n))&
                        *(1.-step_scalar(y(m),equator-ck_equator_gap,ck_gap_step)+&
                            step_scalar(y(m),equator+ck_equator_gap,ck_gap_step))
              else
                forcing_rhs(:,j)=rho1*real(cmplx(coef1(j),-coef2(j)) &
                  *fx(l1:l2)*fy(m)*fz(n))&
                  *(1.-step_scalar(y(m),equator-ck_equator_gap,ck_gap_step)+&
                      step_scalar(y(m),equator+ck_equator_gap,ck_gap_step))
              endif
              if (lhelical_test) then
                f(l1:l2,m,n,jf)=forcing_rhs(:,j)
              else
                f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+forcing_rhs(:,j)
              endif
              if (ltestfield_forcing) then
                do jtest=1,12
                  iaxtest=iaatest+3*(jtest-1)
                  jf=j+iaxtest-1
                  f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+forcing_rhs(:,j)
                enddo
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
      use Mpicomm
      use General
      use Sub
      use EquationOfState, only: cs0
!      use SpecialFunctions
!
      integer, save :: ifirst=0
      real, dimension(3) :: ee
      real, dimension(nx,3) :: capitalT,capitalS,capitalH,psi
      real, dimension(nx,3,3) :: psi_ij,Tij
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: emm,l,j,jf,Legendrel,lmindex,ilread,ilm,&
                 aindex,ckno,ilist
      real :: a_ell,anum,adenom,jlm,ylm,rphase1,fnorm,alphar,Balpha,&
              psilm,RYlm,IYlm
      real :: rz,rindex,ralpha,&
              rmin,rmax,rphase2
      real, dimension(mx) :: Z_psi
!
!========================================================
      if (.not. lspherical_coords) call warning('chandra-kendall forcing:','This forcing works only in spherical coordinates!')
      if (ifirst==0) then
! If this is the first time this function is being called allocate \psi.
! Next read from file "alpha_in.dat" the two values of \ell. If this two
! matches Legendrel_min and Legendrel_max proceed.
        if (lroot) print*,'Helical forcing in spherical polar coordinate'
        if (lroot) print*,'allocating psif ..'
        allocate(psif(mx,my,mz))
! Read the list of values for emm, ell and alpha. This code is designed for 25 such
        allocate(cklist(nlist_ck,5))
        if (lroot) print*, '..done'
        open(unit=76,file="alpha_in.dat",status="old")
        read(76,*) ckno,rmin,rmax
        if (.not. (ckno==nlist_ck)) then
          call stop_it("CK forcing aborting:  The list does not match check entries in alpha_in.dat")
        else
        endif
        if (lroot) then
          if (.not.((helsign==1).or.(helsign==-1))) &
            call stop_it("CK forcing: helsign must be +1 or -1, aborting")
        else
        endif
! ----------
        do ilread=1,nlist_ck
          read(76,*) (cklist(ilread,ilm),ilm=1,5)
        enddo
        close(76)
        if (lfastCK) then
          allocate(Zpsi_list(mx,nlist_ck,3))
          allocate(RYlm_list(my,mz,nlist_ck),IYlm_list(my,mz,nlist_ck))
          do ilist=1,nlist_ck
            emm = cklist(ilist,1)
            Legendrel = cklist(ilist,2)
            do n=n1-nghost,n2+nghost
              do m=m1-nghost,m2+nghost
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
              do l=l1-nghost,l2+nghost
                alphar=Balpha*x(l)
                call sp_besselj_l(jlm,Legendrel,alphar)
                call sp_bessely_l(ylm,Legendrel,alphar)
                Zpsi_list(l,ilist,aindex) = (a_ell*jlm+ylm)
              enddo
            enddo
          enddo
        else
        endif
        ifirst= ifirst+1
        if (lroot) write(*,*) 'dhruba: first time in Chandra-Kendall successful'
      else
      endif
! This is designed from 5 emm values and for each one 5 ell values. Total 25 values
   call random_number_wrapper(rindex)
   lmindex=nint(rindex*(nlist_ck-1))+1
   emm = cklist(lmindex,1)
   Legendrel = cklist(lmindex,2)
   call random_number_wrapper(ralpha)
   aindex=nint(ralpha*2)
   Balpha = cklist(lmindex,3+aindex)
! Now calculate the "potential" for the helical forcing. The expression
! is taken from Chandrasekhar and Kendall.
! Now construct the Z_psi(r)
   call random_number_wrapper(rphase1)
   rphase1=rphase1*2*pi
   if (lfastCK) then
     do n=n1-nghost,n2+nghost
       do m=m1-nghost,m2+nghost
         psilm=0.
         psilm= RYlm_list(m,n,lmindex)*cos(rphase1)- &
           IYlm_list(m,n,lmindex)*sin(rphase1)
         psif(:,m,n) = psilm*Zpsi_list(:,lmindex,aindex+1)
         if (ck_equator_gap/=0)&
           psif(:,m,n)=psif(:,m,n)*(1.-step_scalar(y(m),pi/2.-ck_equator_gap,ck_gap_step)+&
             step_scalar(y(m),pi/2+ck_equator_gap,ck_gap_step))
       enddo
     enddo
   else
     call sp_bessely_l(anum,Legendrel,Balpha*x(l1))
     call sp_besselj_l(adenom,Legendrel,Balpha*x(l1))
     a_ell = -anum/adenom
!        write(*,*) 'dhruba:',anum,adenom,Legendrel,Bessel_alpha,x(l1)
     do l=l1-nghost,l2+nghost
       alphar=Balpha*x(l)
       call sp_besselj_l(jlm,Legendrel,alphar)
       call sp_bessely_l(ylm,Legendrel,alphar)
       Z_psi(l) = (a_ell*jlm+ylm)
     enddo
!-------
     do n=n1-nghost,n2+nghost
       do m=m1-nghost,m2+nghost
         psilm=0.
         call sp_harm_real(RYlm,Legendrel,emm,y(m),z(n))
         call sp_harm_imag(IYlm,Legendrel,emm,y(m),z(n))
         psilm= RYlm*cos(rphase1)-IYlm*sin(rphase1)
         psif(:,m,n) = Z_psi*psilm
         if (ck_equator_gap/=0)&
           psif(:,m,n)=psif(:,m,n)*(1.-step_scalar(y(m),pi/2.-ck_equator_gap,ck_gap_step)+&
           step_scalar(y(m),pi/2+ck_equator_gap,ck_gap_step))
       enddo
     enddo
   endif
! ----- Now calculate the force from the potential and add this to
! velocity
! get a random unit vector with three components ee_r, ee_theta, ee_phi
! psi at present is just Z_{ell}^m. We next do a sum over random coefficients
! get random psi.
!      write(*,*) 'mmin=',mmin
!! ----------now generate and add the force ------------
   call random_number_wrapper(rz)
   ee(3) = rz
   call random_number_wrapper(rphase2)
   rphase2 = pi*rphase2
   ee(1) = sqrt(1-rz*rz)*cos(rphase2)
   ee(2) = sqrt(1-rz*rz)*sin(rphase2)
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
       if ((y(m)<pi/2.).or.(lsamesign)) then
         capitalS = float(helsign)*(1./Balpha)*capitalS
       else
         capitalS = -float(helsign)*(1./Balpha)*capitalS
       endif
       capitalH = capitalT + capitalS
       do j=1,3
         jf = iuu+j-1
         if (r_ff /= 0.) then
            capitalH(:,j)=profx_ampl*capitalH(:,j)
         endif
         if (lhelical_test) then
           if (lwrite_psi) then
             f(l1:l2,m,n,jf) = psif(l1:l2,m,n)
           else
             f(l1:l2,m,n,jf) = fnorm*capitalH(:,j)
           endif
       else
! stochastic euler scheme of integration[sqrt(dt) is already included in fnorm]
           f(l1:l2,m,n,jf) = f(l1:l2,m,n,jf)+ fnorm*capitalH(:,j)
         endif
       enddo
     enddo
   enddo
!
    endsubroutine forcing_chandra_kendall
!***********************************************************************
    subroutine forcing_cktest(f)
!
! Testing the Chandrasekhar-Kendall forcing function
!  22-june-08/dhruba: adapted from forcing_chandrasekhar_kendall
!
      use Mpicomm
      use General
      use Sub
      use EquationOfState, only: cs0
!      use SpecialFunctions
!
      integer, save :: ifirst
      real, dimension(3) :: ee
      real, dimension(nx,3) :: capitalT,capitalS,capitalH,psi
      real, dimension(nx,3,3) :: psi_ij,Tij
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: emm,l,j,jf,Legendrel,lmindex,ilread,ilm,&
                 aindex,ckno,ilist
      real :: a_ell,anum,adenom,jlm,ylm,rphase1,fnorm,alphar,Balpha,&
              psilm,RYlm,IYlm
      real :: rz,ralpha,&
              rmin,rmax,rphase2
      real, dimension(mx) :: Z_psi
!
!========================================================
      if (.not. lspherical_coords) call warning('chandra-kendall forcing:','This forcing works only in spherical coordinates!')
      if (ifirst==0) then
! If this is the first time this function is being called allocate \psi.
! Next read from file "alpha_in.dat" the two values of \ell. If this two
! matches Legendrel_min and Legendrel_max proceed.
        if (lroot) print*,'Testing helical forcing in spherical polar coordinate'
        if (lroot) print*,'allocating psif ..'
        allocate(psif(mx,my,mz))
! Read the list of values for emm, ell and alpha. This code is designed for 25 such
        allocate(cklist(nlist_ck,5))
        if (lroot) print*, '..done'
        open(unit=76,file="alpha_in.dat",status="old")
        read(76,*) ckno,rmin,rmax
        if (.not. (ckno==nlist_ck)) then
          call stop_it("CK forcing aborting:  The list does not match check entries in alpha_in.dat")
        else
        endif
        if (lroot) then
          if (.not.((helsign==1).or.(helsign==-1))) &
            call stop_it("CK forcing: helsign must be +1 or -1, aborting")
        else
        endif
! ----------
        do ilread=1,nlist_ck
          read(76,*) (cklist(ilread,ilm),ilm=1,5)
        enddo
        close(76)
        if (lfastCK) then
          allocate(Zpsi_list(mx,nlist_ck,3))
          allocate(RYlm_list(my,mz,nlist_ck),IYlm_list(my,mz,nlist_ck))
          do ilist=1,nlist_ck
            emm = cklist(ilist,1)
            Legendrel = cklist(ilist,2)
            do n=n1-nghost,n2+nghost
              do m=m1-nghost,m2+nghost
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
              do l=l1-nghost,l2+nghost
                alphar=Balpha*x(l)
                call sp_besselj_l(jlm,Legendrel,alphar)
                call sp_bessely_l(ylm,Legendrel,alphar)
                Zpsi_list(l,ilist,aindex) = (a_ell*jlm+ylm)
              enddo
            enddo
          enddo
        else
        endif
        icklist=0
        ifirst= ifirst+1
        if (lroot) write(*,*) 'dhruba: first time in Chandra-Kendall successful'
      else
      endif
! This is designed from 5 emm values and for each one 5 ell values. Total 25 values
   icklist=icklist+1
   if (icklist==(nlist_ck+1)) &
            call stop_it("CK testing: no more value in list; ending")
   lmindex=icklist
   emm = cklist(lmindex,1)
   Legendrel = cklist(lmindex,2)
   call random_number_wrapper(ralpha)
   aindex=nint(ralpha*2)
   Balpha = cklist(lmindex,3)
! Now calculate the "potential" for the helical forcing. The expression
! is taken from Chandrasekhar and Kendall.
! Now construct the Z_psi(r)
   call random_number_wrapper(rphase1)
   rphase1=rphase1*2*pi
   if (lfastCK) then
     do n=n1-nghost,n2+nghost
       do m=m1-nghost,m2+nghost
         psilm=0.
         psilm= RYlm_list(m,n,lmindex)*cos(rphase1)- &
           IYlm_list(m,n,lmindex)*sin(rphase1)
         psif(:,m,n) = psilm*Zpsi_list(:,lmindex,aindex+1)
       enddo
     enddo
   else
     call sp_bessely_l(anum,Legendrel,Balpha*x(l1))
     call sp_besselj_l(adenom,Legendrel,Balpha*x(l1))
     a_ell = -anum/adenom
!        write(*,*) 'dhruba:',anum,adenom,Legendrel,Bessel_alpha,x(l1)
     do l=l1-nghost,l2+nghost
       alphar=Balpha*x(l)
       call sp_besselj_l(jlm,Legendrel,alphar)
       call sp_bessely_l(ylm,Legendrel,alphar)
       Z_psi(l) = (a_ell*jlm+ylm)
     enddo
!-------
     do n=n1-nghost,n2+nghost
       do m=m1-nghost,m2+nghost
         psilm=0.
         call sp_harm_real(RYlm,Legendrel,emm,y(m),z(n))
         call sp_harm_imag(IYlm,Legendrel,emm,y(m),z(n))
         psilm= RYlm*cos(rphase1)-IYlm*sin(rphase1)
         psif(:,m,n) = Z_psi*psilm
       enddo
     enddo
   endif
! ----- Now calculate the force from the potential and add this to
! velocity
! get a random unit vector with three components ee_r, ee_theta, ee_phi
! psi at present is just Z_{ell}^m. We next do a sum over random coefficients
! get random psi.
!      write(*,*) 'mmin=',mmin
!! ----------now generate and add the force ------------
   call random_number_wrapper(rz)
   ee(3) = rz
   call random_number_wrapper(rphase2)
   rphase2 = pi*rphase2
   ee(1) = sqrt(1-rz*rz)*cos(rphase2)
   ee(2) = sqrt(1-rz*rz)*sin(rphase2)
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
       capitalH = capitalT + capitalS
       do j=1,3
         jf = iuu+j-1
         f(l1:l2,m,n,jf) = fnorm*capitalH(:,j)
       enddo
     enddo
   enddo
!
    endsubroutine forcing_cktest
!***********************************************************************
    subroutine forcing_GP(f)
!
!  Add Galloway-Proctor forcing function.
!
!  24-jul-06/axel: coded
!
      use Mpicomm
      use General
      use Sub
!
      real :: irufm
      real, dimension (1) :: fsum_tmp,fsum
      real, dimension (nx) :: ruf,rho
      real, dimension (nx,3) :: variable_rhs,forcing_rhs,force_all
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx) :: cosx,sinx
      real :: cost,sint,cosym,sinym
      integer :: j,jf
      real :: fact
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
              rho=exp(f(l1:l2,m,n,ilnrho))
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
          fsum_tmp(1)=irufm
          call mpireduce_sum(fsum_tmp,fsum,1)
          irufm=fsum(1)
          call mpibcast_real(irufm,1)
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
    subroutine forcing_TG(f)
!
!  Add Taylor-Green forcing function.
!
!   9-oct-04/axel: coded
!
      use Mpicomm
      use General
      use Sub
!
      real :: irufm
      real, dimension (1) :: fsum_tmp,fsum
      real, dimension (nx) :: ruf,rho
      real, dimension (nx,3) :: variable_rhs,forcing_rhs,force_all
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx), save :: sinx,cosx
      real, dimension (my), save :: siny,cosy
      real, dimension (mz), save :: cosz
      logical, dimension (3), save :: extent
      integer, save :: ifirst
      integer :: j,jf
      real :: fact
!
      if (ifirst==0) then
        if (lroot) print*,'forcing_TG: calculate sinx,cosx,siny,cosy,cosz'
        sinx=sin(k1_ff*x)
        cosx=cos(k1_ff*x)
        siny=sin(k1_ff*y)
        cosy=cos(k1_ff*y)
        cosz=cos(k1_ff*z)
        extent(1)=nx/=1
        extent(2)=ny/=1
        extent(3)=nz/=1
      endif
      ifirst=ifirst+1
!
      if (ip<=6) print*,'forcing_TG: dt, ifirst=',dt,ifirst
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
            if (extent(j)) then
              jf=j+ifff-1
              f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+forcing_rhs(:,j)
            endif
          enddo
          if (lout) then
            if (idiag_rufm/=0) then
              rho=exp(f(l1:l2,m,n,ilnrho))
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
          fsum_tmp(1)=irufm
          call mpireduce_sum(fsum_tmp,fsum,1)
          irufm=fsum(1)
          call mpibcast_real(irufm,1)
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
      use General
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      real :: irufm
      real, dimension (1) :: fsum_tmp,fsum
      real, dimension (nx) :: ruf,rho
      real, dimension (nx,3) :: variable_rhs,forcing_rhs,force_all,bb,fxb
      real, dimension (mx), save :: sinx,cosx
      real, dimension (my), save :: siny,cosy
      real, dimension (mz), save :: sinz,cosz
      integer, save :: ifirst
      integer :: j,jf
      real :: fact
!
!  at the first step, the sin and cos functions are calculated for all
!  x,y,z points and are then saved and used for all subsequent steps
!  and pencils
!
      if (ip<=6) print*,'forcing_ABC: ifirst=',ifirst
      if (ifirst==0) then
        if (lroot) print*,'forcing_ABC: calculate sinx,cosx,siny,cosy,sinz,cosz'
        sinx=sin(k1_ff*x); cosx=cos(k1_ff*x)
        siny=sin(k1_ff*y); cosy=cos(k1_ff*y)
        sinz=sin(k1_ff*z); cosz=cos(k1_ff*z)
      endif
      ifirst=ifirst+1
      if (ip<=6) print*,'forcing_ABC: dt, ifirst=',dt,ifirst
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
              rho=exp(f(l1:l2,m,n,ilnrho))
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
          fsum_tmp(1)=irufm
          call mpireduce_sum(fsum_tmp,fsum,1)
          irufm=fsum(1)
          call mpibcast_real(irufm,1)
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
    subroutine forcing_nocos(f)
!
!  Add no-cosine forcing function.
!
!  27-oct-04/axel: coded
!
      use Mpicomm
      use General
      use Sub
!
      real :: irufm
      real, dimension (1) :: fsum_tmp,fsum
      real, dimension (nx) :: ruf,rho
      real, dimension (nx,3) :: variable_rhs,forcing_rhs,force_all
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx), save :: sinx
      real, dimension (my), save :: siny
      real, dimension (mz), save :: sinz
      logical, dimension (3), save :: extent
      integer, save :: ifirst
      integer :: j,jf
      real :: fact
!
      if (ifirst==0) then
        if (lroot) print*,'forcing_nocos: calculate sinx,siny,sinz'
        sinx=sin(k1_ff*x)
        siny=sin(k1_ff*y)
        sinz=sin(k1_ff*z)
        extent(1)=nx/=1
        extent(2)=ny/=1
        extent(3)=nz/=1
      endif
      ifirst=ifirst+1
!
      if (ip<=6) print*,'forcing_hel: dt, ifirst=',dt,ifirst
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
            if (extent(j)) then
              jf=j+ifff-1
              f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+forcing_rhs(:,j)
            endif
          enddo
          if (lout) then
            if (idiag_rufm/=0) then
              rho=exp(f(l1:l2,m,n,ilnrho))
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
          fsum_tmp(1)=irufm
          call mpireduce_sum(fsum_tmp,fsum,1)
          irufm=fsum(1)
          call mpibcast_real(irufm,1)
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
      use Mpicomm
      use General
      use Sub
      use EquationOfState, only: cs0
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: force_ampl, force_tmp
!
      real, dimension (1) :: fsum_tmp,fsum
      real, dimension (3) :: fran
      real, dimension (nx) :: radius2,gaussian,ruf,rho
      real, dimension (nx,3) :: variable_rhs,force_all,delta
      logical, dimension (3), save :: extent
      integer :: j,jf
      real :: irufm,fact,width_ff21
!
!  check length of time step
!
      if (ip<=6) print*,'forcing_gaussianpot: dt=',dt
!
!  check whether there is any extent in each of the three directions
!
      extent(1)=nx/=1
      extent(2)=ny/=1
      extent(3)=nz/=1
!
!  generate random numbers
!
      if (t>tsforce) then
        if (lrandom_location) then
          call random_number_wrapper(fran)
          location=fran*Lxyz+xyz0
        else
          location=location_fixed
        endif
!
!  It turns out that in the presence of shear, and even for weak shear,
!  vorticitity is being produced. In order to check whether the shearing
!  periodic boundaries are to blame, we can cut the y extent of forcing
!  locations by half.
!
        if (lforce_cuty) location(2)=location(2)*.5
!
!  reset location(i) to x, y, or z
!
        if (.not.extent(1)) location(1)=x(1)
        if (.not.extent(2)) location(2)=y(1)
        if (.not.extent(3)) location(3)=z(1)
!
!  write location to file
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
!  Set forcing amplitude (same value for each location by default)
!
      if (iforce_tprofile=='nothing') then
        force_tmp=force_ampl
      elseif (iforce_tprofile=='sin^2') then
        force_tmp=force_ampl*sin(pi*(tsforce-t)/dtforce)**2
      else
        call  fatal_error('forcing_gaussianpot','iforce_tprofile not good')
      endif
!
!  Possibility of outputting data at each time step (for testing)
!
      if (lroot .and. lwrite_gausspot_to_file_always) then
        open(1,file=trim(datadir)//'/gaussian_pot_forcing_always.dat', &
            status='unknown',position='append')
        write(1,'(5f14.7)') t, location, force_tmp
        close (1)
      endif
!
!  Let explosion last dtforce_duration or, by default, until next explosion.
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
!  Also define width_ff21 = 1/width^2
!
        width_ff21=1./width_ff**2
        fact=2.*width_ff21*force_tmp*dt*sqrt(cs0*width_ff/max(dtforce+.5*dt,dt))
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
              if (.not.extent(j)) delta(:,j)=0.
            enddo
!
!  compute gaussian blob and set to forcing function
!
            radius2=delta(:,1)**2+delta(:,2)**2+delta(:,3)**2
            gaussian=fact*exp(-radius2*width_ff21)
            variable_rhs=f(l1:l2,m,n,iffx:iffz)
            do j=1,3
              if (extent(j)) then
                jf=j+ifff-1
                f(l1:l2,m,n,jf)=f(l1:l2,m,n,jf)+gaussian*delta(:,j)
              endif
            enddo
!
!  test
!
!--         if (icc/=0) f(l1:l2,m,n,icc)=f(l1:l2,m,n,icc)+gaussian
!
!  diagnostics
!
            if (lout) then
              if (idiag_rufm/=0) then
                rho=exp(f(l1:l2,m,n,ilnrho))
                call multsv_mn(rho/dt,spread(gaussian,2,3)*delta,force_all)
                call dot_mn(variable_rhs,force_all,ruf)
                irufm=irufm+sum(ruf)
              endif
            endif
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
          fsum_tmp(1)=irufm
          call mpireduce_sum(fsum_tmp,fsum,1)
          irufm=fsum(1)
          call mpibcast_real(irufm,1)
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
    subroutine calc_force_ampl(f,fx,fy,fz,coef,force_ampl)
!
!  calculates the coefficient for a forcing that satisfies
!  <rho*u*f> = constant.
!
!   7-sep-02/axel: coded
!
      use Sub
      use Mpicomm
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3) :: uu
      real, dimension (nx) :: rho,udotf
      real, dimension (1) :: fsum_tmp,fsum
      complex, dimension (mx) :: fx
      complex, dimension (my) :: fy
      complex, dimension (mz) :: fz
      complex, dimension (3) :: coef
      real :: rho_uu_ff,force_ampl
      integer :: j
!
      rho_uu_ff=0.
      do n=n1,n2
        do m=m1,m2
          rho=exp(f(l1:l2,m,n,ilnrho))
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
      fsum_tmp(1)=rho_uu_ff
      call mpireduce_sum(fsum_tmp,fsum,1)
      if (lroot) rho_uu_ff=fsum(1)/nwgrid
!      if (lroot) rho_uu_ff=rho_uu_ff/nwgrid
      call mpibcast_real(rho_uu_ff,1)
!
!  scale forcing function
!  but do this only when rho_uu_ff>0.; never allow it to change sign
!
!
!print*,fname(idiag_urms)
!
        if (headt) print*,'calc_force_ampl: divide forcing function by rho_uu_ff=',rho_uu_ff
        !      force_ampl=work_ff/(.1+max(0.,rho_uu_ff))
        force_ampl=work_ff/rho_uu_ff
        if (force_ampl > max_force) force_ampl=max_force
        if (force_ampl < -max_force) force_ampl=-max_force
!
    endsubroutine calc_force_ampl
!***********************************************************************
    subroutine forcing_hel_noshear(f)
!
!  add helical forcing function, using a set of precomputed wavevectors
!
!  10-apr-00/axel: coded
!
      use Mpicomm
      use General
      use Sub
      use EquationOfState, only: cs0
!
      real :: phase,ffnorm
      real, save :: kav
      real, dimension (2) :: fran
      real, dimension (nx) :: radius,tmpx
!      real, dimension (mz) :: tmpz
      real, dimension (mx,my,mz,mfarray) :: f
      complex, dimension (mx) :: fx
      complex, dimension (my) :: fy
      complex, dimension (mz) :: fz
      complex, dimension (3) :: coef
      integer, parameter :: mk=3000
      real, dimension(mk), save :: kkx,kky,kkz
      integer, save :: ifirst,nk
      integer :: ik,j,jf,kx,ky,kz,kex,key,kez,kkex,kkey,kkez
      real :: k2,k,ex,ey,ez,kde,sig,fact
      real, dimension(3) :: e1,e2,ee,kk
      real :: norm,phi
!
      if (ifirst==0) then
        if (lroot) print*,'force_hel_noshear: opening k.dat'
        open(9,file='k.dat',status='old')
        read(9,*) nk,kav
        if (lroot) print*,'force_hel_noshear: average k=',kav
        if (nk>mk) then
          if (lroot) print*,'force_hel_noshear: dimension mk in forcing_hel is insufficient'
          print*,'nk=',nk,'mk=',mk
          call mpifinalize
        endif
        read(9,*) (kkx(ik),ik=1,nk)
        read(9,*) (kky(ik),ik=1,nk)
        read(9,*) (kkz(ik),ik=1,nk)
        close(9)
      endif
      ifirst=ifirst+1
!
!  generate random coefficients -1 < fran < 1
!  ff=force*Re(exp(i(kx+phase)))
!  |k_i| < akmax
!
      call random_number_wrapper(fran)
      phase=pi*(2*fran(1)-1.)
      ik=nk*.9999*fran(2)+1
      if (ip<=6) print*,'force_hel_noshear: ik,phase,kk=',ik,phase,kkx(ik),kky(ik),kkz(ik),dt,ifirst
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
        call random_number_wrapper(phi); phi = phi*2*pi
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
      if (ip<=12) print*,'force_hel_noshear: k,kde,ffnorm,kav,dt,cs0=',k,kde,ffnorm,kav,dt,cs0
      if (ip<=12) print*,'force_hel_noshear: k*sqrt(k2-kde**2)=',k*sqrt(k2-kde**2)
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
!-      if (lroot .and. ifirst==1) print*,'forcing_hel_noshear: include z-profile'
!-      tmpz=(z/height_ff)**2
!-      fz=fz*exp(-tmpz**5/max(1.-tmpz,1e-5))
!-    endif
!
!  need to discuss with axel
!
!  possibly multiply forcing by sgn(z) and radial profile
!
!      if (r_ff/=0.) then
!        if (lroot .and. ifirst==1) &
!             print*,'forcing_hel_noshear: applying sgn(z)*xi(r) profile'
!        !
!        ! only z-dependent part can be done here; radial stuff needs to go
!        ! into the loop
!        !
!        tmpz = tanh(z/width_ff)
!        fz = fz*tmpz
!      endif
!
      if (ip<=5) print*,'force_hel_noshear: fx=',fx
      if (ip<=5) print*,'force_hel_noshear: fy=',fy
      if (ip<=5) print*,'force_hel_noshear: fz=',fz
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
!  add differential rotation
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
!  add blobs in entropy every dforce time units
!
!  28-jul-02/axel: coded
!
      !use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, save :: tforce=0.
      integer, save :: ifirst=0
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
      if (ifirst==0) then
        call read_snaptime(trim(file),tforce,nforce,dforce,t)
        ifirst=1
      endif
!
!  Check whether we want to do forcing at this time.
!
      call update_snaptime(file,tforce,nforce,dforce,t,lforce,ch)
      if (lforce) then
        call blob(force,f,iss,radius_ff,0.,0.,.5)
      endif
!
    endsubroutine forcing_blobs
!***********************************************************************
    subroutine forcing_hel_smooth(f)
!
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,3) :: force1,force2,force_vec
      real, dimension (nx) :: ruf,rho
      real, dimension (nx,3) :: variable_rhs,forcing_rhs,force_all
      real :: phase1,phase2,p_weight
      real :: kx01,ky1,kz1,kx02,ky2,kz2
      real :: mulforce_vec,irufm
      real, dimension (1) :: fsum_tmp,fsum
      integer, parameter :: mk=3000
      real, dimension(mk), save :: kkx,kky,kkz
      integer, save :: ifirst,nk
      integer :: ik1,ik2,ik
      real, save :: kav
!
      if (ifirst==0) then
         if (lroot) print*,'forcing_hel_smooth: opening k.dat'
         open(9,file='k.dat',status='old')
         read(9,*) nk,kav
         if (lroot) print*,'forcing_hel_smooth: average k=',kav
         if (nk>mk) then
            if (lroot) print*, &
                 'forcing_hel_smooth: dimension mk in forcing_hel_smooth is insufficient'
            print*,'nk=',nk,'mk=',mk
            call mpifinalize
         endif
         read(9,*) (kkx(ik),ik=1,nk)
         read(9,*) (kky(ik),ik=1,nk)
         read(9,*) (kkz(ik),ik=1,nk)
         close(9)
      endif
      ifirst=ifirst+1
!
!  Re-calculate forcing wave numbers if necessary
!
      !tsforce is set to -10 in cdata.f90. It should also be saved in a file
      !so that it will be read again on restarts.
      if (t > tsforce) then
         if (tsforce < 0) then
            call random_number_wrapper(fran1)
         else
            fran1=fran2
         endif
         call random_number_wrapper(fran2)
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
!
!  Calculate forcing function
!
      call hel_vec(f,kx01,ky1,kz1,phase1,kav,ifirst,force1)
      call hel_vec(f,kx02,ky2,kz2,phase2,kav,ifirst,force2)
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
              rho=exp(f(l1:l2,m,n,ilnrho))
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
        fsum_tmp(1)=irufm
        call mpireduce_sum(fsum_tmp,fsum,1)
        irufm=fsum(1)
        call mpibcast_real(irufm,1)
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
    subroutine hel_vec(f,kx0,ky,kz,phase,kav,ifirst,force1)
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
      use Mpicomm
      use General
      use Sub
      use EquationOfState, only: cs0
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: phase,ffnorm
      real :: kav
      real, dimension (nx) :: radius,tmpx
      real, dimension (mz) :: tmpz
      real, dimension (mx,my,mz,3) :: force1
      complex, dimension (mx) :: fx
      complex, dimension (my) :: fy
      complex, dimension (mz) :: fz
      complex, dimension (3) :: coef
      integer :: j,jf
      integer :: ifirst
      real :: kx0,kx,ky,kz,k2,k,force_ampl
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
        call random_number_wrapper(phi); phi = phi*2*pi
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
      if (ip<=9) print*,'hel_vec: k,kde,ffnorm,kav,dt,cs0=',k,kde, &
                                                        ffnorm,kav,dt,cs0
      if (ip<=9) print*,'hel_vec: k*sqrt(k2-kde**2)=',k*sqrt(k2-kde**2)
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
!-      if (lroot .and. ifirst==1) print*,'hel_vec: include z-profile'
!-      tmpz=(z/height_ff)**2
!-      fz=fz*exp(-tmpz**5/max(1.-tmpz,1e-5))
!-    endif
!
!  possibly multiply forcing by sgn(z) and radial profile
!
      if (r_ff/=0.) then
        if (lroot .and. ifirst==1) &
             print*,'hel_vec: applying sgn(z)*xi(r) profile'
        !
        ! only z-dependent part can be done here; radial stuff needs to go
        ! into the loop
        !
        tmpz = tanh(z/width_ff)
        fz = fz*tmpz
      endif
!
      if (ip<=5) print*,'hel_vec: fx=',fx
      if (ip<=5) print*,'hel_vec: fy=',fy
      if (ip<=5) print*,'hel_vec: fz=',fz
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
        if (lwork_ff) call calc_force_ampl(f,fx,fy,fz,coef,force_ampl)
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
    subroutine calc_fluxring_cylindrical(force)
!
!   4-aug-11/dhruba+axel: adapted from fluxring_cylindrical
!
      real, dimension (nx,3), intent(out) :: force
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
      force = fcont_ampl*force
!
    endsubroutine calc_fluxring_cylindrical
!***********************************************************************
    subroutine calc_lforcing_cont_pars(f)
!
!  precalculate parameters that are new at each timestep,
!  but the same for all pencils
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(in) :: f
!
!  for the AKA effect, calculate auxiliary functions phi1_ff and phi2_ff
!
      if (iforcing_cont=='AKA') then
        phi1_ff=cos(kf_fcont*y+omega_fcont*t)
        phi2_ff=cos(kf_fcont*x-omega_fcont*t)
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_lforcing_cont_pars
!***********************************************************************
    subroutine pencil_criteria_forcing()
!
!  All pencils that the Forcing module depends on are specified here.
!
!  24-mar-08/axel: adapted from density.f90
!
      if (lforcing_cont) lpenc_requested(i_fcont)=.true.
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
      use Debug_IO, only: output_pencil
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(inout) :: f,p
!
!  calculate forcing
!
      if (lpencil(i_fcont)) then
!
        if (headtt .and. lroot) print*,'forcing: add continuous forcing'
        call forcing_cont(p%fcont)
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_pencils_forcing
!***********************************************************************
    subroutine forcing_cont(force)
!
! 9-apr-10/MR: added RobertsFlow_exact forcing, compensates \nu\nabla^2 u and u.grad u for Roberts geometry
! 4-nov-11/MR:                                  now also compensates Coriolis force
!
! Note: It is not enough to set lforcing_cont = T in input parameters of forcing
! one must also set  lforcing_cont_uu = T in hydro for the continious is time
! forcing to be added to velocity.
!
      use Sub, only: quintic_step, quintic_der_step
      use Mpicomm, only: stop_it
      use Gravity, only: gravz
      use SharedVariables, only: get_shared_variable
      use Viscosity, only: getnu
!
      real, dimension (nx,3), intent(out) :: force
!
      real            :: fact, fact2, fpara, dfpara, sqrt21k1, kf, kx, ky, nu, arg
      real, pointer   :: gravx
      integer         :: ierr
      integer :: i2d1=1,i2d2=2,i2d3=3
!
        select case (iforcing_cont)
        case('ABC')
          fact=ampl_ff/sqrt(3.)
          fact2=1.-eps_fcont
          force(:,1)=fact*(sinz(n    )+fact2*cosy(m)    )
          force(:,2)=fact*(sinx(l1:l2)+fact2*cosz(n)    )
          force(:,3)=fact*(siny(m    )+fact2*cosx(l1:l2))
        case ('AKA')
          fact=sqrt(2.)*ampl_ff
          force(:,1)=fact*phi1_ff(m    )
          force(:,2)=fact*phi2_ff(l1:l2)
          force(:,3)=fact*(phi1_ff(m)+phi2_ff(l1:l2))
        case ('grav_z')
          force(:,1)=0
          force(:,2)=0
          force(:,3)=gravz*ampl_ff*cos(omega_ff*t)
        case ('grav_xz')
          call get_shared_variable('gravx', gravx, ierr)
          if (ierr/=0) call stop_it("forcing: "//&
            "there was a problem when getting gravx")
          force(:,1)=gravx*ampl_ff*cos(omega_ff*t)
          force(:,2)=0
          force(:,3)=gravz*ampl_ff*cos(omega_ff*t)
        case('KolmogorovFlow-x')
          fact=ampl_ff
          force(:,1)=0
          force(:,2)=fact*cosx(l1:l2)
          force(:,3)=0
        case('KolmogorovFlow-z')
          fact=ampl_ff
          force(:,1)=fact*cosz(n)
          force(:,2)=0.
          force(:,3)=0.
        case ('nocos')
          fact=ampl_ff
          force(:,1)=fact*sinz(n)
          force(:,2)=fact*sinx(l1:l2)
          force(:,3)=fact*siny(m)
        case('RobertsFlow')
          fact=ampl_ff
          force(:,1)=-fact*cosx(l1:l2)*siny(m)
          force(:,2)=+fact*sinx(l1:l2)*cosy(m)
          force(:,3)=+relhel*fact*cosx(l1:l2)*cosy(m)*sqrt2
        case('RobertsFlow2d')
          fact=ampl_ff
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
          kx=k1_ff; ky=k1_ff
          kf=sqrt(kx*kx+ky*ky)
          call getnu(nu_input=nu)
          fact=ampl_ff*kf*kf*nu
          fact2=ampl_ff*ampl_ff*kx*ky
!
!!print*, 'forcing: kx, ky, kf, fact, fact2=', kx, ky, kf, fact, fact2
          force(:,1)=-fact*ky*cosx(l1:l2)*siny(m) - fact2*ky*sinx(l1:l2)*cosx(l1:l2)
          force(:,2)=+fact*kx*sinx(l1:l2)*cosy(m) - fact2*kx*siny(m)*cosy(m)
          force(:,3)=+fact*kf*cosx(l1:l2)*cosy(m)

          if ( Omega/=0. .and. theta==0. ) then              ! Obs, only implemented for rotation axis in z direction.
            fact = 2.*ampl_ff*Omega
            force(:,1)= force(:,1)-fact*kx*sinx(l1:l2)*cosy(m)
            force(:,2)= force(:,2)-fact*ky*cosx(l1:l2)*siny(m)
          endif
!
        case ('RobertsFlow-zdep')
          if (headtt) print*,'z-dependent Roberts flow; eps_fcont=',eps_fcont
          fpara=quintic_step(z(n),-1.+eps_fcont,eps_fcont) &
               -quintic_step(z(n),+1.-eps_fcont,eps_fcont)
          dfpara=quintic_der_step(z(n),-1.+eps_fcont,eps_fcont)&
                -quintic_der_step(z(n),+1.-eps_fcont,eps_fcont)
!
!  abbreviations
!
          sqrt21k1=1./(sqrt2*k1_ff)
!
!  amplitude factor missing in upper lines
!
          force(:,1)=-ampl_ff*cosx(l1:l2)*siny(m) &
                -dfpara*ampl_ff*sinx(l1:l2)*cosy(m)*sqrt21k1
          force(:,2)=+ampl_ff*sinx(l1:l2)*cosy(m) &
                -dfpara*ampl_ff*cosx(l1:l2)*siny(m)*sqrt21k1
          force(:,3)=+fpara*ampl_ff*cosx(l1:l2)*cosy(m)*sqrt2
!
!  f=(sinx,0,0)
!
        case ('sinx')
          if (lgentle.and.t<tgentle) then
            fact=.5*ampl_ff*(1.-cos(pi*t/tgentle))
          else
            fact=ampl_ff
          endif
          force(:,1)=fact*sinx(l1:l2)
          force(:,2)=0.
          force(:,3)=0.
!
!  f=(0,0,cosx)
!
        case ('(0,0,cosx)')
          force(:,1)=0.
          force(:,2)=0.
          force(:,3)=ampl_ff*cosx(l1:l2)
!
        case ('TG')
          fact=2.*ampl_ff
          force(:,1)=+fact*sinx(l1:l2)*cosy(m)*cosz(n)
          force(:,2)=-fact*cosx(l1:l2)*siny(m)*cosz(n)
          force(:,3)=0.
!
! Continuous emf required in Induction equation for the Mag Buoy Inst
!
        case('mbi_emf')
          fact=2.*ampl_bb*eta_bb/width_bb**2
          arg=(z(n)-z_bb)/width_bb
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
!  fluxring_cylindrical
!
        case('fluxring_cylindrical')
          call calc_fluxring_cylindrical(force)
!
        case default
          call stop_it('forcing: no continuous iforcing_cont specified')
        endselect
!
    endsubroutine forcing_cont
!***********************************************************************
    subroutine forcing_continuous(df,p)
!
!  add a continuous forcing term (used to be in hydro.f90)
!
!  17-sep-06/axel: coded
!
      use Diagnostics
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: fxb
      real, dimension (nx) :: uf
      type (pencil_case) :: p
!
!  diagnostics
!
      if (ldiagnos) then
        if (idiag_rufm/=0) then
          call dot_mn(p%uu,p%fcont,uf)
          call sum_mn_name(p%rho*uf,idiag_rufm)

        endif
        if (lmagnetic) then
          if (idiag_fxbxm/=0.or.idiag_fxbym/=0.or.idiag_fxbzm/=0) then
            call cross(p%fcont,p%bb,fxb)
            call sum_mn_name(fxb(:,1),idiag_fxbxm)
            call sum_mn_name(fxb(:,2),idiag_fxbym)
            call sum_mn_name(fxb(:,3),idiag_fxbzm)
          endif
        endif
      endif
!
      call keep_compiler_quiet(df)
!
    endsubroutine forcing_continuous
!***********************************************************************
    subroutine read_forcing_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      call keep_compiler_quiet(present(iostat))
!
    endsubroutine read_forcing_init_pars
!***********************************************************************
    subroutine write_forcing_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_forcing_init_pars
!***********************************************************************
    subroutine read_forcing_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=forcing_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=forcing_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_forcing_run_pars
!***********************************************************************
    subroutine write_forcing_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=forcing_run_pars)
!
    endsubroutine write_forcing_run_pars
!***********************************************************************
    subroutine input_persistent_forcing(id,done)
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
      endselect
!
      if (lroot) print *, 'input_persistent_forcing: ', location, tsforce
!
    endsubroutine input_persistent_forcing
!***********************************************************************
    logical function output_persistent_forcing()
!
!  Writes out the time of the next SNI
!  This is used, for example, for forcing functions with temporal
!  memory, such as in the paper by Mee & Brandenburg (2006, MNRAS)
!
!  21-dec-05/tony: coded
!  13-Dec-2011/Bourdin.KIS: reworked
!
      use IO, only: write_persist
!
      if (lroot .and. (tsforce>=0.)) &
          print *, 'output_persistent_forcing: ', location, tsforce
!
!  write details
!
      output_persistent_forcing = .false.
!
      if (write_persist ('FORCING_LOCATION', id_record_FORCING_LOCATION, location)) &
          output_persistent_forcing = .true.
      if (write_persist ('FORCING_TSFORCE', id_record_FORCING_TSFORCE, tsforce)) &
          output_persistent_forcing = .true.
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
        idiag_rufm=0; idiag_ufm=0; idiag_ofm=0; idiag_ffm=0
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
        call parse_name(iname,cname(iname),cform(iname),'ruxfxm',idiag_ruxfxm)
        call parse_name(iname,cname(iname),cform(iname),'ruxfym',idiag_ruxfym)
        call parse_name(iname,cname(iname),cform(iname),'ruyfxm',idiag_ruyfxm)
        call parse_name(iname,cname(iname),cform(iname),'ruyfym',idiag_ruyfym)
        call parse_name(iname,cname(iname),cform(iname),'ruzfzm',idiag_ruzfzm)
        call parse_name(iname,cname(iname),cform(iname),'ufm',idiag_ufm)
        call parse_name(iname,cname(iname),cform(iname),'ofm',idiag_ofm)
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
      if (iforce=='chandra-kendall') then
        print*,'Deallocating arrays relevent to Chanrasekhar-Kendall forcing ..'
        deallocate(psif,cklist)
        if (lfastCK)  deallocate(Zpsi_list,RYlm_list)
        print*,'Done.'
      endif
!
    endsubroutine forcing_clean_up
!*******************************************************************
endmodule Forcing
