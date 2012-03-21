! $Id$
!
!  This module calculates a number of outputs and removes a mean
!  (phi-averaged) emf from the simulations with net vertical fields
!  on global accretion disks
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
!
!-------------------------------------------------------------------
!
! HOW TO USE THIS FILE
! --------------------
!
! The rest of this file may be used as a template for your own
! special module.  Lines which are double commented are intended
! as examples of code.  Simply fill out the prototypes for the
! features you want to use.
!
! Save the file with a meaningful name, eg. geo_kws.f90 and place
! it in the $PENCIL_HOME/src/special directory.  This path has
! been created to allow users ot optionally check their contributions
! in to the Pencil-Code CVS repository.  This may be useful if you
! are working on/using the additional physics with somebodyelse or
! may require some assistance from one of the main Pencil-Code team.
!
! To use your additional physics code edit the Makefile.local in
! the src directory under the run directory in which you wish to
! use your additional physics.  Add a line with all the module
! selections to say something like:
!
!    SPECIAL=special/rtime_phiavg
!
! Where nstar it replaced by the filename of your new module
! upto and not including the .f90
!
!--------------------------------------------------------------------

module Special

  use Cdata
  use Cparam
  use EquationOfState
  use Messages
  use Sub, only: keep_compiler_quiet

  implicit none

  include '../special.h'

! input parameters
  real, dimension(nrcylrun)  :: rho_tmp
  real, dimension(nrcylrun,3) :: u_tmp,b_tmp
  integer, dimension(nrcylrun) :: k_tmp
!
  real, dimension (nrcylrun,3) :: bavg_coarse,uavg_coarse
  real, dimension (nrcylrun) :: rcyl_coarse,rhoavg_coarse
  real, dimension (nx,3) :: bavg,uavg
  real, dimension (nx) :: rhoavg
  real :: drc,r1,r2,B_ext=0.,rt_int=0.,rt_ext=impossible
  logical :: llarge_scale_Bz=.false.
  integer :: dummy=0,nd
  logical :: laverage_smooth=.true.,lmedian_smooth=.false.
  logical :: lcalc_density_pars=.true.
!
!  start parameters
!
  namelist /special_init_pars/ dummy
!
!   run parameters
!
  namelist /special_run_pars/ &
       B_ext,llarge_scale_Bz,rt_int,rt_ext,laverage_smooth,lmedian_smooth,&
       lcalc_density_pars
!
! Keep some over used pencils
!
!!
!! Declare any index variables necessary for main or
!!
!!   integer :: iSPECIAL_VARIABLE_INDEX=0
!
!! other variables (needs to be consistent with reset list below)
!
  !integer :: idiag_dtcrad=0
  !integer :: idiag_dtchi=0
!
! hydro diagnostics
!
  integer :: idiag_urm=0  ,idiag_upm=0  ,idiag_uzzm=0 
  integer :: idiag_ur2m=0 ,idiag_up2m=0 ,idiag_uzz2m=0
  integer :: idiag_urupm=0,idiag_uzupm=0,idiag_uruzm=0
!
! magnetic diagnostics
!
  integer :: idiag_brm=0  ,idiag_bpm=0  ,idiag_bzm=0
  integer :: idiag_br2m=0 ,idiag_bp2m=0 ,idiag_bzz2m=0
  integer :: idiag_brbpm=0,idiag_bzbpm=0,idiag_brbzm=0
!
! 1D average diagnostics
! 
  integer :: idiag_brbpmr=0,idiag_urupmr=0,idiag_mdotmr=0
!
  contains

!***********************************************************************
    subroutine register_special()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  6-oct-03/tony: coded
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f,lstarting)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  06-oct-03/tony: coded
!
      use EquationOfState
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
      real :: rloop_int,rloop_ext,tmp,drc1
      integer :: ir
!
!  Initialize any module variables which are parameter dependant
!
      if (r_int.ne.0) rt_int=r_int
      if (r_ext.ne.impossible) rt_ext=r_ext
!
      if (lroot) print*,'rt_int,rt_ext',rt_int,rt_ext
!
      tmp = (rt_ext - rt_int)/nrcylrun
!
      drc=tmp
      drc1=1./drc
!
      do ir=1,nrcylrun
        rloop_int = rt_int + (ir-1)*drc
        rloop_ext = rt_int +  ir   *drc
        rcyl_coarse(ir)= 0.5*(rloop_int + rloop_ext)
      enddo
!  dimensions
      nd=3
      if (nzgrid==1) nd=2
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise special condition; called from start.f90
!  06-oct-2003/tony: coded
!
      use EquationOfState
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_interdep_special(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified here.
!
!  18-07-06/tony: coded                                                         
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_special
!***********************************************************************
    subroutine pencil_criteria_special()
!
!  All pencils that this special module depends on are specified here.
!
!  18-07-06/tony: coded
!
      if (lcalc_density_pars) lpenc_requested(i_rho)=.true.
      if (lhydro)             lpenc_requested(i_uu)=.true.
      if (lmagnetic)          lpenc_requested(i_bb)=.true.
      if (lmagnetic.or.lhydro) then
        lpenc_diagnos(i_phix)=.true.
        lpenc_diagnos(i_phiy)=.true.
        lpenc_diagnos(i_pomx)=.true.
        lpenc_diagnos(i_pomy)=.true.
      endif
!      
      if (llarge_scale_Bz) then 
        lpenc_requested(i_uxb)=.true.
        lpenc_requested(i_phix)=.true.
        lpenc_requested(i_phiy)=.true.
        lpenc_requested(i_pomx)=.true.
        lpenc_requested(i_pomy)=.true.
        lpenc_requested(i_rcyl_mn1)=.true.
      endif
        
!      if (lmagnetic) lpenc_requested(i_bavg)=.true.
!      if (lhydro)    lpenc_requested(i_uavg)=.true.
!      if (ldensity)  lpenc_requested(i_rhoavg)=.true.

    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine dspecial_dt(f,df,p)
!
!  calculate right hand side of ONE OR MORE extra coupled PDEs
!  along the 'current' Pencil, i.e. f(l1:l2,m,n) where
!  m,n are global variables looped over in equ.f90
!
!  Due to the multi-step Runge Kutta timestepping used one MUST always
!  add to the present contents of the df array.  NEVER reset it to zero.
!
!  several precalculated Pencils of information are passed if for
!  efficiency.
!
!   06-oct-03/tony: coded
!
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p

!
      intent(in) :: f,p
      intent(inout) :: df
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dSPECIAL_dt'
!!      if (headtt) call identify_bcs('ss',iss)
!
!!
!! SAMPLE DIAGNOSTIC IMPLEMENTATION
!!
      !if (ldiagnos) then
      !  if (idiag_dtcrad/=0) &
      !    call max_mn_name(sqrt(advec_crad2)/cdt,idiag_dtcrad,l_dt=.true.)
      !  if (idiag_dtchi/=0) &
      !    call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
      !endif

!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine get_slices_special(f,slices)
!
!  Write slices for animation of special variables.
!
!  26-jun-06/tony: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_special
!***********************************************************************
    subroutine read_special_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=special_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=special_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=special_init_pars)

    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=special_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=special_run_pars,ERR=99)
      endif

99    return
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=special_run_pars)

    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  reads and registers print parameters relevant to special
!
!   06-oct-03/tony: coded
!
      use Diagnostics
!
!  define diagnostics variable
!
      integer :: iname,inamer
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
        !hydro
        idiag_urm=0  ;idiag_upm=0  ;idiag_uzzm=0 
        idiag_ur2m=0 ;idiag_up2m=0 ;idiag_uzz2m=0
        idiag_urupm=0;idiag_uzupm=0;idiag_uruzm=0
        !magnetic
        idiag_brm=0  ;idiag_bpm=0  ;idiag_bzm=0
        idiag_br2m=0 ;idiag_bp2m=0 ;idiag_bzz2m=0
        idiag_brbpm=0;idiag_bzbpm=0;idiag_brbzm=0
        !1d
        idiag_brbpmr=0;idiag_urupmr=0
      endif
!
      do iname=1,nname
        !hydro
        call parse_name(iname,cname(iname),cform(iname),'ur2m',idiag_ur2m)
        call parse_name(iname,cname(iname),cform(iname),'up2m',idiag_up2m)
        call parse_name(iname,cname(iname),cform(iname),'uzz2m',idiag_uzz2m)
        call parse_name(iname,cname(iname),cform(iname),'urupm',idiag_urupm)
        call parse_name(iname,cname(iname),cform(iname),'urm',idiag_urm)
        call parse_name(iname,cname(iname),cform(iname),'upm',idiag_upm)
        call parse_name(iname,cname(iname),cform(iname),'uzzm',idiag_uzzm)
        call parse_name(iname,cname(iname),cform(iname),'uzupm',idiag_uzupm)
        call parse_name(iname,cname(iname),cform(iname),'uruzm',idiag_uruzm)
        !magnetic
        call parse_name(iname,cname(iname),cform(iname),'br2m',idiag_br2m)
        call parse_name(iname,cname(iname),cform(iname),'bp2m',idiag_bp2m)
        call parse_name(iname,cname(iname),cform(iname),'bzz2m',idiag_bzz2m)
        call parse_name(iname,cname(iname),cform(iname),'brbpm',idiag_brbpm)
        call parse_name(iname,cname(iname),cform(iname),'brm',idiag_brm)
        call parse_name(iname,cname(iname),cform(iname),'bpm',idiag_bpm)
        call parse_name(iname,cname(iname),cform(iname),'bzm',idiag_bzm)
        call parse_name(iname,cname(iname),cform(iname),'bzbpm',idiag_bzbpm)
        call parse_name(iname,cname(iname),cform(iname),'brbzm',idiag_brbzm)
      enddo
!
      do inamer=1,nnamer
        call parse_name(inamer,cnamer(inamer),cform(inamer),'urupmr',idiag_urupmr)
        call parse_name(inamer,cnamer(inamer),cform(inamer),'brbpmr',idiag_brbpmr)
        call parse_name(inamer,cnamer(inamer),cform(inamer),'mdotmr',idiag_mdotmr)
      enddo
!
!  write column where which special variable is stored
!
      if (lwr) then
        !hydro
        write(3,*) 'i_urm=',idiag_urm
        write(3,*) 'i_upm=',idiag_upm
        write(3,*) 'i_uzzm=',idiag_uzzm
        write(3,*) 'i_ur2m=',idiag_ur2m
        write(3,*) 'i_up2m=',idiag_up2m
        write(3,*) 'i_uzz2m=',idiag_uzz2m
        write(3,*) 'i_urupm=',idiag_urupm
        !magnetic
        write(3,*) 'i_brm=',idiag_brm
        write(3,*) 'i_bpm=',idiag_bpm
        write(3,*) 'i_bzm=',idiag_bzm
        write(3,*) 'i_br2m=',idiag_br2m
        write(3,*) 'i_bp2m=',idiag_bp2m
        write(3,*) 'i_bzz2m=',idiag_bzz2m
        write(3,*) 'i_brbpm=',idiag_brbpm
        !1d
        write(3,*) 'i_brbpmr=',idiag_brbpmr
        write(3,*) 'i_urupmr=',idiag_urupmr
        write(3,*) 'i_mdotmr=',idiag_mdotmr
      endif
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!   06-oct-03/tony: coded
!   called on main loop - the average must be calculated BEFORE
!
      use Mpicomm
      use General, only: spline
!
      real, dimension(mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      integer :: i,j
      logical :: err
!
! expand it onto the pencil with spline interpolation
!
      if (lcalc_density_pars) & 
           call spline(rcyl_coarse,rhoavg_coarse,p%rcyl_mn,rhoavg,nrcylrun,nx,err)
      do j=1,nd
        if (lhydro) &
             call spline(rcyl_coarse,uavg_coarse(:,j),p%rcyl_mn,uavg(:,j),nrcylrun,nx,err)
        if (lmagnetic) &
             call spline(rcyl_coarse,bavg_coarse(:,j),p%rcyl_mn,bavg(:,j),nrcylrun,nx,err)
      enddo
!
! fill in with pencil values the parts of the array that are away from the interpolation
!
      do i=1,nx
        if ((p%rcyl_mn(i).lt.rcyl_coarse(1)).or.(p%rcyl_mn(i).gt.rcyl_coarse(nrcylrun))) then
          if (lcalc_density_pars) rhoavg(i) = p%rho(i)
          if (lhydro) then
            uavg(i,1)=p%uu(i,1)*p%pomx(i)+p%uu(i,2)*p%pomy(i)
            uavg(i,2)=p%uu(i,1)*p%phix(i)+p%uu(i,2)*p%phiy(i)
            uavg(i,3)=p%uu(i,3)
          endif
          if (lmagnetic) then
            bavg(i,1)=p%bb(i,1)*p%pomx(i)+p%bb(i,2)*p%pomy(i)
            bavg(i,2)=p%bb(i,1)*p%phix(i)+p%bb(i,2)*p%phiy(i)
            bavg(i,3)=p%bb(i,3)
          endif
        endif
      enddo
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine calc_lspecial_pars(f)
!
!  dummy routine
!
!  15-jan-08/axel: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_lspecial_pars
!***********************************************************************
    subroutine special_calc_density(f,df,p)
!
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_density
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!
!   16-jul-06/wlyra: coded
!
      use Diagnostics
!
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
      real, dimension (nx) :: ur,up,uz,urad,uphi
      real :: fac
!
      urad=p%uu(:,1)*p%pomx+p%uu(:,2)*p%pomy
      uphi=p%uu(:,1)*p%phix+p%uu(:,2)*p%phiy
      ur=urad      - uavg(:,1)
      up=uphi      - uavg(:,2)
      uz=p%uu(:,3) - uavg(:,3)
      if (nzgrid/=1) then 
        fac=Lxyz(3)
      else
        fac=1.
      endif

      if (ldiagnos) then
        if (idiag_urm/=0)    call sum_lim_mn_name(ur,idiag_urm,p)
        if (idiag_upm/=0)    call sum_lim_mn_name(up,idiag_upm,p)
        if (idiag_uzzm/=0)   call sum_lim_mn_name(uz,idiag_uzzm,p)
        if (idiag_ur2m/=0)   call sum_lim_mn_name(p%rho*ur**2,idiag_ur2m,p)
        if (idiag_up2m/=0)   call sum_lim_mn_name(p%rho*up**2,idiag_up2m,p)
        if (idiag_uzz2m/=0)  call sum_lim_mn_name(p%rho*uz**2,idiag_uzz2m,p)
        if (idiag_urupm/=0)  call sum_lim_mn_name(p%rho*ur*up,idiag_urupm,p)
        if (idiag_uzupm/=0)  call sum_lim_mn_name(p%rho*uz*up,idiag_uzupm,p)
        if (idiag_uruzm/=0)  call sum_lim_mn_name(p%rho*ur*uz,idiag_uruzm,p)
      endif
!
      if (l1davgfirst) then 
        if (idiag_urupmr/=0) &
             call phizsum_mn_name_r(p%rho*ur*up,idiag_urupmr)
        if (idiag_mdotmr/=0) &
             call phizsum_mn_name_r(p%rho*urad*2*pi*p%rcyl_mn*fac,idiag_mdotmr)
      endif
!
      call keep_compiler_quiet(f,df)
!
    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine special_calc_magnetic(f,df,p)
!
!   06-oct-03/tony: coded
!
     use Diagnostics
     use Mpicomm

     real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
     real, dimension (mx,my,mz,mvar), intent(inout) :: df
     type (pencil_case), intent(in) :: p
     real, dimension (nx) :: br,bp,bz
     real, dimension(nx,3) :: puxb
     integer :: i
!
! Remove mean electromotive force from induction equation.
! Activated only when large Bz fields and are present 
! keplerian advection
!
     if (llarge_scale_Bz) then
       if ((lmedian_smooth).and.(laverage_smooth)) then
         call stop_it("can't have both median and average at the same time")
       endif
!smooth a rapidly varying function - median box size=5
       if (lmedian_smooth) then
         do i=3,nx-2 
           puxb(i,:)=1./35*( -3*p%uxb(i-2,:)& 
                +12*p%uxb(i-1,:)&
                +17*p%uxb(i  ,:)&
                +12*p%uxb(i+1,:)&
                -3*p%uxb(i+2,:))
         enddo
         puxb(1,:)=p%uxb(1,:) ;  puxb(nx,:)=p%uxb(nx,:)    
         puxb(2,:)=p%uxb(2,:) ;  puxb(nx-1,:)=p%uxb(nx-1,:)
       elseif (laverage_smooth) then
         !use averages
         puxb(:,1)=uavg(:,2)*bavg(:,3)*x(l1:l2)*p%rcyl_mn1
         puxb(:,2)=uavg(:,2)*bavg(:,3)*y(m)*p%rcyl_mn1
       else
         call stop_it("choose average or median smooth for removing mean emf")
       endif
!
       df(l1:l2,m,n,iax) = df(l1:l2,m,n,iax) - puxb(:,1)
       df(l1:l2,m,n,iay) = df(l1:l2,m,n,iay) - puxb(:,2)
     endif
!
     br=p%bb(:,1)*p%pomx+p%bb(:,2)*p%pomy - bavg(:,1)
     bp=p%bb(:,1)*p%phix+p%bb(:,2)*p%phiy - bavg(:,2)
     bz=p%bb(:,3)                         - bavg(:,3)
!
     if (ldiagnos) then
       if (idiag_brm/=0)    call sum_lim_mn_name(br,idiag_brm,p)
       if (idiag_bpm/=0)    call sum_lim_mn_name(bp,idiag_bpm,p)
       if (idiag_bzm/=0)    call sum_lim_mn_name(bz,idiag_bzm,p)
       if (idiag_br2m/=0)   call sum_lim_mn_name(br**2,idiag_br2m,p)
       if (idiag_bp2m/=0)   call sum_lim_mn_name(bp**2,idiag_bp2m,p)
       if (idiag_bzz2m/=0)  call sum_lim_mn_name(bz**2,idiag_bzz2m,p)
       if (idiag_brbpm/=0)  call sum_lim_mn_name(br*bp,idiag_brbpm,p)
       if (idiag_bzbpm/=0)  call sum_lim_mn_name(bz*bp,idiag_bzbpm,p)
       if (idiag_brbzm/=0)  call sum_lim_mn_name(br*bz,idiag_brbzm,p)
     endif
!
     if (l1davgfirst.and.idiag_brbpmr/=0) &
          call phizsum_mn_name_r(br*bp,idiag_brbpmr)
!
     call keep_compiler_quiet(f,df)
!
    endsubroutine special_calc_magnetic
!***********************************************************************
    subroutine special_calc_entropy(f,df,p)
!
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_entropy
!***********************************************************************
    subroutine special_before_boundary(f)
!
!   Possibility to modify the f array before the boundaries are
!   communicated.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   Called from equ, but before the evolving loop that calls the dynamical equations.
!   Set the 1D coarse average here. This average is LOCAL to this special module   
!
!   calculate the averages
!
!   06-jul-06/tony: coded
!
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(nrcylrun) :: ktot1,s_rho,rho_sum
      real, dimension(nrcylrun,3) :: s_u,s_b,u_sum,b_sum
      real, dimension(nx,3) :: uuf,bbf,pbb
      integer, dimension(nrcylrun) :: k,ktot
      real, dimension(nx) :: prcyl_mn,prcyl_mn1,sin,cos,puu1,puu2,puu3,prho
      real, dimension(nx) :: pbb1,pbb2,pbb3
      integer :: i,j,ir
      logical :: lfp,llp
      real :: rloop_int,rloop_ext
!
      if (lmagnetic)             bavg_coarse=0.
      if (lhydro)                uavg_coarse=0
      if (lcalc_density_pars)  rhoavg_coarse=0.
!
      do m=m1,m2
        do n=n1,n2
!           
          lfp=((m==m1).and.(n==n1))
          llp=((m==m2).and.(n==n2))
          k=0;s_u=0;s_b=0;s_rho=0
!
          prcyl_mn=max(sqrt(x(l1:l2)**2+y(m)**2),tini)
          prcyl_mn1=1./prcyl_mn
          sin=y(m)*prcyl_mn1 ; cos= x(l1:l2)*prcyl_mn1
!
          if (lcalc_density_pars) prho=f(l1:l2,m,n,ilnrho)
          if (lhydro) then
            puu1=f(l1:l2,m,n,iux);puu2=f(l1:l2,m,n,iuy);puu3=f(l1:l2,m,n,iuz)
            uuf(:,1)=  puu1*cos+puu2*sin
            uuf(:,2)= -puu1*sin+puu2*cos
            uuf(:,3)=  puu3
          endif
!
          if (lmagnetic) then
            call curl(f,iaa,pbb)
            pbb1=pbb(:,1);pbb2=pbb(:,2);pbb3=pbb(:,3)+B_ext
            bbf(:,1)=  pbb1*cos+pbb2*sin
            bbf(:,2)= -pbb1*sin+pbb2*cos
            bbf(:,3)=  pbb3
          endif
!
          do ir=1,nrcylrun
            rloop_int = rt_int + (ir-1)*drc
            rloop_ext = rt_int +  ir   *drc
            do i=1,nx
              if ((prcyl_mn(i).le.rloop_ext).and.(prcyl_mn(i).ge.rloop_int)) then
                k(ir)=k(ir)+1
                if (lcalc_density_pars) s_rho(ir) = s_rho(ir) + prho(i)
                do j=1,nd
                  if (lhydro)           s_u(ir,j) = s_u(ir,j) + uuf(i,j)
                  if (lmagnetic)        s_b(ir,j) = s_b(ir,j) + bbf(i,j)
                enddo
              endif
            enddo
          enddo
!
          if (lfp) then
            k_tmp=k
            if (lcalc_density_pars) rho_tmp=s_rho
            if (lhydro)               u_tmp=s_u
            if (lmagnetic)            b_tmp=s_b
          else
            k_tmp = k_tmp + k
            if (lcalc_density_pars) rho_tmp = rho_tmp+s_rho
            if (lhydro)             u_tmp   =   u_tmp+s_u
            if (lmagnetic)          b_tmp   =   b_tmp+s_b
          endif
!
          if (llp) then
            call mpireduce_sum_int(k_tmp,ktot,nrcylrun)
            if (lcalc_density_pars) call mpireduce_sum(rho_tmp,   rho_sum     ,nrcylrun)
            do j=1,nd
              if (lhydro)           call mpireduce_sum(  u_tmp(:,j),u_sum(:,j),nrcylrun)
              if (lmagnetic)        call mpireduce_sum(  b_tmp(:,j),b_sum(:,j),nrcylrun)
            enddo
!
            call mpibcast_int(ktot,nrcylrun)
            if (lcalc_density_pars) call mpibcast_real(rho_sum     ,nrcylrun)
            do j=1,nd
              if (lhydro)           call mpibcast_real(  u_sum(:,j),nrcylrun)
              if (lmagnetic)        call mpibcast_real(  b_sum(:,j),nrcylrun)
            enddo
!
            if (any(ktot == 0)) &
                 call error("set_new_average","ktot=0")
            ktot1=1./ktot
            if (lcalc_density_pars) rhoavg_coarse   =rho_sum     *ktot1
            do j=1,nd
              if (lhydro)           uavg_coarse(:,j)=  u_sum(:,j)*ktot1
              if (lmagnetic)        bavg_coarse(:,j)=  b_sum(:,j)*ktot1
            enddo
          endif
!
        enddo
      enddo
!
      call keep_compiler_quiet(f)
!
    endsubroutine special_before_boundary
!**********************************************************************
!**********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************

endmodule Special

