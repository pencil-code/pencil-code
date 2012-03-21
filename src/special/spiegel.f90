! $Id$
!
!  This module is currently fairly much obsolete.
!  This module incorporates all the modules used for Natalia's
!  neutron star -- disk coupling simulations (referred to as nstar)
!
!  This sample modules solves a special set of problems related
!  to computing the accretion through a thin disk onto a rigid surface
!  and the spreading of accreted material along the neutron star's surface.
!  One-dimensional problems along the disc and perpendicular to it
!  can also be considered.
!
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
!    SPECIAL=special/nstar
!
! Where nstar it replaced by the filename of your new module
! upto and not including the .f90
!
!--------------------------------------------------------------------
!
module Special
!
  use Cdata
  use Cparam
  use EquationOfState
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include '../special.h'
!
  ! input parameters
!
!
  character (len=labellen) :: initnstar='default'
  real :: rho_bot=1.,rho_top=1.
!
  real :: uu_left=0.
  real :: uy_left=0.,uy_right=0.
 real :: R_star=0.
  real :: M_star=0.
  real :: T_bot=0.
!
  logical :: ldecelerat_zone=.false.
  logical :: lsurface_zone=.false.
  logical :: lnstar_T_const=.false.
  logical :: lnstar_entropy=.false.
!
  integer :: ac_dc_size=5
!
  logical :: lraddif_local=.false.
!
! Keep some over used pencils
!
  real, dimension(nx) :: z_2
!
! start parameters
  namelist /neutron_star_init_pars/ &
      initnstar, ldecelerat_zone, lsurface_zone, &
      rho_bot,rho_top, &
       lraddif_local, R_star, M_star, &
      T_bot,&
      uu_left, uy_left, uy_right, &
      lnstar_entropy
!
!
! run parameters
  namelist /neutron_star_run_pars/ &
      rho_bot,rho_top, &
     ldecelerat_zone, lsurface_zone,lraddif_local, &
     R_star, M_star, T_bot, &
       lnstar_entropy
!
! Declare any index variables necessary for main or
!
!   integer :: iSPECIAL_VARIABLE_INDEX=0

! other variables (needs to be consistent with reset list below)
!
  integer :: idiag_dtcrad=0
  integer :: idiag_dtchi=0
!
  contains
!
!***********************************************************************
    subroutine register_special()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  6-oct-03/tony: coded
!
!  Identify CVS/SVN version information.
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
      real, dimension (mx,my,mz,mvar+maux) :: f
      logical :: lstarting
!
!  Initialize any module variables which are parameter dependent
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
!!
      select case (initnstar)
        case ('default')
          if (lroot) print*,'init_special: Default neutron star setup'
          call density_init(f)
          call entropy_init(f)
          call velocity_init(f)
!
        case default
          !
          !  Catch unknown values
          !
          if (lroot) print*,'init_special: No such value for initnstar: ', trim(initnstar)
          call stop_it("")
      endselect
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special()
!
!  All pencils that this special module depends on are specified here.
!
!  18-07-06/tony: coded
!
!
!
!
   if (lnstar_entropy) then
         lpenc_requested(i_TT)=.true.
          lpenc_requested(i_lnTT)=.true.
         lpenc_requested(i_cs2)=.true.
         lpenc_requested(i_ss)=.true.
         lpenc_requested(i_rho)=.true.
      endif
!
     if (lraddif_local) then
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_TT)=.true.
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_gss)=.true.
        lpenc_requested(i_del2lnrho)=.true.
        lpenc_requested(i_del2ss)=.true.
        lpenc_requested(i_rho)=.true.
      endif
!
!
!
!
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
      use Diagnostics
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
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
      if (ldiagnos) then
        if (idiag_dtcrad/=0) &
          call max_mn_name(sqrt(advec_crad2)/cdt,idiag_dtcrad,l_dt=.true.)
        if (idiag_dtchi/=0) &
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
      endif
!
! Keep compiler quiet by ensuring every parameter is used
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=neutron_star_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=neutron_star_init_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
      integer, intent(in) :: unit
!
      write(unit,NML=neutron_star_init_pars)
!
    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=neutron_star_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=neutron_star_run_pars,ERR=99)
      endif
!
99    return
endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
      integer, intent(in) :: unit
!
      write(unit,NML=neutron_star_run_pars)
!
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
      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_dtcrad=0
        idiag_dtchi=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'dtcrad',idiag_dtcrad)
        call parse_name(iname,cname(iname),cform(iname),'dtchi',idiag_dtchi)
      enddo
!
!  write column where which magnetic variable is stored
      if (lwr) then
        write(3,*) 'i_dtcrad=',idiag_dtcrad
        write(3,*) 'i_dtchi=',idiag_dtchi
      endif
!
    endsubroutine rprint_special
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
!   calculate a additional 'special' term on the right hand side of the
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      use EquationOfState
!
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
      integer :: i, l_sz, tmp_int
      real :: cs2_star
!
   !    call eoscalc(ilnrho_lnTT,log(rho_star),log(T_bot), cs2=cs2_star)
!
!  mass sources and sinks for the boundary layer on NS in 1D approximation
!
!
      if (lsurface_zone) then
!
!
          if ( dt .GT.0.) then
            l_sz=l2-5
!
          !  df(l_sz:l2,m,n,ilnrho)=df(l_sz:l2,m,n,ilnrho)&
          !       -1./(5.*dt)*(1.-rho_surf/exp(f(l_sz:l2,m,n,ilnrho)))
!
!
!
!
         endif
      endif
!
!
!
! Keep compiler quiet by ensuring every parameter is used
!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_density
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   16-jul-06/natalia: coded
!
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      integer :: j,l_sz
!
! add effective gravity term = -Fgrav+Fcentrifugal
! Natalia
!
     l_sz=l2-5
!
       df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)-M_star/z(n)**3*x(l1:l2)
    !   df(l1:l2,m,n,iuz)=df(l_sz:l2,m,n,iuz)&
    !               -1./(2.*dt)*(f(l_sz:l2,m,n,iuz)-0.)
    !   df(l1:l2,m,n,iuy)=df(l_sz:l2,m,n,iuy)&
    !               -1./(2.*dt)*(f(l_sz:l2,m,n,iuy)-0.)
!
!
! surface zone in a case of a Keplerian disk
!
      if (lsurface_zone) then
        if ( dt .gt.0.) then
!
           df(l_sz:l2,m,n,iux)=df(l_sz:l2,m,n,iux)&
                   -1./(2.*dt)*(f(l_sz:l2,m,n,iux)-0.)
!
         endif
!
      endif
!
    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine special_calc_magnetic(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!
!!
!!  df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) + SOME NEW TERM
!!  df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) + SOME NEW TERM
!!  df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) + SOME NEW TERM
!!
!!
!
! Keep compiler quiet by ensuring every parameter is used
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_magnetic
!!***********************************************************************
    subroutine special_calc_entropy(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
      integer :: j, l_sz, l_sz_1
!
        if (lraddif_local) call raddif_local(f,df,p)
!
!
       if (lsurface_zone) then
          if ( dt .GT.0.) then
            l_sz=l2-5
            l_sz_1=nxgrid-5
!
        !  if (lnstar_1D) then
     !       df(l_sz:l2,m,n,iss)=df(l_sz:l2,m,n,iss) &
     !       -1./(5.*dt)*(f(l_sz:l2,m,n,iss)-log(T_disk)/gamma) &
     !       /p%rho(l_sz_1:nxgrid)/p%TT(l_sz_1:nxgrid)
        !  else
!
        !    do j=l_sz,l2
        !     df(j,m,n,iss)=df(j,m,n,iss)&
        !       -1./(5.*dt)*(f(j,m,n,iss)-f(j-1,m,n,iss))
        !    enddo
!
         ! endif
!
!
         endif
      endif
!
! Keep compiler quiet by ensuring every parameter is used
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_entropy
!***********************************************************************
    subroutine special_boundconds(f,bc)
!
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      type (boundary_condition) :: bc
!
      select case (bc%bcname)
       case ('stp')
         select case (bc%location)
         case (iBC_X_TOP)
           call bc_BL_x(f,-1, bc)
         case (iBC_X_BOT)
           call bc_BL_x(f,-1, bc)
         endselect
         bc%done=.true.
      endselect
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(bc%bcname)
!
    endsubroutine special_boundconds
!***********************************************************************
!
!  PRIVATE UTITLITY ROUTINES
!
!***********************************************************************
!
 subroutine raddif_local(f,df,p)
!
!  heat conduction
!  Natalia (NS)
!   12-apr-06/axel: adapted from Wolfgang's more complex version
!
      use Sub, only: max_mn_name,dot
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      type (pencil_case) :: p
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: glnT,glnThcond !,glhc
      real, dimension (nx) :: chix,diffus_chi1
      real, dimension (nx) :: thdiff,g2,thdiff_1D
      real, dimension (nx) :: hcond
      real ::  beta
      integer :: l_sz, l_sz_1
!
      intent(in) :: f,p
      intent(out) :: df
!
!
!  Heat conduction
!
      chix = p%rho1*p%rho1*p%TT**3*16./3.*sigmaSB/kappa_es!hcond
      glnT = gamma*p%gss + gamma_m1*p%glnrho ! grad ln(T)
      glnThcond = glnT !... + glhc/spread(hcond,2,3)    ! grad ln(T*hcond)
      call dot(glnT,glnThcond,g2)
!
!AB:  derivs of chix missing??
!
      thdiff = chix * (gamma*p%del2ss+gamma_m1*p%del2lnrho + g2)
!print*,' p%del2ss(10)'   ,p%del2ss(10)
!
   !  add heat conduction to entropy equation
    !
!
!
         df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff
         if (headtt) print*,'calc_heatcond_diffusion: added thdiff'
!
        df(l1:l2,m,n,iuz) = &
         df(l1:l2,m,n,iuz)-p%rho1*16./3.*sigmaSB/c_light*p%TT**4*glnT(:,3)
!
        df(l1:l2,m,n,iuy) = &
          df(l1:l2,m,n,iuy)-p%rho1*16./3.*sigmaSB/c_light*p%TT**4*glnT(:,2)
!
        df(l1:l2,m,n,iux) = &
         df(l1:l2,m,n,iux)-p%rho1*16./3.*sigmaSB/c_light*p%TT**4*glnT(:,1)
!
!  include constraint from radiative time step
!
      if (lfirst.and.ldt) then
        advec_crad2=p%rho1*16./3.*sigmaSB/c_light*p%TT**4
      endif
!
     if (headtt) print*,'calc_radiation_pressure: added to z-component'
!
!  check maximum diffusion from thermal diffusion
!  With heat conduction, the second-order term for entropy is
!  gamma*chix*del2ss
!
      if (lfirst.and.ldt) then
! Calculate timestep limitation
        diffus_chi=max(diffus_chi,gamma*chix*dxyz_2)
!
     !   diffus_chi1=min(gamma*chix*dxyz_2, &
      !              real(sigmaSB*kappa_es*p%TT**3*4.*p%cp1tilde))
     !   diffus_chi=max(diffus_chi,diffus_chi1)
      endif
    endsubroutine raddif_local
!*************************************************************************
    subroutine density_init(f)
!
!Natalia
!Initialization of density in a case of the step-like distribution
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (my,mz) :: lnrho_2d
!
      real ::   ln_ro_l, ln_ro_r, ln_ro_u
      real :: cs2_star
!
        ln_ro_r=log(rho_bot)
        ln_ro_u=log(rho_top)
!
        do n=n1,n2; do m=m1,m2
          f(l1:l2,m,n,ilnrho)=(x(l1:l2)-0.)/Lxyz(1)*(ln_ro_u-ln_ro_r)+ln_ro_r
        enddo; enddo
!
    endsubroutine density_init
!***************************************************************
    subroutine entropy_init(f)
!Natalia
!Initialization of entropy in a case of the step-like distribution
 use EquationOfState
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (nx) ::  lnrho, lnTT,ss
      integer ::  mi,ni
!
!
!
!
 !     lnTT=log(T_bot)!  log(T0)
! if (T_disk.EQ.0) then
 !    T_disk=cs0**2/gamma_m1
 !  endif
!
!
!
      do ni=n1,n2;
       do mi=m1,m2;
!
      ! lnrho=f(l1:l2,mi,ni,ilnrho)
      ! const_tmp=M_star/sigmaSB*c_light*3./4.
!
     !  call eoscalc(4,lnrho,lnTT,ss=ss)
!
      ! f(l1:l2,mi,ni,iss)=ss
!
!
        f(l1:l2,mi,ni,iss)=-f(l1:l2,mi,ni,ilnrho)*gamma_m1/gamma
!
!
       end do
     end do
!
!
    endsubroutine entropy_init
!**********************************************************************
    subroutine velocity_init(f)
!Natalia
!Initialization of velocity in a case of the step-like distribution
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
!
        f(:,:,:,iuz)=uu_left
        f(:,:,:,iuy)=uu_left
        f(:,:,:,iux)=uu_left
!
!
!
    endsubroutine  velocity_init
!***********************************************************************
!
    subroutine bc_BL_x(f,sgn,bc)
!
! Natalia
!  11-may-06
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer :: sgn
      type (boundary_condition) :: bc
      integer :: i,j,ni,mi
        real :: lnrho,lnTT,ss
!
      j=bc%ivar
      if (bc%location==iBC_X_BOT) then
      ! bottom boundary
!
    !         f(l1,m1:m2,n1:n2,j)=bc%value1
!
       if (j==5) then
!
      do ni=n1,n2;
       do mi=m1,m2;
!
       lnrho=f(l1,mi,ni,ilnrho)
       lnTT=log(T_bot)
!
         call eoscalc(4,lnrho,lnTT,ss=ss)
!
        f(l1,mi,ni,iss)=ss
!
       enddo
      enddo
          else
            f(l1,m1:m2,n1:n2,j)=bc%value1
          endif
!
!
           do i=1,nghost; f(l2+i,:,:,j)=2*f(l2,:,:,j)+sgn*f(l2-i,:,:,j); enddo
!
      elseif (bc%location==iBC_X_TOP) then
      ! top boundary
!
          if (j == 1) then
            f(l2,m1:m2,n1:n2,j)=bc%value1
            do i=1,nghost; f(l2+i,:,:,j)=2*f(l2,:,:,j)+sgn*f(l2-i,:,:,j); enddo
          else
            f(l2+1,:,:,j)=0.25*(  9*f(l2,:,:,j)- 3*f(l2-1,:,:,j)- 5*f(l2-2,:,:,j)+ 3*f(l2-3,:,:,j))
            f(l2+2,:,:,j)=0.05*( 81*f(l2,:,:,j)-43*f(l2-1,:,:,j)-57*f(l2-2,:,:,j)+39*f(l2-3,:,:,j))
            f(l2+3,:,:,j)=0.05*(127*f(l2,:,:,j)-81*f(l2-1,:,:,j)-99*f(l2-2,:,:,j)+73*f(l2-3,:,:,j))
          endif
!
      else
        print*, "bc_BL_x: ", bc%location, " should be `top(", &
                        iBC_X_TOP,")' or `bot(",iBC_X_BOT,")'"
      endif
!
    endsubroutine bc_BL_x
 !***********************************************************************
!
!***********************************************************************
    subroutine special_before_boundary(f)
!
!   Possibility to modify the f array before the boundaries are
!   communicated.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-jul-06/tony: coded
!
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine special_before_boundary
!
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************
endmodule Special
