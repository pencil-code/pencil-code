! $Id$

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
!  This file adds the terms to the shearing box equations that 
!  come from an underlying large scale pressure gradient. If the pressure
!  falls like p=p0/(r/R)**beta, then a Taylor expansion around R leads to
!  p=p0*(1-beta/R0*x), since r=R0+x. This extra term enters on the equations
!  as 
!
!     d(ux)/dt = usual terms + beta*p0/(rho*R0) - beta*p0/(rho0*R0) 
!     d(EE)/dt = usual terms + beta*E0*ux/R0
!
!  The last term on the RHS of the momentum equation is because 
!  the underlying pressure gradient also leads to a further reduction 
!  of the rotational velocity. Having only the second term implemented 
!  would lead to a large scale wind as the momentum equation tries to 
!  adjust to the new rotational speed. 
!  
!  10-apr-09/wlad: coded
!
module Special
!
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include '../special.h'
! Global arrays
  real :: dummy,rho01,cs201,gammam11,Bshear=0.,p0,p01,nw1
  real :: cv1
  logical :: lunstratified=.false.,lstratification=.true.
  logical :: lstatic_stratification=.false.
  real, dimension (nx) :: strat
  real, dimension (nz) :: rtime_strat
!
  namelist /special_init_pars/ dummy
!   
  namelist /special_run_pars/ Bshear,lunstratified,&
      lstatic_stratification
!
  integer :: idiag_pstratm=0,idiag_pstratmax=0,idiag_pstratmin=0
!
  contains

!***********************************************************************
    subroutine register_special()
!
!  Configure pre-initialised (i.e. before parameter read) variables 
!  which should be know to be able to evaluate
!
!  14-jul-09/wlad: coded
!
      use Cdata
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
!  14-jul-09/wlad: coded
!
      use Cdata 
      use Mpicomm
      use EquationOfState, only: rho0,gamma_m1,cs20,gamma1,get_cp1,gamma
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real :: cp1
      logical :: lstarting
!
!  Define if we want stratification
!
      if (nzgrid==1.or.&
         (nzgrid/=1.and.lunstratified)) lstratification=.false.
!
!  Shortcuts for inverse variables
!
      rho01=1./rho0
      cs201=1./cs20
      gammam11=1./gamma_m1
!
!  Shortcut for the initial pressure
!
      p0 =rho0*cs20*gamma1
      p01=1./p0
!
!  Shortcut for normalization nxgrid*nygrid 
!
      nw1=1./(nxgrid*nygrid)
!
!  Get cv1
!
      call get_cp1(cp1)
      cv1=gamma*cp1
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine pencil_criteria_special()
! 
!  All pencils that this special module depends on are specified here.
! 
!  18-07-06/wlad: coded
!
      lpenc_requested(i_rho1)=.true.
      lpenc_requested(i_uu)=.true.
!
      if (lentropy.or.&
          (ltemperature.and.(.not.ltemperature_nolog))) &
          lpenc_requested(i_TT1)=.true.
!
      if (ltemperature.or.&
          (lentropy.and.pretend_lnTT)) &
          lpenc_requested(i_cv1)=.true.
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!   14-jul-09/wlad: coded
!
      use Mpicomm
      use Gravity, only: potential
      use EquationOfState, only: gamma,cs20
!
      real, dimension(mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      real, dimension (nx) :: pot
      integer :: i
!
!  Calculate the stratification pencil
!
      if (lstratification) then 
!
!  Choose between static (exp{-z^2/2H^2}) or 
!  time-varying stratification (calculated in 
!  special_before_boundary. 
!
        if (lstatic_stratification) then 
          call potential(x(l1:l2),y(m),z(n),pot=pot)
          strat=exp(-gamma*pot*cs201)
        else
          do i=1,nx;strat(i)=rtime_strat(n-n1+1);enddo
        endif
      else !no stratification  
        strat=1.
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine dspecial_dt(f,df,p)
!
      use Diagnostics
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
!
      if (ldiagnos) then
        if (idiag_pstratm/=0) &
             call sum_mn_name(strat,idiag_pstratm)
        if (idiag_pstratmax/=0) &
             call max_mn_name(strat,idiag_pstratmax)
        if (idiag_pstratmin/=0) &
             call max_mn_name(-strat,idiag_pstratmin,lneg=.true.)
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine special_before_boundary(f)
!
!  This subroutine calculates the mean stratification of the 
!  pressure. This is used on the baroclinic terms on the 
!  radial velocity and the entropy equation. 
!
!   dux/dt = ... -1/rho* dp/dx
!
!  where p = p0 (1-beta*x/R) * f(z), with z being the 
!  stratification function. The initial case is of course 
!  f(z)=exp[-z^2/(2H^2)], and that is good for most cases. 
!  (since the RHS of the vertical equation doesn't change in 
!  time, Omega^2*z). But just for the sake of completeness, 
!  this is a way to have a stratification function f(z,t) that
!  also changes in time. The stratification is this case is 
!  computed as the mean pressure at each level xy.
!
      use Mpicomm
      use Sub
      use EquationOfState, only: gamma1,gamma_m1,&
          cs20,lnrho0,rho0
      use Gravity
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (nx) :: cs2,rho,lnrho,ss
      real, dimension (nx,nz) :: pp_tmp,pp_sumy
      real, dimension (nz) :: pp_sum
      integer :: nn,i
!
      if (lstratification.and.&
          (.not.lstatic_stratification)) then 
!
!  Calculate mean pressure stratification at each xy-level
!        
        pp_tmp=0.
!
        do n=n1,n2
          do m=m1,m2
!
!  Get density and sound speed
!
            if (ldensity_nolog) then 
              rho=f(l1:l2,m,n,ilnrho)
              lnrho=log(rho)
            else
              lnrho=f(l1:l2,m,n,ilnrho)
              rho=exp(lnrho)
            endif
!
            if (ltemperature.or.pretend_lnTT) then 
              call fatal_error("lazy wlad didn't want to code entropy","")
            else
              ss=f(l1:l2,m,n,iss)
            endif
            cs2=cs20*exp(cv1*ss+gamma_m1*(lnrho-lnrho0))
!            
!  Average (sum) the pressure of all y-grid cells. The result
!  is xz dependent.
!
            nn=n-n1+1
            pp_tmp(:,nn) = pp_tmp(:,nn)+rho*cs2*gamma1
!
          enddo
        enddo
!
!  Now perform the y-sum over all y-processors (that's
!  why lsumy=T). The result is still xz dependent.
!
        !call mpiallreduce_sum(pp_tmp,pp_sumy,(/nx,nz/),LSUMY=.true.)
        call mpiallreduce_sum(pp_tmp,pp_sumy,(/nx,nz/),idir=2)
!
!  Now sum in x. The result is z-dependent.
!
        pp_sum=0.
        do n=1,nz;do i=1,nx
          pp_sum(n)=pp_sum(n)+pp_sumy(i,n)
        enddo;enddo
!
!  Now divide the sum by the number of grid cells to get 
!  the average (pp_sum*nw1). And normalize by the midplane
!  initial pressure p0=rho0*cs0^2/gamma
!
        rtime_strat=pp_sum*nw1*p01
!
      endif

    endsubroutine special_before_boundary
!***********************************************************************
    subroutine read_special_init_pars(unit,iostat)
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
!
      if (present(iostat)) then
        read(unit,NML=special_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=special_run_pars,ERR=99)
      endif
!
99    return
endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
      integer, intent(in) :: unit
!
      write(unit,NML=special_run_pars)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
      use Diagnostics
!
!  reads and registers print parameters relevant to special
!
!   14-jul-09/wlad: coded
!
      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Write information to index.pro
!
      if (lreset) then
        idiag_pstratm=0;idiag_pstratmax=0;idiag_pstratmin=0 
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'pstratm',idiag_pstratm)
        call parse_name(iname,cname(iname),cform(iname),'pstratmax',idiag_pstratmax)
        call parse_name(iname,cname(iname),cform(iname),'pstratmin',idiag_pstratmin)
      enddo
!
      if (lwr) then
        write(3,*) 'i_pstratm=',idiag_pstratm
        write(3,*) 'i_pstratmax=',idiag_pstratmax
        write(3,*) 'i_pstratmin=',idiag_pstratmin
      endif
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the 
!   momentum equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
      use Cdata
      use Messages, only: fatal_error
      use EquationOfState, only : rho0
!      
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
      real, dimension (nx) :: rho1
!
!  Modified momentum equation
!
      rho1=p%rho1*strat
!      
      df(l1:l2,m,n,iux)= df(l1:l2,m,n,iux)+Bshear*p0*(rho1-rho01)
!
      call keep_compiler_quiet(f)
!
    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine special_calc_entropy(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the 
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   14-jul-09/wlad: coded
!
      use Cdata
      use Messages, only: fatal_error
!
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
      real, dimension (nx) :: rhs
!
!  Right hand side on the energy equation - background energy gradient
!
      rhs=Bshear*p0*p%uu(:,1)*gammam11*strat
!
      if (lentropy) then 
        if (pretend_lnTT) then 
          df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss) + p%cv1*p%rho1*p%TT1*rhs
        else
          df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss) + p%rho1*p%TT1*rhs
        endif
      else if (ltemperature) then 
        if (ltemperature_nolog) then 
          df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT) + p%cv1*p%rho1*rhs
        else
          df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT) + p%cv1*p%rho1*p%TT1*rhs
        endif
      else
        print*,"You want to use a global baroclinic term but    "
        print*,"you are NOT solving the energy equation. Better "
        print*,"stop and check."
        call fatal_error("global_baroclinic","")
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine special_calc_entropy
!***********************************************************************
!********************************************************************
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

