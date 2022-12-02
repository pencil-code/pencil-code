! $Id$
!
!  Solve the lucky droplet model for many realizations.
!  The different realizations correspond to "meshpoints".
!  To add the contributions for each step, we use the usual
!  time step in the Pencil Code, so t is just the step, and
!  the accumulated collision times (after 125 steps or so)
!  for all realizations at the same time are the values in
!  the f-array.
!
!  16-apr-20/axel: adapted from nospecial.f90
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of special_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
!***************************************************************
!
module Special
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include '../special.h'
!
! Declare index of variables
!
  integer :: iaxionSU2back_Q=0, iaxionSU2back_Qdot=0, iaxionSU2back_chi=0, iaxionSU2back_chidot=0
  integer :: iaxionSU2back_psi=0, iaxionSU2back_psidot=0, iaxionSU2back_TR=0, iaxionSU2back_TRdot=0
!
  ! input parameters
  real :: k=1e-2, fdecay=.003, g=1.11e-2, lam=500., mu=1.5e-4
  real :: Q=3e-4, Qdot=0., chi_prefactor=.49, chidot=0., H=1.04e-6
  real :: chi, psi, psidot, a, TR, TRdot
  character(len=50) :: init_axionSU2back='zero'
  namelist /special_init_pars/ &
    k, fdecay, g, lam, mu, Q, Qdot, chi_prefactor, chidot, H
!
  ! run parameters
  namelist /special_run_pars/ &
    k, fdecay, g, lam, mu, H
!
! other variables (needs to be consistent with reset list below)
!
  integer :: idiag_Q =0 ! DIAG_DOC: $Q$
!
  contains
!
!***********************************************************************
    subroutine register_special()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  19-feb-2019/axel: coded
!
      use FArrayManager
      use SharedVariables, only: put_shared_variable
!
!  Identify CVS/SVN version information.
!
      if (lroot) call svn_id( &
           "$Id$")
!
!  Set iaxionSU2back to consecutive numbers
!
      call farray_register_pde('axionSU2back_Q',iaxionSU2back_Q)
      call farray_register_pde('axionSU2back_Qdot',iaxionSU2back_Qdot)
      call farray_register_pde('axionSU2back_chi',iaxionSU2back_chi)
      call farray_register_pde('axionSU2back_chidot',iaxionSU2back_chidot)
      call farray_register_pde('axionSU2back_psi',iaxionSU2back_psi)
      call farray_register_pde('axionSU2back_psidot',iaxionSU2back_psidot)
      call farray_register_pde('axionSU2back_TR',iaxionSU2back_TR)
      call farray_register_pde('axionSU2back_TRdot',iaxionSU2back_TRdot)
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  19-feb-2019/axel: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  Initialize any module variables which are parameter dependent
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise special condition; called from start.f90
!   2-dec-2022/axel: coded
!
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      intent(inout) :: f
!
!  Initial condition; same for every population.
!
      select case (init_axionSU2back)
        case ('nothing'); if (lroot) print*,'nothing'
        case ('standard')
          a=exp(H*t)
          psi=(a/sqrt(2.*k))
          psidot=psi*k
          TR=(a/sqrt(2.*k))
          TRdot=TR*k
          chi=chi_prefactor*pi*fdecay
          f(l1:l2,m,n,iaxionSU2back_Q)=Q
          f(l1:l2,m,n,iaxionSU2back_Qdot)=Qdot
          f(l1:l2,m,n,iaxionSU2back_chi)=chi
          f(l1:l2,m,n,iaxionSU2back_chidot)=chidot
          f(l1:l2,m,n,iaxionSU2back_psi)=psi
          f(l1:l2,m,n,iaxionSU2back_psidot)=psidot
          f(l1:l2,m,n,iaxionSU2back_TR)=TR
          f(l1:l2,m,n,iaxionSU2back_TRdot)=TRdot
!
        case default
          !
          !  Catch unknown values
          !
          if (lroot) print*,'init_axionSU2back: No such value for init_axionSU2back: ', trim(init_axionSU2back)
          call stop_it("")
      endselect
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special()
!
!  All pencils that this special module depends on are specified here.
!
!  18-07-06/tony: coded
!
    endsubroutine pencil_criteria_special
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
    subroutine calc_pencils_special(f,p)
!
!  Calculate Hydro pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!   24-nov-04/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_special
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
!   2-dec-2022/axel: coded
!
      use General, only: random_number_wrapper
      use Diagnostics
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: Q, Qdot, chi, chidot
      real, dimension (nx) :: psi, psidot, TR, TRdot
      real, dimension (nx) :: Uprime, lamf, mQ, xi, a, epsQE, epsQB
      real :: Mpl2=1., Hdot=0.
      type (pencil_case) :: p
!
      intent(in) :: f,p
      intent(inout) :: df
!
      call keep_compiler_quiet(p)
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dSPECIAL_dt'
!
!  Set the 8 variable
!
      Q=f(l1:l2,m,n,iaxionSU2back_Q)
      Qdot=f(l1:l2,m,n,iaxionSU2back_Qdot)
      chi=f(l1:l2,m,n,iaxionSU2back_chi)
      chidot=f(l1:l2,m,n,iaxionSU2back_chidot)
      psi=f(l1:l2,m,n,iaxionSU2back_psi)
      psidot=f(l1:l2,m,n,iaxionSU2back_psidot)
      TR=f(l1:l2,m,n,iaxionSU2back_TR)
      TRdot=f(l1:l2,m,n,iaxionSU2back_TRdot)
!
!  Set parameters
!
      Uprime=-mu**4/fdecay*sin(chi/fdecay)
      lamf=lam/fdecay
      mQ=g*Q/H
      xi=lamf*chidot/(2.*H)
      a=exp(H*t)
      epsQE=(Qdot+H*Q)**2/(Mpl2*H**2)
      epsQB=g**2*Q**4/(Mpl2*H**2)
!
!  background
!
      df(l1:l2,m,n,iaxionSU2back_Q)=df(l1:l2,m,n,iaxionSU2back_Q)+Qdot
      df(l1:l2,m,n,iaxionSU2back_Qdot)=df(l1:l2,m,n,iaxionSU2back_Qdot) &
        +g*lamf*chidot*Q**2-3.*H*Qdot-(Hdot+2*H**2)*Q-2.*g**2*Q**3
      df(l1:l2,m,n,iaxionSU2back_chi)=df(l1:l2,m,n,iaxionSU2back_chi)+chidot
      df(l1:l2,m,n,iaxionSU2back_chidot)=df(l1:l2,m,n,iaxionSU2back_chidot) &
        -3.*g*lamf*Q**2*(Qdot+H*Q)-3.*H*chidot-Uprime
!
!  perturbation
!
      df(l1:l2,m,n,iaxionSU2back_psi)=df(l1:l2,m,n,iaxionSU2back_psi)+psidot
      df(l1:l2,m,n,iaxionSU2back_psidot)=df(l1:l2,m,n,iaxionSU2back_psidot) &
        -H*psidot-(k**2/a**2-2.*H**2)*psi-2.*H*sqrt(epsQE)*TRdot+2.*H**2*sqrt(epsQB)*(mQ-k/(a*H))*TR
      df(l1:l2,m,n,iaxionSU2back_TR)=df(l1:l2,m,n,iaxionSU2back_TR)+TRdot
      df(l1:l2,m,n,iaxionSU2back_TRdot)=df(l1:l2,m,n,iaxionSU2back_TRdot) &
        -H*TRdot-(k**2/a**2+2.*H**2*(mQ*xi-k/(a*H)*(mQ+xi)))*TR+2.*H*sqrt(epsQE)*psidot &
        +2.*H**2*(sqrt(epsQB)*(mQ-k/(a*H))+sqrt(epsQE))*psi

!  diagnostics
!
      if (ldiagnos) then
        if (idiag_Q/=0) call sum_mn_name(Q,idiag_Q)
      endif
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=special_init_pars, IOSTAT=iostat)
!
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=special_init_pars)
!
    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=special_run_pars, IOSTAT=iostat)
!
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=special_run_pars)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  reads and registers print parameters relevant to special
!
!  19-feb-2019/axel: coded
!
      use Diagnostics
      use FArrayManager, only: farray_index_append
      use Sub
!
!   SAMPLE IMPLEMENTATION
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
        idiag_Q=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'Q' ,idiag_Q)
      enddo
!
    endsubroutine rprint_special
!***********************************************************************
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
