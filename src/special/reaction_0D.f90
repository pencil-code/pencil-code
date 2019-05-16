! $Id$
!
!  Solve for a set of two ODEs, used to test time step
!
!  19-apr-09/axel: adapted from nospecial.f90
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of special_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 3
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
   integer :: ispecial=0,ispecial1=0,ispecial2=0
!
  ! input parameters
  integer :: AA0=10, DD0=0, LL0=0
  character(len=50) :: init_qq='set'
  namelist /special_init_pars/ &
    init_qq, AA0, DD0, LL0
!
  ! run parameters
  real :: kp=.0, km=.0, kC=.9, kX=.1
  real :: fidelity=1., betaD=.0, betaL=.0
  real :: fidelity_factor1, fidelity_factor2
  real, dimension (nx) :: Ntot
  namelist /special_run_pars/ &
    kp, km, kC, kX, fidelity, betaD, betaL
!
! other variables (needs to be consistent with reset list below)
!
  integer :: idiag_eem=0  ! DIAG_DOC: $\langle \eta \rangle$
  integer :: idiag_ee1m=0 ! DIAG_DOC: $\langle |\eta| \rangle$
  integer :: idiag_ee2m=0 ! DIAG_DOC: $\langle \eta^2 \rangle$
  integer :: idiag_ee3m=0 ! DIAG_DOC: $\langle \eta^3 \rangle$
  integer :: idiag_ee4m=0 ! DIAG_DOC: $\langle \eta^4 \rangle$
  integer :: idiag_ee10=0 ! DIAG_DOC: $\langle \eta_{10\%} \rangle$
  integer :: idiag_ee50=0 ! DIAG_DOC: $\langle \eta_{50\%} \rangle$
  integer :: idiag_ee90=0 ! DIAG_DOC: $\langle \eta_{90\%} \rangle$
  integer :: idiag_ee99=0 ! DIAG_DOC: $\langle \eta_{99\%} \rangle$
  integer :: idiag_AAm=0  ! DIAG_DOC: $\langle [A] \rangle$
  integer :: idiag_DDm=0  ! DIAG_DOC: $\langle [D] \rangle$
  integer :: idiag_LLm=0  ! DIAG_DOC: $\langle [L] \rangle$
  integer :: idiag_DLm=0  ! DIAG_DOC: $\langle [D]+[L] \rangle$
  integer :: idiag_kC=0   ! DIAG_DOC: $k_C$
  integer :: idiag_A1=0   ! DIAG_DOC: $A_1$
  integer :: idiag_A2=0   ! DIAG_DOC: $A_2$
  integer :: idiag_A3=0   ! DIAG_DOC: $A_3$
  integer :: idiag_A4=0   ! DIAG_DOC: $A_4$
  integer :: idiag_A5=0   ! DIAG_DOC: $A_5$
  integer :: idiag_D1=0   ! DIAG_DOC: $D_1$
  integer :: idiag_D2=0   ! DIAG_DOC: $D_2$
  integer :: idiag_D3=0   ! DIAG_DOC: $D_3$
  integer :: idiag_D4=0   ! DIAG_DOC: $D_4$
  integer :: idiag_D5=0   ! DIAG_DOC: $D_5$
  integer :: idiag_L1=0   ! DIAG_DOC: $L_1$
  integer :: idiag_L2=0   ! DIAG_DOC: $L_2$
  integer :: idiag_L3=0   ! DIAG_DOC: $L_3$
  integer :: idiag_L4=0   ! DIAG_DOC: $L_4$
  integer :: idiag_L5=0   ! DIAG_DOC: $L_5$
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
!
!  Identify CVS/SVN version information.
!
      if (lroot) call svn_id( &
           "$Id$")
!
!  Set ichemistry to consecutive numbers nvar+1, nvar+2, ..., nvar+nchemspec.
!
      call farray_register_pde('qq',iqq,vector=3)
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
!  Cannot have kp+km+kX > 1.
!
      if ((kp+km+kC+kX)>1.) call fatal_error('initialize_special','kp+km+kC+kX > 1')
!
!  Compute fidelity factors.
!
      fidelity_factor1=.5*(1.+fidelity)
      fidelity_factor2=.5*(1.-fidelity)
!
!  Compute total number of molecules of all populations.
!
      Ntot=f(l1:l2,m1,n1,iqq)+f(l1:l2,m1,n1,iqq+1)+f(l1:l2,m1,n1,iqq+2)
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
!  19-feb-2019/axel: coded
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
      select case (init_qq)
        case ('nothing'); if (lroot) print*,'init_qq: nothing'
        case ('racemic'); f(:,:,:,iqq)=AA0
        case ('set')
          f(:,:,:,iqq  )=AA0
          f(:,:,:,iqq+1)=DD0
          f(:,:,:,iqq+2)=LL0
!
        case default
          !
          !  Catch unknown values
          !
          if (lroot) print*,'init_qq: No such value for init_qq: ', trim(init_qq)
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
!  19-feb-2019/axel: coded
!
      use General, only: random_number_wrapper
      use Diagnostics
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: rrr
      real, dimension (nx) ::     d1, d2, d3, d4, d5, d6, d7
      real, dimension (nx) :: r0, r1, r2, r3, r4, r5, r6, r7
      real, dimension (nx) :: qA, qD, qL
      real, dimension (nx) :: dqA, dqD, dqL
      real, dimension (nx) :: ee, ee1, ee10, ee50, ee90, ee99
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
!  Define state vector.
!  The positions l1-l2 refer to different populations.
!
      qA=f(l1:l2,m,n,1)
      qD=f(l1:l2,m,n,2)
      qL=f(l1:l2,m,n,3)
!
!  Determine boundary points
!
      d1=kp*.5
      d2=kp*.5
      d3=km*.5
      d4=km*.5
      d5=kC*(fidelity_factor1*qD+fidelity_factor2*qL+betaD*qA)/Ntot
      d6=kC*(fidelity_factor1*qL+fidelity_factor2*qD+betaL*qA)/Ntot
      d7=kX
!
      r0=0.
      r1=r0+d1
      r2=r1+d2
      r3=r2+d3
      r4=r3+d4
      r5=r4+d5
      r6=r5+d6
      r7=r6+d7  !(Enantiomeric cross inhibition)
!
!  Produce random numbers
!
      call random_number_wrapper(rrr)
!
!  Initialize to zero
!
        dqA=0
        dqD=0
        dqL=0
!
!  Spontaneous formation of chiral molecule *or*
!  autocatalysis of D. 
!
      where (( (rrr >= r0 .and. rrr < r1) .or. &
               (rrr >= r4 .and. rrr < r5) ) .and. qA /= 0)
        dqA=-1
        dqD=+1
        dqL= 0
      endwhere
!
!  Spontaneous formation of chiral molecule *or*
!  autocatalysis of L. 
!
      where (( (rrr >= r1 .and. rrr < r2) .or. &
               (rrr >= r5 .and. rrr < r6) ) .and. qA /= 0)
        dqA=-1
        dqD= 0
        dqL=+1
      endwhere
!
!  Spontaneous racemizatio of D.
!
      where ((rrr >= r2 .and. rrr < r3) .and. qD > 0)
        dqA= 1
        dqD=-1
        dqL= 0
      endwhere
!
!  Spontaneous racemizatio of L.
!
      where ((rrr >= r3 .and. rrr < r4) .and. qL > 0)
        dqA= 1
        dqD= 0
        dqL=-1
      endwhere
!
!  Enantiomeric cross inhibition
!
      where (rrr >= r6 .and. rrr < r7 .and. qD > 0 .and. qL > 0)
        dqA=+2
        dqD=-1
        dqL=-1
      endwhere
!
      df(l1:l2,m,n,1)=df(l1:l2,m,n,1)+dqA
      df(l1:l2,m,n,2)=df(l1:l2,m,n,2)+dqD
      df(l1:l2,m,n,3)=df(l1:l2,m,n,3)+dqL
!
!  diagnostics
!
      if (ldiagnos) then
        where (qD+qL/=0.)
          ee=(qD-qL)/(qD+qL)
        elsewhere
          ee=0.
        endwhere
        ee1=abs(ee)
        if (idiag_eem /=0) call sum_mn_name(ee   ,idiag_eem)
        if (idiag_ee1m/=0) call sum_mn_name(ee1  ,idiag_ee1m)
        if (idiag_ee2m/=0) call sum_mn_name(ee**2,idiag_ee2m)
        if (idiag_ee3m/=0) call sum_mn_name(ee**3,idiag_ee3m)
        if (idiag_ee4m/=0) call sum_mn_name(ee**4,idiag_ee4m)
        if (idiag_AAm/=0) call sum_mn_name(qA/Ntot,idiag_AAm)
        if (idiag_DDm/=0) call sum_mn_name(qD/Ntot,idiag_DDm)
        if (idiag_LLm/=0) call sum_mn_name(qL/Ntot,idiag_LLm)
        if (idiag_DLm/=0) call sum_mn_name((qD+qL)/Ntot,idiag_DLm)
        if (idiag_kC /=0) call   save_name(kC            ,idiag_kC)
        if (idiag_A1 /=0) call   save_name(f(l1+0 ,m,n,1),idiag_A1)
        if (idiag_A2 /=0) call   save_name(f(l1+1 ,m,n,1),idiag_A2)
        if (idiag_A3 /=0) call   save_name(f(l1+2 ,m,n,1),idiag_A3)
        if (idiag_A4 /=0) call   save_name(f(l1+3 ,m,n,1),idiag_A4)
        if (idiag_A5 /=0) call   save_name(f(l1+4 ,m,n,1),idiag_A5)
        if (idiag_D1 /=0) call   save_name(f(l1+0 ,m,n,2),idiag_D1)
        if (idiag_D2 /=0) call   save_name(f(l1+1 ,m,n,2),idiag_D2)
        if (idiag_D3 /=0) call   save_name(f(l1+2 ,m,n,2),idiag_D3)
        if (idiag_D4 /=0) call   save_name(f(l1+3 ,m,n,2),idiag_D4)
        if (idiag_D5 /=0) call   save_name(f(l1+4 ,m,n,2),idiag_D5)
        if (idiag_L1 /=0) call   save_name(f(l1+0 ,m,n,3),idiag_L1)
        if (idiag_L2 /=0) call   save_name(f(l1+1 ,m,n,3),idiag_L2)
        if (idiag_L3 /=0) call   save_name(f(l1+2 ,m,n,3),idiag_L3)
        if (idiag_L4 /=0) call   save_name(f(l1+3 ,m,n,3),idiag_L4)
        if (idiag_L5 /=0) call   save_name(f(l1+4 ,m,n,3),idiag_L5)
        if (idiag_ee10/=0) then
          where (ee1>.1); ee10=1.; elsewhere; ee10=0.; endwhere
          call sum_mn_name(ee10,idiag_ee10)
        endif
        if (idiag_ee50/=0) then
          where (ee1>.5); ee50=1.; elsewhere; ee50=0.; endwhere
          call sum_mn_name(ee50,idiag_ee50)
        endif
        if (idiag_ee90/=0) then
          where (ee1>.9); ee90=1.; elsewhere; ee90=0.; endwhere
          call sum_mn_name(ee90,idiag_ee90)
        endif
        if (idiag_ee99/=0) then
          where (ee1>.99); ee99=1.; elsewhere; ee99=0.; endwhere
          call sum_mn_name(ee99,idiag_ee99)
        endif
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
        idiag_AAm=0; idiag_DDm=0; idiag_LLm=0; idiag_DLm=0; idiag_kC=0
        idiag_A1=0; idiag_A2=0; idiag_A3=0; idiag_A4=0; idiag_A5=0
        idiag_D1=0; idiag_D2=0; idiag_D3=0; idiag_D4=0; idiag_D5=0
        idiag_L1=0; idiag_L2=0; idiag_L3=0; idiag_L4=0; idiag_L5=0
        idiag_eem=0; idiag_ee1m=0; idiag_ee2m=0; idiag_ee3m=0; idiag_ee4m=0
        idiag_ee10=0; idiag_ee50=0; idiag_ee90=0; idiag_ee99=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'eem',idiag_eem)
        call parse_name(iname,cname(iname),cform(iname),'ee1m',idiag_ee1m)
        call parse_name(iname,cname(iname),cform(iname),'ee2m',idiag_ee2m)
        call parse_name(iname,cname(iname),cform(iname),'ee3m',idiag_ee3m)
        call parse_name(iname,cname(iname),cform(iname),'ee4m',idiag_ee4m)
        call parse_name(iname,cname(iname),cform(iname),'ee10',idiag_ee10)
        call parse_name(iname,cname(iname),cform(iname),'ee50',idiag_ee50)
        call parse_name(iname,cname(iname),cform(iname),'ee90',idiag_ee90)
        call parse_name(iname,cname(iname),cform(iname),'ee99',idiag_ee99)
        call parse_name(iname,cname(iname),cform(iname),'AAm',idiag_AAm)
        call parse_name(iname,cname(iname),cform(iname),'DDm',idiag_DDm)
        call parse_name(iname,cname(iname),cform(iname),'LLm',idiag_LLm)
        call parse_name(iname,cname(iname),cform(iname),'DLm',idiag_DLm)
        call parse_name(iname,cname(iname),cform(iname),'kC',idiag_kC)
        call parse_name(iname,cname(iname),cform(iname),'A1',idiag_A1)
        call parse_name(iname,cname(iname),cform(iname),'A2',idiag_A2)
        call parse_name(iname,cname(iname),cform(iname),'A3',idiag_A3)
        call parse_name(iname,cname(iname),cform(iname),'A4',idiag_A4)
        call parse_name(iname,cname(iname),cform(iname),'A5',idiag_A5)
        call parse_name(iname,cname(iname),cform(iname),'D1',idiag_D1)
        call parse_name(iname,cname(iname),cform(iname),'D2',idiag_D2)
        call parse_name(iname,cname(iname),cform(iname),'D3',idiag_D3)
        call parse_name(iname,cname(iname),cform(iname),'D4',idiag_D4)
        call parse_name(iname,cname(iname),cform(iname),'D5',idiag_D5)
        call parse_name(iname,cname(iname),cform(iname),'L1',idiag_L1)
        call parse_name(iname,cname(iname),cform(iname),'L2',idiag_L2)
        call parse_name(iname,cname(iname),cform(iname),'L3',idiag_L3)
        call parse_name(iname,cname(iname),cform(iname),'L4',idiag_L4)
        call parse_name(iname,cname(iname),cform(iname),'L5',idiag_L5)
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
