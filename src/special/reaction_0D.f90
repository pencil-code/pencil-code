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
  integer :: AA0=0
  character(len=50) :: init_qq='racemic'
  namelist /special_init_pars/ &
    init_qq, AA0
!
  ! run parameters
  real :: kp=.1, km=.1, kC=0., kX=.1
  namelist /special_run_pars/ &
    kp, km, kC, kX
!
! other variables (needs to be consistent with reset list below)
!
  integer :: idiag_AAm=0  ! DIAG_DOC: $\langle [A] \rangle$
  integer :: idiag_DDm=0  ! DIAG_DOC: $\langle [D] \rangle$
  integer :: idiag_LLm=0  ! DIAG_DOC: $\langle [L] \rangle$
  integer :: idiag_kC=0   ! DIAG_DOC: $k_C$
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
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  Cannot have kp+km+kX > 1.
!
      if ((kp+km+kX)>1.) call fatal_error('reaction_0D','initialize_special')
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
!  06-oct-2003/tony: coded
!
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      intent(inout) :: f
!
!  initial condition
!
      select case (init_qq)
        case ('nothing'); if (lroot) print*,'init_qq: nothing'
        case ('racemic'); f(:,:,:,iqq)=AA0
        case ('all_equal')
          f(:,:,:,iqq+1)=nx*AA0/3
          f(:,:,:,iqq+2)=nx*AA0/3
          f(:,:,:,iqq)=nx-f(:,:,:,iqq+1)-f(:,:,:,iqq+2)
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
!   06-oct-03/tony: coded
!
      use General, only: random_number_wrapper
      use Diagnostics
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: rrr
      real ::     d1, d2, d3, d4, d5, d6, d7
      real :: r0, r1, r2, r3, r4, r5, r6, r7
      real :: Dm, Lm, Dm_per_proc,Lm_per_proc
      integer, dimension (nx) :: dqA, dqD, dqL
      integer, dimension (nx) :: qA, qD, qL
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
!  Define state vector
!
      qA=f(l1:l2,m,n,1)
      qD=f(l1:l2,m,n,2)
      qL=f(l1:l2,m,n,3)
!
!  For autocatalysis, need to know the mean number of molecules.
!
      if (kC/=0.) then
        Dm_per_proc=sum(qD)
        Lm_per_proc=sum(qL)
        call mpireduce_sum(Dm_per_proc,Dm,1)
        call mpireduce_sum(Lm_per_proc,Lm,1)
        call mpibcast_real(Dm)
        call mpibcast_real(Lm)
      endif
!
!  Recompute kC
!
      if (kC/=0..and.(Dm+Lm)/=0.) kC=(1.-(kp+km+kX))*nx/(Dm+Lm)
!
!  Determine boundary points
!
      d1=.5*kp
      d2=   kp
      d3=.5*km
      d4=   km
      d5=   kC*Dm/nx
      d6=   kC*Lm/nx
      d7=   kX
!
      r0=0.
      r1=   d1
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
print*,'AXEL-: ',r4, r5, r6, r7
print*,'AXEL0: ',rrr
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
      where ((rrr >= r2 .and. rrr < r3) .and. qD > 0)
        dqA= 1
        dqD=-1
        dqL= 0
      endwhere
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
print*,'AXEL1: ',qA
print*,'AXEL2: ',dqA
print*
!
!  diagnostics
!
      if (ldiagnos) then
        if (idiag_AAm/=0) call sum_mn_name(f(l1:l2,m,n,1),idiag_AAm)
        if (idiag_DDm/=0) call sum_mn_name(f(l1:l2,m,n,2),idiag_DDm)
        if (idiag_LLm/=0) call sum_mn_name(f(l1:l2,m,n,3),idiag_LLm)
        if (idiag_kC /=0) call   save_name(kC            ,idiag_kC)
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
!   06-oct-03/tony: coded
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
        idiag_AAm=0; idiag_DDm=0; idiag_LLm=0; idiag_kC=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'AAm',idiag_AAm)
        call parse_name(iname,cname(iname),cform(iname),'DDm',idiag_DDm)
        call parse_name(iname,cname(iname),cform(iname),'LLm',idiag_LLm)
        call parse_name(iname,cname(iname),cform(iname),'kC',idiag_kC)
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
    include 'special_dummies.inc'
!********************************************************************
endmodule Special
