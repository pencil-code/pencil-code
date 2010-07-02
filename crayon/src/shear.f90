! $Id: shear.f90 13341 2010-02-23 13:20:53Z AxelBrandenburg $
!
!  This modules deals with all aspects of shear; if no
!  shear is invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  shear relevant subroutines listed in here.
!  Shear can either be given relative to Omega (using qshear),
!  or in absolute fashion via the parameters Sshear.
!
module Shear
!
  use Cdata
  use Messages
  use Sub
!
  implicit none
!
  real :: x0_shear=0.0
  logical :: lshearadvection_as_shift=.false.
  logical :: lmagnetic_stretching=.true.,lrandomx0=.false.
!
  include 'shear.h'
!
  namelist /shear_init_pars/ &
      qshear,Sshear,deltay,Omega,lshearadvection_as_shift, &
      lmagnetic_stretching,lrandomx0,x0_shear
!
  namelist /shear_run_pars/ &
      qshear,Sshear,deltay,Omega,lshearadvection_as_shift, &
      lmagnetic_stretching,lrandomx0,x0_shear
!
  integer :: idiag_dtshear=0    ! DIAG_DOC: advec\_shear/cdt
  integer :: idiag_deltay=0     ! DIAG_DOC: deltay
!
  contains
!***********************************************************************
    subroutine register_shear()
!
!  Initialise variables.
!
!  2-july-02/nils: coded
!
      lshear=.true.
!
!  Identify version number.
!
      if (lroot) call svn_id( &
           "$Id: shear.f90 13341 2010-02-23 13:20:53Z AxelBrandenburg $")
!
    endsubroutine register_shear
!***********************************************************************
    subroutine initialize_shear()
!
!  21-nov-02/tony: coded
!  08-jul-04/anders: Sshear calculated whenever qshear /= 0
!
!  Calculate shear flow velocity; if qshear is given then Sshear=-qshear*Omega
!  is calculated. Otherwise Sshear keeps its value from the input list.
!
      use SharedVariables, only: put_shared_variable
!
      if (qshear/=0.0) Sshear=-qshear*Omega
      if (lroot .and. ip<=12) &
          print*,'initialize_shear: Sshear,qshear=',Sshear,qshear
!
    endsubroutine initialize_shear
!***********************************************************************
    subroutine read_shear_init_pars(unit,iostat)
!
!  Read initial shear parameters.
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=shear_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=shear_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_shear_init_pars
!***********************************************************************
    subroutine write_shear_init_pars(unit)
!
!  Write initial shear parameters.
!
      integer, intent(in) :: unit
!
      write(unit,NML=shear_init_pars)
!
    endsubroutine write_shear_init_pars
!***********************************************************************
    subroutine read_shear_run_pars(unit,iostat)
!
!  Read run shear parameters.
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=shear_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=shear_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_shear_run_pars
!***********************************************************************
    subroutine write_shear_run_pars(unit)
!
!  Write run shear parameters.
!
      integer, intent(in) :: unit
!
      write(unit,NML=shear_run_pars)
!
    endsubroutine write_shear_run_pars
!***********************************************************************
    subroutine shear_before_boundary(f)
!
!  Actions to take before boundary conditions are set.
!
!   1-may-08/anders: coded
!
      use General
      use Mpicomm
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  Possible to shear around a random position in x, to let all points
!  be subjected to shear in a statistically equal way.
!
      if (itsub==1) then
        if (lrandomx0) then
          if (lroot) then
            call random_number_wrapper(x0_shear)
            x0_shear=x0_shear*Lxyz(1)+xyz0(1)
          endif
          call mpibcast_real(x0_shear,1,0)
        endif
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine shear_before_boundary
!***********************************************************************
    subroutine pencil_criteria_shear()
!
!  All pencils that the Shear module depends on are specified here.
!
!  01-may-09/wlad: coded
!
      if (lhydro)    lpenc_requested(i_uu)=.true.
      if (lmagnetic) lpenc_requested(i_aa)=.true.
!
    endsubroutine pencil_criteria_shear
!***********************************************************************
    subroutine pencil_interdep_shear(lpencil_in)
!
!  Interdependency among pencils from the Shear module is specified here.
!
!  01-may-09/wlad: coded
!
      use Sub, only: keep_compiler_quiet
      logical, dimension(npencils) :: lpencil_in
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_shear
!***********************************************************************
    subroutine calc_pencils_shear(f,p)
!
!  Calculate Shear pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  01-may-09/wlad: coded
!
      use Sub, only: keep_compiler_quiet
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f,p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_shear
!***********************************************************************
    subroutine shearing(f,df,p)
!
!  Calculates the shear terms -uy0*df/dy (shearing sheat approximation).
!
!  2-jul-02/nils: coded
!  6-jul-02/axel: runs through all nvar variables; added timestep check
! 16-aug-02/axel: use now Sshear which is calculated in param_io.f90
! 20-aug-02/axel: added magnetic stretching term
!
      use Deriv
      use Diagnostics
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: uy0,dfdy
      integer :: j
!
      intent(in)  :: f
!
!  Print identifier.
!
      if (headtt.or.ldebug) print*, 'shearing: Sshear,qshear=', Sshear, qshear
!
!  Add shear term, -uy0*df/dy, for all variables.
!
      uy0=Sshear*(x(l1:l2)-x0_shear)
!
!  Advection of all variables by shear flow.
!
      do j=1,nvar
        call der(f,j,dfdy,2)
        df(l1:l2,m,n,j)=df(l1:l2,m,n,j)-uy0*dfdy
      enddo
!
!  Advection of background velocity profile. Appears like a correction
!  to the Coriolis force, but is actually not related to the Coriolis
!  force.
!
      if (lhydro) df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-Sshear*p%uu(:,1)
!
!  Magnetic stretching term (can be turned off for debugging purposes).
!
      if (lmagnetic .and. lmagnetic_stretching) then
        df(l1:l2,m,n,iax)=df(l1:l2,m,n,iax)-Sshear*p%aa(:,2)
      endif
!
!  Take shear into account for calculating time step.
!
      if (lfirst.and.ldt) advec_shear=abs(uy0*dy_1(m))
!
!  Calculate shearing related diagnostics.
!
      if (ldiagnos) then
        if (idiag_dtshear/=0) &
            call max_mn_name(advec_shear/cdt,idiag_dtshear,l_dt=.true.)
      endif
!
    endsubroutine shearing
!***********************************************************************
    subroutine advance_shear(dt_shear)
!
!  Advance shear distance, deltay, using dt. Using t instead introduces
!  significant errors when nt = t/dt exceeds ~100,000 steps.
!  This formulation works also when Sshear is changed during the run.
!
!  18-aug-02/axel: incorporated from nompicomm.f90
!
      use Diagnostics
      use Mpicomm, only: stop_it
!
      real :: dt_shear
!
!  Must currently use lshearadvection_as_shift=T when Sshear is positive.
!
      if (Sshear>0. .and. ncpus/=1 .and. headt) then
        print*
        print*, 'NOTE: for Sshear > 0, MPI is not completely correct.'
        print*, 'It is better to use lshearadvection_as_shift=T.'
        print*
      endif
!
!  Make sure deltay is in the range 0 <= deltay < Ly (assuming Sshear<0).
!
      deltay=deltay-Sshear*Lx*dt_shear
      deltay=deltay-int(deltay/Ly)*Ly
!
!  Print identifier.
!
      if (headtt.or.ldebug) print*, 'advance_shear: deltay=',deltay
!
!  Calculate shearing related diagnostics
!
      if (ldiagnos) then
        if (idiag_deltay/=0) &
            call save_name(deltay,idiag_deltay)
      endif
!
    endsubroutine advance_shear
!***********************************************************************
    subroutine boundcond_shear(f,ivar1,ivar2)
!
!  Shearing boundary conditions, called from the Boundconds module.
!
!  02-oct-07/anders: coded
!
      use Mpicomm, only: initiate_shearing, finalize_shearing
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: ivar1, ivar2
!
      if (ip<12.and.headtt) print*, &
          'boundconds_x: use shearing sheet boundary condition'
!
      call initiate_shearing(f,ivar1,ivar2)
      call finalize_shearing(f,ivar1,ivar2)
!
    endsubroutine boundcond_shear
!***********************************************************************
    subroutine rprint_shear(lreset,lwrite)
!
!  Reads and registers print parameters relevant to shearing.
!
!   2-jul-04/tobi: adapted from entropy
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
!  Reset everything in case of reset.
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_dtshear=0
        idiag_deltay=0
      endif
!
!  iname runs through all possible names that may be listed in print.in.
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'dtshear',idiag_dtshear)
        call parse_name(iname,cname(iname),cform(iname),'deltay',idiag_deltay)
      enddo
!
!  Write column where which shear variable is stored.
!
      if (lwr) then
        write(3,*) 'i_dtshear=',idiag_dtshear
        write(3,*) 'i_deltay=',idiag_deltay
      endif
!
    endsubroutine rprint_shear
!***********************************************************************
endmodule Shear
