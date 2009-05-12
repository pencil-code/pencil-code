! $Id$
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
  real :: x0_shear=0.0, Bshear=0.01
  logical :: lshearadvection_as_shift=.false.
  logical :: lmagnetic_stretching=.true.,lrandomx0=.false.
  logical, target :: lcoriolis_force=.true., lcentrifugal_force=.false.
  logical :: lglobal_baroclinic=.false.
!
  include 'shear.h'
!
  namelist /shear_init_pars/ &
      qshear,Sshear,deltay,Omega,lshearadvection_as_shift, &
      lmagnetic_stretching,lrandomx0,x0_shear,lcoriolis_force,lcentrifugal_force
!
  namelist /shear_run_pars/ &
      qshear,Sshear,deltay,Omega,lshearadvection_as_shift, &
      lmagnetic_stretching,lrandomx0,x0_shear, &
      lcoriolis_force,lcentrifugal_force,&
      lglobal_baroclinic,Bshear
!
  integer :: idiag_dtshear=0    ! DIAG_DOC: advec\_shear/cdt
  integer :: idiag_deltay=0     ! DIAG_DOC: deltay
!
  contains
!***********************************************************************
    subroutine register_shear()
!
!  Initialise variables
!
!  2-july-02/nils: coded
!
      use Mpicomm
!
      lshear = .true.
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id$")
!
    endsubroutine register_shear
!***********************************************************************
    subroutine initialize_shear()
!
!  21-nov-02/tony: coded
!  08-jul-04/anders: Sshear calculated whenever qshear /= 0
!
!  calculate shear flow velocity; if qshear is given then Sshear=-qshear*Omega
!  is calculated. Otherwise Sshear keeps its value from the input list.
!
      use SharedVariables, only: put_shared_variable
!
      integer :: ierr=0
!
      if (qshear/=0.0) Sshear=-qshear*Omega
      if (lroot .and. ip<=12) &
          print*,'initialize_shear: Sshear,qshear=',Sshear,qshear
!
!  Share lcoriolis_force and lcentrifugal_force so the Particles module knows
!  whether to apply them or not.
!
      if (lparticles.and.Omega/=0.0.and.(.not.lhydro)) then
        call put_shared_variable('lcoriolis_force',&
            lcoriolis_force,ierr)
        if (ierr/=0) call fatal_error('register_shear',&
            'there was a problem sharing lcoriolis_force')
        call put_shared_variable('lcentrifugal_force',&
            lcentrifugal_force,ierr)
        if (ierr/=0) call fatal_error('register_shear',&
            'there was a problem sharing lcentrifugal_force')
      endif
!
    endsubroutine initialize_shear
!***********************************************************************
    subroutine read_shear_init_pars(unit,iostat)
!
!  read initial shear parameters
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
!  write initial shear parameters
!
      integer, intent(in) :: unit
!
      write(unit,NML=shear_init_pars)
!
    endsubroutine write_shear_init_pars
!***********************************************************************
    subroutine read_shear_run_pars(unit,iostat)
!
!  read run shear parameters
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
!  write run shear parameters
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
      if (lglobal_baroclinic) then 
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
      endif
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
!  Calculates the shear terms, -uy0*df/dy (shearing sheat approximation)
!
!  2-jul-02/nils: coded
!  6-jul-02/axel: runs through all nvar variables; added timestep check
! 16-aug-02/axel: use now Sshear which is calculated in param_io.f90
! 20-aug-02/axel: added magnetic stretching term
!
      use Deriv
      use Diagnostics
      use Fourier, only: fourier_shift_yz_y
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (ny,nz) :: f_tmp_yz
      real, dimension (nx) :: uy0,dfdy
      integer :: j,k
!
      intent(in)  :: f
!
!  print identifier
!
      if (headtt.or.ldebug) print*, 'shearing: Sshear,qshear=', Sshear, qshear
!
!  add shear term, -uy0*df/dy, for all variables
!
      uy0=Sshear*(x(l1:l2)-x0_shear)
!
!  Advection of all variables by shear flow.
!
      if (.not. lshearadvection_as_shift) then
        do j=1,nvar
          call der(f,j,dfdy,2)
          df(l1:l2,m,n,j)=df(l1:l2,m,n,j)-uy0*dfdy
        enddo
      endif
!
! Taking care of the fact that the Coriolis force changes when
! we have got shear. The rest of the Coriolis force is calculated
! in hydro.
!
      if (lhydro) df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-Sshear*p%uu(:,1)
!
!  Loop over dust species
!
      if (ldustvelocity) then
        do k=1,ndustspec
!
!  Correct Coriolis force term for all dust species
!
           df(l1:l2,m,n,iudy(k)) = df(l1:l2,m,n,iudy(k)) &
              - Sshear*f(l1:l2,m,n,iudx(k))
!
!  End loop over dust species
!
        enddo
      endif
!
!  Magnetic stretching term (can be turned off for debugging purposes).
!
      if (lmagnetic .and. lmagnetic_stretching) then
        df(l1:l2,m,n,iax)=df(l1:l2,m,n,iax)-Sshear*p%aa(:,2)
      endif
!
!  Testfield stretching term
!  Loop through all the dax/dt equations and add -S*ay contribution
!
      if (ltestfield) then
        do j=iaatest,iaxtestpq,3
          df(l1:l2,m,n,j)=df(l1:l2,m,n,j)-Sshear*f(l1:l2,m,n,j+1)
        enddo
      endif
!
!  Testflow stretching term
!  Loop through all the duy/dt equations and add -S*ux contribution
!
      if (ltestflow) then
        do j=iuutest+1,iuxtestpq,3
          df(l1:l2,m,n,j)=df(l1:l2,m,n,j)-Sshear*f(l1:l2,m,n,j-1)
        enddo
      endif
!
!  Meanfield stretching term
!  Loop through all the dax/dt equations and add -S*ay contribution
!
      if (iam/=0) then
        df(l1:l2,m,n,iamx)=df(l1:l2,m,n,iamx)-Sshear*f(l1:l2,m,n,iamy)
      endif
!
!  Global baroclinic term 
!
      if (lglobal_baroclinic) call global_baroclinic(f,df,p)
!
!  Take shear into account for calculating time step
!
      if (lfirst.and.ldt.and.(.not.lshearadvection_as_shift)) &
          advec_shear=abs(uy0*dy_1(m))
!
!  Calculate shearing related diagnostics
!
      if (ldiagnos) then
        if (idiag_dtshear/=0) &
            call max_mn_name(advec_shear/cdt,idiag_dtshear,l_dt=.true.)
      endif
!
    endsubroutine shearing
!***********************************************************************
    subroutine advance_shear(f,df,dt_shear)
!
!  Advance shear distance, deltay, using dt. Using t instead introduces
!  significant errors when nt = t/dt exceeds ~100,000 steps.
!  This formulation works also when Sshear is changed during the run.
!
!  18-aug-02/axel: incorporated from nompicomm.f90
!
      use Diagnostics
      use Fourier, only: fourier_shift_y
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real :: dt_shear
!
      real, dimension (nx,ny,nz) :: tmp
      real, dimension (nx) :: uy0
      integer :: l, ivar
!
!  Must currently use lshearadvection_as_shift=T when Sshear is positive.
!
      if (Sshear>0. .and. .not. lshearadvection_as_shift) then
        if (lroot) print*, 'Note: must use lshearadvection_as_shift=T with positive values of Sshear'
        call stop_it('')
      endif
!
!  Make sure deltay is in the range 0 <= deltay < Ly (assuming Sshear<0).
!
      deltay=deltay-Sshear*Lx*dt_shear
      deltay=deltay-int(deltay/Ly)*Ly
!
!  Solve for advection by shear motion by shifting all variables and their
!  time derivative (following Gammie 2001). Removes time-step constraint
!  from shear motion.
!
      if (lshearadvection_as_shift) then
        uy0=Sshear*(x(l1:l2)-x0_shear)
        do ivar=1,mvar
          tmp=f(l1:l2,m1:m2,n1:n2,ivar)
          call fourier_shift_y(tmp,uy0*dt_shear)
          f(l1:l2,m1:m2,n1:n2,ivar)=tmp
        enddo
        if (itsub/=itorder) then
          do ivar=1,mvar
            tmp=df(l1:l2,m1:m2,n1:n2,ivar)
            call fourier_shift_y(tmp,uy0*dt_shear)
            df(l1:l2,m1:m2,n1:n2,ivar)=tmp
          enddo
        endif
      endif
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
      if (lshearadvection_as_shift) then
        call fourier_shift_ghostzones(f,ivar1,ivar2)
      else
        call initiate_shearing(f,ivar1,ivar2)
        if (nprocy>1 .or. (.not. lmpicomm)) call finalize_shearing(f,ivar1,ivar2)
      endif
!
    endsubroutine boundcond_shear
!***********************************************************************
    subroutine fourier_shift_ghostzones(f,ivar1,ivar2)
!
!  Shearing boundary conditions by Fourier interpolation.
!
!  02-oct-07/anders: coded
!
      use Fourier, only: fourier_shift_yz_y
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: ivar1, ivar2
!
      real, dimension (ny,nz) :: f_tmp_yz
      integer :: i, ivar
!
      if (nxgrid/=1) then
        f(l2+1:mx,m1:m2,n1:n2,ivar1:ivar2)=f(l1:l1+2,m1:m2,n1:n2,ivar1:ivar2)
        f( 1:l1-1,m1:m2,n1:n2,ivar1:ivar2)=f(l2-2:l2,m1:m2,n1:n2,ivar1:ivar2)
      endif
!
      if (nygrid/=1) then
        do ivar=ivar1,ivar2
          do i=1,3
            f_tmp_yz=f(l1-i,m1:m2,n1:n2,ivar)
            call fourier_shift_yz_y(f_tmp_yz,+deltay)
            f(l1-i,m1:m2,n1:n2,ivar)=f_tmp_yz
            f_tmp_yz=f(l2+i,m1:m2,n1:n2,ivar)
            call fourier_shift_yz_y(f_tmp_yz,-deltay)
            f(l2+i,m1:m2,n1:n2,ivar)=f_tmp_yz
          enddo
        enddo
      endif
!
    endsubroutine fourier_shift_ghostzones
!***********************************************************************
    subroutine global_baroclinic(f,df,p)
!
      use EquationOfState, only: rho0,cs20,gamma11,gamma1
      use Messages, only: fatal_error
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: rhs
      type (pencil_case) :: p     
      real :: P0
!
!  x-momentum
!      
      P0=rho0*cs20*gamma11
      df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)-Bshear*P0/rho0*(p%rho1*rho0-1.)
!
!  Right hand side on the energy equation
!
      rhs=Bshear*P0*p%uu(:,1)/gamma1 
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
      endsubroutine global_baroclinic
!***********************************************************************
    subroutine rprint_shear(lreset,lwrite)
!
!  reads and registers print parameters relevant to shearing
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
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_dtshear=0
        idiag_deltay=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'dtshear',idiag_dtshear)
        call parse_name(iname,cname(iname),cform(iname),'deltay',idiag_deltay)
      enddo
!
!  write column where which shear variable is stored
!
      if (lwr) then
        write(3,*) 'i_dtshear=',idiag_dtshear
        write(3,*) 'i_deltay=',idiag_deltay
      endif
!
    endsubroutine rprint_shear
!***********************************************************************
  endmodule Shear
