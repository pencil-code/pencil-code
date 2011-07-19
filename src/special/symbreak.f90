
! $Id$
!
!  Solve for a set of two ODEs, used to test time step
!
!  19-apr-09/axel: adapted from nospecial.f90
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 4
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Special

  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet

  implicit none

  include '../special.h'

!
! Declare index of variables
! 
   integer :: ispecial=0,ispecial1=0,ispecial2=0, & 
              ispecial3=0,ispecial4=0

  ! input parameters
  real :: beta_real,beta_imag,mu_real,mu_imag,gam,gam_imag,rho,Lreal_ini,Limag_ini,Rreal_ini,&
      Rimag_ini
  real, dimension (ninit) :: ampl=0.
  character (len=labellen), dimension(ninit) :: init='nothing'
  namelist /special_init_pars/ &
    init,ampl,beta_real,beta_imag,gam,gam_imag,mu_real,mu_imag,rho,Lreal_ini,&
    Limag_ini,Rreal_ini,Rimag_ini

  ! run parameters
  namelist /special_run_pars/ &
    beta_real,beta_imag,gam,gam_imag,rho,mu_real,mu_imag
!
! other variables (needs to be consistent with reset list below)
!
  integer :: idiag_Lreal=0,idiag_Limag=0,idiag_Rreal=0,idiag_Rimag
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
      use FArrayManager
!
!  Identify CVS/SVN version information.
!
      if (lroot) call svn_id( &
           "$Id$")
!
!      call farray_register_pde('spec_3vec',ispecial,vector=3)
      call farray_register_pde('spec_4vec',ispecial,vector=4)
!
      ispecial1=ispecial
      ispecial2=ispecial+1
      ispecial3=ispecial+2
      ispecial4=ispecial+3
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f,lstarting)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
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
      use Initcond
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
!
      intent(inout) :: f
!
!  initial condition
!
      do j=1,ninit
        select case (init(j))
        case ('zero'); f(:,:,:,ispecial1:ispecial4)=0.
        case ('nothing'); if (lroot) print*,'init_special: nothing'
        case ('gaussian-noise'); call gaunoise(ampl(j),f,ispecial1,ispecial3)
        case ('set')
          f(:,:,:,ispecial1)=f(:,:,:,ispecial1)+Lreal_ini
          f(:,:,:,ispecial2)=f(:,:,:,ispecial2)+Limag_ini
          f(:,:,:,ispecial3)=f(:,:,:,ispecial3)+Rreal_ini
          f(:,:,:,ispecial4)=f(:,:,:,ispecial4)+Rimag_ini
        case default
         if (lroot) print*,'init_special: No such value: ',trim(init(j))
          call stop_it("")
        endselect
      enddo
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
      use Diagnostics
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: Lreal,Limag,Rreal,Rimag,ddphi_phi4 
      type (pencil_case) :: p

!
      intent(in) :: f,p
      intent(inout) :: df
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dSPECIAL_dt'
!
!  Solve the equations for symmetry breaking
!
        Lreal=f(l1:l2,m,n,ispecial1)
        Limag=f(l1:l2,m,n,ispecial2)
        Rreal=f(l1:l2,m,n,ispecial3)
        Rimag=f(l1:l2,m,n,ispecial4)
!
!        write(*,*) 'DM',ispecial1,gam*Lreal,Lreal
        df(l1:l2,m,n,ispecial1)=df(l1:l2,m,n,ispecial1)&
             +gam*Lreal-gam_imag*Limag& 
!             -(beta_real*Lreal-beta_imag*Limag)*((Lreal*Lreal+Limag*Limag)&
!                  +(Rreal*Rreal+Rimag*Rimag)) 
             -(beta_real*Lreal-beta_imag*Limag)*(Rreal*Rreal+Rimag*Rimag) &
             -(mu_real*Lreal-mu_imag*Limag)*(Lreal*Lreal+Limag*Limag) 
        df(l1:l2,m,n,ispecial2)=df(l1:l2,m,n,ispecial2)&
             +gam*Limag+gam_imag*Lreal&
!             -(beta_real*Limag+beta_imag*Lreal)*((Lreal*Lreal+Limag*Limag)&
!                  + (Rreal*Rreal+Rimag*Rimag))
             -(beta_real*Limag+beta_imag*Lreal)*(Rreal*Rreal+Rimag*Rimag) &
             -(mu_real*Limag+mu_imag*Lreal)*(Lreal*Lreal+Limag*Limag)
        df(l1:l2,m,n,ispecial3)=df(l1:l2,m,n,ispecial3)&
              +gam*Rreal-gam_imag*Rimag&
!             -(beta_real*Rreal-beta_imag*Rimag)*((Lreal*Lreal+Limag*Limag)&
!                  + (Rreal*Rreal+Rimag*Rimag))
             -(beta_real*Rreal-beta_imag*Rimag)*(Lreal*Lreal+Limag*Limag) &
             -(mu_real*Rreal-mu_imag*Rimag)*(Rreal*Rreal+Rimag*Rimag) 
        df(l1:l2,m,n,ispecial4)=df(l1:l2,m,n,ispecial4)&
              +gam*Rimag+gam_imag*Rreal&
!             -(beta_real*Rimag+beta_imag*Rreal)*((Lreal*Lreal+Limag*Limag)&
!                  + (Rreal*Rreal+Rimag*Rimag))
             -(beta_real*Rimag+beta_imag*Rreal)*(Lreal*Lreal+Limag*Limag) &
             -(mu_real*Rimag+mu_imag*Rreal)*(Rreal*Rreal+Rimag*Rimag)
!
!  diagnostics
!
      if (ldiagnos) then
        if (lroot.and.m==mpoint.and.n==npoint) then
          if (idiag_Lreal/=0) &
             call save_name(Lreal(lpoint-nghost),idiag_Lreal)
          if (idiag_Limag/=0) &
             call save_name(Limag(lpoint-nghost),idiag_Limag)
          if (idiag_Rreal/=0) & 
             call save_name(Rreal(lpoint-nghost),idiag_Rreal)
          if (idiag_Rimag/=0) & 
             call save_name(Rimag(lpoint-nghost),idiag_Rimag)
        endif
      endif
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
!  read namelist
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
!
      integer, intent(in) :: unit
!
!  write name list
!
      write(unit,NML=special_init_pars)
!
    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
!  write name list
!
      if (present(iostat)) then
        read(unit,NML=special_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=special_run_pars,ERR=99)
      endif
!
99  endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
      integer, intent(in) :: unit
!
!  write name list
!
      write(unit,NML=special_run_pars)
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
        idiag_Lreal=0; idiag_Limag=0 
        idiag_Rreal=0; idiag_Rimag=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'Lreal',idiag_Lreal)
        call parse_name(iname,cname(iname),cform(iname),'Limag',idiag_Limag)
        call parse_name(iname,cname(iname),cform(iname),'Rreal',idiag_Rreal)
        call parse_name(iname,cname(iname),cform(iname),'Rimag',idiag_Rimag)
      enddo
!
!  write column where which variable is stored
!
      if (lwr) then
        write(3,*) 'i_Lreal=',idiag_Lreal
        write(3,*) 'i_Limag=',idiag_Limag
        write(3,*) 'i_Rreal=',idiag_Rreal
        write(3,*) 'i_Rimag=',idiag_Rimag
      endif
!
    endsubroutine rprint_special
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
!   continuity equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p

!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!
!!
!!  df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + SOME NEW TERM
!!
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_density
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the
!   momentum equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p

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
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine special_calc_magnetic(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the
!   induction equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p

!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!
!!
!!  df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) + SOME NEW TERM
!!  df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) + SOME NEW TERM
!!  df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) + SOME NEW TERM
!!
      call keep_compiler_quiet(f,df)
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
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p

!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!
!!
!!  df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) + SOME NEW TERM
!!
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_entropy
!***********************************************************************
    subroutine special_after_timestep(f,df,dt_)
!
!   Possibility to modify the f and df after df is updated
!   Used for the fargo shift, for instance.
!
!   27-nov-08/wlad: coded
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real :: dt_
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(dt_)
!
    endsubroutine  special_after_timestep
!********************************************************************
    subroutine special_calc_particles(fp)
!
!   Called before the loop, in case some particle value is needed 
!   for the special density/hydro/magnetic/entropy
!
!   20-nov-08/wlad: coded
!
      real, dimension (:,:), intent(in) :: fp
!
      call keep_compiler_quiet(fp)
!
    endsubroutine special_calc_particles
!***********************************************************************
    subroutine special_calc_particles_nbody(fsp)
!
!   Called before the loop, in case some massive particles value 
!   is needed for the special density/hydro/magnetic/entropy
!
!   20-nov-08/wlad: coded
!
      real, dimension (:,:), intent(in) :: fsp
!
      call keep_compiler_quiet(fsp)
!
    endsubroutine special_calc_particles_nbody
!***********************************************************************
    subroutine special_boundconds(f,bc)
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      type (boundary_condition) :: bc
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(bc)
!
    endsubroutine special_boundconds
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
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine special_before_boundary
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

