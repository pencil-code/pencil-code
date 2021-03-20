! $Id$
!
!  Electric field, dE/dt = curlB, currently only for the special case
!  of no fluid induction.
!
!  25-feb-07/axel: adapted from nospecial.f90
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of disp_current_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 3
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED EEE2; EEE(3)
!***************************************************************
!
module disp_current
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
  ! input parameters
  real :: ampl=1e-3, constf=1. !,etaphi=0.,kx=1., ky=0., kz=0.
  real :: ampl_Ex=0.0, ampl_Ey=0.0, ampl_Ez=0.0
  real :: kx_Ex=0.0, kx_Ey=0.0, kx_Ez=0.0
  real :: ky_Ex=0.0, ky_Ey=0.0, ky_Ez=0.0
  real :: kz_Ex=0.0, kz_Ey=0.0, kz_Ez=0.0
  real :: phase_Ex=0.0, phase_Ey=0.0, phase_Ez=0.0
  character(len=50) :: initEE='zero'
  namelist /disp_current_init_pars/ &
    initEE, constf, &
    ampl_Ex, ampl_Ey, ampl_Ez, &
    kx_Ex, kx_Ey, kx_Ez, &
    ky_Ex, ky_Ey, ky_Ez, &
    kz_Ex, kz_Ey, kz_Ez, &
    phase_Ex, phase_Ey, phase_Ez
!
  ! run parameters
  namelist /disp_current_run_pars/ &
    constf
!
! Declare any index variables necessary for main or
!
  real :: c_light2
!
! other variables (needs to be consistent with reset list below)
!
  integer :: idiag_EErms=0       ! DIAG_DOC: $\left<\Ev^2\right>^{1/2}$
  integer :: idiag_EEmax=0       ! DIAG_DOC: $\max(|\Ev|)$
!
! xy averaged diagnostics given in xyaver.in
!
  integer :: idiag_EExmz=0       ! XYAVG_DOC: $\left<{\cal E}_x\right>_{xy}$
  integer :: idiag_EEymz=0       ! XYAVG_DOC: $\left<{\cal E}_y\right>_{xy}$
  integer :: idiag_EEzmz=0       ! XYAVG_DOC: $\left<{\cal E}_z\right>_{xy}$
!
  contains
!
!***********************************************************************
    subroutine register_special()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  18-mar-21/axel: coded Faraday displacement current
!
      use FArrayManager
!
!  It would have been more consistent to call the indices to the
!  three components iex, iey, and iez, but the names
!  iEEx, iEEy, and iEEz were already in use.
!
      call farray_register_pde('ee',iEE,vector=3)
      iEEx=iEE; iEEy=iEEx+1; iEEz=iEEx+2
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  20-mar-21/axel: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  Initialize module variables which are parameter dependent
!  If one really wants to work with c_light /= 1,
!  then one needs to override this.
!
      if (c_light/=1.) call fatal_error('disp_current', "use unit_system='set'")
      c_light2=c_light**2
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
      use Initcond
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      select case (initEE)
        case ('nothing'); if (lroot) print*,'initEE: nothing'
        case ('zero'); f(:,:,:,iEEx:iEEz)=0.
        case ('coswave-phase')
          call coswave_phase(f,iEEx,ampl_Ex,kx_Ex,ky_Ex,kz_Ex,phase_Ex)
          call coswave_phase(f,iEEy,ampl_Ey,kx_Ey,ky_Ey,kz_Ey,phase_Ey)
          call coswave_phase(f,iEEz,ampl_Ez,kx_Ez,ky_Ez,kz_Ez,phase_Ez)
!
        case default
          !
          !  Catch unknown values
          !
          if (lroot) print*,'initEE: No such value for initEE: ', trim(initEE)
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
!  25-feb-07/axel: adapted
!
      lpenc_requested(i_aa)=.true.
      lpenc_requested(i_EEE)=.true.
      lpenc_requested(i_curlB)=.true.
      lpenc_requested(i_del2a)=.true.

      if (idiag_EErms/=0 .or. idiag_EEmax/=0) lpenc_diagnos(i_EEE2)=.true.
      if (idiag_EExmz/=0 .or. idiag_EEymz/=0 .or. idiag_EEzmz/=0 ) lpenc_diagnos(i_EEE)=.true.
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
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!      logical, dimension(:),              intent(in)   :: lpenc_loc
!
      intent(in) :: f
      intent(inout) :: p
!
! EEE
      p%EEE=f(l1:l2,m,n,iEEx:iEEz)
! EEE2
      call dot2_mn(p%EEE,p%EEE2)
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
!   18-mar-21/axel: coded Faraday displacement current
!
      use Diagnostics
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx,3) :: gphi
      real, dimension (nx) :: phi,del2phi
!
      intent(in) :: f,p
      intent(inout) :: df
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dSPECIAL_dt'
      if (headtt) call identify_bcs('EE',iEE)
!
!  solve: dE/dt = curlB - (const/t^2)*aa
!  Calculate curlB as -del2a, because curlB leads to instability.
!
      if (lmagnetic) then
        if (t==0.) call fatal_error('disp_current: dspecial_dt', 't=0 not allowed')
        df(l1:l2,m,n,iEEx:iEEz)=df(l1:l2,m,n,iEEx:iEEz)-c_light2*p%del2a-constf/t**2*p%aa
        df(l1:l2,m,n,iax:iaz)=df(l1:l2,m,n,iax:iaz)-p%EEE
      endif
!
!  timestep constraint
!
      if (lfirst.and.ldt) advec_cs2=max(advec_cs2,c_light2*dxyz_2)
!
!  diagnostics
!
      if (ldiagnos) then
        call sum_mn_name(p%EEE2,idiag_EErms,lsqrt=.true.)
        call max_mn_name(p%EEE2,idiag_EEmax,lsqrt=.true.)
!
        call xysum_mn_name_z(p%EEE(:,1),idiag_EExmz)
        call xysum_mn_name_z(p%EEE(:,2),idiag_EEymz)
        call xysum_mn_name_z(p%EEE(:,3),idiag_EEzmz)
!
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
      read(parallel_unit, NML=disp_current_init_pars, IOSTAT=iostat)
!
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=disp_current_init_pars)
!
    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=disp_current_run_pars, IOSTAT=iostat)
!
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=disp_current_run_pars)
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
!  define counters
!
      integer :: iname,inamez
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
        idiag_EErms=0; idiag_EEmax=0
        cformv=''
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'EErms',idiag_EErms)
        call parse_name(iname,cname(iname),cform(iname),'EEmax',idiag_EEmax)
      enddo
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'EExmz',idiag_EExmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'EEymz',idiag_EEymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'EEzmz',idiag_EEzmz)
      enddo
!
!  check for those quantities for which we want video slices
!
      if (lwrite_slices) then
        where(cnamev=='EE') cformv='DEFINED'
      endif
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
      endif
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine get_slices_special(f,slices)
!
!  Write slices for animation of electric potential
!
!  26-feb-07/axel: adapted from gross_pitaevskii
!
      use Slices_methods, only: assign_slices_vec
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      type (slice_data) :: slices
!
      integer :: inamev
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  Electric field.
!
      case ('EE'); call assign_slices_vec(slices,f,iEE)
!
      endselect
!
    endsubroutine get_slices_special
!***********************************************************************
!
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../disp_current_dummies.inc'
!********************************************************************
!
endmodule disp_current
