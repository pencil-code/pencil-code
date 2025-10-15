! $Id$
!
!  Lorenz gauge, dphi/dt = -cphi2*divA
!
!  25-feb-07/axel: adapted from nospecial.f90
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
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include '../special.h'
!
!  square of wave speed for gauge field
!
  real :: cphi2
!
  ! input parameters
  real :: cphi=1.,etaphi=0.,ampl=1e-3,kx=1.,ky=0.,kz=0.
  character(len=50) :: init='zero'
  namelist /special_init_pars/ &
    cphi,etaphi,init,ampl,kx,ky,kz
!
  ! run parameters
  namelist /special_run_pars/ &
    cphi,etaphi
!
! Declare any index variables necessary for main or
!
   integer :: iphi=0
!
! other variables (needs to be consistent with reset list below)
!
  integer :: idiag_phim=0       ! DIAG_DOC: $\left<\phi\right>$
  integer :: idiag_phipt=0      ! DIAG_DOC: $\phi(x1,y1,z1)$
  integer :: idiag_phip2=0      ! DIAG_DOC: $\phi(x2,y2,z2)$
  integer :: idiag_phibzm=0     ! DIAG_DOC: $\left<\phi B_z\right>$
  integer :: idiag_phibzmz=0    ! DIAG_DOC: $\left<\phi B_z\right>_{xy}$
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
      call farray_register_pde('phi',iphi)
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
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  Initialize module variables which are parameter dependent
!  wave speed of gauge potential
!
      cphi2=cphi**2
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
      select case (init)
        case ('nothing'); if (lroot) print*,'init_special: nothing'
        case ('zero'); f(:,:,:,iphi)=0.
        case ('sinwave-x'); call sinwave(ampl,f,iphi,kx=kx)
        case ('sinwave-y'); call sinwave(ampl,f,iphi,ky=ky)
        case ('sinwave-z'); call sinwave(ampl,f,iphi,kz=kz)
!
        case default
          !
          !  Catch unknown values
          !
          if (lroot) print*,'init_special: No such value for init: ', trim(init)
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
      if (cphi/=0.) then
        lpenc_requested(i_diva)=.true.
      endif
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
!!      if (headtt) call identify_bcs('ss',iss)
!
!  solve gauge condition
!
      if (lmagnetic) then
        call grad(f,iphi,gphi)
        if (etaphi/=0.) then
          call del2(f,iphi,del2phi)
          df(l1:l2,m,n,iphi)=df(l1:l2,m,n,iphi)-cphi2*p%diva+etaphi*del2phi
        else
          df(l1:l2,m,n,iphi)=df(l1:l2,m,n,iphi)-cphi2*p%diva
        endif
        df(l1:l2,m,n,iax:iaz)=df(l1:l2,m,n,iax:iaz)-gphi
      endif
!
!  diagnostics
!
      if (ldiagnos) then
        phi=f(l1:l2,m,n,iphi)
        if (idiag_phim/=0) call sum_mn_name(phi,idiag_phim)
        if (idiag_phibzm/=0) call sum_mn_name(phi*p%bb(:,3),idiag_phibzm)
!
!  check for point 1
!
        if (lroot.and.m==mpoint.and.n==npoint) then
          if (idiag_phipt/=0) call save_name(phi(lpoint-nghost),idiag_phipt)
        endif
!
!  check for point 2
!
        if (lroot.and.m==mpoint2.and.n==npoint2) then
          if (idiag_phip2/=0) call save_name(phi(lpoint2-nghost),idiag_phip2)
        endif
!
      endif
!
      if (l1davgfirst .or. (ldiagnos .and. ldiagnos_need_zaverages)) then
        if (idiag_phibzmz/=0)   call xysum_mn_name_z(p%bb(:,3),idiag_phibzmz)
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
        idiag_phim=0; idiag_phipt=0; idiag_phip2=0
        idiag_phibzm=0; idiag_phibzmz=0
        cformv=''
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'phim',idiag_phim)
        call parse_name(iname,cname(iname),cform(iname),'phibzm',idiag_phibzm)
        call parse_name(iname,cname(iname),cform(iname),'phipt',idiag_phipt)
        call parse_name(iname,cname(iname),cform(iname),'phip2',idiag_phip2)
      enddo
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'phibzmz',idiag_phibzmz)
      enddo
!
!  check for those quantities for which we want video slices
!
      if (lwrite_slices) then
        where(cnamev=='potturb') cformv='DEFINED'
      endif
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
        call farray_index_append('i_phim',idiag_phim)
        call farray_index_append('i_phibzm',idiag_phibzm)
        call farray_index_append('i_phibzmz',idiag_phibzmz)
        call farray_index_append('i_phipt',idiag_phipt)
        call farray_index_append('i_phip2',idiag_phip2)
        call farray_index_append('iphi',iphi)
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
!  phi
!
        case ('phi'); call assign_slices_vec(slices,f,iphi)
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
    include '../special_dummies.inc'
!********************************************************************
!
endmodule Special
