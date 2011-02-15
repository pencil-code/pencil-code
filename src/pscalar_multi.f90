! $Id$
!
!  This module solves multiple passive scalar advection equations.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lpscalar = .true.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
! PENCILS PROVIDED cc(npscalar); gcc(3,npscalar); ugcc(npscalar)
!
!***************************************************************
module Pscalar
!
  use Cdata
  use Cparam
  use Messages
!
  implicit none
!
  include 'pscalar.h'
!
!  Init parameters
!
  real, dimension(3) :: u0_advec = 0.
  logical :: lper_unit_volume = .false.
!
  namelist /pscalar_init_pars/ u0_advec, lper_unit_volume
!
!  Run parameters
!
  real :: diffcc_shock = 0.
!
  namelist /pscalar_run_pars/ u0_advec, diffcc_shock
!
!  Diagnostic variables
!
  integer :: idiag_ccmax = 0, idiag_ccmin = 0, idiag_ccm = 0
  integer :: idiag_rhoccm = 0
!
!  Module specific variables
!
  logical :: lconst_advection = .false.
  logical :: lshock_diffusion = .false.
!
!  Dummy variables needed by other modules
!
  integer :: idiag_gcc2m = 0, idiag_cc2m = 0
  real :: rhoccm = 0., cc2m = 0., gcc2m = 0.
!
  contains
!***********************************************************************
    subroutine register_pscalar()
!
!  Initialise variables which should know that we solve for passive
!  scalar: icc; increase nvar accordingly
!
!  09-feb-11/ccyang: coded
!
      use FArrayManager
!
!  This module only implements nolog version.
!
      lpscalar_nolog=.true.
!
      call farray_register_pde('cc', icc, vector=npscalar)
      ilncc=0        ! needed for idl
!
      if (lroot) print*, 'Number of passive scalars = ', npscalar
!
!  Identify version number.
!
      if (lroot) call svn_id("$Id")
!
    endsubroutine register_pscalar
!***********************************************************************
    subroutine initialize_pscalar(f)
!
!  Perform any necessary post-parameter read initialization.
!
!  02-feb-11/ccyang: coded
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
!  Check if there exists constant background advection.
!
      if (any(u0_advec /= 0.)) then
        lconst_advection = .true.
        if (lroot) print*, 'initialize_pscalar: constant background advection = ', u0_advec
      endif
!
!  Shock Diffusion
!
      if (diffcc_shock /= 0.) then
        lshock_diffusion = .true.
        if (lroot) print*, 'initialize_pscalar: shock diffusion, diffcc_shock = ', diffcc_shock
      endif
!
!      call init_lncc(f)
!
    endsubroutine initialize_pscalar
!***********************************************************************
    subroutine init_lncc(f)
!
!  Initialise passive scalar field.
!  It does not make sense to initialize multiple scalar fields here; use the
!  initial_conditions facility instead.
!
!  10-feb-11/ccyang: coded
!
      use InitialCondition, only: initial_condition_lncc
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      call initial_condition_lncc(f)
!
    endsubroutine init_lncc
!***********************************************************************
    subroutine pencil_criteria_pscalar()
!
!  All pencils that the Pscalar module depends on are specified here.
!
!  04-feb-11/ccyang: coded
!
!  Dynamics
!
      lpenc_requested(i_ugcc) = .true.
!
      if (lper_unit_volume) then
        lpenc_requested(i_cc) = .true.
        lpenc_requested(i_divu) = .true.
      endif
!
      if (lshock_diffusion .and. lper_unit_volume) then
        lpenc_requested(i_gcc) = .true.
        lpenc_requested(i_shock) = .true.
        lpenc_requested(i_gshock) = .true.
      endif
!
!  Diagnostics
!
      if (idiag_ccmax /= 0 .or. idiag_ccmin /= 0 .or. idiag_ccm /= 0) lpenc_diagnos(i_cc) = .true.
      if (idiag_rhoccm /= 0) then
        lpenc_diagnos(i_rho) = .true.
        lpenc_diagnos(i_cc) = .true.
      endif
!
    endsubroutine pencil_criteria_pscalar
!***********************************************************************
    subroutine pencil_interdep_pscalar(lpencil_in)
!
!  Interdependency among pencils provided by the Pscalar module
!  is specified here.
!
!  31-jan-11/ccyang: coded
!
      logical, dimension(npencils), intent(inout) :: lpencil_in
!
      if (lpencil_in(i_ugcc)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_gcc)=.true.
      endif
!
    endsubroutine pencil_interdep_pscalar
!***********************************************************************
    subroutine calc_pencils_pscalar(f,p)
!
!  Calculates Pscalar pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  04-feb-11/ccyang: coded
!
      use Sub, only: grad, u_dot_grad
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      type(pencil_case), intent(inout) :: p
!
      real, dimension(nx,3) :: uu
!
      integer :: i, icc1
!
      species: do i = 1, npscalar
        icc1 = icc + i - 1
        if (lpencil(i_cc)) p%cc(:,i) = f(l1:l2,m,n,icc1)
        if (lpencil(i_gcc)) call grad(f, icc1, p%gcc(:,:,i))
        if (lpencil(i_ugcc)) then
          if (i == 1) then
            uu = p%uu
            if (lconst_advection) uu = uu + spread(u0_advec, 1, nx)
          endif
          call u_dot_grad(f, icc1, p%gcc(:,:,i), uu, p%ugcc(:,i))
        endif
      enddo species
!
    endsubroutine calc_pencils_pscalar
!***********************************************************************
    subroutine dlncc_dt(f,df,p)
!
!  Passive scalar evolution
!
!  02-feb-11/ccyang: coded
!
      use Diagnostics
      use Sub, only: del2, dot_mn
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      type(pencil_case), intent(in) :: p
!
      real, dimension(nx) :: del2cc, gshockgcc
!
      integer :: i, icc1
!
      species: do i = 1, npscalar
        icc1 = icc + i - 1
!
        df(l1:l2,m,n,icc1) = df(l1:l2,m,n,icc1) - p%ugcc(:,i)
        if (lper_unit_volume) df(l1:l2,m,n,icc1) = df(l1:l2,m,n,icc1) - p%cc(:,i) * p%divu
!
!  Shock Diffusion
!
        if (lshock_diffusion) then
          if (lper_unit_volume) then
            call del2(f, icc1, del2cc)
            call dot_mn(p%gshock, p%gcc(:,:,i), gshockgcc)
            df(l1:l2,m,n,icc1) = df(l1:l2,m,n,icc1) + diffcc_shock * (p%shock * del2cc + gshockgcc)
          else
!
!  TO DO: correction for shock mass diffusion if lper_unit_volume=F
!
            call fatal_error('dlncc_dt', 'shock diffusion is not implemented for lper_unit_volume = F')
          endif
        endif
!
!  Diagnostics
!
        if (ldiagnos) then
          if (idiag_ccmax /= 0) call max_mn_name(p%cc(:,i), idiag_ccmax)
          if (idiag_ccmin /= 0) call max_mn_name(-p%cc(:,i), idiag_ccmin, lneg=.true.)
          if (idiag_ccm /= 0) call sum_mn_name(p%cc(:,i)/real(npscalar), idiag_ccm)
          if (idiag_rhoccm /= 0) call sum_mn_name(p%rho*p%cc(:,i)/real(npscalar), idiag_rhoccm)
        endif
!
      enddo species
!
    endsubroutine dlncc_dt
!***********************************************************************
    subroutine read_pscalar_init_pars(unit,iostat)
!
!  Reads the namelist pscalar_init_pars.
!
!  31-jan-11/ccyang: coded
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      integer :: status
!
      read(unit,NML=pscalar_init_pars,IOSTAT=status)
      if (present(iostat)) then
        iostat=status
      elseif (status/=0) then
        call fatal_error('read_pscalar_init_pars', 'failed to read the pscalar_init_pars')
      endif
!
    endsubroutine read_pscalar_init_pars
!***********************************************************************
    subroutine write_pscalar_init_pars(unit)
!
!  Writes the namelist pscalar_init_pars.
!
      integer, intent(in) :: unit
!
      integer :: status
!
      write(unit,NML=pscalar_init_pars,IOSTAT=status)
      if (status/=0) call fatal_error('write_pscalar_init_pars', 'failed to write the pscalar_init_pars')
!
    endsubroutine write_pscalar_init_pars
!***********************************************************************
    subroutine read_pscalar_run_pars(unit,iostat)
!
!  Reads the nameplist pscalar_run_pars
!
!  31-jan-11/ccyang: coded
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      integer :: status
!
      read(unit,NML=pscalar_run_pars,IOSTAT=status)
      if (present(iostat)) then
        iostat=status
      elseif (status/=0) then
        call fatal_error('read_pscalar_run_pars', 'failed to read the pscalar_run_pars')
      endif
!
    endsubroutine read_pscalar_run_pars
!***********************************************************************
    subroutine write_pscalar_run_pars(unit)
!
!  Writes the namelist pscalar_run_pars
!
!  31-jan-11/ccyang: coded
!
      integer, intent(in) :: unit
!
      integer :: status
!
      write(unit,NML=pscalar_run_pars,IOSTAT=status)
      if (status/=0) call fatal_error('write_pscalar_run_pars', 'failed to write the pscalar_run_pars')
!
    endsubroutine write_pscalar_run_pars
!***********************************************************************
    subroutine rprint_pscalar(lreset,lwrite)
!
!  Reads and registers print parameters relevant for passive scalar.
!
!   31-jan-11/ccyang: coded.
!
      use Diagnostics
!
      logical, intent(in) :: lreset
      logical, intent(in), optional :: lwrite
!
      character(len=80) :: fmt
      integer :: iname
      integer :: i
!
!  Reset everything in case of reset.
!
      if (lreset) then
        idiag_ccmax=0
        idiag_ccmin=0
        idiag_ccm=0
        idiag_rhoccm=0
      endif
!
!  Turn on requested diagnositc variables.
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'ccmax',idiag_ccmax)
        call parse_name(iname,cname(iname),cform(iname),'ccmin',idiag_ccmin)
        call parse_name(iname,cname(iname),cform(iname),'ccm',idiag_ccm)
        call parse_name(iname,cname(iname),cform(iname),'rhoccm',idiag_ccm)
      enddo
!
!  Write column where which passive scalar variable is stored.
!
      if (present(lwrite)) then
        if (lwrite) then
          write(3,*) 'ilncc = 0'
          if (npscalar > 1) then
            write(fmt,'(I2)') npscalar - 1
            fmt = '(1X, A, ' // trim(fmt) // '(I2, A2), I2, A)'
            write(3,fmt) 'icc = [', (i, ', ', i = icc, icc+npscalar-2), icc+npscalar-1, ']'
          else
            write(3,*) 'icc = ', icc
          endif
        endif
      endif
!
    endsubroutine rprint_pscalar
!***********************************************************************
    subroutine get_slices_pscalar(f,slices)
!
!  Write slices for animation of pscalar variables.
!
!  31-jan-11/ccyang: adapted from pscalar.f90
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      type(slice_data), intent(inout) :: slices
!
      integer :: icc1
!
!  Loop over slices.
!
      select case (trim(slices%name))
!
!  Passive scalar.
!
        case ('cc')
          if (0 <= slices%index .and. slices%index < npscalar) then
            icc1 = icc + slices%index
            slices%index = slices%index + 1
            slices%yz = f(slices%ix,m1:m2,n1:n2,icc1)
            slices%xz = f(l1:l2,slices%iy,n1:n2,icc1)
            slices%xy = f(l1:l2,m1:m2,slices%iz,icc1)
            slices%xy2 = f(l1:l2,m1:m2,slices%iz2,icc1)
            if (lwrite_slice_xy3) slices%xy3 = f(l1:l2,m1:m2,slices%iz3,icc1)
            if (lwrite_slice_xy4) slices%xy4 = f(l1:l2,m1:m2,slices%iz4,icc1)
            slices%ready = .true.
          else
            slices%ready = .false.
          endif
!
!  Logarithmic passive scalar.
!
        case ('lncc')
          if (0 <= slices%index .and. slices%index < npscalar) then
            icc1 = icc + slices%index
            slices%index = slices%index + 1
            slices%yz = alog(f(slices%ix,m1:m2,n1:n2,icc1))
            slices%xz = alog(f(l1:l2,slices%iy,n1:n2,icc1))
            slices%xy = alog(f(l1:l2,m1:m2,slices%iz,icc1))
            slices%xy2 = alog(f(l1:l2,m1:m2,slices%iz2,icc1))
            if (lwrite_slice_xy3) slices%xy3 = alog(f(l1:l2,m1:m2,slices%iz3,icc1))
            if (lwrite_slice_xy4) slices%xy4 = alog(f(l1:l2,m1:m2,slices%iz4,icc1))
            slices%ready = .true.
          else
            slices%ready = .false.
          endif
!
      endselect
!
    endsubroutine get_slices_pscalar
!***********************************************************************
    subroutine calc_mpscalar
!
!  Stub
!
    endsubroutine calc_mpscalar
!***********************************************************************
endmodule Pscalar
