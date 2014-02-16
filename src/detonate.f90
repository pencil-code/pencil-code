! $Id$
!
!  This module searches for gravitationally collapsing sites and
!  detonates them.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: ldetonate = .true.
!
! MAUX CONTRIBUTION 1
! COMMUNICATED AUXILIARIES 1
!
!***************************************************************
module Detonate
!
  use Cparam
  use Cdata
  use Messages, only: svn_id, fatal_error, warning
!
  implicit none
!
  include 'detonate.h'
!
!  Runtime parameters
!
  character(len=labellen) :: det_scale = 'linear'    ! Scaling of the energy input at detonation with density
  integer :: det_njeans = 4    ! Minimum number of cells per Jeans length for stability
  integer :: det_radius = 2    ! Radius of a local region in number of cells
  real :: det_factor = 1.0     ! Scale factor of the energy input
!
  namelist /detonate_run_pars/ det_scale, det_njeans, det_radius, det_factor
!
!  Diagnostic variables
!
  integer :: idiag_detn = 0      ! DIAG_DOC: Number of detonated sites
  integer :: idiag_dettot = 0    ! DIAG_DOC: Total energy input
!
!  Module variables
!
  logical, dimension(-nghost:nghost,-nghost:nghost,-nghost:nghost) :: mask_sphere = .false.
  integer :: nxg = 0, nyg = 0, nzg = 0
  integer :: power = 0
  real :: deth0 = 0.0
  real :: jeans = 0.0
  real :: tgentle = 0.0
!
  contains
!***********************************************************************
    subroutine register_detonate()
!
!  Set up indices for variables.
!
!  06-feb-14/ccyang: coded
!
      use FArrayManager, only: farray_register_auxiliary
!
      if (lroot) call svn_id("$Id$")
!
      call farray_register_auxiliary('detonate', idet, communicated=.true.)
!
    endsubroutine register_detonate
!***********************************************************************
    subroutine initialize_detonate(f, lstarting)
!
!  Initializes module specific variables for later use.
!
!  14-feb-14/ccyang: coded
!
      use General, only: keep_compiler_quiet
      use EquationOfState, only: cs20, gamma, gamma_m1
      use SharedVariables, only: get_shared_variable
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      logical, intent(in) :: lstarting
!
      real, pointer :: tsg
      integer :: istat
      integer :: i, j, k
      real :: volume
!
      call keep_compiler_quiet(f)
!
!  This module currently only works with thermal_energy.
!
      if (.not. lthermal_energy) call fatal_error('initialize_detonate', 'currently only works with thermal_energy. ')
!
!  If called by start.f90: nothing to do.
!
      if (lstarting) return
!
!  Calculate the conversion factor between density and energy for Jeans
!  stability.
!
      if (nzgrid /= 1) then
        call fatal_error('initialize_detonate', '3D Jeans stability criterion is under construction. ')
      else
        jeans = real(det_njeans) * G_Newton * dxmax / (gamma * gamma_m1)
      endif
!
!  Make a spherical mask.
!
      nxg = merge(nghost, 0, nxgrid > 1)
      nyg = merge(nghost, 0, nygrid > 1)
      nzg = merge(nghost, 0, nzgrid > 1)
      forall (i=-nxg:nxg, j=-nyg:nyg, k=-nzg:nzg, i*i+j*j+k*k <= det_radius**2) mask_sphere(i,j,k) = .true.
!
!  Check the scaling of the energy input with density.
!
      scale: select case (det_scale)
!
      case ('linear') scale
!       Calculate the volume of the mask.
        ndim: select case (count((/ nxgrid > 1, nygrid > 1, nzgrid > 1 /)))
        case (0) ndim
          volume = 0.0
        case (1) ndim
          volume = real(2 * det_radius)
        case (2) ndim
          volume = pi * real(det_radius**2)
        case (3) ndim
          volume = pi / 3.0 * real(4 * det_radius**3)
        endselect ndim
        power = 1
        deth0 = det_factor * volume * cs20
!
      case ('quadratic', 'binding') scale
        power = 2
        if (nzgrid > 1) then
          call fatal_error('initialize_detonate', 'quadratic scaling for 3D is under construction. ')
        else
          deth0 = det_factor * (8.0 * pi / 3.0) * G_Newton * (det_radius * dxmax)**3
        endif
!
      case default scale
        call fatal_error('initialize_detonate', 'unknown det_scale')
      endselect scale
!
!  Get tselfgrav_gentle.
!
      selfgrav: if (lselfgravity) then
        call get_shared_variable('tselfgrav_gentle', tsg, istat)
        if (istat /= 0) call warning('initialize_detonate', 'unable to get tselfgrav_gentle')
        tgentle = tsg
      endif selfgrav
!
    endsubroutine initialize_detonate
!***********************************************************************
    subroutine read_detonate_run_pars(unit,iostat)
!
!  Reads the runtime parameters.
!
!  06-feb-14/ccyang: coded
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      integer :: stat
!
      read(unit, NML=detonate_run_pars, IOSTAT=stat)
      if (present(iostat)) then
        iostat = stat
      else if (stat /= 0) then
        call fatal_error('read_detonate_run_pars', 'cannot read detonate_run_pars. ')
      endif
!
    endsubroutine read_detonate_run_pars
!***********************************************************************
    subroutine write_detonate_run_pars(unit)
!
!  Writes the runtime parameters.
!
!  06-feb-14/ccyang: coded
!
      integer, intent(in) :: unit
!
      integer :: stat
!
      write(unit, NML=detonate_run_pars, IOSTAT=stat)
      if (stat /= 0) call warning('write_detonate_run_pars', 'cannot write detonate_run_pars. ')
!
    endsubroutine write_detonate_run_pars
!***********************************************************************
    subroutine rprint_detonate(lreset, lwrite)
!
!  Reads and registers print parameters.
!
!  06-feb-14/ccyang: dummy
!
      use Diagnostics, only: parse_name
      use General, only: keep_compiler_quiet
!
      logical, intent(in) :: lreset
      logical, intent(in), optional :: lwrite
!
      integer :: iname
!
!  Reset diagnostic indices if requested.
!
      reset: if (lreset) then
        idiag_detn = 0
        idiag_dettot = 0
      endif reset
!
!  Parse names in print.in.
!
      diagnostics: do iname = 1, nname
        call parse_name(iname, cname(iname), cform(iname), 'detn', idiag_detn)
        call parse_name(iname, cname(iname), cform(iname), 'dettot', idiag_dettot)
      enddo diagnostics
!
!  lwrite should be phased out.
!
      if (present(lwrite)) call keep_compiler_quiet(lwrite)
!
    endsubroutine rprint_detonate
!***********************************************************************
    subroutine detonate_before_boundary(f)
!
!  Possibility to modify the f before boundary conmmunications.
!
!  14-feb-14/ccyang: coded
!
      use Boundcond, only: zero_ghosts, update_ghosts
      use Mpicomm, only: mpiallreduce_or, mpireduce_sum_int, mpireduce_sum
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      logical, dimension(mx,my,mz) :: mask
      logical :: flag
      integer :: ndet
      real :: esum
!
!  Initialize the detonation energy to tiny positive value.
!
      f(l1:l2,m1:m2,n1:n2,idet) = tiny(1.0)
!
!  Detonate only before one full time-step.
!
      first: if (lfirst) then
!
!  Check Jeans stability for each cell.
!
        call zero_ghosts(f, ieth)
        call update_ghosts(f, ieth)
        jeans: if (ldensity_nolog) then
          call zero_ghosts(f, irho)
          call update_ghosts(f, irho)
          mask = jeans_unstable(f(:,:,:,irho), f(:,:,:,ieth))
        else jeans
          call zero_ghosts(f, ilnrho)
          call update_ghosts(f, ilnrho)
          mask = jeans_unstable(exp(f(:,:,:,ilnrho)), f(:,:,:,ieth))
        endif jeans
!
!  Sift out collapsing sites.
!
        if (any(mask)) call collapsing(f, mask)
!
!  Detonate them.
!
        call mpiallreduce_or(any(mask), flag)
        if (flag) call set_detonations(f, mask)
!
!  Diagnostics
!
        diagnos: if (ldiagnos .and. flag) then
!
          detn: if (idiag_detn /= 0) then
            call mpireduce_sum_int(count(mask), ndet)
            if (lroot) fname(idiag_detn) = real(ndet)
          endif detn
!
          dettot: if (idiag_dettot /= 0) then
            call mpireduce_sum(total_energy(f), esum)
            if (lroot) fname(idiag_dettot) = esum
          endif dettot
!
        endif diagnos
!
      endif first
!
    endsubroutine detonate_before_boundary
!***********************************************************************
!***********************************************************************
!  LOCAL PROCEDURES GO BELOW HERE.
!***********************************************************************
!***********************************************************************
    elemental logical function jeans_unstable(rho, eth)
!
!  Returns true if a cell is considered Jeans unstable given its density
!  and internal energy.
!
!  08-feb-14/ccyang: coded.
!
      real, intent(in) :: rho, eth
!
      real :: eth_min
!
!  Find the minimum internal energy for Jeans stability.
!
      if (nzgrid /= 1) then
!       3D run: rho is volume mass density
      else
!       2D run: rho is surface mass density
        eth_min = jeans * rho**2
      endif
!
!  Compare it with the actual energy.
!
      jeans_unstable = eth <= eth_min
!
    endfunction jeans_unstable
!***********************************************************************
    subroutine collapsing(f, unstable)
!
!  Searches for collapsing sites in Jeans unstable regions.
!
!  12-feb-14/ccyang: coded
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      logical, dimension(mx,my,mz), intent(inout) :: unstable
!
      logical, dimension(mx,my,mz) :: work
      integer, dimension(3) :: center
      integer :: ll1, ll2, mm1, mm2, nn1, nn2
      integer :: i, j, k
      real :: divu
!
!  Loop over each cell and check if it is a local maximum and converging.
!
      work = .false.
      center = (/ nxg + 1, nyg + 1, nzg + 1 /)
      zscan: do k = n1, n2
        nn1 = k - nzg
        nn2 = k + nzg
        yscan: do j = m1, m2
          mm1 = j - nyg
          mm2 = j + nyg
          xscan: do i = l1, l2
            ll1 = i - nxg
            ll2 = i + nxg
            localmax: if (all(unstable(ll1:ll2,mm1:mm2,nn1:nn2) .and. mask_sphere(-nxg:nxg,-nyg:nyg,-nzg:nzg) &
                                                                .eqv. mask_sphere(-nxg:nxg,-nyg:nyg,-nzg:nzg)) .and. &
                          all(maxloc(f(ll1:ll2,mm1:mm2,nn1:nn2,irho), &
                                     mask=mask_sphere(-nxg:nxg,-nyg:nyg,-nzg:nzg)) == center)) then
              divu = 0.0
              if (nxgrid > 1) divu = divu + derivative(f(ll1:ll2,j,k,iux), dx_1(i))
              if (nygrid > 1) divu = divu + derivative(f(i,mm1:mm2,k,iuy), dy_1(j))
              if (nzgrid > 1) divu = divu + derivative(f(i,j,nn1:nn2,iuz), dz_1(k))
              if (divu < 0.0) work(i,j,k) = .true.
            endif localmax
          enddo xscan
        enddo yscan
      enddo zscan
      unstable = work
!
    endsubroutine collapsing
!***********************************************************************
    subroutine set_detonations(f, mask)
!
!  Detonates cells specified by mask.
!
!  14-feb-14/ccyang: dummy
!
      use Boundcond, only: zero_ghosts, update_ghosts
      use Sub, only: smooth
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      logical, dimension(mx,my,mz), intent(in) :: mask
!
      real :: deth
!
      if (t < tgentle) then
        deth = 0.5 * deth0 * (1.0 - cos(pi * t / tgentle))
      else
        deth = deth0
      endif
!
!  Put in point energy where mask is .true.
!
      if (ldensity_nolog) then
        where(mask(l1:l2,m1:m2,n1:n2)) f(l1:l2,m1:m2,n1:n2,idet) = deth * f(l1:l2,m1:m2,n1:n2,irho)**power
      else
        where(mask(l1:l2,m1:m2,n1:n2)) f(l1:l2,m1:m2,n1:n2,idet) = deth * exp(power * f(l1:l2,m1:m2,n1:n2,ilnrho))
      endif
!
!  Smooth it.
!
      call zero_ghosts(f, idet)
      call update_ghosts(f, idet)
      call smooth(f, idet)
!
!  Add it to the energy field.
!
      f(l1:l2,m1:m2,n1:n2,ieth) = f(l1:l2,m1:m2,n1:n2,ieth) + f(l1:l2,m1:m2,n1:n2,idet)
      call zero_ghosts(f, ieth)
      call update_ghosts(f, ieth)
!
    endsubroutine set_detonations
!***********************************************************************
    real function derivative(a, dx1)
!
!  Returns the first derivative at the center of a seven-point data
!  array.
!
!  19-jan-13/ccyang: coded.
!
      real, dimension(-3:3), intent(in) :: a
      real, intent(in) :: dx1
!
      real, parameter :: sixtieth = 1.0 / 60.0
!
      derivative = sixtieth * dx1 * ((a(3) - a(-3)) - 9.0 * (a(2) - a(-2)) + 45.0 * (a(1) - a(-1)))
!
    endfunction derivative
!***********************************************************************
    real function total_energy(f)
!
!  Returns total energy input by detonation.
!
!  14-feb-14/ccyang: coded
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      real, dimension(nx) :: dx
      integer :: j, k
      real :: dy, dz
!
!  Initialization
!
      total_energy = 0.0
      if (nxgrid > 1) then
        dx = xprim(l1:l2)
      else
        dx = 1.0
      endif
!
!  Integrate the det field.
!
      zdir: do k = n1, n2
        dz = merge(zprim(k), 1.0, nzgrid > 1)
        ydir: do j = m1, m2
          dy = merge(yprim(j), 1.0, nygrid > 1)
          total_energy = total_energy + sum(dz * dy * dx * f(l1:l2,j,k,idet), mask=f(l1:l2,j,k,idet)>0.0)
        enddo ydir
      enddo zdir
!
    endfunction total_energy
!***********************************************************************
endmodule Detonate
