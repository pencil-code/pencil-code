! $Id$
!
!  This module sets up nonlinear streaming instability with
!  multiple particle species.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!***************************************************************
module InitialCondition
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
  use Particles_sub, only: assign_species, dragforce_equi_multispecies 
!
  implicit none
!
  include '../initial_condition.h'
!
! Input Parameters
!
  complex, dimension(4*(npar_species+1)) :: si_ev = (0.0, 0.0)
  logical :: lsi_random = .false.
  logical :: ltaus_log_center = .true.
  real :: logtausmin = -4.0, logtausmax = -1.0
  real :: dlnndlntaus = -4.0
  real :: dlnrhodlnr = -0.1
  real :: si_kx = 0.0, si_kz = 0.0, si_amp = 1E-6
!
  namelist /initial_condition_pars/ &
    lsi_random, ltaus_log_center, logtausmin, logtausmax, dlnndlntaus, dlnrhodlnr, si_kx, si_kz, si_ev, si_amp
!
! Module Variables
!
  real, dimension(npar_species) :: taus = 0.0, eps0 = 0.0, rhopj = 0.0, vpx0 = 0.0, vpy0 = 0.0
  real :: eta_vK = 0.0, ux0 = 0.0, uy0 = 0.0
  real :: amp_scale = 0.0
!
  contains
!***********************************************************************
    subroutine initialize_initial_condition(f)
!
! Initialize any module variables which are parameter dependent.
!
! 25-jul-20/ccyang: coded
!
      use EquationOfState, only: cs0, rho0
      use Mpicomm, only: mpibcast
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
!
      integer :: i
      real :: dlogtaus
!
      if (lrun) return
      call keep_compiler_quiet(f)
!
! Assemble stopping times.
!
      dlogtaus = (logtausmax - logtausmin) / real(npar_species)
      gettaus: if (ltaus_log_center) then
        taus = logtausmin + real((/ (i, i = 1, npar_species) /) - 0.5) * dlogtaus
        if (lroot) print *, "initialize_initial_condition: log(taus) = ", taus
        taus = 10.0**taus
      else gettaus
        taus = 0.5 * 10.0**logtausmin * (10.0**(real((/ (i, i = 0, npar_species - 1) /)) * dlogtaus) + &
                                         10.0**(real((/ (i, i = 1, npar_species) /)) * dlogtaus))
        if (lroot) print *, "initialize_initial_condition: taus = ", taus
      endif gettaus
!
! Evaluate the radial pressure gradient support.
!
      eta_vK = -0.5 * dlnrhodlnr * cs0
      if (lroot) print*, 'initialize_initial_condition: eta * v_K = ', eta_vK
!
! Find the density ratio for each species.
!
      eps0 = taus**(4.0 + dlnndlntaus)
      eps0 = eps_dtog / sum(eps0) * eps0
!
! Find the mass of each particle.
!
      rhopj = rho0 / real(npar / (npar_species * nxgrid * nygrid * nzgrid)) * eps0
      if (lroot) print *, "initialize_initial_condition: rhopj = ", rhopj
!
! Find the equilibrium velocities.
!
      eqvel: if (lroot) then
        call dragforce_equi_multispecies(npar_species, taus, eps0, eta_vK, vpx0, vpy0, ux0, uy0)
        print *, "initialize_initial_condition: ux0, uy0 = ", ux0, uy0
        print *, "initialize_initial_condition: vpx0 = ", vpx0
        print *, "initialize_initial_condition: vpy0 = ", vpy0
      endif eqvel
!
      call mpibcast(ux0)
      call mpibcast(uy0)
      call mpibcast(vpx0, npar_species)
      call mpibcast(vpy0, npar_species)
!
! Save the equilibrium velocities to a file.
!
      record: if (lroot) then
        open(10, file="data/multisp_drag_eq.dat", form="unformatted", action="write")
        write(10) ux0, uy0, vpx0, vpy0
        close(10)
      endif record
!
! Note the method of perturbations.
!
      perturb: if (lsi_random) then
        if (lroot) print *, "initialize_initial_condition: randomly perturb particle positions, amplitude = ", si_amp
      else perturb
        amp_scale = si_amp * eps_dtog / sum(abs(si_ev(8::4)))
        if (lroot) print *, "initialize_initial_condition: exact wave mode, scaled amplitude = ", amp_scale
      endif perturb
!
    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine initial_condition_uu(f)
!
! Initialize the velocity field.
!
! 21-jul-20/ccyang: coded
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      real, dimension(nx) :: argx, sinkx, coskx, dux, duy, duz
      real :: argz, sinkz, coskz
!
! Assign the equilibrium velocities to the gas
!
      f(l1:l2,m1:m2,n1:n2,iux) = f(l1:l2,m1:m2,n1:n2,iux) + ux0
      f(l1:l2,m1:m2,n1:n2,iuy) = f(l1:l2,m1:m2,n1:n2,iuy) + uy0
!
      if (lsi_random) return
!
! Perturb the gas.
!
      argx = si_kx * x(l1:l2)
      sinkx = sin(argx)
      coskx = cos(argx)
!
      duz = amp_scale * eta_vK
      dux = duz * (real(si_ev(1)) * coskx - aimag(si_ev(1)) * sinkx)
      duy = duz * (real(si_ev(2)) * coskx - aimag(si_ev(2)) * sinkx)
      duz = duz * (real(si_ev(3)) * sinkx + aimag(si_ev(3)) * coskx)
!
      nloop: do n = n1, n2
        argz = si_kz * z(n)
        sinkz = sin(argz)
        coskz = cos(argz)
!
        mloop: do m = m1, m2
          f(l1:l2,m,n,iux) = f(l1:l2,m,n,iux) + dux * coskz
          f(l1:l2,m,n,iuy) = f(l1:l2,m,n,iuy) + duy * coskz
          f(l1:l2,m,n,iuz) = f(l1:l2,m,n,iuz) - duz * sinkz
        enddo mloop
      enddo nloop
!
    endsubroutine initial_condition_uu
!***********************************************************************
    subroutine initial_condition_lnrho(f)
!
! Initialize logarithmic density.
!
! 21-jul-20/ccyang: coded
!
      use EquationOfState, only: rho0
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      real, dimension(nx) :: argx, drho
      real :: coskz
!
      random: if (lsi_random) then
!
! Uniform gas density,
!
        f(l1:l2,m1:m2,n1:n2,irho) = rho0
!
      else random
!
! Or perturb the gas density.
!
        argx = si_kx * x(l1:l2)
        drho = amp_scale * rho0 * (real(si_ev(4)) * cos(argx) - aimag(si_ev(4)) * sin(argx))
!
        nloop: do n = n1, n2
          coskz = cos(si_kz * z(n))
          do m = m1, m2
            f(l1:l2,m,n,irho) = rho0 + drho * coskz
          enddo
        enddo nloop
      endif random
!
! Convert to logarithmic density.
!
      f(l1:l2,m1:m2,n1:n2,ilnrho) = log(f(l1:l2,m1:m2,n1:n2,irho))
!
    endsubroutine initial_condition_lnrho
!***********************************************************************
    subroutine initial_condition_xxp(f, fp)
!
! Initialize particles' positions.
!
! 25-jul-20/ccyang: coded
!
      use General, only: random_number_wrapper
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(:,:), intent(inout) :: fp
!
      integer, dimension(npar_species) :: ip
      real, dimension(npar_species) :: ar, ai, a1, a2, a3
      integer :: npps, npx, npz, ix, iz, is, k
      real :: dxp, dzp, xp, yp, zp, dxp1, dzp1
      real :: argx, argz, ampl
      real :: sinp, sinm, cosp, cosm
      real :: sinp2, sinm2, cosp2, cosm2
      real :: sin2kz, cos2kx, sin2kx
      real :: c1x, c1z, c2x, c2z
!
      call keep_compiler_quiet(f)
!
! Sanity checks
!
      proc: if (mod(npar, ncpus) /= 0 .or. npar_loc /= npar / ncpus) then
        if (lroot) print *, "initial_condition_xxp: npar, ncpus, npar_loc = ", npar, ncpus, npar_loc
        call fatal_error("initial_condition_xxp", "particles are not evenly distributed among processors")
      endif proc
!
      species: if (mod(npar_loc, npar_species) /= 0) then
        if (lroot) print *, "initial_condition_xxp: npar_loc, npar_species = ", npar_loc, npar_species
        call fatal_error("initial_condition_xxp", "particle species are not evenly divided")
      endif species
!
! Find the initial ID for each species.
!
      npps = npar_loc / npar_species
      ip = (/ ((k - 1) * (npar / npar_species) + iproc * npps, k = 1, npar_species) /)
!
! Find the spacing between particles.
!
      getnp: if (nzgrid > 1) then
        npx = nint(sqrt(Lxyz_loc(1) * real(npps) / Lxyz_loc(3)))
        npz = npps / npx
      else getnp
        npx = npps
        npz = 1
      endif getnp
!
      grid: if (npx * npz /= npps) then
        if (lroot) print *, "initial_condition_xxp: Lx_loc, Lz_loc = ", Lxyz_loc(1), Lxyz_loc(3)
        if (lroot) print *, "initial_condition_xxp: npps, npx, npz = ", npps, npx, npz
        call fatal_error("initial_condition_xxp", "cannot find equal spacing between particles")
      endif grid
!
      dxp = Lxyz_loc(1) / real(npx)
      dzp = Lxyz_loc(3) / real(npz)
      if (lroot) print *, "initial_condition_xxp: npx, npz = ", npx, npz
      if (lroot) print *, "initial_condition_xxp: dxp, dzp = ", dxp, dzp
!
      random: if (lsi_random) then
!
! Uniform distribution plus random perturbations:
!
        ampl = 3.2 * sqrt(real(npar / nxzgrid)) * si_amp * sqrt(dx * dz) / pi
        k = 0
        yp = xyz0(2) + 0.5 * Lxyz(2)
        zloop1: do iz = 1, npz
          zp = xyz0_loc(3) + (real(iz) - 0.5) * dzp
          xloop1: do ix = 1, npx
            xp = xyz0_loc(1) + (real(ix) - 0.5) * dxp
            sloop1: do is = 1, npar_species
              call random_number_wrapper(argx)
              call random_number_wrapper(argz)
              argx = ampl * sqrt(-2.0 * log(argx)) * sqrt(rhop_swarm / rhopj(is))
              argz = 2.0 * pi * argz
              dxp1 = argx * sin(argz)
              dzp1 = argx * cos(argz)
!
              k = k + 1
              fp(k,ixp) = xp + dxp1
              fp(k,iyp) = yp
              fp(k,izp) = zp + dzp1
              ip(is) = ip(is) + 1
              ipar(k) = ip(is)
            enddo sloop1
          enddo xloop1
        enddo zloop1
!
      else random
!
! Exact wave mode:
!
        c1x = si_kx**2 + si_kz**2
        c2x = c1x**2
        coeff: if (c1x > 0.0) then
          c1x = 0.5 / c1x
          c2x = 1.0 / c2x
        endif coeff
        c1z = c1x * si_kz
        c2z = c2x * si_kz**3
        c1x = c1x * si_kx
        c2x = c2x * si_kx**3
!
        ar = amp_scale * real(si_ev(8::4)) / eps0
        ai = amp_scale * aimag(si_ev(8::4)) / eps0
        a1 = 0.25 * (ar**2 - ai**2)
        a2 = 0.5 * ar * ai
        a3 = 0.25 * (ar**2 + ai**2)
!
        k = 0
        yp = xyz0(2) + 0.5 * Lxyz(2)
!
        zloop2: do iz = 1, npz
          zp = xyz0_loc(3) + (real(iz) - 0.5) * dzp
          argz = si_kz * zp
          sin2kz = sin(2.0 * argz)
!
          xloop2: do ix = 1, npx
            xp = xyz0_loc(1) + (real(ix) - 0.5) * dxp
            argx = si_kx * xp
            cos2kx = cos(2.0 * argx)
            sin2kx = sin(2.0 * argx)
!
            sinp = sin(argx + argz)
            sinm = sin(argx - argz)
            cosp = cos(argx + argz)
            cosm = cos(argx - argz)
!
            sinp2 = sin(2.0 * (argx + argz))
            sinm2 = sin(2.0 * (argx - argz))
            cosp2 = cos(2.0 * (argx + argz))
            cosm2 = cos(2.0 * (argx - argz))
!
            sloop2: do is = 1, npar_species
              dxp1 = -c1x * (ar(is) * (sinp + sinm) + ai(is) * (cosp + cosm) &
                           - a1(is) * (sinp2 + sinm2) - a2(is) * (cosp2 + cosm2)) &
                     + c2x * (a2(is) * cos2kx + a1(is) * sin2kx)
              dzp1 = -c1z * (ar(is) * (sinp - sinm) + ai(is) * (cosp - cosm) &
                           - a1(is) * (sinp2 - sinm2) - a2(is) * (cosp2 - cosm2)) &
                     + c2z * a3(is) * sin2kz
!
              k = k + 1
              fp(k,ixp) = xp + dxp1
              fp(k,iyp) = yp
              fp(k,izp) = zp + dzp1
              ip(is) = ip(is) + 1
              ipar(k) = ip(is)
            enddo sloop2
          enddo xloop2
        enddo zloop2
      endif random
!
    endsubroutine initial_condition_xxp
!***********************************************************************
    subroutine initial_condition_vvp(f, fp)
!
! Initialize particles' mass and velocity.
!
! 25-jul-20/ccyang: coded
!
      use SharedVariables, only: get_shared_variable
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(:,:), intent(inout) :: fp
!
      real, dimension(:), pointer :: tausp_species, tausp1_species
!
      real, dimension(npar_species) :: vpx, vpy
      real :: argx, argz, sinkx, coskx, sinkz, coskz
      real :: dv, dvpx, dvpy, dvpz
!
      integer :: k, p, i
      real :: ux, uy, hgas, zp
!
      call keep_compiler_quiet(f)
!
! Override the stopping times in particles_dust.
!
      multisp: if (npar_species > 1) then
        call get_shared_variable("tausp_species", tausp_species)
        call get_shared_variable("tausp1_species", tausp1_species)
        tausp_species = taus / omega
        tausp1_species = omega / taus
        if (lroot) print *, "initial_condition_vvp: override tausp_species = ", tausp_species
      endif multisp
!
! Assign the mass and the equilibrium velocity of each particle.
!
      equil: do k = 1, npar_loc
        p = assign_species(ipar(k))
        fp(k,irhopswarm) = rhopj(p)
        fp(k,ivpx) = fp(k,ivpx) + vpx0(p)
        fp(k,ivpy) = fp(k,ivpy) + vpy0(p)
      enddo equil
!
      if (lsi_random) return
!
! Perturb the velocity.
!
      dv = amp_scale * eta_vK
      perturb: do k = 1, npar_loc
        argx = si_kx * fp(k,ixp)
        argz = si_kz * fp(k,izp)
        sinkx = sin(argx)
        coskx = cos(argx)
        sinkz = sin(argz)
        coskz = cos(argz)
!
        p = assign_species(ipar(k))
        i = 4 * p + 1
        dvpx = dv * (real(si_ev(i)) * coskx - aimag(si_ev(i)) * sinkx) * coskz
        dvpy = dv * (real(si_ev(i+1)) * coskx - aimag(si_ev(i+1)) * sinkx) * coskz
        dvpz = -dv * (real(si_ev(i+2)) * sinkx + aimag(si_ev(i+2)) * coskx) * sinkz
!
        fp(k,ivpx) = fp(k,ivpx) + dvpx
        fp(k,ivpy) = fp(k,ivpy) + dvpy
        fp(k,ivpz) = fp(k,ivpz) + dvpz
      enddo perturb
!
    endsubroutine initial_condition_vvp
!***********************************************************************
    subroutine read_initial_condition_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=initial_condition_pars, IOSTAT=iostat)
!
    endsubroutine read_initial_condition_pars
!***********************************************************************
    subroutine write_initial_condition_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=initial_condition_pars)
!
    endsubroutine write_initial_condition_pars
!***********************************************************************
!
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from noinitial_condition.f90 for any    **
!**  InitialCondition routines not implemented in this file        **
!**                                                                **
    include '../initial_condition_dummies.inc'
!********************************************************************
endmodule InitialCondition
