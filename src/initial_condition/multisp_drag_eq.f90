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
  use Particles_sub, only: dragforce_equi_multispecies 
!
  implicit none
!
  include '../initial_condition.h'
!
! Module Variables
!
  real, dimension(npar_species) :: taus = 0.0
  real :: eta_vK = 0.0
!
! Input Parameters
!
  real :: logtausmin = -4.0, logtausmax = -1.0
  real :: dlnndlntaus = -4.0
  real :: dlnrhodlnr = -0.1
  real :: zp0 = 0.0
!
  namelist /initial_condition_pars/ &
    logtausmin, logtausmax, dlnndlntaus, dlnrhodlnr, zp0
!
  contains
!***********************************************************************
    subroutine initialize_initial_condition(f)
!
!  Initialize any module variables which are parameter dependent.
!
!  02-oct-19/ccyang: coded
!
      use EquationOfState, only: cs0
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
!
      integer :: i
      real :: dlogtaus
!
      call keep_compiler_quiet(f)
!
!  Assemble stopping times.
!
      if (npar_species > 1) then
        dlogtaus = (logtausmax - logtausmin) / real(npar_species - 1)
      else
        dlogtaus = 0.5 * (logtausmax + logtausmin)
      endif
!
      taus = logtausmin + real((/ (i, i = 0, npar_species - 1) /)) * dlogtaus
      if (lroot) print *, "initialize_initial_condition: log(taus) = ", taus
      taus = 10.0**taus
!
!  Evaluate the radial pressure gradient support.
!
      eta_vK = -0.5 * dlnrhodlnr * cs0
      if (lroot) print*, 'initialize_initial_condition: eta * v_K = ', eta_vK
!
    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine initial_condition_uu(f)
!
!  Initialize the velocity field.
!
!  03-oct-19/ccyang: coded
!
      use EquationOfState, only: cs0, rho0
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      real, dimension(npar_species) :: eps, eps0, vpx, vpy
!
      integer :: n
      real :: ux, uy, hgas
!
! Find the mid-plane density for each species.
!
      eps0 = taus**(4.0 + dlnndlntaus)
      eps0 = eps_dtog / sum(eps0) * eps0
      if (zp0 > 0.0) eps0 = eps0 / zp0
!
! Assign the equilibrium velocities to the gas
!
      ueq: if (zp0 > 0.0) then
!       Stratified case
        hgas = cs0 / omega
        zloop: do n = n1, n2
          eps = eps0 * exp(0.5 * ((z(n) / hgas)**2 - (z(n) / zp0)**2))
          call dragforce_equi_multispecies(npar_species, taus, eps, eta_vK, vpx, vpy, ux, uy)
          f(l1:l2,m1:m2,n,iux) = f(l1:l2,m1:m2,n,iux) + ux
          f(l1:l2,m1:m2,n,iuy) = f(l1:l2,m1:m2,n,iuy) + uy
        enddo zloop
      else ueq
!       Unstratified case
        call dragforce_equi_multispecies(npar_species, taus, eps0, eta_vK, vpx, vpy, ux, uy)
        f(l1:l2,m1:m2,n1:n2,iux) = f(l1:l2,m1:m2,n1:n2,iux) + ux
        f(l1:l2,m1:m2,n1:n2,iuy) = f(l1:l2,m1:m2,n1:n2,iuy) + uy
!       Save the equilibrium velocities to a file.
        open(10, file="data/multisp_drag_eq.dat", form="unformatted", action="write")
        write(10) ux, uy, vpx, vpy
        close(10)
      endif ueq
!
    endsubroutine initial_condition_uu
!***********************************************************************
    subroutine initial_condition_vvp(f, fp)
!
!  Initialize particles' mass and velocity.
!
!  01-oct-19/ccyang: coded
!
      use EquationOfState, only: cs0, rho0
      use SharedVariables, only: get_shared_variable
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(:,:), intent(inout) :: fp
!
      real, dimension(:), pointer :: tausp_species, tausp1_species
!
      real, dimension(npar_species) :: rhopj, eps0, eps, vpx, vpy
!
      integer :: k, p
      real :: ux, uy, hgas, zp
!
      call keep_compiler_quiet(f)
!
!  Find the solid-to-gas mass ratio for each species.
!
      eps0 = taus**(4.0 + dlnndlntaus)
      eps0 = eps_dtog / sum(eps0) * eps0
!
!  Find the mass of each particle.
!
      rhop: if (zp0 > 0.0) then
!      Stratified case
        hgas = cs0 / omega
        rhopj = sqrt(2.0 * pi) * real(nxgrid * nygrid) * rho0 * hgas / (real(npar / npar_species) * dz) * eps0
        eps0 = eps0 / zp0
      else rhop
!      Unstratified case
        rhopj = rho0 / real(npar / (npar_species * nxgrid * nygrid * nzgrid)) * eps0
        call dragforce_equi_multispecies(npar_species, taus, eps0, eta_vK, vpx, vpy, ux, uy)
      endif rhop
      if (lroot) print *, "initial_condition_vvp: rhopj = ", rhopj
!
!  Assign the mass and velocity of each particle.
!
      ploop: do k = 1, npar_loc
!
        p = npar_species * (ipar(k) - 1) / npar + 1
        fp(k,irhopswarm) = rhopj(p)
!
        strat: if (zp > 0.0) then
          zp = fp(k,izp)
          eps = eps0 * exp(0.5 * ((zp / hgas)**2 - (zp / zp0)**2))
          call dragforce_equi_multispecies(npar_species, taus, eps, eta_vK, vpx, vpy, ux, uy)
        endif strat
!
        fp(k,ivpx) = fp(k,ivpx) + vpx(p) 
        fp(k,ivpy) = fp(k,ivpy) + vpy(p)
      enddo ploop
!
!  Override the stopping times in particles_dust.
!
      call get_shared_variable("tausp_species", tausp_species)
      call get_shared_variable("tausp1_species", tausp1_species)
      tausp_species = taus / omega
      tausp1_species = omega / taus
      if (lroot) print *, "initial_condition_vvp: override tausp_species = ", tausp_species
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
