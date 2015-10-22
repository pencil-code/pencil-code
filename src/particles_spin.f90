! $Id$
!
!  This module takes care of everything related to particle spin
!  including lifting forces. The module maintains a full f-array
!  vorticity field, to be able to interpolate on the flow vorticity.
!
!  The module should be considered experimental as it is virtually
!  untested (as of aug-08).
!
!  NOTE: all code relating to particle spin or the magnus force
!        have been commented out.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MPVAR CONTRIBUTION 3
! CPARAM logical, parameter :: lparticles_spin = .true.
!
!***************************************************************
module Particles_spin
!
  use Cdata
  use Cparam
  use Messages
  use Particles_cdata
!
  implicit none
!
  include 'particles_spin.h'
!
  character(len=labellen), dimension(ninit) :: initsp = 'nothing'
  real, dimension(ninit) :: amplsp = 0.0
  logical :: lsaffman_lift = .false.
  logical :: lmagnus_lift = .false.
!
  namelist /particles_spin_init_pars/ initsp, amplsp, lsaffman_lift, lmagnus_lift
!
  namelist /particles_spin_run_pars/ lsaffman_lift, lmagnus_lift
!
  integer :: idiag_psxm = 0, idiag_psym = 0, idiag_pszm = 0
!
  contains
!***********************************************************************
    subroutine register_particles_spin()
!
!  Set up indices for access to the fp and dfp arrays.
!
!  21-jul-08/kapelrud: coded
!
!      use FArrayManager
!
      if (lroot) call svn_id("$Id$")
!
!  Indices for particle spin
!
      ipsx = npvar + 1
      pvarname(ipsx) = 'ipsx'
      ipsy = npvar + 2
      pvarname(ipsy) = 'ipsy'
      ipsz = npvar + 3
      pvarname(ipsz) = 'ipsz'
!
      npvar = npvar + 3
!
!  Check that the fp and dfp arrays are big enough.
!
      bound: if (npvar > mpvar) then
        if (lroot) print *, 'npvar = ', npvar, ', mpvar = ', mpvar
        call fatal_error('register_particles_spin', 'npvar > mpvar')
      endif bound
!
    endsubroutine register_particles_spin
!***********************************************************************
    subroutine initialize_particles_spin(f)
!
!  Perform any post-parameter-read initialization, i.e., calculate
!  derived parameters.
!
!  21-jul-08/kapelrud: coded
!  22-oct-15/ccyang: continued.
!
      use General, only: keep_compiler_quiet
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
      if (lstart) return
!
!  Sanity check.
!
      if (lsaffman_lift) call fatal_error('initialize_particles_spin', 'Saffman lift is currently not supported. ')
!
      magnus: if (lmagnus_lift) then
        if (.not. lparticles_radius .and. particle_radius <= 0.0) &
            call fatal_error('initialize_particles_spin', 'Magnus lift requires the radius of each constituent particle. ')
        if (mpmat <= 0.0) &
            call fatal_error('initialize_particles_spin', 'Magnus lift requires the mass of each constituent particle. ')
      endif magnus
!
!  Request interpolation of variables:
!
      interp%luu = interp%luu .or. lsaffman_lift .or. lmagnus_lift
      interp%loo = interp%loo .or. lsaffman_lift .or. lmagnus_lift
      interp%lrho = interp%lrho .or. lsaffman_lift .or. lmagnus_lift
!
    endsubroutine initialize_particles_spin
!***********************************************************************
    subroutine init_particles_spin(f, fp)
!
!  Initial spin of particles.
!
!  21-jul-08/kapelrud: coded.
!  07-oct-15/ccyang: continued.
!
      use General, only: keep_compiler_quiet
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mpar_loc,mparray), intent(inout) :: fp
!
      integer :: j
!
      call keep_compiler_quiet(f)
!
      loop: do j = 1, ninit
        init: select case (initsp(j))
!
!  Do nothing.
!
        case ('nothing') init
          if (lroot .and. j == 1) print *, 'init_particles_spin: no initial condition for particle spin'
!
!  Zero out all spins.
!
        case ('zero') init
          if (lroot) print *, 'init_particles_spin: zero particle spin'
          fp(1:npar_loc,ipsx:ipsz) = 0.0
!
!  Random magnitude and orientation.
!
        case ('random') init
          if (lroot) print *, 'init_particles_spin: spins of random magnitude and orientation; amplsp = ', amplsp(j)
          call gaunoise_vect(amplsp(j), fp, ipsx, ipsz)
!
!  Unknown initial condition.
!
        case default init
          call fatal_error('init_particles_spin', 'unknown initsp = ' // initsp(j))
!
        endselect init
      enddo loop
!
    endsubroutine init_particles_spin
!***********************************************************************
    subroutine pencil_criteria_par_spin()
!
!  All pencils that the Particles_spin module depends on are specified
!  here.
!
!  06-oct-15/ccyang: stub.
!
    endsubroutine pencil_criteria_par_spin
!***********************************************************************
    subroutine dps_dt_pencil(f, df, fp, dfp, p, ineargrid)
!
!  Evolution of particle spin (called in the pencil loop.)
!
!  06-oct-15/ccyang: stub.
!
      use General, only: keep_compiler_quiet
!      use Viscosity, only: getnu
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mx,my,mz,mvar), intent(in) :: df
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      real, dimension(mpar_loc,mpvar), intent(inout) :: dfp
      type(pencil_case), intent(in) :: p
      integer, dimension(mpar_loc,3), intent(in) :: ineargrid
!
      logical :: lfirstcall = .true.
!      real, dimension(3) :: tau
      logical :: lheader
!      integer :: k
!      real :: ip_tilde, nu
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
!
!      call getnu(nu_input=nu)
!
!  Print out header information in first time step.
!
      lheader = lfirstcall .and. lroot
      lfirstcall = .false.
!
!  Identify module and boundary conditions.
!
      if (lheader) print *, 'dps_dt_pencil: Calculate dps_dt (currently do nothing)'
!!
!!  Calculate torque on particle due to the shear flow, and
!!  update the particles' spin.
!!
!     if (lmagnus_lift) then
!       do k=k1_imn(imn),k2_imn(imn)
!!
!!  Calculate angular momentum
!!
!         ip_tilde=0.4*mpmat*fp(k,iap)**2
!!
!         tau=8.0*pi*interp_rho(k)*nu*fp(k,iap)**3* &
!             (0.5*interp_oo(k,:)-fp(k,ipsx:ipsz))
!         dfp(k,ipsx:ipsz)=dfp(k,ipsx:ipsz)+tau/ip_tilde
!       enddo
!     endif
!
    endsubroutine dps_dt_pencil
!***********************************************************************
    subroutine dps_dt(f, df, fp, dfp, ineargrid)
!
!  Evolution of particle spin (called after the pencil loop.)
!
!  25-jul-08/kapelrud: coded
!
      use Particles_sub, only: sum_par_name
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mx,my,mz,mvar), intent(in) :: df
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      real, dimension(mpar_loc,mpvar), intent(in) :: dfp
      integer, dimension(mpar_loc,3), intent(in) :: ineargrid
!
!  Diagnostics
!
      diag: if (ldiagnos) then
        if (idiag_psxm /= 0) call sum_par_name(fp(1:npar_loc,ipsx), idiag_psxm)
        if (idiag_psym /= 0) call sum_par_name(fp(1:npar_loc,ipsy), idiag_psym)
        if (idiag_pszm /= 0) call sum_par_name(fp(1:npar_loc,ipsz), idiag_pszm)
      endif diag
!
    endsubroutine dps_dt
!***********************************************************************
    subroutine read_particles_spin_init_pars(iostat)
!
!  Read initialization parameters from namelist particles_spin_init_pars.
!
!  06-oct-15/ccyang: coded.
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_spin_init_pars, IOSTAT=iostat)
!
    endsubroutine read_particles_spin_init_pars
!***********************************************************************
    subroutine write_particles_spin_init_pars(unit)
!
!  Write initialization parameters from namelist particles_spin_init_pars.
!
!  06-oct-15/ccyang: coded.
!
      integer, intent(in) :: unit
!
      write(unit, NML=particles_spin_init_pars)
!
    endsubroutine write_particles_spin_init_pars
!***********************************************************************
    subroutine read_particles_spin_run_pars(iostat)
!
!  Read runtime parameters from namelist particles_spin_run_pars.
!
!  06-oct-15/ccyang: coded.
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_spin_run_pars, IOSTAT=iostat)
!
    endsubroutine read_particles_spin_run_pars
!***********************************************************************
    subroutine write_particles_spin_run_pars(unit)
!
!  Write runtime parameters from namelist particles_spin_run_pars.
!
!  06-oct-15/ccyang: coded.
!
      integer, intent(in) :: unit
!
      write(unit, NML=particles_spin_run_pars)
!
    endsubroutine write_particles_spin_run_pars
!***********************************************************************
    subroutine rprint_particles_spin(lreset, lwrite)
!
!  Read and register print parameters relevant for particles spin.
!
!  21-jul-08/kapelrud: coded.
!  06-oct-15/ccyang: continued.
!
      use Diagnostics
!
      logical, intent(in) :: lreset
      logical, intent(in), optional :: lwrite
!
      logical :: lwr
      integer :: iname
!
!  Write information to index.pro
!
      lwr = .false.
      if (present(lwrite)) lwr = lwrite
!
      indices: if (lwr) then
        write(3,*) "ipsx = ", ipsx
        write(3,*) "ipsy = ", ipsy
        write(3,*) "ipsz = ", ipsz
      endif indices
!
!  Reset everything in case of reset
!
      reset: if (lreset) then
        idiag_psxm = 0
        idiag_psym = 0
        idiag_pszm = 0
      endif reset
!
!  Run through all possible names that may be listed in print.in
!
      if (lroot .and. ip < 14) print *, 'rprint_particles_spin: run through parse list'
      diag: do iname = 1, nname
        call parse_name(iname, cname(iname), cform(iname), 'psxm', idiag_psxm)
        call parse_name(iname, cname(iname), cform(iname), 'psym', idiag_psym)
        call parse_name(iname, cname(iname), cform(iname), 'pszm', idiag_pszm)
      enddo diag
!
    endsubroutine rprint_particles_spin
!***********************************************************************
    subroutine calc_liftforce(fp, k, rep, liftforce)
!
!  Calculate lifting forces for a given particle. It should be possible
!  to make this a routine operating on pencils.
!
!  22-jul-08/kapelrud: coded
!
      real, dimension(mparray), intent(in) :: fp
      integer, intent(in) :: k
      real, intent(in) :: rep
      real, dimension(3), intent(out) :: liftforce
!
      real,dimension(3) :: dlift
!
!  Initialization
!
      liftforce = 0.0
!
!  Find Saffman lift.
!
!      if (lsaffman_lift) then
!        call calc_saffman_liftforce(fp,k,rep,dlift)
!        liftforce=liftforce+dlift
!      endif
!
!  Find Magnus list.
!
     magnus: if (lmagnus_lift) then
       call calc_magnus_liftforce(fp, k, rep, dlift)
       liftforce = liftforce + dlift
     endif magnus
!
    endsubroutine calc_liftforce
!***********************************************************************
    subroutine calc_saffman_liftforce(fp,k,rep,dlift)
!
!  Calculate the Saffman lifting force for a given particles.
!
!  16-jul-08/kapelrud: coded
!
      use Particles_cdata
      use Sub, only: cross
      use Viscosity, only: getnu
!
      real,dimension(mparray) :: fp
      integer :: k
      real,dimension(3) :: dlift
      real :: rep
!
      intent(in) :: fp, k, rep
      intent(out) :: dlift
!
      real :: csaff,diameter,beta,oo,nu
!
      call getnu(nu_input=nu)
!
      if (.not.lparticles_radius) then
        if (lroot) print*,'calc_saffman_liftforce: '//&
             'Particle_radius module must be enabled!'
        call fatal_error('calc_saffman_liftforce','')
      endif
!
      diameter=2*fp(iap)
      oo=sqrt(sum(interp_oo(k,:)**2))
!
      beta=diameter**2*oo/(2.0*rep*nu)
      if (beta<0.005) then
        beta=0.005
      elseif (beta>0.4) then
        beta=0.4
      endif
!
      if (rep<=40) then
        csaff=(1-0.3314*beta**0.5)*exp(-rep/10.0)+0.3314*beta**0.5
      else
        csaff=0.0524*(beta*rep)**0.5
      endif
!
      call cross(interp_uu(k,:)-fp(ivpx:ivpz),interp_oo(k,:),dlift)
      dlift=1.61*csaff*diameter**2*nu**0.5*&
                 interp_rho(k)*oo**(-0.5)*dlift/mpmat
!
    endsubroutine calc_saffman_liftforce
!***********************************************************************
    subroutine calc_magnus_liftforce(fp, k, rep, dlift)
!
!  Calculate the Magnus liftforce for a given spinning particle.
!
!  22-jul-08/kapelrud: coded
!  22-oct-15/ccyang: continued.
!
      use Sub, only: cross
      use Viscosity, only: getnu
!
      real, dimension(mparray), intent(in) :: fp
      integer, intent(in) :: k
      real, intent(in) :: rep
      real, dimension(3), intent(out) :: dlift
!
      real, dimension(3) :: ps_rel, uu_rel
      real :: const_lr, spin_omega, area, nu, ap
      real :: ps2, uu2
!
!  Find the relative velocity and spin.
!
      uu_rel = interp_uu(k,:) - fp(ivpx:ivpz)
      ps_rel = fp(ipsx:ipsz) - 0.5 * interp_oo(k,:)
      uu2 = sum(uu_rel**2)
      ps2 = sum(ps_rel**2)
!
      lift: if (uu2 > 0.0 .and. ps2 > 0.0) then
!
!  Get the radius of the constituent particle.
!
        if (lparticles_radius) then
          ap = fp(iap)
        else
          ap = particle_radius
        endif
!
!  Projected area of the particle
!
        area = pi * ap**2
!
!  Calculate the Magnus lift coefficent
!
        uu2 = sqrt(uu2)
        spin_omega = ap * sqrt(sum(fp(ipsx:ipsz)**2)) / uu2
        const_lr = min(0.5, 0.5 / spin_omega)
        call cross(uu_rel, ps_rel, dlift)
        call getnu(nu_input=nu)
        dlift = 0.25 * interp_rho(k) * (rep * nu / uu2) * const_lr * area / mpmat / sqrt(ps2) * dlift
!
      else lift
        dlift = 0.0
!
      endif lift
!
    endsubroutine calc_magnus_liftforce
!***********************************************************************
    subroutine gaunoise_vect(ampl, fp, ivar1, ivar2)
!
!  Add Gaussian noise for variables ivar1:ivar2
!
!  07-oct-15/ccyang: adapted from gaunoise_vect in Initcond.
!
      use General, only: random_number_wrapper
!
      real, intent(in) :: ampl
      real, dimension(mpar_loc,mparray), intent(inout) :: fp
      integer, intent(in) :: ivar1, ivar2
!
      real, dimension(npar_loc) :: r, p, tmp
      integer :: i
!
      if (lroot .and. ip <= 8) print *, 'gaunoise_vect: ampl, ivar1, ivar2=', ampl, ivar1, ivar2
      comp: do i = ivar1, ivar2
        random: if (modulo(i - ivar1, 2) == 0) then
          call random_number_wrapper(r)
          call random_number_wrapper(p)
          tmp = sqrt(-2.0 * log(r)) * sin(twopi * p)
        else random
          tmp = sqrt(-2.0 * log(r)) * cos(twopi * p)
        endif random
        fp(1:npar_loc,i) = fp(1:npar_loc,i) + ampl * tmp
      enddo comp
!
    endsubroutine gaunoise_vect
!***********************************************************************
endmodule Particles_spin
