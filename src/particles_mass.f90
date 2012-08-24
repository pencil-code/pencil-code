! $Id$
!
!  This module takes care of everything related to the mass represented by
!  each (super)particle.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MPVAR CONTRIBUTION 1
! CPARAM logical, parameter :: lparticles_mass=.true.
!
!***************************************************************
module Particles_mass
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
  use Particles_sub
!
  implicit none
!
  include 'particles_mass.h'
!
  real :: rhop_swarm0=1.0, rhop_swarm1=1.0
  real :: gravr_swarm0=1.0, gravr_swarm1=1.0
  real :: eps_dtog=0.01
  real, pointer :: rhs_poisson_const
  character (len=labellen), dimension(ninit) :: initrhopswarm='nothing'
!
  namelist /particles_mass_init_pars/ &
      initrhopswarm, rhop_swarm0, rhop_swarm1, gravr_swarm0, gravr_swarm1, &
      eps_dtog
!
  namelist /particles_mass_run_pars/ &
      initrhopswarm, rhop_swarm0, rhop_swarm1, gravr_swarm0, gravr_swarm1
!
  contains
!***********************************************************************
    subroutine register_particles_mass()
!
!  Set up indices for access to the fp and dfp arrays.
!
!  22-nov-10/anders+michiel: adapted
!
      if (lroot) call svn_id( &
          "$Id$")
!
!  Index for particle mass.
!
      irhopswarm=npvar+1
      pvarname(npvar+1)='irhopswarm'
!
!  Increase npvar accordingly.
!
      npvar=npvar+1
!
!  Check that the fp and dfp arrays are big enough.
!
      if (npvar > mpvar) then
        if (lroot) write(0,*) 'npvar = ', npvar, ', mpvar = ', mpvar
        call fatal_error('register_particles_mass: npvar > mpvar','')
      endif
!
    endsubroutine register_particles_mass
!***********************************************************************
    subroutine initialize_particles_mass(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  22-nov-10/anders+michiel: adapted
!
      use SharedVariables, only: get_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      if (lselfgravity) then
        call get_shared_variable('rhs_poisson_const',rhs_poisson_const)
      endif
!
      if (.not.(lcartesian_coords.and.(all(lequidist)))) call fatal_error( &
           'initialize_particles_mass', 'particles_mass only implemented '// &
           'for Cartesian equidistant grids.')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_particles_mass
!***********************************************************************
    subroutine init_particles_mass(f,fp)
!
!  Initial particle mass.
!
!  22-nov-10/anders+michiel: adapted
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mpvar) :: fp
!
      real :: rhom
      integer :: j, k
!
      do j=1,ninit
!
        select case (initrhopswarm(j))
!
        case ('nothing')
          if (lroot.and.j==1) print*, 'init_particles_mass: nothing'
!
        case ('constant')
          if (lroot) then
            print*, 'init_particles_mass: constant particle mass density'
            print*, 'init_particles_mass: rhop_swarm0=', rhop_swarm0
          endif
          fp(1:npar_loc,irhopswarm)=rhop_swarm0
!
        case ('constant-rhop')
          fp(1:npar_loc,irhopswarm)=rhop_swarm0/(float(npar)/nwgrid)
!
        case ('constant-1')
          if (lroot) then
            print*, 'init_particles_mass: set particle 1 mass density'
            print*, 'init_particles_mass: rhop_swarm1=', rhop_swarm1
          endif
          do k=1,npar_loc
            if (ipar(k)==1) fp(k,irhopswarm)=rhop_swarm1
          enddo
!
        case ('constant-grav')
          if (lroot) then
            print*, 'init_particles_mass: constant particle gravity'
            print*, 'init_particles_mass: gravr_swarm=', gravr_swarm0
          endif
          if (.not. lselfgravity) then
            if (lroot) print*, 'init_particles_mass: need selfgravity '// &
                'module for this initial condition'
            call fatal_error('init_particles_mass','')
          endif
          fp(1:npar_loc,irhopswarm)=gravr_swarm0/ &
              (rhs_poisson_const/(4*pi)*dx**3)
!
        case ('constant-grav-1')
          if (lroot) then
            print*, 'init_particles_mass: set particle 1 gravity'
            print*, 'init_particles_mass: gravr_swarm1=', gravr_swarm1
          endif
          if (.not. lselfgravity) then
            if (lroot) print*, 'init_particles_mass: need selfgravity '// &
                'module for this initial condition'
            call fatal_error('init_particles_mass','')
          endif
          do k=1,npar_loc
            if (ipar(k)==1) &
                fp(k,irhopswarm)=gravr_swarm1/(rhs_poisson_const/(4*pi)*dx**3)
          enddo
!
        case ('rhopswarm-from-particles-module','particles-to-gas-ratio')
          fp(1:npar_loc,irhopswarm)=rhop_swarm
!
        case ('particles-to-gas-ratio-thin-disk')
          rhom=sqrt(2*pi)*1.0*1.0/Lz  ! rhom = Sigma/Lz, Sigma=sqrt(2*pi)*H*rho1
          fp(1:npar_loc,irhopswarm)=eps_dtog*rhom/(real(npar)/nwgrid)
!
        case default
          if (lroot) print*, 'init_particles_mass: '// &
              'No such such value for initrhopswarm: ', trim(initrhopswarm(j))
          call fatal_error('init_particles_mass','')
        endselect
!
      enddo
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_particles_mass
!***********************************************************************
    subroutine pencil_criteria_par_mass()
!
!  All pencils that the Particles_mass module depends on are specified
!  here.
!
!  22-nov-10/anders+michiel: adapted
!
    endsubroutine pencil_criteria_par_mass
!***********************************************************************
    subroutine drhopswarm_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of particle mass.
!
!  22-nov-10/anders+michiel: adapted
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
!
      logical :: lheader, lfirstcall=.true.
!
      intent (in) :: f, fp
      intent (out) :: dfp
!
!  Print out header information in first time step.
!
      lheader=lfirstcall .and. lroot
!
!  Identify module and boundary conditions.
!
      if (lheader) print*, 'drhopswarm_dt_pencil: Calculate drhopswarm_dt'
!
      lfirstcall=.false.
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(p)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine drhopswarm_dt_pencil
!***********************************************************************
    subroutine drhopswarm_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of internal particle number.
!
!  22-nov-10/anders+michiel: adapted
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine drhopswarm_dt
!***********************************************************************
    subroutine read_particles_mass_init_pars(unit,iostat)
!
!  22-nov-10/anders+michiel: adapted
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_mass_init_pars,ERR=99,IOSTAT=iostat)
      else
        read(unit,NML=particles_mass_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_mass_init_pars
!***********************************************************************
    subroutine write_particles_mass_init_pars(unit)
!
!  22-nov-10/anders+michiel: adapted
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_mass_init_pars)
!
    endsubroutine write_particles_mass_init_pars
!***********************************************************************
    subroutine read_particles_mass_run_pars(unit,iostat)
!
!  22-nov-10/anders+michiel: adapted
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_mass_run_pars,ERR=99,IOSTAT=iostat)
      else
        read(unit,NML=particles_mass_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_mass_run_pars
!***********************************************************************
    subroutine write_particles_mass_run_pars(unit)
!
!  22-nov-10/anders+michiel: adapted
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_mass_run_pars)
!
    endsubroutine write_particles_mass_run_pars
!***********************************************************************
    subroutine rprint_particles_mass(lreset,lwrite)
!
!  Read and register print parameters relevant for particle mass.
!
!  22-nov-10/anders+michiel: adapted
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname
      logical :: lwr
!
!  Write information to index.pro.
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
      if (lwr) write(3,*) 'irhopswarm=', irhopswarm
!
      call keep_compiler_quiet(lreset)
!
    endsubroutine rprint_particles_mass
!***********************************************************************
endmodule Particles_mass
