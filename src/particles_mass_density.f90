! $Id: particles_number.f90 14421 2010-07-23 23:55:02Z Bourdin.KIS $
!
!  This module takes care of everything related to the mass density
!  represented by each (super)particle.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MPVAR CONTRIBUTION 1
! CPARAM logical, parameter :: lparticles_mass_density=.true.
!
!***************************************************************
module Particles_mass_density
!
  use Cdata
  use Messages
  use Particles_cdata
  use Particles_sub
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'particles_mass_density.h'
!
  real :: rhop_swarm0=1.0, rhop_swarm1=1.0
  character (len=labellen), dimension(ninit) :: initrhopswarm='nothing'
!
  namelist /particles_mass_density_init_pars/ &
      initrhopswarm, rhop_swarm0, rhop_swarm1
!
  namelist /particles_mass_density_run_pars/ &
      initrhopswarm, rhop_swarm0, rhop_swarm1
!
  contains
!***********************************************************************
    subroutine register_particles_mass_density()
!
!  Set up indices for access to the fp and dfp arrays.
!
!  22-nov-10/anders+michiel: adapted
!
      if (lroot) call svn_id( &
          "$Id: particles_number.f90 14421 2010-07-23 23:55:02Z Bourdin.KIS $")
!
!  Index for particle mass density.
!
      irhopswarm=npvar+1
!
!  Increase npvar accordingly.
!
      npvar=npvar+1
!
!  Check that the fp and dfp arrays are big enough.
!
      if (npvar > mpvar) then
        if (lroot) write(0,*) 'npvar = ', npvar, ', mpvar = ', mpvar
        call fatal_error('register_particles_mass_density: npvar > mpvar','')
      endif
!
    endsubroutine register_particles_mass_density
!***********************************************************************
    subroutine initialize_particles_mass_density(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  22-nov-10/anders+michiel: adapted
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_particles_mass_density
!***********************************************************************
    subroutine init_particles_mass_density(f,fp)
!
!  Initial particle mass density.
!
!  22-nov-10/anders+michiel: adapted
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mpvar) :: fp
!
      integer :: j, k
!
      do j=1,ninit
!
        select case (initrhopswarm(j))
!
        case ('nothing')
          if (lroot.and.j==1) print*, 'init_particles_mass_density: nothing'
!
        case ('constant')
          if (lroot) then
            print*, 'init_particles_mass_density: constant particle mass density'
            print*, 'init_particles_mass_density: rhop_swarm0=', rhop_swarm0
          endif
          fp(1:npar_loc,irhopswarm)=rhop_swarm0
!
        case ('constant-1')
          if (lroot) then
            print*, 'init_particles_mass_density: set particle 1 mass density'
            print*, 'init_particles_mass_density: rhop_swarm1=', rhop_swarm1
          endif
          do k=1,npar_loc
            if (ipar(k)==1) fp(k,irhopswarm)=rhop_swarm1
          enddo
!
        case ('particles-to-gas-ratio')
          fp(1:npar_loc,irhopswarm)=rhop_swarm
!
        endselect
!
      enddo
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_particles_mass_density
!***********************************************************************
    subroutine pencil_criteria_par_mass_density()
!
!  All pencils that the Particles_mass_density module depends on are specified
!  here.
!
!  22-nov-10/anders+michiel: adapted
!
    endsubroutine pencil_criteria_par_mass_density
!***********************************************************************
    subroutine drhopswarm_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of particle mass density.
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
    subroutine read_particles_dens_init_pars(unit,iostat)
!
!  22-nov-10/anders+michiel: adapted
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_mass_density_init_pars,ERR=99,IOSTAT=iostat)
      else
        read(unit,NML=particles_mass_density_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_dens_init_pars
!***********************************************************************
    subroutine write_particles_dens_init_pars(unit)
!
!  22-nov-10/anders+michiel: adapted
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_mass_density_init_pars)
!
    endsubroutine write_particles_dens_init_pars
!***********************************************************************
    subroutine read_particles_dens_run_pars(unit,iostat)
!
!  22-nov-10/anders+michiel: adapted
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_mass_density_run_pars,ERR=99,IOSTAT=iostat)
      else
        read(unit,NML=particles_mass_density_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_dens_run_pars
!***********************************************************************
    subroutine write_particles_dens_run_pars(unit)
!
!  22-nov-10/anders+michiel: adapted
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_mass_density_run_pars)
!
    endsubroutine write_particles_dens_run_pars
!***********************************************************************
    subroutine rprint_particles_mass_density(lreset,lwrite)
!
!  Read and register print parameters relevant for particle mass density.
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
    endsubroutine rprint_particles_mass_density
!***********************************************************************
endmodule Particles_mass_density
