! $Id: particles_dust.f90 21950 2014-07-08 08:53:00Z michiel.lambrechts $
!
!  This module takes care of everything related to inertial particles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MPVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
! CPARAM logical, parameter :: lparticles_temperature=.true.
!
!! PENCILS PROVIDED TTp
!
!***************************************************************
module Particles_temperature
!
  use Cdata
  use Cparam
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
  use Particles_map
  use Particles_mpicomm
  use Particles_sub
  use Particles_radius
 !
  implicit none
!
  include 'particles_temperature.h'
!
  real :: init_part_temp, emissivity
  character (len=labellen), dimension (ninit) :: init_particle_temperature='nothing'
!
  namelist /particles_temp_init_pars/ &
      init_particle_temperature, init_part_temp, emissivity
!
  namelist /particles_temp_run_pars/ &
      emissivity
!
  contains
!***********************************************************************
    subroutine register_particles_temp()
!
!  Set up indices for access to the fp and dfp arrays
!
!  27-aug-14/jonas+nils: coded
!
      use FArrayManager, only: farray_register_auxiliary
!
      if (lroot) call svn_id( &
           "$Id: particles_dust.f90 21950 2014-07-08 08:53:00Z michiel.lambrechts $")
!
!  Indices for particle position.
!
      iTp=npvar+1
      pvarname(npvar+1)='iTp'
!
!  Increase npvar accordingly.
!
      npvar=npvar+1
!
!  Check that the fp and dfp arrays are big enough.
!
      if (npvar > mpvar) then
        if (lroot) write(0,*) 'npvar = ', npvar, ', mpvar = ', mpvar
        call fatal_error('register_particles_temp','npvar > mpvar')
      endif
!
    endsubroutine register_particles_temp
!***********************************************************************
    subroutine initialize_particles_temp(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  28-aug-14/jonas+nils: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      
!
    end subroutine initialize_particles_temp
!***********************************************************************
    subroutine init_particles_temp(f,fp)
!
!  Initial particle temperature
!
!  28-aug-14/jonas+nils: coded
!
      use General, only: random_number_wrapper
      use Mpicomm, only: mpireduce_sum, mpibcast_real
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mpvar) :: fp
      integer :: j
!

!
      intent (out) :: f, fp
!
!  Initial particle position.
!
      fp(1:npar_loc,iTp)=0.
      do j=1,ninit
!
        select case (init_particle_temperature(j))
!
        case ('nothing')
          if (lroot .and. j==1) print*, 'init_particles: nothing'
        case ('constant')
          if (lroot) print*, 'init_particles_temp: Constant temperature'
          fp(1:npar_loc,iTp)=fp(1:npar_loc,iTp)+init_part_temp
        case default
          if (lroot) &
              print*, 'init_particles_temp: No such such value for init_particle_temperature: ', &
              trim(init_particle_temperature(j))
          call fatal_error('init_particles_temp','')
!
        endselect
!
      enddo
!
    endsubroutine init_particles_temp
!***********************************************************************
    subroutine dpTT_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of particle temperature.
!
!  28-aug-14/jonas+nils: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dpTT_dt
!***********************************************************************
    subroutine dpTT_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of particle temperature.
!
!  28-aug-14/jonas+nils: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
!
      intent (in) :: f, df, fp, ineargrid
      intent (inout) :: dfp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(p)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dpTT_dt_pencil
!***********************************************************************
    subroutine read_particles_temp_init_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_temp_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_temp_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_temp_init_pars
!***********************************************************************
    subroutine write_particles_temp_init_pars(unit)
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_temp_init_pars)
!
    endsubroutine write_particles_temp_init_pars
!***********************************************************************
    subroutine read_particles_temp_run_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_temp_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_temp_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_temp_run_pars
!***********************************************************************
    subroutine write_particles_temp_run_pars(unit)
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_temp_run_pars)
!
    endsubroutine write_particles_temp_run_pars
!***********************************************************************
    subroutine rprint_particles_temp(lreset,lwrite)
!
!  Read and register print parameters relevant for particles temperature.
!
!  28-aug-14/jonas+nils: coded
!
      logical :: lreset
      logical, optional :: lwrite
!
      logical :: lwr
!
!  Write information to index.pro
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
      if (lwr) write(3,*) 'iox=', iox
!
      call keep_compiler_quiet(lreset)
!
    endsubroutine rprint_particles_temp
!***********************************************************************
    subroutine particles_temp_prepencil_calc(f)
!
!  28-aug-14/jonas+nils: coded
!
      real,dimension(mx,my,mz,mfarray),intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine particles_temp_prepencil_calc
!***********************************************************************
    subroutine pencil_criteria_par_temperature()
!
!  All pencils that the Particles_temperature module depends on are specified here.
!
!  28-aug-14/jonas+nils: coded
!
    endsubroutine pencil_criteria_par_temperature
!***********************************************************************
  end module Particles_temperature
