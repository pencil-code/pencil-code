! $Id: particles_temperature.f90 21950 2014-07-08 08:53:00Z michiel.lambrechts $
!
!  This module takes care of everything related to the mass of the particles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MPVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
! CPARAM logical, parameter :: lparticles_mass=.true.
!
!! PENCILS PROVIDED TTp
!
!***************************************************************
module Particles_mass
!
  use Cdata
  use Cparam
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
!  use Particles_map
  use Particles_mpicomm
  use Particles_sub
  use Particles_chemistry
!
  implicit none
!
  include 'particles_mass.h'
!
  logical :: lpart_mass_backreac=.true.
  real :: mass_const
  character (len=labellen), dimension (ninit) :: init_particle_mass='nothing'
!
  namelist /particles_mass_init_pars/ &
      init_particle_mass, mass_const
!
  namelist /particles_mass_run_pars/ &
      lpart_mass_backreac
!
  integer :: idiag_mpm=0
!
  contains
!***********************************************************************
    subroutine register_particles_mass()
!
!  Set up indices for access to the fp and dfp arrays.
!
!  23-sep-14/Nils: adapted
!
      if (lroot) call svn_id( &
          "$Id: particles_mass.f90 20849 2013-08-06 18:45:43Z anders@astro.lu.se $")
!
!  Index for particle mass.
!
      imp=npvar+1
      pvarname(npvar+1)='imp'
print*,'imp,npvar=',imp,npvar

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
!  23-sep-14/Nils: adapted
!
      use SharedVariables, only: get_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
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
!  23-sep-14/Nils: adapted
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mpvar) :: fp
!
      real :: rhom
      integer :: j, k
!
      do j=1,ninit
!
        select case (init_particle_mass(j))
!
        case ('nothing')
          if (lroot.and.j==1) print*, 'init_particles_mass: nothing'
!
        case ('constant')
          if (lroot) then
            print*, 'init_particles_mass: constant particle mass'
            print*, 'init_particles_mass: mass_const=', mass_const
          endif
          fp(1:npar_loc,imp)=mass_const
!
        endselect
!
      enddo
!
    endsubroutine init_particles_mass
!***********************************************************************
    subroutine pencil_criteria_par_mass()
!
!  All pencils that the Particles_mass module depends on are specified
!  here.
!
!  23-sep-14/Nils: adapted
!
    endsubroutine pencil_criteria_par_mass
!***********************************************************************
    subroutine dpmass_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of particle temperature.
!
!  23-sep-14/Nils: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
!  Diagnostic output
!
      if (ldiagnos) then
        if (idiag_mpm/=0)  call sum_par_name(fp(1:npar_loc,imp),idiag_mpm)
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dpmass_dt
!***********************************************************************
    subroutine dpmass_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of particle temperature.
!
!  23-sep-14/Nils: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
      real, dimension(nx) :: feed_back, volume_pencil
      real :: volume_cell
      real, dimension (mpar_loc) :: St, Rc_hat
      integer :: k
!
      intent (in) :: f, fp, ineargrid
      intent (inout) :: dfp, df
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
      call keep_compiler_quiet(ineargrid)
!
      feed_back=0.
!
!  Check if particles chemistry is turned on
!
      if (lparticles_chemistry) then
!
!  Get total surface area and molar reaction rate of carbon
!
        call get_St(St,fp)
        call get_R_c_hat(Rc_hat,fp)

!  Loop over all particles in current pencil.
!
        do k=k1_imn(imn),k2_imn(imn)
!
!  Calculate the change in particle mass
!
          dfp(k,imp)=-St(k)*Rc_hat(k)*mol_mass_carbon
!
        enddo
!
      else
        do k=k1_imn(imn),k2_imn(imn)
          dfp(k,imp)=-1e-3
        enddo
      endif
!
    endsubroutine dpmass_dt_pencil
!***********************************************************************
    subroutine read_particles_mass_init_pars(unit,iostat)
!
!  23-sep-14/Nils: adapted
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
!  23-sep-14/Nils: adapted
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_mass_init_pars)
!
    endsubroutine write_particles_mass_init_pars
!***********************************************************************
    subroutine read_particles_mass_run_pars(unit,iostat)
!
!  23-sep-14/Nils: adapted
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
!  23-sep-14/Nils: adapted
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
!  23-sep-14/Nils: adapted
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
      if (lwr) write(3,*) 'imp=', imp
!
      call keep_compiler_quiet(lreset)
!
    endsubroutine rprint_particles_mass
!***********************************************************************
endmodule Particles_mass
