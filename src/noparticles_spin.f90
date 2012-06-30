! $Id$
!
!  This module takes care of everything related to particle spin.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_spin=.false.
!
!***************************************************************
module Particles_spin
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Particles_cdata
!
  implicit none
!
  include 'particles_spin.h'
!
  contains
!***********************************************************************
    subroutine register_particles_spin()
!
!  Set up indices for access to the fp and dfp arrays
!
!  23-jul-08/kapelrud: dummy
!
    endsubroutine register_particles_spin
!***********************************************************************
    subroutine initialize_particles_spin(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  23-jul-08/kapelrud: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_particles_spin
!***********************************************************************
    subroutine init_particles_spin(f,fp)
!
!  Initial spin of particles.
!
!  21-jul-08/kapelrud: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mpvar) :: fp
!
      intent(inout) :: f
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
!
    endsubroutine init_particles_spin
!***********************************************************************
    subroutine prepare_curl_vectorfield(f)
!
!  Prepare the curl(uu) field here so that ghost zones can be communicated
!  between processors before the spin is calculated in dpl_dt_pencil.
!
!  22-jul-08/kapelrud: coded
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine prepare_curl_vectorfield
!***********************************************************************
    subroutine pencil_criteria_par_spin()
!
!  All pencils that the Particles_spin module depends on are specified here.
!
!  21-nov-06/anders: coded
!
    endsubroutine pencil_criteria_par_spin
!***********************************************************************
    subroutine dps_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of particle spin.
!
!  22-aug-05/anders: coded
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
    endsubroutine dps_dt_pencil
!***********************************************************************
    subroutine dps_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of particle spin.
!
!  21-nov-06/anders: coded
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
    endsubroutine dps_dt
!***********************************************************************
    subroutine read_particles_spin_init_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_particles_spin_init_pars
!***********************************************************************
    subroutine write_particles_spin_init_pars(unit)
!
      integer, intent (in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_particles_spin_init_pars
!***********************************************************************
    subroutine read_particles_spin_run_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_particles_spin_run_pars
!***********************************************************************
    subroutine write_particles_spin_run_pars(unit)
!
      integer, intent (in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_particles_spin_run_pars
!***********************************************************************
    subroutine rprint_particles_spin(lreset,lwrite)
!
!  Read and register print parameters relevant for particles spin.
!
!  21-jul-08/kapelrud: adapted from particles_radius
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
    endsubroutine rprint_particles_spin
!***********************************************************************
    subroutine calc_liftforce(fp,k,rep,liftforce)
!
!  22-jul-08/kapelrud: dummy coded
!
      real,dimension(mpvar) :: fp
      integer :: k
      real,dimension(3) :: liftforce
      real :: rep
!
      intent(in) :: fp, k, rep
      intent(out) :: liftforce
!
      liftforce=0.0
!
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(k)
      call keep_compiler_quiet(rep)
!
    endsubroutine calc_liftforce
!***********************************************************************
    subroutine particles_spin_prepencil_calc(f)
!
!  12-aug-08/kapelrud: dummy coded
!
      real,dimension(mx,my,mz,mfarray),intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine particles_spin_prepencil_calc
!***********************************************************************

endmodule Particles_spin
