! $Id$
!
!  This module takes care of everything related to massive particles
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_nbody=.false.
!
!***************************************************************
module Particles_nbody

  use Cdata
  use Messages
  use Particles_cdata
  use Particles_sub

  implicit none

  include 'particles_nbody.h'

  contains

!***********************************************************************
    subroutine register_particles_nbody()
!
!  Set up indices for access to the f, fsp and dfsp arrays.
!
!  14-jun-06/anders: adapted
!
!
    endsubroutine register_particles_nbody
!***********************************************************************
    subroutine initialize_particles_nbody(lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  14-jun-06/anders: adapted
!
      logical :: lstarting
!
      if (NO_WARN) print*, lstarting
!
    endsubroutine initialize_particles_nbody
!***********************************************************************
    subroutine init_particles_nbody(f,fp)
!
! Initial positions and velocities of nbody particles
! Overwrite the position asserted by the dust module
!
! 28-aug-06/wlad: dummy
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mpar_loc,mpvar) :: fp
!
      if (NO_WARN) print*,f,fp
!
    endsubroutine init_particles_nbody
!***********************************************************************
    subroutine pencil_criteria_par_nbody()
!
!  All pencils that the Particles_nbody module depends on are specified here.
!
!  02-jul-06/anders: adapted
!
    endsubroutine pencil_criteria_par_nbody
!***********************************************************************
    subroutine pencil_interdep_par_nbody(lpencil_in)
!
!  Interdependency among pencils provided by the Particles_nbody module
!  is specified here.
!
!  02-jul-06/anders: adapted
!
      logical, dimension(npencils) :: lpencil_in
!
      if (NO_WARN) print*, lpencil_in
!
    endsubroutine pencil_interdep_par_nbody
!***********************************************************************
    subroutine calc_pencils_par_nbody(f,p)
!
!  Calculate particle pencils.
!
!  02-jul-06/anders: adapted
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      if (NO_WARN) print*, f, p
!
    endsubroutine calc_pencils_par_nbody
!***********************************************************************
    subroutine dvvp_dt_nbody_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of nbody particles' velocities (called from main pencil loop)
!
!  27-aug-06/wlad: coded
!
      use Messages, only: fatal_error
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
!
      intent (in) :: f, df, p, fp, dfp, ineargrid
!
      if (NO_WARN) print*, f, df, fp, dfp, p, ineargrid
!
    endsubroutine dvvp_dt_nbody_pencil
!***********************************************************************
    subroutine dxxp_dt_nbody(dfp)
!
      real, dimension(mpar_loc,mpvar) :: dfp
!
      if (NO_WARN) print*,dfp
!
    endsubroutine dxxp_dt_nbody
!***********************************************************************
    subroutine dvvp_dt_nbody(f,df,fp,dfp,ineargrid)
!
!  Evolution of nbody particles velocities
!
!  27-aug-06/wlad: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      if (NO_WARN) print*, f, df, fp, dfp, ineargrid
!
    endsubroutine dvvp_dt_nbody
!***********************************************************************
    subroutine remove_particles_sink_nbody(f,fp,dfp,ineargrid)
!
      real,    dimension (mx,my,mz,mfarray) :: f
      real,    dimension (mpar_loc,mpvar)   :: fp, dfp
      integer, dimension (mpar_loc,3)       :: ineargrid
!
      if (NO_WARN) print*, f, fp, dfp, ineargrid
!
    endsubroutine remove_particles_sink_nbody
!***********************************************************************
    subroutine create_sink_particles_nbody(f,fp,dfp,ineargrid)
!
      real,    dimension (mx,my,mz,mfarray) :: f
      real,    dimension (mpar_loc,mpvar)   :: fp, dfp
      integer, dimension (mpar_loc,3)       :: ineargrid
!
      if (NO_WARN) print*, f, fp, dfp, ineargrid
!
    endsubroutine create_sink_particles_nbody
!***********************************************************************
    subroutine calc_nbodygravity_particles(f)
!
      real, dimension(mx,my,mz,mfarray) :: f
      if (NO_WARN) print*, f
!      
    endsubroutine calc_nbodygravity_particles
!***********************************************************************
    subroutine read_particles_nbody_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
      if (present(iostat) .and. (NO_WARN)) print*,iostat
      if (NO_WARN) print*,unit
    endsubroutine read_particles_nbody_init_pars
!***********************************************************************
    subroutine write_particles_nbody_init_pars(unit)
      integer, intent(in) :: unit
      if (NO_WARN) print*,unit
    endsubroutine write_particles_nbody_init_pars
!***********************************************************************
    subroutine read_particles_nbody_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
      if (present(iostat) .and. (NO_WARN)) print*,iostat
      if (NO_WARN) print*,unit
    endsubroutine read_particles_nbody_run_pars
!***********************************************************************
    subroutine write_particles_nbody_run_pars(unit)
      integer, intent(in) :: unit
      if (NO_WARN) print*,unit
    endsubroutine write_particles_nbody_run_pars
!***********************************************************************
    subroutine get_totalmass(tmass)
      real :: tmass
      if (NO_WARN) print*,tmass
    endsubroutine get_totalmass
!***********************************************************************
    subroutine bcast_nbodyarray(fp)
      real, dimension(mpar_loc,mpvar) :: fp
      if (NO_WARN) print*, fp
    endsubroutine bcast_nbodyarray
!************************************************************************
    subroutine particles_nbody_special
      real :: dummy_=0.
      if (NO_WARN) print*,dummy_
    endsubroutine particles_nbody_special
!************************************************************************
    subroutine particles_nbody_read_snapshot(filename)
!
! Read nbody particle info
!
! 01-apr-08/wlad: dummy
!
      character (len=*) :: filename
      if (NO_WARN) print*, filename
!
    endsubroutine particles_nbody_read_snapshot
!************************************************************************
    subroutine particles_nbody_write_snapshot(snapbase,enum,flist)
!
! Input and output of information about the massive particles
!
! 01-apr-08/wlad: dummy
!
      character (len=*) :: snapbase,flist
      logical :: enum
      optional :: flist
!          
      if (NO_WARN) print*, snapbase, enum, flist
!
    endsubroutine particles_nbody_write_snapshot
!************************************************************************
    subroutine particles_nbody_write_spdim(filename)
!
!  Write nspar, mspvar and mspar to file.
!
!  01-apr-08/wlad: dummy
!
      character (len=*) :: filename
!
      if (NO_WARN) print*, filename
!
    endsubroutine particles_nbody_write_spdim
!***********************************************************************
    subroutine rprint_particles_nbody(lreset,lwrite)
!
!  Read and register print parameters relevant for particle self-gravity.
!
!  14-jun-06/anders: adapted
!
      use Cdata
      use Sub, only: parse_name
!
      logical :: lreset
      logical, optional :: lwrite
      logical :: lwr
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
      if (NO_WARN) print*,lreset
!
    endsubroutine rprint_particles_nbody
!***********************************************************************
endmodule Particles_nbody
