! $Id: noparticles_number.f90,v 1.1 2005-11-24 15:37:08 ajohan Exp $
!
!  This module takes care of everything related to particle number.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MPVAR CONTRIBUTION 1
! CPARAM logical, parameter :: lparticles_number=.false.
!
!***************************************************************
module Particles_number

  use Cdata
  use Particles_cdata
  use Particles_sub

  implicit none

  include 'particles_number.h'

  contains

!***********************************************************************
    subroutine register_particles_number()
!
!  Set up indices for access to the fp and dfp arrays.
!
!  24-nov-05/anders: dummy
!
    endsubroutine register_particles_number
!***********************************************************************
    subroutine initialize_particles_number(lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  24-nov-05/anders: dummy
!
      logical :: lstarting
!
      if (NO_WARN) print*, lstarting
!
    endsubroutine initialize_particles_number
!***********************************************************************
    subroutine init_particles_number(f,fp)
!
!  Initial internal particle number.
!
!  24-nov-05/anders: dummy
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mpar_loc,mpvar) :: fp
!
      if (NO_WARN) print*, f, fp
!
    endsubroutine init_particles_number
!***********************************************************************
    subroutine dnptilde_dt(f,df,fp,dfp)
!
!  Evolution of internal particle number.
!
!  24-oct-05/anders: dummy
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
!
      if (NO_WARN) print*, f, df, fp, dfp
!
    endsubroutine dnptilde_dt
!***********************************************************************
    subroutine read_particles_num_init_pars(unit,iostat)
!    
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (NO_WARN) print*, unit, iostat
!
    endsubroutine read_particles_num_init_pars
!***********************************************************************
    subroutine write_particles_num_init_pars(unit)
!    
      integer, intent (in) :: unit
!
      if (NO_WARN) print*, unit
!
    endsubroutine write_particles_num_init_pars
!***********************************************************************
    subroutine read_particles_num_run_pars(unit,iostat)
!    
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (NO_WARN) print*, unit, iostat
!
    endsubroutine read_particles_num_run_pars
!***********************************************************************
    subroutine write_particles_num_run_pars(unit)
!    
      integer, intent (in) :: unit
!
      if (NO_WARN) print*, unit
!
    endsubroutine write_particles_num_run_pars
!***********************************************************************
    subroutine rprint_particles_number(lreset,lwrite)
!   
!  Read and register print parameters relevant for internal particle number.
!
!  24-nov-05/anders: dummy
!
      use Cdata
      use Sub, only: parse_name
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname
      logical :: lwr
!
!  Write information to index.pro
! 
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!        
      if (lwr) write(3,*) 'inptilde=', inptilde
!
    endsubroutine rprint_particles_number
!***********************************************************************

endmodule Particles_number
