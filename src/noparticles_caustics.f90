! $Id$
!
!  This module writes information about the local state of the gas at
!  the positions of a selected number of particles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_caustics=.false.
!
!***************************************************************
module Particles_caustics
!
  use Cdata
  use Messages
  use General
  use Particles_cdata
  use Particles_map
  use Particles_mpicomm
!
  implicit none
!
  include 'particles_caustics.h'
!
  contains
!***********************************************************************
    subroutine register_particles_caustics
!
!  Set up indices for access to the fp and dfp arrays
!
!  May-16/dhruba: coded
!
      use FArrayManager, only: farray_register_auxiliary
!
      if (lroot) call svn_id( &
          "$Id$")
!
    endsubroutine register_particles_caustics
!***********************************************************************
    subroutine initialize_particles_caustics(f)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  May-16/dhruba+nishant+akshay: coded
!
      use General, only: keep_compiler_quiet
      use FArrayManager
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_particles_caustics
!***********************************************************************
    subroutine init_particles_caustics(f,fp,ineargrid)
!
      real, dimension (mx,my,mz,mfarray), intent (in) :: f
      real, dimension (mpar_loc,mparray), intent (out) :: fp
      integer, dimension (mpar_loc,3), intent (in) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine init_particles_caustics
!***********************************************************************
    subroutine dcaustics_dt(f,df,fp,dfp,ineargrid)
!
      use Diagnostics
!
      real, dimension (mx,my,mz,mfarray), intent (in) :: f
      real, dimension (mx,my,mz,mvar), intent (inout) :: df
      real, dimension (mpar_loc,mparray), intent (in) :: fp
      real, dimension (mpar_loc,mpvar), intent (inout) :: dfp
      integer, dimension (mpar_loc,3), intent (in) :: ineargrid
      
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)     
!
    endsubroutine dcaustics_dt
!***********************************************************************
    subroutine dcaustics_dt_pencil(f,df,fp,dfp,p,ineargrid,k,taup1)
!
      use Sub, only: linarray2matrix,matrix2linarray
!
      real, dimension (mx,my,mz,mfarray),intent(in) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (mpar_loc,mparray), intent(in) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      integer, dimension (mpar_loc,3) :: ineargrid
      integer :: k
      real :: taup1
      intent (inout) :: df, dfp,ineargrid
      intent (in) :: k,taup1
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(fp,dfp)
      call keep_compiler_quiet(p)
      call keep_compiler_quiet(ineargrid)     
      call keep_compiler_quiet(k)
      call keep_compiler_quiet(taup1)
!
    endsubroutine dcaustics_dt_pencil
!***********************************************************************
    subroutine read_pcaustics_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat

      iostat=0
!
    endsubroutine read_pcaustics_init_pars
!***********************************************************************
    subroutine write_pcaustics_init_pars(unit)
!
      integer, intent(in) :: unit
      call keep_compiler_quiet(unit)
!
    endsubroutine write_pcaustics_init_pars
!***********************************************************************
    subroutine read_pcaustics_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat

      iostat=0
!
    endsubroutine read_pcaustics_run_pars
!***********************************************************************
    subroutine write_pcaustics_run_pars(unit)
!
      integer, intent(in) :: unit
      call keep_compiler_quiet(unit)
!
    endsubroutine write_pcaustics_run_pars
!***********************************************************************
    subroutine rprint_particles_caustics(lreset,lwrite)
!
!  Read and register print parameters relevant for particles.
!
!  may-2016/dhruba+akshay: coded
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset,lwrite)
!    
    endsubroutine rprint_particles_caustics
!***********************************************************************
    subroutine reset_caustics(fp)
      real, dimension (mpar_loc,mparray), intent (out) :: fp
      call keep_compiler_quiet(fp)
    endsubroutine reset_caustics
!***********************************************************************
    subroutine reinitialize_caustics(fp)
      real, dimension (mpar_loc,mparray), intent (out) :: fp
      call keep_compiler_quiet(fp)

   endsubroutine reinitialize_caustics
!***********************************************************************
  endmodule Particles_caustics
