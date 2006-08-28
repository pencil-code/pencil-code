! $Id: noparticles_nbody.f90,v 1.3 2006-08-28 20:35:03 wlyra Exp $
!
!  This module takes care of everything related to particle self-gravity.
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
    subroutine init_particles_nbody(f,fp,fsp)
!
! Initial positions and velocities of sink particles
! Overwrite the position asserted by the dust module
!
! 28-aug-06/wlad: dummy
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mpar_loc,mpvar) :: fp
      real, dimension(nspar,mpvar) :: fsp
!
      if (NO_WARN) print*,f,fp,fsp
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
!  Evolution of sink particles' velocities (called from main pencil loop)
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
    subroutine dvvp_dt_nbody(f,df,fp,dfp,fsp,dfsp,ineargrid)
!
!  Evolution of sink particles velocities
!
!  27-aug-06/wlad: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      real, dimension (nspar,mpvar) :: fsp,dfsp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      if (NO_WARN) print*, f, df, fp, dfp, fsp, dfsp, ineargrid
!
    endsubroutine dvvp_dt_nbody
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
    subroutine get_particles_interdistances(fsp,rp_mn,rpcyl_mn)
!
! 18-jul-06/wlad: dummy subroutine
!
      real, dimension (nspar,mpvar) :: fsp
      real, dimension (nx,nspar) :: rp_mn,rpcyl_mn
      integer :: k
!
      intent(out) :: rp_mn,rpcyl_mn
!
       do k=1,nspar
          rp_mn(:,k)    = 0.
          rpcyl_mn(:,k) = 0.
       enddo
!
       if (NO_WARN) print*, fsp
!
     endsubroutine get_particles_interdistances
!************************************************************************
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
