! $Id: particles_selfgravity.f90,v 1.1 2006-06-15 19:34:43 ajohan Exp $
!
!  This module takes care of everything related to particle self-gravity.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_selfgravity=.true.
!
! MAUX CONTRIBUTION 4
! COMMUNICATED AUXILIARIES 3
!
!***************************************************************
module particles_selfgravity

  use Cdata
  use Messages
  use Particles_cdata
  use Particles_sub

  implicit none

  include 'particles_selfgravity.h'
  
  real :: dummy
  logical :: lselfgravity_particles=.true.

  namelist /particles_selfgrav_init_pars/ &
      lselfgravity_particles

  namelist /particles_selfgrav_run_pars/ &
      lselfgravity_particles

  contains

!***********************************************************************
    subroutine register_particles_selfgravity()
!
!  Set up indices for access to the f, fp and dfp arrays.
!
!  14-jun-06/anders: adapted
!
      use Messages, only: fatal_error, cvs_id
!
      logical, save :: first=.true.
!
      if (.not. first) &
          call fatal_error('register_particles_selfgravity: called twice','')
      first = .false.
!
      if (lroot) call cvs_id( &
           "$Id: particles_selfgravity.f90,v 1.1 2006-06-15 19:34:43 ajohan Exp $")
!
!  Index for gradient for the self-potential and for the smooth particle
!  density field.
!
      igpotselfx = mvar + naux_com + 1; naux = naux + 1; naux_com = naux_com + 1
      igpotselfy = mvar + naux_com + 1; naux = naux + 1; naux_com = naux_com + 1
      igpotselfz = mvar + naux_com + 1; naux = naux + 1; naux_com = naux_com + 1
      irhop      = mvar + naux + 1 + (maux_com - naux_com); naux = naux + 1
!
!  Check that we aren't registering too many auxiliary variables
!
      if (naux > maux) then
        if (lroot) write(0,*) 'naux = ', naux, ', maux= ', maux
        call fatal_error('register_shock','naux > maux')
      endif
!
    endsubroutine register_particles_selfgravity
!***********************************************************************
    subroutine initialize_particles_selfgravity(lstarting)
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
    endsubroutine initialize_particles_selfgravity
!***********************************************************************
    subroutine calc_selfpotential_particles(f,rhs_poisson,rhs_poisson_const,    lcontinued)
!
!  Calculate the gravity potential of the dust particles.
!
!  13-jun-06/anders: coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (nx,ny,nz) :: rhs_poisson
      real :: rhs_poisson_const
      logical :: lcontinued
!
      if (lselfgravity_particles) then
        if (lcontinued) then  ! Potential has already been zeroed by the gas.
          rhs_poisson=rhs_poisson+ &
              rhs_poisson_const*f(l1:l2,m1:m2,n1:n2,irhop)
        else                  ! Must zero potential from last time-step.
          rhs_poisson= &
              rhs_poisson_const*f(l1:l2,m1:m2,n1:n2,irhop)
        endif
      endif
!
    endsubroutine calc_selfpotential_particles
!***********************************************************************
    subroutine dvvp_dt_selfgravity(f,df,fp,dfp,ineargrid)
!
!  Add self-gravity to particle equation of motion.
!
!  14-jun-06/anders: coded
!
      use Messages, only: fatal_error
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (3) :: gpotself
      integer :: k
      logical :: lheader, lfirstcall=.true.
!
      intent (in) :: f, fp
      intent (out) :: dfp
!
!  Print out header information in first time step.
!
      lheader=lfirstcall .and. lroot
!
      lfirstcall=.false.
!
      if (lheader) print*, 'dvvp_dt_selfgravity: add self-gravity'
!
!  Interpolate the gradient of the potential to the location of the particle.
!
      if (lselfgravity_particles) then
        do k=1,npar_loc
          call interpolate_quadratic( &
              f,igpotselfx,igpotselfz,fp(k,ixp:izp),gpotself, &
              ineargrid(k,:),ipar(k) )
          dfp(k,ivpx:ivpz)=dfp(k,ivpx:ivpz)-gpotself
        enddo
      endif
!
    endsubroutine dvvp_dt_selfgravity
!***********************************************************************
    subroutine read_particles_selfg_init_pars(unit,iostat)
!    
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_selfgrav_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_selfgrav_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_selfg_init_pars
!***********************************************************************
    subroutine write_particles_selfg_init_pars(unit)
!    
      integer, intent (in) :: unit
!
      write(unit,NML=particles_selfgrav_init_pars)
!
    endsubroutine write_particles_selfg_init_pars
!***********************************************************************
    subroutine read_particles_selfg_run_pars(unit,iostat)
!    
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_selfgrav_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_selfgrav_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_selfg_run_pars
!***********************************************************************
    subroutine write_particles_selfg_run_pars(unit)
!    
      integer, intent (in) :: unit
!
      write(unit,NML=particles_selfgrav_run_pars)
!
    endsubroutine write_particles_selfg_run_pars
!***********************************************************************
    subroutine rprint_particles_selfgravity(lreset,lwrite)
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
!
      logical :: lwr
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Write information to index.pro
! 
      if (lwr) then
        write(3,*) 'igpotselfx=', igpotselfx
        write(3,*) 'igpotselfy=', igpotselfy
        write(3,*) 'igpotselfz=', igpotselfz
        write(3,*) 'irhop     =', irhop
      endif
!
    endsubroutine rprint_particles_selfgravity
!***********************************************************************

endmodule particles_selfgravity
