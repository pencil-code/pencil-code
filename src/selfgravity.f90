! $Id: selfgravity.f90,v 1.5 2006-06-14 23:57:00 ajohan Exp $

!
!  This module takes care of self gravity by solving the Poisson equation
!    (d^2/dx^2 + d^2/dy^2 + d^2/dz^2)phi = 4*pi*G*rho
!  for the potential phi.
!

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 1
! COMMUNICATED AUXILIARIES 1
!
! PENCILS PROVIDED gpotself
!
!***************************************************************

module Selfgravity

  use Cdata
  use Cparam
  use Messages

  implicit none

  include 'selfgravity.h'

  real :: rhs_poisson_const=1.0
  logical :: lselfgravity_gas=.true., lselfgravity_dust=.false.
  
  namelist /selfgrav_init_pars/ &
      rhs_poisson_const, lselfgravity_gas, lselfgravity_dust
      
  namelist /selfgrav_run_pars/ &
      rhs_poisson_const, lselfgravity_gas, lselfgravity_dust

  contains

!***********************************************************************
    subroutine register_selfgravity()
!
!  Initialise self gravity variables.
!
!  15-may-06/anders+jeff: adapted
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_selfgravity: called twice')
      first = .false.
!
!  Set indices for auxiliary variables
!
      ipotself = mvar + naux_com + 1
      naux = naux + 1
      naux_com = naux_com + 1
!
!  Identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: selfgravity.f90,v 1.5 2006-06-14 23:57:00 ajohan Exp $")
!
!  Put variable name in array
!
      varname(ipotself) = 'potself'
!
!  Set lselfgravity.
!
      lselfgravity=.true.
!
    endsubroutine register_selfgravity
!***********************************************************************
    subroutine initialize_selfgravity()
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  15-may-06/anders+jeff: adapted
!
      if (.not.lpoisson) then
        if (lroot) print*, 'initialize_selfgravity: must choose a Poisson '// &
            'solver in Makefile.local for self-gravity'
        call fatal_error('initialize_selfgravity','')
      endif
!
      if (lselfgravity_gas.and..not.(lhydro.and.ldensity)) then
        if (lroot) print*, 'initialize_selfgravity: must choose a hydro '// &
            'and a density module in Makefile.local for self-gravity'
        call fatal_error('initialize_selfgravity','')
      endif
!
      if (lselfgravity_dust.and..not.(ldustvelocity.and.ldustdensity)) then
        if (lroot) then
          print*, 'initialize_selfgravity: must choose a dust velocity '// &
              'and a dust density module in Makefile.local for '// &
              'self-gravity on the dust fluid'
        endif
        call fatal_error('initialize_selfgravity','')
      endif
!
    endsubroutine initialize_selfgravity
!***********************************************************************
    subroutine pencil_criteria_selfgravity()
! 
!  All pencils that the Selfgravity module depends on are specified here.
! 
!  15-may-06/anders+jeff: adapted
!
      lpenc_requested(i_gpotself)=.true.
!
    endsubroutine pencil_criteria_selfgravity
!***********************************************************************
    subroutine pencil_interdep_selfgravity(lpencil_in)
!
!  Interdependency among pencils from the Selfgravity module is specified here.
!
!  15-may-06/anders+jeff: adapted
!
      logical, dimension(npencils) :: lpencil_in
!
      if (NO_WARN) print*, lpencil_in !(keep compiler quiet)
!
    endsubroutine pencil_interdep_selfgravity
!***********************************************************************
    subroutine calc_pencils_selfgravity(f,p)
!
!  Calculate Selfgravity pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  15-may-06/anders+jeff: coded
!
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      type (pencil_case) :: p
!      
      intent(in) :: f
      intent(inout) :: p
!
      if (lpencil(i_gpotself)) call grad(f,ipotself,p%gpotself)
!
    endsubroutine calc_pencils_selfgravity
!***********************************************************************
    subroutine calc_selfpotential(f)
!
!  Calculate the potential of the self gravity.
!
!  15-may-06/anders+jeff: coded
!
      use Particles_main, only: particles_calc_selfpotential
      use Poisson
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      real, dimension (nx,ny,nz) :: rhs_poisson
!
!  Consider self-gravity from gas and dust density or from either one.
!
      if (lhydro.and.ldensity.and.lselfgravity_gas) then
        if (ldensity_nolog) then
          rhs_poisson=rhs_poisson_const*f(l1:l2,m1:m2,n1:n2,ilnrho)
        else
          rhs_poisson=rhs_poisson_const*exp(f(l1:l2,m1:m2,n1:n2,ilnrho))
        endif
      endif
!  Dust.
      if (ldustdensity.and.ldustvelocity.and.lselfgravity_dust) then
        if (ldustdensity_log) then
          if (lselfgravity_gas) then  ! No need to zero rhs.
            rhs_poisson = rhs_poisson + &
                rhs_poisson_const*exp(f(l1:l2,m1:m2,n1:n2,ind(1)))
          else                        ! Must zero rhs.
            rhs_poisson = rhs_poisson_const*exp(f(l1:l2,m1:m2,n1:n2,ind(1)))
          endif
        else
          if (lselfgravity_gas) then  ! No need to zero rhs.
            rhs_poisson = rhs_poisson + &
                rhs_poisson_const*f(l1:l2,m1:m2,n1:n2,ind(1))
          else                        ! Must zero rhs.
            rhs_poisson = rhs_poisson_const*f(l1:l2,m1:m2,n1:n2,ind(1))
          endif
        endif
      endif
!
!  Contribution from particles is taken care of by the particle modules.
!
      if (lparticles) &
          call particles_calc_selfpotential(f,rhs_poisson,rhs_poisson_const, &
          lselfgravity_gas.or.lselfgravity_dust)
!
!  Send the right-hand-side of the Poisson equation to the Poisson solver and
!  receive the self-gravity potential back.
!
      call poisson_solver_fft(rhs_poisson)
!
!  Put potential into f array.
!
      f(l1:l2,m1:m2,n1:n2,ipotself) = rhs_poisson
!
    endsubroutine calc_selfpotential
!***********************************************************************
    subroutine duu_dt_selfgrav(f,df,p)
!
!  Add self gravity acceleration on gas.
!
!  15-may-06/anders+jeff: coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      integer :: k
!
      intent(in) :: f,p
      intent(out) :: df
!
!  Add self-gravity acceleration on the gas and on the dust.
!
      if (lhydro) df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) - p%gpotself
      if (ldustvelocity) df(l1:l2,m,n,iudx(1):iudz(1)) = &
          df(l1:l2,m,n,iudx(1):iudz(1)) - p%gpotself
!
      if (NO_WARN) print*, f, p !(keep compiler quiet)
!        
    endsubroutine duu_dt_selfgrav
!***********************************************************************
    subroutine read_selfgravity_init_pars(unit,iostat)
!
!  Read self gravity init parameters
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=selfgrav_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=selfgrav_init_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_selfgravity_init_pars
!***********************************************************************
    subroutine write_selfgravity_init_pars(unit)
!
!  Write self gravity init parameters
!
      integer, intent(in) :: unit
!
      write(unit,NML=selfgrav_init_pars)
!
    endsubroutine write_selfgravity_init_pars
!***********************************************************************
    subroutine read_selfgravity_run_pars(unit,iostat)
!
!  Read self gravity run parameters
!    
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=selfgrav_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=selfgrav_run_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_selfgravity_run_pars
!***********************************************************************
    subroutine write_selfgravity_run_pars(unit)
!
!  Write self gravity run parameters
!
      integer, intent(in) :: unit
!
      write(unit,NML=selfgrav_run_pars)
!
    endsubroutine write_selfgravity_run_pars
!***********************************************************************
    subroutine rprint_selfgravity(lreset,lwrite)
!
!  reads and registers print parameters relevant for gravity advance
!  dummy routine
!
!  16-may-06/anders+jeff: adapted
!
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
      if(NO_WARN) print*, lreset  !(to keep compiler quiet)
!
!  write column where which variable is stored
!
      if (lwr) then
        write(3,*) 'ipotself=', ipotself
      endif
!
    endsubroutine rprint_selfgravity
!***********************************************************************


endmodule Selfgravity
