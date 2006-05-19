! $Id: selfgravity.f90,v 1.3 2006-05-19 20:40:03 joishi Exp $

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

  real :: rhs_const=1.0
  
  namelist /selfgrav_init_pars/ &
      rhs_const
      
  namelist /selfgrav_run_pars/ &
      rhs_const

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
           "$Id: selfgravity.f90,v 1.3 2006-05-19 20:40:03 joishi Exp $")
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
      use Poisson
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      call poisson_solver_fft(f,ilnrho,ipotself,rhs_const)
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
!  Add self gravity acceleration on the gas.
!
      if (lhydro) df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) - p%gpotself
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
