! $Id: interstellar.f90,v 1.1 2002-11-19 14:12:22 mee Exp $

!  This modules solves contains ISM and SNe 

module Interstellar

  use Cparam
  use Cdata

  implicit none

  ! input parameters
  integer :: dummy
  namelist /interstellar_init_pars/ dummy

  ! run parameters
  namelist /interstellar_run_pars/ dummy

 
  contains

!***********************************************************************
    subroutine register_interstellar()
!
!  19-nov-02/tony: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_interstellar called twice')
      first = .false.
!
      linterstellar = .true.
!
      if ((ip<=8) .and. lroot) then
        print*, 'Register_lncc'
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: interstellar.f90,v 1.1 2002-11-19 14:12:22 mee Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('Register_interstellar: nvar > mvar')
      endif
!
    endsubroutine register_interstellar


endmodule interstellar



