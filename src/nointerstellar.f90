! $Id: nointerstellar.f90,v 1.7 2003-06-16 04:41:10 brandenb Exp $

!  This modules solves contains ISM and SNe 

module Interstellar

  use Cparam
  use Cdata

  implicit none

  integer, parameter :: ninterstellarsave=1
  real, dimension(ninterstellarsave) :: interstellarsave

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
      if (.not. first) call stop_it('register_nointerstellar called twice')
      first = .false.
!
      linterstellar = .false.
!
!      if ((ip<=8) .and. lroot) then
!        print*, 'Register_nointerstellar'
!      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: nointerstellar.f90,v 1.7 2003-06-16 04:41:10 brandenb Exp $")
!
!      if (nvar > mvar) then
!        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
!        call stop_it('Register_nointerstellar: nvar > mvar')
!      endif
!
    endsubroutine register_interstellar
!***********************************************************************
    subroutine initialize_interstellar(lstart)
!
!  Perform any post-parameter-read initialization eg. set derived 
!  parameters
!
!  24-nov-02/tony: coded - dummy
!
      logical :: lstart
!
! (to keep compiler quiet)
      if (ip==0) print*,lstart
!
    endsubroutine initialize_interstellar
!***********************************************************************
    subroutine calc_heat_cool_interstellar(df,rho1,TT1)
!
!  adapted from calc_heat_cool
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: rho1,TT1
!
! (to keep compiler quiet)
      if (ip==0) print*,'calc_heat_cool_interstellar', &
                  df(1,1,1,1),rho1(1),TT1(1)
!
    endsubroutine calc_heat_cool_interstellar
!***********************************************************************
    subroutine check_SN(f)
!
!  dummy routine for checking for SNe (interstellar)
!
    use Cdata
!
    real, dimension(mx,my,mz,mvar+maux) :: f
!
! (to keep compiler quiet)
      if (ip==0) print*,'SN check',f(1,1,1,1)
!
    endsubroutine check_SN

endmodule interstellar
