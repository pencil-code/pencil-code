! $Id: hydro.f90,v 1.4 2002-05-11 12:18:48 dobler Exp $

module Hydro

  use Cparam

  implicit none

  integer :: i_t=0,i_urms=0,i_umax=0

  contains

!***********************************************************************
    subroutine register_hydro
!
!  Initialise variables which should know that we solve the hydro
!  equations: iuu, etc; increase nvar accordingly.
!
!  6-nov-01/wolf: coded
!
      use Cdata
      use Mpicomm, only: lroot,stop_it
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_hydro called twice')
      first = .false.
!
      iuu = nvar+1             ! indices to access uu and lam
      iux = iuu
      iuy = iuu+1
      iuz = iuu+2
      ilnrho = iuu+3
      nvar = nvar+4             ! added 4 variables
!
      if ((ip<=8) .and. lroot) then
        print*, 'Register_hydro:  nvar = ', nvar
        print*, 'iuu,ilnrho = ', iuu,ilnrho
        print*, 'iux,iuy,iuz = ', iux,iuy,iuz
      endif
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$RCSfile: hydro.f90,v $", &
           "$Revision: 1.4 $", &
           "$Date: 2002-05-11 12:18:48 $")
!
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('Register_hydro: nvar > mvar')
      endif
!
    endsubroutine register_hydro
!***********************************************************************
    subroutine rprint_hydro
!
!  reads and registers print parameters relevant for hydro part
!
!   3-may-02/axel: coded
!
      use Cdata
      use Sub
!
      integer :: iname
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'t',i_t)
        call parse_name(iname,cname(iname),cform(iname),'urms',i_urms)
        call parse_name(iname,cname(iname),cform(iname),'umax',i_umax)
      enddo
!
    endsubroutine rprint_hydro
!***********************************************************************

endmodule Hydro
