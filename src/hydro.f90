! $Id: hydro.f90,v 1.3 2002-05-04 15:08:29 dobler Exp $

module Hydro

  use Cparam

  implicit none

  integer :: i_t=0,i_urms=0,i_umax=0

  contains

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
