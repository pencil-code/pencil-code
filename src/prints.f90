! $Id: prints.f90,v 1.3 2002-05-04 15:13:42 dobler Exp $

module Print

  use Cdata
  use Hydro
  use Magnetic

  implicit none

  contains

!***********************************************************************
    subroutine prints
!
!  reads and registers print parameters gathered from the different
!  modules and marked in `print.in'
!
!   3-may-02/axel: coded
!
      use Cdata
      use Sub
!
      logical,save :: first=.true.
!
!  this needs to be made more sophisticated of course...
!
      if(lroot) then
        if (first) write( 6,*) 'it,t,dt,',cname(1:nname)
        write( 6,*) it-1,t_diag,dt,fname(1:nname)
      endif
      first = .false.
!
    endsubroutine Prints
!***********************************************************************

endmodule Print
