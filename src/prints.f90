! $Id: prints.f90,v 1.2 2002-05-04 12:40:13 dobler Exp $

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
!  this needs to be made more sophisticated of course...
!
      if(lroot) then
        write( 6,*) 'it,t,dt,',cname(1:nname)
        write( 6,*) it-1,t_diag,dt,fname(1:nname)
        write( 6,*) 
      endif
!
    endsubroutine Prints
!***********************************************************************

endmodule Print
