! $Id: prints.f90,v 1.1 2002-05-04 09:11:59 brandenb Exp $

module Print

  use Cdata
  use Hydro
  use Magnetic

  implicit none

  contains

!***********************************************************************
    subroutine prints
!
!  reads and registers print parameters relevant for hydro part
!
!   3-may-02/axel: coded
!
      use Cdata
      use Sub
!
!  this needs to be made more sophisticated of course...
!
      if(lroot) then
        write( 6,*) 'it,t,',cname(1:nname)
        write( 6,*) it-1,t_diag,fname(1:nname)
        write( 6,*) 
      endif
!
    endsubroutine Prints
!***********************************************************************

endmodule Print
