! $Id: prints.f90,v 1.4 2002-05-04 15:49:26 brandenb Exp $

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
      character (len=320) :: fform
      character (len=1) :: comma=','
      integer :: iname
!
!  produce the format
!
      fform='(i10,f10.3,'//cform(1)
      do iname=2,nname
        fform=trim(fform)//comma//cform(iname)
      enddo
      fform=trim(fform)//')'
print*,fform
!
!  this needs to be made more sophisticated of course...
!
      if(lroot) then
        write(6,fform) it-1,t_diag,dt,fname(1:nname)
      endif
      first = .false.
!
    endsubroutine Prints
!***********************************************************************

endmodule Print
