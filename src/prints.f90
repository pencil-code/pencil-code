! $Id: prints.f90,v 1.9 2002-05-19 22:45:59 brandenb Exp $

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
      fform='(i10,f10.3,1pg10.3,0p,'//cform(1)
      do iname=2,nname
        fform=trim(fform)//comma//cform(iname)
      enddo
      fform=trim(fform)//')'
!print*,'PRINTS: form = ',trim(fform)
!print*,'PRINTS: args = ',it-1,t_diag,dt,fname(1:nname)
!
!  this needs to be made more sophisticated of course...
!
      if(lroot) then
        open(3,file='check')
        write(3,trim(fform))  it-1,t_diag,dt,fname(1:nname)
        write(6,trim(fform))  it-1,t_diag,dt,fname(1:nname)
        write(20,trim(fform)) it-1,t_diag,dt,fname(1:nname)
        close(3)
      endif
      first = .false.
!
    endsubroutine Prints
!***********************************************************************

endmodule Print
