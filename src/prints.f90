! $Id: prints.f90,v 1.12 2002-05-29 07:09:06 brandenb Exp $

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
!  27-may-02/axel: it,t,dt added as extra save parameters
!
      use Cdata
      use Sub
      use Hydro
!
      logical,save :: first=.true.
      character (len=320) :: fform,legend
      character (len=1) :: comma=','
      integer :: iname
!
!  If the timestep (=dt) is to be outputted, it is known only after
!  rk_2n, so the best place to enter it into the save list is here
!
      if (i_t/=0) call save_name(tdiagnos,i_t)
      if (i_dt/=0) call save_name(dt,i_dt)
      if (i_it/=0) call save_name(float(it-1),i_it)
!
!  produce the format
!  must set cform(1) explicitly, and then do iname>=2 in loop
!
      fform='('//cform(1)
      legend=cname(1)
      do iname=2,nname
        fform=trim(fform)//comma//cform(iname)
        legend=trim(legend)//comma//cname(iname)
      enddo
      fform=trim(fform)//')'
!
!! print*,'prints: form = ',trim(fform)
!! print*,'prints: args = ',fname(1:nname)
!
!  This treats all numbers as floating point numbers.
!  Only those numbers are given (and computed) that are
!  also listed in print.in.
!
      if(lroot) then
        open(3,file='check')
        write(3,trim(fform))  fname(1:nname)  ! write to `check'
        write(6,trim(fform))  fname(1:nname)  ! write to standard output
        write(20,trim(fform)) fname(1:nname)  ! write to `fort.20'
        close(3)
      endif
      first = .false.
!
    endsubroutine Prints
!***********************************************************************

endmodule Print
