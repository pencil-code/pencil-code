! $Id: prints.f90,v 1.19 2002-06-15 11:29:35 brandenb Exp $

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
!   7-jun-02/axel: dtc (=dt/cdt) added as extra save parameter
!
      use Cdata
      use Sub
      use Hydro
      use Magnetic
!
      logical,save :: first=.true.
      character (len=320) :: fform,legend
      character (len=1) :: comma=','
      integer :: iname
      real :: bmz
!
!  If the timestep (=dt) is to be outputted, it is known only after
!  rk_2n, so the best place to enter it into the save list is here
!
      if (i_t/=0) call save_name(tdiagnos,i_t)
      if (i_dt/=0) call save_name(dt,i_dt)
      if (i_it/=0) call save_name(float(it-1),i_it)
      if (i_dtc/=0) call save_name(dt/(dxmin*cs0),i_dtc)
!
!  Magnetic energy in horizontally averaged field
!  The bxmz and bymz must have been calculated,
!  so they are present on the root processor.
!
      if (lmagnetic.and.i_bmz/=0.and.lroot) then
        if(i_bxmz==0.or.i_bymz==0) then
          if(first) print*
          if(first) print*,"NOTE: to get bmz, bxmz and bymz must also be set in xyaver"
          if(first) print*,"      This may be because we renamed zaver.in into xyaver.in"
          if(first) print*,"      We proceed, but you'll get bmz=0"
        endif
        bmz=sqrt(sum(fnamez(:,i_bxmz)**2+fnamez(:,i_bymz)**2)/nz)
        call save_name(bmz,i_bmz)
      endif
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
        if(first) print*
!
!  write legend to extra file
!  (might want to do only once after each lreset)
!
        open(1,file='tmp/legend.dat')
        write(1,*) legend
        close(1)
!
!  append to diagnostics file
!
        open(1,file='tmp/n.dat',position='append')
        write(1,trim(fform)) fname(1:nname)  ! write to `n.dat'
        write(6,trim(fform)) fname(1:nname)  ! write to standard output
        close(1)
!
      endif
      first = .false.
!
    endsubroutine Prints
!***********************************************************************
    subroutine write_xyaverages
!
!  reads and registers print parameters gathered from the different
!  modules and marked in `print.in'
!
!   6-jun-02/axel: coded
!
      use Cdata
      use Sub
      use Hydro
!
      logical,save :: first=.true.
!
      if(lroot.and.nnamez>0) then
        open(1,file='tmp/xyaverages.dat',position='append')
        write(1,'(1pe12.5)') t
        write(1,'(1p,8e10.3)') fnamez(:,1:nnamez)
        close(1)
      endif
      first = .false.
!
    endsubroutine write_xyaverages
!***********************************************************************

endmodule Print
