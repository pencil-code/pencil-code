! $Id: prints.f90,v 1.34 2003-01-23 10:52:16 nilshau Exp $

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
!  08-oct-02/tony: added safe_character_assign when appending to fform
!
      use Cdata
      use Sub
      use Hydro
      use Magnetic
      use General, only: safe_character_assign
!
      logical,save :: first=.true.
      character (len=320) :: fform,legend,line
      character (len=1) :: comma=','
      integer :: iname,index_d,index_a
!
!  If the timestep (=dt) is to be written, it is known only after
!  rk_2n, so the best place to enter it into the save list is here
!  Use 1.*(it-1) to have floating point or double prrecision.
!
      if (lroot) then
        if (i_t/=0)   call save_name(tdiagnos,i_t)
        if (i_dt/=0)  call save_name(dt,i_dt)
        if (i_it/=0)  call save_name(1.*(it-1),i_it)
        if (i_dtc/=0) call save_name(dt/(dxmin*cs0),i_dtc)
        if (lmagnetic) call calc_mfield
        if (lhydro)    call calc_mflow
!
!  produce the format
!  must set cform(1) explicitly, and then do iname>=2 in loop
!
        fform='('//cform(1)
        legend=noform(cname(1))
        do iname=2,nname
          call safe_character_assign(fform,  trim(fform)//comma//cform(iname))
          call safe_character_assign(legend, trim(legend)//noform(cname(iname)))
        enddo
        call safe_character_assign(fform, trim(fform)//')')
!
!! print*,'prints: form = ',trim(fform)
!! print*,'prints: args = ',fname(1:nname)
!
!  This treats all numbers as floating point numbers.
!  Only those numbers are given (and computed) that are
!  also listed in print.in.
!
        if(first) print*
        if(first) write(*,'(" ",A)') trim(legend)
!
!  write legend to extra file
!  (might want to do only once after each lreset)
!
        if(first) then
          open(1,file=trim(datadir)//'/legend.dat')
          write(1,'(" ",A)') trim(legend)
          close(1)
        endif
!
!  put output line into a string and remove spurious dots
!
        if(ldebug) print*,'bef. writing prints'
        write(line,trim(fform)) fname(1:nname)
        index_d=index(line,'. ')
        if (index_d >= 1) then
          line(index_d:index_d)=' '
        endif
!
!  if the line contains multiple asterisks then comment out line
!
        index_a=index(line,'***')
        if (index_a > 0) then
          line(1:1)='#'
        endif
!
!  append to diagnostics file
!
        open(1,file=trim(datadir)//'/time_series.dat',position='append')
        if(first) write(1,'("#",A)') trim(legend)
        !write(1,trim(fform)) fname(1:nname)  ! write to `time_series.dat'
        !write(6,trim(fform)) fname(1:nname)  ! write to standard output
        write(1,'(a)') trim(line)
        write(6,'(a)') trim(line)
        close(1)
!
      endif
!
      if(ldebug) print*,'exit prints'
      first = .false.
    endsubroutine prints
!***********************************************************************
    subroutine write_xyaverages
!
!  reads and registers print parameters gathered from the different
!  modules and marked in `xyaver.in'
!
!   6-jun-02/axel: coded
!
      logical,save :: first=.true.
!
      if(lroot.and.nnamez>0) then
        open(1,file=trim(datadir)//'/xyaverages.dat',position='append')
        write(1,'(1pe12.5)') t
        write(1,'(1p,8e10.3)') fnamez(:,:,1:nnamez)
        close(1)
      endif
      first = .false.
!
    endsubroutine write_xyaverages
!***********************************************************************
    subroutine write_zaverages
!
!  reads and registers print parameters gathered from the different
!  modules and marked in `zaver.in'
!
!  19-jun-02/axel: adapted from write_xyaverages
!
      logical,save :: first=.true.
!
      if(lroot.and.nnamexy>0) then
        open(1,file=trim(datadir)//'/zaverages.dat',position='append')
        write(1,'(1pe12.5)') t
        write(1,'(1p,8e10.3)') fnamexy(:,:,:,1:nnamexy)
        close(1)
      endif
      first = .false.
!
    endsubroutine write_zaverages
!***********************************************************************

endmodule Print
