! $Id: prints.f90,v 1.51 2003-10-12 05:45:16 brandenb Exp $

module Print

  use Cdata
  use Hydro
  use Magnetic

  implicit none

  contains

!***********************************************************************
    subroutine initialize_prints()
!
!  setup variables needed for output of diagnostic quantities and
!  averages.
!
!  14-aug-03/axel: added dxy, dyz, and dxz
!
      use Cdata
!
      integer :: i
      real :: dxeff,dyeff,dzeff
      logical, save :: first=.true.
!
      if (first) then
        !
        !  Initialize rcyl for the phi-averages grid. Does not need to be
        !  done after each reload of run.in, but this is the easiest way
        !  of doing it.
        !
        drcyl = dx
        rcyl = (/ ((i-0.5)*drcyl, i=1,nrcyl) /)
        !
        !  Calculate the three surface elements.
        !  Take care of degenerate dimensions.
        !
        if(nxgrid==1) then; dxeff=1.; else; dxeff=dx; endif
        if(nygrid==1) then; dyeff=1.; else; dyeff=dy; endif
        if(nzgrid==1) then; dzeff=1.; else; dzeff=dz; endif
        !
        dsurfxy=dxeff*dyeff
        dsurfyz=dyeff*dzeff
        dsurfzx=dzeff*dxeff
        !
        !  Calculate the volume element.
        !
        dvol=dxeff*dyeff*dzeff
      endif
!
      first = .false.
!
    endsubroutine initialize_prints
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
!   8-may-03/tarek: changed .or. to +. Could not compile with NAGf95 with .or.
!
      use Cdata
      use Sub
      use Hydro
      use Magnetic
      use Pscalar
      use Mpicomm
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
        if (lpscalar)  call calc_mpscalar
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
!  if the line contains unreadible characters, then comment out line
!
        index_a=(index(line,'***') +  index(line,'???'))
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
!  calculate brms (this requires that brms is set in print.in)
!  broadcast result to other processors
!
      if (i_brms/=0) then
        if (iproc==0) brms=fname(i_brms)
        call mpibcast_real(brms,1)
      endif
!
      if(ldebug) print*,'exit prints'
      first = .false.
    endsubroutine prints
!***********************************************************************
    subroutine write_1daverages()
!
!  Write 1d averages (z-averages, .., i.e. quantities that are only  1d
!  after averaging). These are written every it1 timesteps (like the
!  diagnostic line) and appended to their individual files.
!
!   7-aug-03/wolf: coded
!
      if (lout) call write_xyaverages()
!
    endsubroutine write_1daverages
!***********************************************************************
    subroutine write_2daverages()
!
!  Write 2d averages (z-averages, phi-averages, .., i.e. quantities that
!  are still 2d after averaging) if it is time.
!  In analogy to 3d output to VARN, the time interval between two writes
!  is determined by a parameter (t2davg) in run.in.
!
!   7-aug-03/wolf: adapted from wsnap
!
      use Param_IO
!
      real, save :: t2davg
      integer, save :: n2davg
      logical, save :: first=.true.
      logical :: lnow
      character (len=4) :: ch
      character (len=135) :: file
!
      file=trim(datadir)//'/t2davg.dat'
      if (first) then
        call read_snaptime(trim(file),t2davg,n2davg,d2davg,t)
      endif
      first = .false.
      !
      call update_snaptime(file,t2davg,n2davg,d2davg,t,lnow,ch,ENUM=.true.)
      if (lnow) then
        if (lwrite_zaverages)   call write_zaverages(ch)
        if (lwrite_phiaverages) call write_phiaverages(ch)
        !
        if(ip<=10.and.lroot) print*, 'write_2daverages: wrote phi(etc.)avgs'//ch
      endif
!
    endsubroutine write_2daverages
!***********************************************************************
    subroutine write_xyaverages()
!
!  Write xy-averages (which are 1d data) that have been requested via
!  `xyaver.in'
!
!   6-jun-02/axel: coded
!
      if(lroot.and.nnamez>0) then
        open(1,file=trim(datadir)//'/xyaverages.dat',position='append')
        write(1,'(1pe12.5)') t
        write(1,'(1p,8e10.3)') fnamez(:,:,1:nnamez)
        close(1)
      endif
!
    endsubroutine write_xyaverages
!***********************************************************************
    subroutine write_zaverages(ch)
!
!  Write z-averages (which are 2d data) that have been requested via
!  `zaver.in'
!
!  19-jun-02/axel: adapted from write_xyaverages
!
      character (len=4) :: ch
!
      if(lroot.and.nnamexy>0) then
        open(1,file=trim(datadir)//'/zaverages.dat',position='append')
        write(1,'(1pe12.5)') t
        write(1,'(1p,8e10.3)') fnamexy(:,:,:,1:nnamexy)
        close(1)
      endif
!
      if(ip==0) print*,ch       ! (keep compiler quiet)
!
    endsubroutine write_zaverages
!***********************************************************************
    subroutine write_phiaverages(ch)
!
!  Write azimuthal averages (which are 2d data) that have been requested
!  via `phiaver.in'
!  File format:
!     1. data
!     2. t, r_phiavg, z_phiavg, dr, dz
!     3. nr_phiavg, nz_phiavg, nvars
!     4. labels
!
!   2-jan-03/wolf: adapted from write_zaverages
!
      use General
!
      integer :: i
      character (len=4) :: ch
      character (len=80) :: fname
      character (len=1024) :: labels
 !
      if(lroot.and.nnamerz>0) then
        call safe_character_assign(fname, &
                                   trim(datadir)//'/averages/PHIAVG'//trim(ch))
        open(1,FILE=fname,FORM='unformatted')
        write(1) nrcyl,n2-n1+1,nnamerz ! sizes (just in case) 
        write(1) t,rcyl,z(n1:n2),drcyl,dz
        write(1) fnamerz(:,1:,:,1:nnamerz) / spread(fnamerz(:,0,:,1:1),2,nz)
        labels = trim(cnamerz(1))
        do i=2,nnamerz
          call safe_character_append(labels,",",trim(cnamerz(i)))
        enddo
        write(1) len(labels),labels
        close(1)
      endif
!
    endsubroutine write_phiaverages
!***********************************************************************

endmodule Print
