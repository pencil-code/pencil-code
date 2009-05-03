! $Id$

module Print

  use Cdata
  use Hydro
  use Magnetic

  implicit none

  private

  public :: initialize_prints, prints
  public :: write_1daverages, write_2daverages
  public :: write_2daverages_prepare
  public :: write_zaverages

  character (len=5) :: ch2davg

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
!drcyl = dx   !this just happens for nrcyl=nx/2
        if (nrcyl/=0) then 
           drcyl = xyz1(1)/nrcyl
        else
           drcyl = 0. 
        endif
        rcyl = (/ ((i-0.5)*drcyl, i=1,nrcyl) /)
!
!  Calculate the three surface elements.
!  Take care of degenerate dimensions.
!
        if (nxgrid==1) then; dxeff=1.; else; dxeff=dx; endif
        if (nygrid==1) then; dyeff=1.; else; dyeff=dy; endif
        if (nzgrid==1) then; dzeff=1.; else; dzeff=dz; endif
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
!  08-oct-02/tony: added safe_character_assign when appending to fform
!   8-may-03/tarek: changed .or. to +. Could not compile with NAGf95 with .or.
!
      use Cdata
      use Sub
      use Hydro
      use Magnetic
      use Pscalar
      use Mpicomm
      use General, only: safe_character_append
!
      logical,save :: first=.true.
      character (len=640) :: fform,legend,line
      character (len=1) :: comma=','
      integer :: iname,index_d,index_a
!
!  Add general (not module-specific) quantities for diagnostic output.
!  If the timestep (=dt) is to be written, it is known only after
!  rk_2n, so the best place to enter it into the save list is here
!  Use 1.*(it-1) to have floating point or double precision.
!
      if (lroot) then
        if (idiag_t/=0)   call save_name(tdiagnos,idiag_t)
        if (idiag_dt/=0)  call save_name(dt,idiag_dt)
        if (idiag_it/=0)  call save_name(1.*(it-1),idiag_it)
      endif
!
!  calculate mean fields
!
      if (lmagnetic)  call calc_mfield
      if (lhydro)     call calc_mflow
      if (lpscalar)   call calc_mpscalar
!
      if (lroot) then
!
!  whenever itype_name=ilabel_max_dt, scale result by dt (for printing
!  Courant time). This trick is necessary, because dt is not known at the
!  time when the corresponding contribution to UUmax is known.
!
        do iname=1,nname
          if (itype_name(iname)==ilabel_max_dt) fname(iname)=dt*fname(iname)
        enddo
!
!  produce the format
!  must set cform(1) explicitly, and then do iname>=2 in loop
!
        fform = '(' // cform(1)
        legend=noform(cname(1))
        do iname=2,nname
          call safe_character_append(fform,  comma // cform(iname))
          call safe_character_append(legend, noform(cname(iname)))
        enddo
        call safe_character_append(fform, ')')

        if (ldebug) then
          write(0,*) 'PRINTS.prints: format = ', trim(fform)
          write(0,*) 'PRINTS.prints: args   = ', fname(1:nname)
        endif
!
!  This treats all numbers as floating point numbers.
!  Only those numbers are given (and computed) that are
!  also listed in print.in.
!
        if (first) write(*,*)
        if (first) write(*,'(" ",A)') trim(legend)
!
!  write legend to extra file
!  (might want to do only once after each lreset)
!
        if (first) then
          open(1,file=trim(datadir)//'/legend.dat')
          write(1,'(" ",A)') trim(legend)
          close(1)
        endif
!
!  put output line into a string and remove spurious dots
!
        if (ldebug) write(*,*) 'bef. writing prints'
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
          line(1:1)=comment_char
        endif
!
!  append to diagnostics file
!
        open(1,file=trim(datadir)//'/time_series.dat',position='append')
        if (first) write(1,"('"//comment_char//"',a)") trim(legend)
        !write(1,trim(fform)) fname(1:nname)  ! write to `time_series.dat'
        !write(6,trim(fform)) fname(1:nname)  ! write to standard output
        write(1,'(a)') trim(line)
        write(6,'(a)') trim(line)
        close(1)
!
      endif                     ! (lroot)
!
!  calculate mean emf (this requires i_... to be set in xyaver.in)
!
!     if (i_alp11z/=0) then
!       if (iproc==0) orms=fname(i_orms)
!
!     endif
!
!  calculate rhoccm and cc2m (this requires that these are set in print.in)
!  broadcast result to other processors. This is needed for calculating PDFs.
!
      if (idiag_rhoccm/=0) then
        if (iproc==0) rhoccm=fname(idiag_rhoccm)
        call mpibcast_real(rhoccm,1)
      endif
!
      if (idiag_cc2m/=0) then
        if (iproc==0) cc2m=fname(idiag_cc2m)
        call mpibcast_real(cc2m,1)
      endif
!
      if (idiag_gcc2m/=0) then
        if (iproc==0) gcc2m=fname(idiag_gcc2m)
        call mpibcast_real(gcc2m,1)
      endif
!
      if (ldebug) write(*,*) 'exit prints'
      first = .false.
!
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
      if (l1dout) then
        call write_xyaverages()
        call write_xzaverages()
        call write_yzaverages()
        call write_phizaverages()
      endif
!
    endsubroutine write_1daverages
!***********************************************************************
    subroutine write_2daverages_prepare()
!
!  Prepare l2davg for writing 2D averages.
!  This needs to be done in the beginning of each time step, so
!  the various routines know that they need to calculate averages.
!
!  23-nov-03/axel: adapted from write_2daverages and wvid_prepare
!
      use Param_IO
      use Sub, only: update_snaptime, read_snaptime
!
      real, save :: t2davg
      integer, save :: n2davg
      logical, save :: first=.true.
      character (len=135) :: file
!
      file=trim(datadir)//'/t2davg.dat'
      if (first) then
        call read_snaptime(trim(file),t2davg,n2davg,d2davg,t)
      endif
      first = .false.
!
!  This routine sets l2davg=T whenever its time to write 2D averages
!
      call update_snaptime(file,t2davg,n2davg,d2davg,t,l2davg,ch2davg,ENUM=.true.)
!
    endsubroutine write_2daverages_prepare
!***********************************************************************
    subroutine write_2daverages()
!
!  Write 2d averages (z-averages, phi-averages, .., i.e. quantities that
!  are still 2d after averaging) if it is time.
!  In analogy to 3d output to VARN, the time interval between two writes
!  is determined by a parameter (t2davg) in run.in.
!  Note: this routine should only be called by the root processor
!
!   7-aug-03/wolf: adapted from wsnap
!
      use Param_IO
!
      if (l2davg.and.lroot) then
        if (lwrite_yaverages)   call write_yaverages(ch2davg)
        if (lwrite_zaverages)   call write_zaverages(ch2davg)
        if (lwrite_phiaverages) call write_phiaverages(ch2davg)
        !
        if (ip<=10) write(*,*) 'write_2daverages: wrote phi(etc.)avgs'//ch2davg
      endif
!
!  Note: zaverages_xy are also needed if bmx and bmy are to be calculated
!  (Of course, yaverages_xz does not need to be calculated for that.)
!  AB: the following lines are still kept in equ.f90
!
!     if (.not.l2davgfirst.and.ldiagnos.and.ldiagnos_need_zaverages) then
!       if (lwrite_zaverages) call zaverages_xy
!     endif
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
      if (lroot.and.nnamez>0) then
        open(1,file=trim(datadir)//'/xyaverages.dat',position='append')
        write(1,'(1pe12.5)') t1ddiagnos
        write(1,'(1p,8e14.5e3)') fnamez(:,:,1:nnamez)
        close(1)
      endif
!
    endsubroutine write_xyaverages
!***********************************************************************
    subroutine write_xzaverages()
!
!  Write xz-averages (which are 1d data) that have been requested via
!  `xzaver.in'
!
!  12-oct-05/anders: adapted from write_xyaverages
!
      if (lroot.and.nnamey>0) then
        open(1,file=trim(datadir)//'/xzaverages.dat',position='append')
        write(1,'(1pe12.5)') t1ddiagnos
        write(1,'(1p,8e14.5e3)') fnamey(:,:,1:nnamey)
        close(1)
      endif
!
    endsubroutine write_xzaverages
!***********************************************************************
    subroutine write_yzaverages()
!
!  Write yz-averages (which are 1d data) that have been requested via
!  `yzaver.in'
!
!   2-oct-05/anders: adapted from write_xyaverages
!
      if (lroot.and.nnamex>0) then
        open(1,file=trim(datadir)//'/yzaverages.dat',position='append')
        write(1,'(1pe12.5)') t1ddiagnos
        write(1,'(1p,8e14.5e3)') fnamex(:,:,1:nnamex)
        close(1)
      endif
!
    endsubroutine write_yzaverages
!***********************************************************************
    subroutine write_phizaverages()
!
!  Write phiz-averages (which are 1d data) that have been requested via
!  `phizaver.in'
!
!  Also write rcyl to the output. It is needed just once, since it will
!  not change over the simulation. The condition "if (it==1)" is not the
!  best one, since it is reset upon restarting of a simulation and 
!  therefore one has to manually remove the extra(s) rcyl from 
!  the phizaverages.dat file
!
!  29-jan-07/wlad: adapted from write_yzaverages
!
      if (lroot.and.nnamer>0) then
        open(1,file=trim(datadir)//'/phizaverages.dat',position='append')
        if (it==1) write(1,'(1p,8e14.5e3)') rcyl
        write(1,'(1pe12.5)') t1ddiagnos
        write(1,'(1p,8e14.5e3)') fnamer(:,1:nnamer)
        close(1)
      endif
!
    endsubroutine write_phizaverages
!***********************************************************************
    subroutine write_yaverages(ch)
!
!  Write y-averages (which are 2d data) that have been requested via
!  `yaver.in'
!
!   7-jun-05/axel: adapted from write_zaverages
!
      character (len=4) :: ch
!
      if (lroot.and.nnamexz>0) then
        open(1, file=trim(datadir)//'/yaverages.dat', form='unformatted', &
            position='append')
        write(1) t2davgfirst
        write(1) fnamexz(:,:,:,1:nnamexz)
        close(1)
      endif
!
      if (NO_WARN) write(*,*) ch       ! (keep compiler quiet)
!
    endsubroutine write_yaverages
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
      if (lroot.and.nnamexy>0) then
        open(1, file=trim(datadir)//'/zaverages.dat', form='unformatted', &
            position='append')
        write(1) t2davgfirst
        write(1) fnamexy(:,:,:,1:nnamexy)
        close(1)
      endif
!
      if (NO_WARN) write(*,*) ch       ! (keep compiler quiet)
!
    endsubroutine write_zaverages
!***********************************************************************
    subroutine write_phiaverages(ch)
!
!  Write azimuthal averages (which are 2d data) that have been requested
!  via `phiaver.in'.
!  Note: fnamerz still has a third dimension indicating ipz, but the way
!  we are writing we automatically end up with the full z-direction
!  written contiguously.
!
!  File format:
!    1. nr_phiavg, nz_phiavg, nvars, nprocz
!    2. t, r_phiavg, z_phiavg, dr, dz
!    3. data
!    4. len(labels),labels
!
!   2-jan-03/wolf: adapted from write_zaverages
!
      use General, only: safe_character_assign,safe_character_append
!
      integer :: i
      character (len=4) :: ch
      character (len=80) :: avgdir,sname,fname
      character (len=1024) :: labels
!
!  write result; normalization is already done in phiaverages_rz
!
      if (lroot.and.nnamerz>0) then
        call safe_character_assign(avgdir, trim(datadir)//'/averages')
        call safe_character_assign(sname, 'PHIAVG'//trim(ch))
        call safe_character_assign(fname, trim(avgdir)//'/'//trim(sname))
        open(1,FILE=fname,FORM='unformatted')
        write(1) nrcyl,nzgrid,nnamerz,nprocz
        write(1) t2davgfirst,rcyl, &
                 z(n1)+(/(i*dz, i=0,nzgrid-1)/), &
                 drcyl,dz
        !ngrs: use pack to explicitly order the array before writing
        !     (write was messing up on copson without this...)
        write(1) pack(fnamerz(:,1:nz,:,1:nnamerz),.true.)
!
!  write labels at the end of file
!
        labels = trim(cnamerz(1))
        do i=2,nnamerz
          call safe_character_append(labels,",",trim(cnamerz(i)))
        enddo
        write(1) len(labels),labels
        close(1)
!
!  write file name to file list
!
        open(1,FILE=trim(avgdir)//'/phiavg.files',POSITION='append')
        write(1,'(A)') trim(sname)
        close(1)
!
      endif
!
    endsubroutine write_phiaverages
!***********************************************************************

endmodule Print
