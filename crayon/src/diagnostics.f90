! $Id: equ.f90 10533 2009-03-26 11:01:45Z ajohan@strw.leidenuniv.nl $
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
!***************************************************************
module Diagnostics
!
  use Cdata
  use Messages
  use Mpicomm
!
  implicit none
!
  private
!
  public :: prints
  public :: diagnostic
  public :: xyaverages_z, xzaverages_y, yzaverages_x
  public :: yaverages_xz, zaverages_xy
  public :: write_1daverages, write_2daverages
  public :: write_2daverages_prepare, write_zaverages
  public :: parse_name, fparse_name, save_name, save_name_halfz
  public :: max_name, sum_name
  public :: max_mn_name, sum_mn_name, integrate_mn_name, sum_weighted_name
  public :: integrate_mn
  public :: sum_mn_name_halfy, surf_mn_name
  public :: sum_mn_name_halfz
  public :: xysum_mn_name_z, xzsum_mn_name_y, yzsum_mn_name_x
  public :: ysum_mn_name_xz, zsum_mn_name_xy
  public :: yzintegrate_mn_name_x, xzintegrate_mn_name_y, xyintegrate_mn_name_z
  public :: allocate_xyaverages, allocate_xzaverages, allocate_yzaverages
  public :: allocate_yaverages, allocate_zaverages
  public :: xyaverages_clean_up, xzaverages_clean_up, yzaverages_clean_up
  public :: yaverages_clean_up, zaverages_clean_up
  public :: get_from_fname
  public :: init_xaver
!
  integer :: mnamer
  character (len=5) :: ch2davg
!
  contains
!***********************************************************************
    subroutine prints
!
!  Reads and registers print parameters gathered from the different
!  modules and marked in `print.in'.
!
!   3-may-02/axel: coded
!
      use General, only: safe_character_append
      use Sub, only: noform
!
      logical,save :: first=.true.
      character (len=640) :: fform,legend,line
      character (len=1), parameter :: comma=','
      integer :: iname,index_d,index_a
!
!  Add general (not module-specific) quantities for diagnostic output. If the
!  timestep (=dt) is to be written, it is known only after rk_2n, so the best
!  place to enter it into the save list is here. Use 1.0*(it-1) to have floating
!  point or double precision.
!
      if (lroot) then
        if (idiag_t/=0)   call save_name(tdiagnos,idiag_t)
        if (idiag_dt/=0)  call save_name(dt,idiag_dt)
        if (idiag_it/=0)  call save_name(one_real*(it-1),idiag_it)
      endif
!
      if (lroot) then
!
!  Whenever itype_name=ilabel_max_dt, scale result by dt (for printing Courant
!  time).
!
        do iname=1,nname
          if (itype_name(iname)==ilabel_max_dt) fname(iname)=dt*fname(iname)
        enddo
!
!  Produce the format.
!  Must set cform(1) explicitly, and then do iname>=2 in loop.
!
        fform = '(' // cform(1)
        legend=noform(cname(1))
        do iname=2,nname
          call safe_character_append(fform,  comma // cform(iname))
          call safe_character_append(legend, noform(cname(iname)))
        enddo
        call safe_character_append(fform, ')')
!
        if (ldebug) then
          write(0,*) 'PRINTS.prints: format = ', trim(fform)
          write(0,*) 'PRINTS.prints: args   = ', fname(1:nname)
        endif
!
!  This treats all numbers as floating point numbers.  Only those numbers are
!  given (and computed) that are also listed in print.in.
!
        if (first) write(*,*)
        if (first) write(*,'(" ",A)') trim(legend)
!
!  Write legend to extra file (might want to do only once after each lreset)
!
        if (first) then
          open(1,file=trim(datadir)//'/legend.dat')
          write(1,'(" ",A)') trim(legend)
          close(1)
        endif
!
!  Put output line into a string and remove spurious dots.
!
        if (ldebug) write(*,*) 'bef. writing prints'
        write(line,trim(fform)) fname(1:nname)
        index_d=index(line,'. ')
        if (index_d >= 1) then
          line(index_d:index_d)=' '
        endif
!
!  If the line contains unreadable characters, then comment out line.
!
        index_a=(index(line,'***') +  index(line,'???'))
        if (index_a > 0) then
          line(1:1)=comment_char
        endif
!
!  Append to diagnostics file.
!
        open(1,file=trim(datadir)//'/time_series.dat',position='append')
        if (first) write(1,"('"//comment_char//"',a)") trim(legend)
        write(1,'(a)') trim(line)
        write(6,'(a)') trim(line)
        close(1)
!
      endif                     ! (lroot)
!
      if (ldebug) write(*,*) 'exit prints'
      first = .false.
!
      fname(1:nname)=0.0
!
    endsubroutine prints
!***********************************************************************
    subroutine diagnostic
!
!  Finalize calculation of diagnostic quantities (0-D).
!
!   2-sep-01/axel: coded
!  14-aug-03/axel: began adding surface integrals
!
      real, dimension (mname) :: fmax_tmp, fsum_tmp, fmax, fsum, fweight_tmp
      real :: vol
      integer :: iname,imax_count,isum_count,nmax_count,nsum_count
      logical :: lweight_comm
      logical, save :: first=.true.
      real, save :: dVol_rel1
!
!  Calculate relative volume integral.
!
      if (first) then
        dVol_rel1=1./(nw*ncpus)
        first=.false.
        if (lroot.and.ip<=10) print*,'dVol_rel1=',dVol_rel1
      endif
!
!  Go through all print names, and sort into communicators
!  corresponding to their type.
!
      imax_count=0
      isum_count=0
      lweight_comm=.false.
      do iname=1,nname
        if (itype_name(iname)<0) then
          imax_count=imax_count+1
          fmax_tmp(imax_count)=fname(iname)
        elseif (itype_name(iname)>0) then
          isum_count=isum_count+1
          fsum_tmp(isum_count)=fname(iname)
          if (itype_name(iname)==ilabel_sum_weighted .or. &
              itype_name(iname)==ilabel_sum_weighted_sqrt .or. &
              itype_name(iname)==ilabel_sum_par) then
            fweight_tmp(isum_count)=fweight(iname)
            lweight_comm=.true.
          endif
        endif
      enddo
      nmax_count=imax_count
      nsum_count=isum_count
!
!  Communicate over all processors.
!
      call mpireduce_max(fmax_tmp,fmax,nmax_count)
      call mpireduce_sum(fsum_tmp,fsum,nsum_count)
      if (lweight_comm) call mpireduce_sum(fweight_tmp,fweight,nsum_count)
!
!  The result is present only on the root processor.
!
      if (lroot) then
!
!  Sort back into original array.
!
        imax_count=0
        isum_count=0
        do iname=1,nname
          if (itype_name(iname)<0) then ! max
            imax_count=imax_count+1
!
            if (itype_name(iname)==ilabel_max)            &
                fname(iname)=fmax(imax_count)
!
            if (itype_name(iname)==ilabel_max_sqrt)       &
                fname(iname)=sqrt(fmax(imax_count))
!
            if (itype_name(iname)==ilabel_max_dt)         &
                fname(iname)=fmax(imax_count)
!
            if (itype_name(iname)==ilabel_max_neg)        &
                fname(iname)=-fmax(imax_count)
!
            if (itype_name(iname)==ilabel_max_reciprocal) &
                fname(iname)=1./fmax(imax_count)
!
          elseif (itype_name(iname)>0) then 
            isum_count=isum_count+1
!
            if (itype_name(iname)==ilabel_sum)            &
                fname(iname)=fsum(isum_count)*dVol_rel1
!
            if (itype_name(iname)==ilabel_sum_sqrt)       &
                fname(iname)=sqrt(fsum(isum_count)*dVol_rel1)
!
            if (itype_name(iname)==ilabel_sum_par)        &
                fname(iname)=fsum(isum_count)/fweight(isum_count)
!
            if (itype_name(iname)==ilabel_integrate)            &
                fname(iname)=fsum(isum_count)
!
             if (itype_name(iname)==ilabel_surf)          &
                 fname(iname)=fsum(isum_count)
!
             if (itype_name(iname)==ilabel_sum_lim) then
                vol=1.
                if (nzgrid/=1)           vol=vol*Lz
                fname(iname)=fsum(isum_count)/vol
             endif
!
            if (itype_name(iname)==ilabel_sum_weighted) then
              if (fweight(isum_count)/=0.0) then
                fname(iname)=fsum(isum_count)/fweight(isum_count)
              else
                fname(iname)=0.0
              endif
            endif
!
            if (itype_name(iname)==ilabel_sum_weighted_sqrt) then
              if (fweight(isum_count)/=0.0) then
                fname(iname)=sqrt(fsum(isum_count)/fweight(isum_count))
              else
                fname(iname)=0.0
              endif
            endif
!
          endif
!
        enddo
!
      endif
!
    endsubroutine diagnostic
!***********************************************************************
    subroutine xyaverages_z()
!
!  Calculate xy-averages (still depending on z)
!  NOTE: these averages depend on z, so after summation in x and y they
!  are still distributed over nprocz CPUs; hence the dimensions of fsumz
!  (and fnamez).
!  In other words: the whole xy-average is present in one and the same fsumz,
!  but the result is not complete on any of the processors before
!  mpireduce_sum has been called. This is simpler than collecting results
!  first in processors with the same ipz and different ipy, and then
!  assemble result from the subset of ipz processors which have ipy=0
!  back on the root processor.
!
!   6-jun-02/axel: coded
!
      real, dimension (nz,nprocz,nnamez) :: fsumz
!
!  Communicate over all processors.
!  The result is only present on the root processor
!
      if (nnamez>0) then
        call mpireduce_sum(fnamez,fsumz,(/nz,nprocz,nnamez/))
        if (lroot) &
            fnamez(:,:,1:nnamez)=fsumz(:,:,1:nnamez)/(nx*ny*nprocx*nprocy)
      endif
!
    endsubroutine xyaverages_z
!***********************************************************************
    subroutine xzaverages_y()
!
!  Calculate xz-averages (still depending on y).
!
!  12-oct-05/anders: adapted from xyaverages_z
!
      real, dimension (ny,nprocy,nnamey) :: fsumy
!
!  Communicate over all processors.
!  The result is only present on the root processor.
!
      if (nnamey>0) then
        call mpireduce_sum(fnamey,fsumy,(/ny,nprocy,nnamey/))
        if (lroot) &
            fnamey(:,:,1:nnamey)=fsumy(:,:,1:nnamey)/(nx*nz*nprocx*nprocz)
      endif
!
    endsubroutine xzaverages_y
!***********************************************************************
    subroutine yzaverages_x()
!
!  Calculate yz-averages (still depending on x).
!
!   2-oct-05/anders: adapted from xyaverages_z
!
      real, dimension (nx,nprocx,nnamex) :: fsumx
!
!  Communicate over all processors.
!  The result is only present on the root processor.
!
      if (nnamex>0) then
        call mpireduce_sum(fnamex,fsumx,(/nx,nprocx,nnamex/))
        if (lroot) &
            fnamex(:,:,1:nnamex)=fsumx(:,:,1:nnamex)/(ny*nz*nprocy*nprocz)
      endif
!
    endsubroutine yzaverages_x
!***********************************************************************
    subroutine yaverages_xz()
!
!  Calculate y-averages (still depending on x and z).
!
!   7-jun-05/axel: adapted from zaverages_xy
!
      real, dimension (nx,nz,nnamexz) :: fsumxz
!
!  Communicate over all processors along y beams.
!  The result is only present on the y-root processors.
!
      if (nnamexz>0) then
        call mpireduce_sum(fnamexz,fsumxz,(/nx,nz,nnamexz/),idir=2)
        if (ipy==0) fnamexz(:,:,1:nnamexz)=fsumxz(:,:,1:nnamexz)/nygrid
      endif
!
    endsubroutine yaverages_xz
!***********************************************************************
    subroutine zaverages_xy()
!
!  Calculate z-averages (still depending on x and y).
!
!  19-jun-02/axel: coded
!
      real, dimension (nx,ny,nnamexy) :: fsumxy
!
!  Communicate over all processors along z beams.
!  The result is only present on the z-root processors.
!
      if (nnamexy>0) then
        call mpireduce_sum(fnamexy,fsumxy,(/nx,ny,nnamexy/),idir=3)
        if (ipz==0) fnamexy(:,:,1:nnamexy)=fsumxy(:,:,1:nnamexy)/nzgrid
      endif
!
    endsubroutine zaverages_xy
!***********************************************************************
    subroutine write_1daverages()
!
!  Write 1d averages (z-averages, .., i.e. quantities that are only  1d
!  after averaging). These are written every it1 timesteps (like the
!  diagnostic line) and appended to their individual files.
!
!   7-aug-03/wolf: coded
!
      call write_xyaverages()
      call write_xzaverages()
      call write_yzaverages()
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
!
!   7-aug-03/wolf: adapted from wsnap
!
      if (lwrite_yaverages)   call write_yaverages()
      if (lwrite_zaverages)   call write_zaverages()
!
      if (ip<=10) write(*,*) 'write_2daverages: wrote phi(etc.)avgs'//ch2davg
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
    subroutine write_yaverages()
!
!  Write y-averages (which are 2d data) that have been requested via yaver.in.
!
!   7-jun-05/axel: adapted from write_zaverages
!
      if (ipy==0.and.nnamexz>0) then
        open(1, file=trim(directory_snap)//'/yaverages.dat', &
            form='unformatted', position='append')
        write(1) t2davgfirst
        write(1) fnamexz(:,:,1:nnamexz)
        close(1)
      endif
!
    endsubroutine write_yaverages
!***********************************************************************
    subroutine write_zaverages()
!
!  Write z-averages (which are 2d data) that have been requested via zaver.in.
!
!  19-jun-02/axel: adapted from write_xyaverages
!
      if (ipz==0.and.nnamexy>0) then
        open(1, file=trim(directory_snap)//'/zaverages.dat', &
            form='unformatted', position='append')
        write(1) t2davgfirst
        write(1) fnamexy(:,:,1:nnamexy)
        close(1)
      endif
!
    endsubroutine write_zaverages
!***********************************************************************
    integer function fparse_name(iname,cname,cform,ctest,itest)
!
!  Parse name and format of scalar print variable
!  On output, ITEST is set to INAME if CNAME matches CTEST
!  and CFORM is set to the format given as default.
!  E.g. if CTEST='bmax' *i.e. we are testing input line CNAME for 'bmax',
!  CNAME='bmax' will be parsed to ITEST=INAME, CFORM='(1pe10.2)',
!  CNAME='bmax(G5.1)' to ITEST=INAME, CFORM='G5.1',
!  CNAME='brms' to ITEST=<unchanged, normally 0>, CFORM='(1pe10.2)'
!
!   4-may-02/axel: coded
!   6-apr-04/wolf: more liberate format reading
!
      use General, only: safe_character_assign
!
      character (len=*) :: cname,cform
      character (len=*) :: ctest
      integer :: iname,itest,iform0,iform1,iform2,length,index_i
!
      intent(in)    :: iname,cname,ctest
      intent(inout) :: itest,cform
!      intent(out)   :: cform
!
!  Check whether format is given.
!
      iform0=index(cname,' ')
      iform1=index(cname,'(')
      iform2=index(cname,')')
!
!  Set format; use default if not given.
!
      if (iform1>0) then
        cform=cname(iform1+1:iform2-1)
        length=iform1-1
      else
        cform='1pE10.2'  !!(the nag-f95 compiler requires a comma after
                         !! 1p [does it?])
        length=iform0-1
      endif
!
!  Fix annoying Fortran 0p/1p stuff (Ew.d --> 1pEw.d, Fw.d --> 0pFw.d).
!
      if ((cform(1:1) == 'e') .or. (cform(1:1) == 'E') &
          .or. (cform(1:1) == 'g') .or. (cform(1:1) == 'G')) then
        call safe_character_assign(cform, '1p'//trim(cform))
      endif
      if ((cform(1:1) == 'f') .or. (cform(1:1) == 'F')) then
        call safe_character_assign(cform, '0p'//trim(cform))
      endif
!
!  If the name matches, we keep the name and can strip off the format.
!  The remaining name can then be used for the legend.
!
      if (cname(1:length)==ctest .and. itest==0) then
        itest=iname
        fparse_name=iname
      else 
        fparse_name=0
      endif
!
!  Integer formats are turned into floating point numbers.
!
      index_i=index(cform,'i')
      if (index_i/=0) then
        cform(index_i:index_i)='f'
        cform=trim(cform)//'.0'
      endif
!
    endfunction fparse_name
!***********************************************************************
    subroutine parse_name(iname,cname,cform,ctest,itest)
!
!   wrapper around fparse_name function: ignores return value 

!   26-nov-09/MR: coded
!
      character (len=*) :: cname,cform
      character (len=*) :: ctest
      integer :: iname,itest
!
      intent(in)    :: iname,cname,ctest
      intent(inout) :: itest,cform

      integer iret;

      iret = fparse_name(iname,cname,cform,ctest,itest)

    endsubroutine parse_name
!***********************************************************************
    subroutine save_name(a,iname)
!
!  Lists the value of a (must be treated as real) in fname array
!
!  26-may-02/axel: adapted from max_mn_name
!
      real :: a
      integer :: iname
!
!  Set corresponding entry in itype_name
!  This routine is to be called only once per step
!
      fname(iname)=a
      itype_name(iname)=ilabel_save
!
   endsubroutine save_name
!***********************************************************************
    subroutine save_name_halfz(a,iname)
!
!  Lists the value of a (must be treated as real) in fname array
!
!  16-may-09/axel: adapted from save_name
!
      real, dimension(2) :: a
      integer :: iname
!
!  Set corresponding entry in itype_name
!  This routine is to be called only once per step
!
      fname_half(iname,1)=a(1)
      fname_half(iname,2)=a(2)
!
   endsubroutine save_name_halfz
!***********************************************************************
    subroutine max_name(a,iname,lneg)
!
!  Successively calculate maximum of a, which is supplied at each call.
!
!  29-aug-05/anders: adapted from save_name
!
      integer :: a, iname
      logical, optional :: lneg
!
      if (present(lneg)) then
        if (lneg) then
          if (a<fname(iname)) fname(iname)=a
        else
          if (a>fname(iname)) fname(iname)=a
        endif
      else
        if (a>fname(iname)) fname(iname)=a
      endif
!
!  Set corresponding entry in itype_name.
!
      if (present(lneg)) then
        itype_name(iname)=ilabel_max_neg
      else
        itype_name(iname)=ilabel_max
      endif
!
    endsubroutine max_name
!***********************************************************************
    subroutine sum_name(a,iname)
!
!  Calculate the summation of a, which is supplied at each call.
!
!  17-jun-09/ccyang: adapted from max_name
!  03-sep-09/MR: corrected to real sum
!
      real, intent(in) :: a
      integer, intent(in) :: iname
!
     if (lfirstpoint) then

       fname(iname)=a
!
!  Set corresponding entry in itype_name.
!
       itype_name(iname)=ilabel_sum

     else
       fname(iname)=fname(iname)+a
     endif
!
    endsubroutine sum_name
!***********************************************************************
    subroutine max_mn_name(a,iname,lsqrt,l_dt,lneg,lreciprocal)
!
!  Successively calculate maximum of a, which is supplied at each call.
!  Start from zero if lfirstpoint=.true.
!
!   1-apr-01/axel+wolf: coded
!   4-may-02/axel: adapted for fname array
!  23-jun-02/axel: allows for taking square root in the end
!
      real, dimension (nx) :: a
      integer :: iname
      logical, optional :: lsqrt,l_dt,lneg,lreciprocal
!
      if (lfirstpoint) then
        fname(iname)=maxval(a)
      else
        fname(iname)=max(fname(iname),maxval(a))
      endif
!
!  Set corresponding entry in itype_name.
!
      if (present(lsqrt)) then
        itype_name(iname)=ilabel_max_sqrt
      elseif (present(l_dt)) then
        itype_name(iname)=ilabel_max_dt
      elseif (present(lneg)) then
        itype_name(iname)=ilabel_max_neg
      elseif (present(lreciprocal)) then
        itype_name(iname)=ilabel_max_reciprocal
      else
        itype_name(iname)=ilabel_max
      endif
!
    endsubroutine max_mn_name
!***********************************************************************
    subroutine sum_mn_name(a,iname,lsqrt,lint,ipart)
!
!  successively calculate sum of a, which is supplied at each call.
!  Start from zero if lfirstpoint=.true.
!  TODO: for nonperiodic arrays we want to multiply boundary data by 1/2.
!
!   1-apr-01/axel+wolf: coded
!   4-may-02/axel: adapted for fname array
!  23-jun-02/axel: allows for taking square root in the end
!  20-jun-07/dhruba: adapted for spherical polar coordinate system
!  30-aug-07/wlad: adapted for cylindrical coordinates
!  22-mar-08/axel: added ladd option, to add to previous values
!
!  Note [24-may-2004, wd]:
!    This routine should incorporate a test for iname /= 0, so instead of
!         if (idiag_b2m/=0)    call sum_mn_name(b2,idiag_b2m)
!    we can just use
!         call sum_mn_name(b2,idiag_b2m)
!  Same holds for similar routines.
!  Update [28-Sep-2004 wd]:
!    Done here, but not yet in all other routines
!
      real, dimension (nx) :: a,a_scaled
      real :: ppart,qpart
      integer :: iname
      integer, optional :: ipart
      logical, optional :: lsqrt, lint
!
!  Only do something if iname is not zero.
!
      if (iname/=0) then
!
!  Set corresponding entry in itype_name.
!
        if (present(lsqrt)) then
          itype_name(iname)=ilabel_sum_sqrt
        elseif (present(lint)) then
          itype_name(iname)=ilabel_integrate
        else
          itype_name(iname)=ilabel_sum
        endif
!
!  Set fraction if old and new stuff.
!
        if (present(ipart)) then
          ppart=1./float(ipart)
          if (lfirstpoint) then
            qpart=1.-ppart
          else
            qpart=1.
          endif
!
!  Add up contributions, taking coordinate system into acount (particles).
!
          fname(iname)=qpart*fname(iname)+ppart*sum(a)
!
!  Normal method.
!
        else
!
!  Scale "a" with volume differential if integration option is set.
!
          if (present(lint)) then
            a_scaled=a*xprim(l1:l2)*yprim(m)*zprim(n)
          else
            a_scaled=a
          endif
!
!  Add up contributions, taking coordinate system into acount (fluid).
!
          if (lfirstpoint) then
            fname(iname)=sum(a_scaled)
          else
            fname(iname)=fname(iname)+sum(a_scaled)
          endif
        endif
!
      endif
!
    endsubroutine sum_mn_name
!***********************************************************************
    subroutine sum_mn_name_halfy(a,iname)
!
!  To calculate averages over half the size of box, useful for simulations
!  which includes equator (either cartesian or spherical).  
!
!  dhruba : aped from sum_mn_name
!
      real, dimension (nx) :: a
      real :: sum_name
      integer :: iname

!
      if (iname /= 0) then
        sum_name=0
!
        if (y(m).ge.yequator)then
          sum_name=fname_half(iname,2)
        else
          sum_name=fname_half(iname,1)
        endif
!
        if (lfirstpoint) then
          fname_half(iname,1)=0
          fname_half(iname,2)=0
          sum_name=sum(a)
        else
          sum_name=sum_name+sum(a)
        endif
!
        if (y(m).ge.yequator)then
          fname_half(iname,2)=sum_name
        else
          fname_half(iname,1)=sum_name
        endif
!
      endif
!
    endsubroutine sum_mn_name_halfy
!***********************************************************************
    subroutine sum_mn_name_halfz(a,iname)
!
! To calculate averages over half the size of box (this time divided along the z
! direction), useful for simulations which includes equator (either cartesian 
! or spherical).  
!
!  7-may-09/dhruba: aped from sum_mn_name_halfy
!
      real, dimension (nx) :: a
      real :: sum_name
      integer :: iname
!
      if (iname /= 0) then
        sum_name=0
!
!  Note: north means 1 and south means 2.
!  However, north corresponds to z < zequator,
!  in analogy to theta < !pi/2 for north.
!
        if (z(n).ge.zequator) then
          sum_name=fname_half(iname,2)
        else
          sum_name=fname_half(iname,1)
        endif
!
        if (lfirstpoint) then
          fname_half(iname,1)=0
          fname_half(iname,2)=0
          sum_name=sum(a)
        else
          sum_name=sum_name+sum(a)
        endif
!
!  North means 1 and south means 2.
!  However, north corresponds to z < zequator,
!  in analogy to theta < !pi/2 for north.
!
        if (z(n).ge.zequator) then
          fname_half(iname,2)=sum_name
        else
          fname_half(iname,1)=sum_name
        endif
!
      endif
!
    endsubroutine sum_mn_name_halfz
!***********************************************************************
    subroutine sum_weighted_name(a,weight,iname,lsqrt)
!
!  Succesively calculate the weighted sum of a. The result is divided by the
!  total weight in the diagnostics subroutine.
!
!  17-apr-06/anders : coded
!
      real, dimension (:) :: a, weight
      integer :: iname
      logical, optional :: lsqrt
!
      integer, save :: it_save=-1, itsub_save=-1
!
      if (iname/=0) then
!
        if (it/=it_save .or. itsub/=itsub_save) then
          fname(iname)=0.0
          fweight(iname)=0.0
          it_save=it
          itsub_save=itsub
        endif
!
        fname(iname)  =fname(iname)  +sum(weight*a)
        fweight(iname)=fweight(iname)+sum(weight)
!
!  Set corresponding entry in itype_name.
!
        if (present(lsqrt)) then
          itype_name(iname)=ilabel_sum_weighted_sqrt
        else
          itype_name(iname)=ilabel_sum_weighted
        endif
!
      endif
!
    endsubroutine sum_weighted_name
!***********************************************************************
    subroutine surf_mn_name(a,iname)
!
!  Successively calculate surface integral. This routine assumes
!  that "a" contains the partial result for each pencil, so here
!  we just need to add up the contributions from all processors.
!  Start from zero if lfirstpoint=.true.
!
!  14-aug-03/axel: adapted from sum_mn_name
!
      real, intent(in) :: a
      integer, intent(in) :: iname
!
      if (lfirstpoint) then
        fname(iname)=a
      else
        fname(iname)=fname(iname)+a
      endif
!
!  Set corresponding entry in itype_name.
!
      itype_name(iname)=ilabel_surf
!
    endsubroutine surf_mn_name
!***********************************************************************
    subroutine integrate_mn(a,inta)
!
!  Successively calculate sum of a, which is supplied at each call.
!  Start from zero if lfirstpoint=.true. ultimately multiply by dv
!  to get the integral.  This differs from sum_mn_name by the
!  setting of ilabel_integrate and hence in the behaviour in the final
!  step.
!
!  Note, for regular integration (uniform meshes) it is better
!  to use the usual sum_mn_name routine with lint=.true.
!
!   1-dec-09/dhruba+piyali: 
!
      real, dimension (nx) :: a,fac
      real :: inta
!
!  initialize by the volume element (which is different for different m and n)
!
      fac=dVol1(l1:l2)*dVol2(m)*dVol3(n)
!
!  initialize if one is on the first point, or add up otherwise
!
      if (lfirstpoint) then
        inta=sum(a*fac)
      else
        inta=inta+sum(a*fac)
      endif
!
    endsubroutine integrate_mn

!***********************************************************************
    subroutine integrate_mn_name(a,iname)
!
!  Successively calculate sum of a, which is supplied at each call.
!  Start from zero if lfirstpoint=.true. ultimately multiply by dv
!  to get the integral.  This differs from sum_mn_name by the
!  setting of ilabel_integrate and hence in the behaviour in the final
!  step.
!
!  Note, for regular integration (uniform meshes) it is better
!  to use the usual sum_mn_name routine with lint=.true.
!
!   30-may-03/tony: adapted from sum_mn_name
!   13-nov-06/tony: modified to handle stretched mesh
!
      real, dimension (nx) :: a,fac
      integer :: iname
!
!  initialize by the volume element (which is different for different m and n)
!
      fac=dVol1(l1:l2)*dVol2(m)*dVol3(n)
!
!  initialize if one is on the first point, or add up otherwise
!
      if (lfirstpoint) then
        fname(iname)=sum(a*fac)
      else
        fname(iname)=fname(iname)+sum(a*fac)
      endif
!
!  set corresponding entry in itype_name
!
      itype_name(iname)=ilabel_integrate
!
    endsubroutine integrate_mn_name
!***********************************************************************
    subroutine xysum_mn_name_z(a,iname)
!
!  Successively calculate sum over x,y of a, which is supplied at each call.
!  The result fnamez is z-dependent.
!  Start from zero if lfirstpoint=.true.
!
!   5-jun-02/axel: adapted from sum_mn_name
!
      real, dimension (nx) :: a
      integer :: iname,n_nghost,isum,lmax
!
!  Initialize to zero, including other parts of the z-array
!  which are later merged with an mpi reduce command.
!
      if (lfirstpoint) fnamez(:,:,iname)=0.0
!
!  n starts with nghost+1=4, so the correct index is n-nghost
!
      lmax=l2
      n_nghost=n-nghost
      if (lav_smallx)  lmax=ixav_max
      if (.not.loutside_avg) then
        do isum=l1,lmax
          fnamez(n_nghost,ipz+1,iname)=fnamez(n_nghost,ipz+1,iname)+ & 
               a(isum-nghost)
        enddo
      endif
!
    endsubroutine xysum_mn_name_z
!***********************************************************************
    subroutine xzsum_mn_name_y(a,iname)
!
!  Successively calculate sum over x,z of a, which is supplied at each call.
!  The result fnamey is y-dependent.
!  Start from zero if lfirstpoint=.true.
!
!  12-oct-05/anders: adapted from xysum_mn_name_z
!
      real, dimension (nx) :: a
      integer :: iname,m_nghost,isum,lmax
!
!  Initialize to zero, including other parts of the z-array
!  which are later merged with an mpi reduce command.
!
      if (lfirstpoint) fnamey(:,:,iname)=0.0
!
!  m starts with mghost+1=4, so the correct index is m-nghost.
!
      m_nghost=m-nghost
      lmax=l2
      if (lav_smallx) lmax=ixav_max
      if (.not.loutside_avg) then
        do isum=l1,lmax
          fnamey(m_nghost,ipy+1,iname)=fnamey(m_nghost,ipy+1,iname)+ &
               a(isum-nghost)
        enddo
      endif
!
    endsubroutine xzsum_mn_name_y
!***********************************************************************
    subroutine yzsum_mn_name_x(a,iname)
!
!  Successively calculate sum over y,z of a, which is supplied at each call.
!  The result fnamex is x-dependent.
!  Start from zero if lfirstpoint=.true.
!
!   2-oct-05/anders: adapted from xysum_mn_name_z
!
      real, dimension (nx) :: a
      integer :: iname
!
!  Initialize to zero.
!
      if (lfirstpoint) fnamex(:,:,iname)=0.0
!
      fnamex(:,ipx+1,iname)=fnamex(:,ipx+1,iname)+a
!
    endsubroutine yzsum_mn_name_x
!***********************************************************************
    subroutine xyintegrate_mn_name_z(a,iname)
!
!   Integrate over x and y. Apply trapezoidal rule properly in the case
!   of non-periodic boundaries.
!
!   18-jun-07/tobi: adapted from xysum_mn_name_z
!
      real, dimension (nx) :: a
      integer :: iname
      real :: fac,suma
!
!  Initialize to zero, including other parts of the z-array
!  which are later merged with an mpi reduce command.
!
      if (lfirstpoint) fnamez(:,:,iname) = 0.0
!
      fac=1.0
!
      if ((m==m1.and.ipy==0).or.(m==m2.and.ipy==nprocy-1)) then
        if (.not.lperi(2)) fac = .5*fac
      endif
!
      if (lperi(1)) then
        suma = fac*sum(a)
      else
        suma = fac*(sum(a(2:nx-1))+.5*(a(1)+a(nx)))
      endif
!
!  n starts with nghost+1=4, so the correct index is n-nghost.
!
      fnamez(n-nghost,ipz+1,iname) = fnamez(n-nghost,ipz+1,iname) + suma
!
    endsubroutine xyintegrate_mn_name_z
!***********************************************************************
    subroutine xzintegrate_mn_name_y(a,iname)
!
!   Integrate over x and z. Apply trapezoidal rule properly in the case
!   of non-periodic boundaries.
!
!   18-jun-07/tobi: adapted from xzsum_mn_name_y
!
      real, dimension (nx) :: a
      integer :: iname
      real :: fac,suma
!
!  Initialize to zero, including other parts of the z-array
!  which are later merged with an mpi reduce command.
!
      if (lfirstpoint) fnamey(:,:,iname) = 0.

      fac = 1.

      if ((n==n1.and.ipz==0).or.(n==n2.and.ipz==nprocz-1)) then
        if (.not.lperi(3)) fac = .5*fac
      endif

      if (lperi(1)) then
        suma = fac*sum(a)
      else
        suma = fac*(sum(a(2:nx-1))+.5*(a(1)+a(nx)))
      endif
!
!  m starts with mghost+1=4, so the correct index is m-nghost.
!
      fnamey(m-nghost,ipy+1,iname) = fnamey(m-nghost,ipy+1,iname) + suma

    endsubroutine xzintegrate_mn_name_y
!***********************************************************************
    subroutine yzintegrate_mn_name_x(a,iname)
!
!   Integrate over y and z. Apply trapezoidal rule properly in the case
!   of non-periodic boundaries.
!
!   18-jun-07/tobi: adapted from yzsum_mn_name_x
!
      real, dimension (nx) :: a
      integer :: iname
      real :: fac
!
!  Initialize to zero.
!
      if (lfirstpoint) fnamex(:,:,iname) = 0.0
!
      fac=1.0
!
      if ((m==m1.and.ipy==0).or.(m==m2.and.ipy==nprocy-1)) then
        if (.not.lperi(2)) fac = .5*fac
      endif
!
      if ((n==n1.and.ipz==0).or.(n==n2.and.ipz==nprocz-1)) then
        if (.not.lperi(3)) fac = .5*fac
      endif
!
      fnamex(:,ipx+1,iname) = fnamex(:,ipx+1,iname) + fac*a
!
    endsubroutine yzintegrate_mn_name_x
!***********************************************************************
    subroutine ysum_mn_name_xz(a,iname)
!
!  Successively calculate sum over y of a, which is supplied at each call.
!  The result fnamexz is xz-dependent.
!  Start from zero if lfirstpoint=.true.
!
!   7-jun-05/axel: adapted from zsum_mn_name_xy
!
      real, dimension (nx) :: a
      integer :: iname,n_nghost
!
!  Initialize to zero, including other parts of the xz-array
!  which are later merged with an mpi reduce command.
!
      if (lfirstpoint) fnamexz(:,:,iname)=0.
!
!  n starts with nghost+1=4, so the correct index is n-nghost.
!
      n_nghost=n-nghost
!
      fnamexz(:,n_nghost,iname)=fnamexz(:,n_nghost,iname)+a
!
    endsubroutine ysum_mn_name_xz
!***********************************************************************
    subroutine zsum_mn_name_xy(a,iname)
!
!  Successively calculate sum over z of a, which is supplied at each call.
!  The result fnamexy is xy-dependent.
!  Start from zero if lfirstpoint=.true.
!
!  19-jun-02/axel: adapted from xysum_mn_name
!
      real, dimension (nx) :: a
      integer :: iname,m_nghost
!
!  Initialize to zero, including other parts of the xy-array
!  which are later merged with an mpi reduce command.
!
      if (lfirstpoint) fnamexy(:,:,iname)=0.
!
!  m starts with nghost+1=4, so the correct index is m-nghost.
!
      m_nghost=m-nghost
      fnamexy(:,m_nghost,iname)=fnamexy(:,m_nghost,iname)+a
!
    endsubroutine zsum_mn_name_xy
!***********************************************************************
    real function get_from_fname(iname)
!
!   Gets value from fname.
!
!   30-oct-09/MR: coded
!
    integer, intent(in) :: iname
!
    if (iname<1.or.iname>nname) then
      call fatal_error('get_from_fname', 'index not in legal range')
      get_from_fname = 0
    endif
!
    get_from_fname = fname(iname)
!
    endfunction get_from_fname
!***********************************************************************
    subroutine allocate_xyaverages
!
!  Allocate the variables needed for xy-averages.
!
!   24-nov-09/anders: copied from allocate_yaverages
!
      integer :: stat
!
!  Allocate and initialize to zero. Setting it to zero is only
!  necessary because of the pencil test, which doesn't compute these
!  averages, and only evaluates its output for special purposes
!  such as computing mean field energies in calc_bmz, for example,
!
      allocate(fnamez(nz,nprocz,nnamez),stat=stat)
      fnamez=0.
!
      if (stat>0) then
        call fatal_error('allocate_xyaverages', &
            'Could not allocate memory for fnamez')
      else
        if (lroot) print*, 'allocate_xyaverages: allocated memory for '// &
            'fnamez  with nnamez  =', nnamez
      endif
!
      allocate(cnamez(nnamez),cformz(nnamez)) 
!
    endsubroutine allocate_xyaverages
!***********************************************************************
    subroutine allocate_xzaverages
!
!  Allocate the variables needed for xz-averages.
!
!   24-nov-09/anders: copied from allocate_yaverages
!
      integer :: stat
!
      allocate(fnamey(ny,nprocy,nnamey),stat=stat)
!
      if (stat>0) then
        call fatal_error('allocate_xzaverages', &
            'Could not allocate memory for fnamey')
      else
        if (lroot) print*, 'allocate_xzaverages: allocated memory for '// &
            'fnamey  with nnamey  =', nnamey
      endif
!
      allocate(cnamey(nnamey),cformy(nnamey)) 
!
    endsubroutine allocate_xzaverages
!***********************************************************************
    subroutine allocate_yzaverages
!
!  Allocate the variables needed for yz-averages.
!
!   24-nov-09/anders: copied from allocate_yaverages
!
      integer :: stat
!
      allocate(fnamex(nx,nprocx,nnamex),stat=stat)
!
      if (stat>0) then
        call fatal_error('allocate_yzaverages', &
            'Could not allocate memory for fnamex')
      else
        if (lroot) print*, 'allocate_yzaverages: allocated memory for '// &
            'fnamex  with nnamex  =', nnamex
      endif
!
      allocate(cnamex(nnamex),cformx(nnamex)) 
!
    endsubroutine allocate_yzaverages
!***********************************************************************
    subroutine allocate_yaverages
!
!  Allocate the variables needed for y-averages.
!
!   12-aug-09/dhruba: coded
!
      integer :: stat
!
      allocate(fnamexz(nx,nz,nnamexz),stat=stat)
!
      if (stat>0) then
        call fatal_error('allocate_yaverages', &
            'Could not allocate memory for fnamexz')
      else
        if (lroot) print*, 'allocate_yaverages : allocated memory for '// &
            'fnamexz with nnamexz =', nnamexz
      endif
!
      allocate(cnamexz(nnamexz),cformxz(nnamexz)) 
!
    endsubroutine allocate_yaverages
!*******************************************************************
    subroutine allocate_zaverages
!
!  Allocate the variables needed for z-averages.
!
!   12-aug-09/dhruba: coded
!
      integer :: stat
!
      allocate(fnamexy(nx,ny,nnamexy),stat=stat)
!
      if (stat>0) then
        call fatal_error('allocate_zaverages', &
            'Could not allocate memory for fnamexy')
      else
        if (lroot) print*, 'allocate_zaverages : allocated memory for '// &
            'fnamexy with nnamexy =', nnamexy
      endif
!
      allocate(cnamexy(nnamexy),cformxy(nnamexy))
!
    endsubroutine allocate_zaverages
!*******************************************************************
    subroutine xyaverages_clean_up
!
!  Allocate the variables needed for xy-averages.
!
!   24-nov-09/anders: copied from yaverages_clean_up
!
      if (allocated(fnamez)) deallocate(fnamez)
      if (allocated(cnamez)) deallocate(cnamez)
      if (allocated(cformz)) deallocate(cformz)
!
    endsubroutine xyaverages_clean_up
!***********************************************************************
    subroutine xzaverages_clean_up
!
!  Allocate the variables needed for xz-averages.
!
!   24-nov-09/anders: copied from yaverages_clean_up
!
      if (allocated(fnamey)) deallocate(fnamey)
      if (allocated(cnamey)) deallocate(cnamey)
      if (allocated(cformy)) deallocate(cformy)
!
    endsubroutine xzaverages_clean_up
!***********************************************************************
    subroutine yzaverages_clean_up
!
!  Allocate the variables needed for yz-averages.
!
!   24-nov-09/anders: copied from yaverages_clean_up
!
      if (allocated(fnamex)) deallocate(fnamex)
      if (allocated(cnamex)) deallocate(cnamex)
      if (allocated(cformx)) deallocate(cformx)
!
    endsubroutine yzaverages_clean_up
!***********************************************************************
    subroutine yaverages_clean_up
!
!  Allocate the variables needed for y-averages.
!
!   12-aug-09/dhruba: coded
!
      if (allocated(fnamexz)) deallocate(fnamexz)
      if (allocated(cnamexz)) deallocate(cnamexz)
      if (allocated(cformxz)) deallocate(cformxz)
!
    endsubroutine yaverages_clean_up
!*******************************************************************
    subroutine zaverages_clean_up
!
!  Allocate the variables needed for z-averages.
!
!   12-aug-09/dhruba: coded
!
      if (allocated(fnamexy)) deallocate(fnamexy)
      if (allocated(cnamexy)) deallocate(cnamexy)
      if (allocated(cformxy)) deallocate(cformxy)
!
    endsubroutine zaverages_clean_up
!*******************************************************************
    subroutine init_xaver
!
!  Initialize variables for x-averaging.
!
!  26-oct-09/dhruba: coded
!
     integer :: lx
!
     if (xav_max.lt.x(l1)) loutside_avg=.true.  
!
     ixav_max=l1   
!
     do lx=l1,l2
       if (x(lx).lt.xav_max) ixav_max=lx
     enddo
!
   endsubroutine init_xaver
!*******************************************************************
endmodule Diagnostics
