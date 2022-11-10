! $Id$
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
  use General, only: safe_sum
!
  implicit none
!
  public :: initialize_diagnostics, prints
  public :: diagnostic, initialize_time_integrals,get_average_density
  public :: xyaverages_z, xzaverages_y, yzaverages_x
  public :: phizaverages_r, yaverages_xz, zaverages_xy
  public :: phiaverages_rz
  public :: write_1daverages, write_2daverages
  public :: write_sound
  public :: write_1daverages_prepare, write_2daverages_prepare
  public :: name_is_present
  public :: expand_cname, parse_name, fparse_name
! GPU-START
  public :: set_type
  public :: save_name
  public :: save_name_halfz, save_name_sound
  public :: max_name, sum_name
  public :: max_mn_name, sum_mn_name, integrate_mn_name, sum_weighted_name
  public :: integrate_mn
  public :: sum_mn_name_halfy, surf_mn_name, sum_lim_mn_name
  public :: sum_mn_name_halfz
  public :: xysum_mn_name_z, xzsum_mn_name_y, yzsum_mn_name_x
  public :: phizsum_mn_name_r, ysum_mn_name_xz, zsum_mn_name_xy
  public :: phisum_mn_name_rz, calc_phiavg_profile
  public :: yzintegrate_mn_name_x, xzintegrate_mn_name_y, xyintegrate_mn_name_z
  public :: ysum_mn_name_xz_npar, xysum_mn_name_z_npar, yzsum_mn_name_x_mpar
  public :: zsum_mn_name_xy_mpar_scal, zsum_mn_name_xy_mpar, &
            zsum_mn_name_xy_arr, zsum_mn_name_xy_arr2
! GPU-END
  public :: allocate_fnames,allocate_vnames,allocate_sound
  public :: allocate_xyaverages, allocate_xzaverages, allocate_yzaverages
  public :: allocate_phizaverages
  public :: allocate_yaverages, allocate_zaverages, allocate_phiaverages, allocate_zaverages_data
  public :: trim_averages
  public :: fnames_clean_up, vnames_clean_up, diagnostics_clean_up
  public :: get_from_fname
  public :: init_xaver
  public :: gen_form_legend
  public :: sign_masked_xyaver
  public :: report_undefined_diagnostics
!
  interface max_name
    module procedure max_name_int
    module procedure max_name_real
  endinterface max_name
!
  interface sum_name
    module procedure sum_name_int
    module procedure sum_name_real
  endinterface sum_name
!
  interface expand_cname
    module procedure expand_cname_short
    module procedure expand_cname_full
  endinterface expand_cname
!
  interface parse_name
    module procedure parse_name_s
    module procedure parse_name_sa
    module procedure parse_name_v
  endinterface parse_name
!
  interface sum_mn_name
    module procedure sum_mn_name_real
    module procedure sum_mn_name_arr2
    module procedure sum_mn_name_std
  endinterface sum_mn_name

  interface zsum_mn_name_xy
    module procedure zsum_mn_name_xy_scal
    module procedure zsum_mn_name_xy_vec
    module procedure zsum_mn_name_xy_arr
    module procedure zsum_mn_name_xy_arr2
  endinterface zsum_mn_name_xy

  interface zsum_mn_name_xy_mpar
    module procedure zsum_mn_name_xy_mpar_scal
    module procedure zsum_mn_name_xy_mpar_vec
  endinterface zsum_mn_name_xy_mpar
!
  private
!
  real, dimension (nrcyl,nx) :: phiavg_profile=0.0
  real :: dVol_rel1

  integer :: mnamer
  character (len=intlen) :: ch1davg, ch2davg
!
! Variables for Yin-Yang grid: z-averages.
!
  real, dimension(:,:,:), allocatable :: fnamexy_cap

  contains
!***********************************************************************
    subroutine initialize_diagnostics
!
!  Setup variables needed for output of diagnostic quantities and
!  averages.
!
!  14-aug-03/axel: added dxy, dyz, and dxz
!  26-aug-13/MR: removed switch first; moved calculation of dVol_rel1 from diagnostic
!  31-mar-15/MR: added switch for proper volume averaging
!
      integer :: i
      real :: dxeff,dyeff,dzeff
      real :: intdr_rel, intdtheta_rel, intdphi_rel, intdz_rel
! 
!  Initialize rcyl for the phi-averages grid. Does not need to be
!  done after each reload of run.in, but this is the easiest way
!  of doing it.
!
      if (nrcyl/=0) then
        drcyl=xyz1(1)/nrcyl
        rcyl=(/ ((i-0.5)*drcyl, i=1,nrcyl) /)
      else
        drcyl=0.0
      endif
!
!  Calculate the three surface elements. Take care of degenerate dimensions.
!
      if (nxgrid==1) then; dxeff=1.; else; dxeff=dx; endif
      if (nygrid==1) then; dyeff=1.; else; dyeff=dy; endif
      if (nzgrid==1) then; dzeff=1.; else; dzeff=dz; endif
!
      dsurfxy=dxeff*dyeff
      dsurfyz=dyeff*dzeff
      dsurfzx=dzeff*dxeff
!
!  Calculate relative volume integral.
!
      if (lproper_averages) then
        dVol_rel1=1./box_volume
      elseif (lspherical_coords) then
!
!  Prevent zeros from less than 3-dimensional runs
!  (maybe this should be 2pi, but maybe not).
!
        if (nxgrid/=1) then
          intdr_rel = (xyz1(1)**3-xyz0(1)**3)/(3.*dx)
        else
          intdr_rel = 1.
        endif
!
        if (nygrid/=1) then
          intdtheta_rel = -(cos(xyz1(2))-cos(xyz0(2)))/dy
        else
          intdtheta_rel = 1.
        endif
!
        if (nzgrid/=1) then
          intdphi_rel = (xyz1(3) - xyz0(3))/dz
        else
          intdphi_rel = 1.
        endif
!
        dVol_rel1=1./(intdr_rel*intdtheta_rel*intdphi_rel)
!
      elseif (lcylindrical_coords) then
!
!  Prevent zeros from less than 3-dimensional runs
!
        if (nxgrid/=1) then
          intdr_rel = (xyz1(1)**2 - xyz0(1)**2)/(2.*dx)
        else
          intdr_rel = 1.
        endif
!
        if (nygrid/=1) then
          intdphi_rel = (xyz1(2) - xyz0(2))/dy
        else
          intdphi_rel = 1.
        endif
!
        if (nzgrid/=1) then
          intdz_rel = (xyz1(3) - xyz0(3))/dz
        else
          intdz_rel = 1.
        endif
!
        dVol_rel1=1./(intdr_rel*intdphi_rel*intdz_rel)
!
      else
        dVol_rel1=1./nwgrid
      endif
!
      if (lroot.and.ip<=10) print*,'dVol_rel1=',dVol_rel1

    endsubroutine initialize_diagnostics
!***********************************************************************
    subroutine prints
!
!  Reads and registers print parameters gathered from the different
!  modules and marked in `print.in'.
!
!   3-may-02/axel: coded
!  20-aug-13/MR: changes for accumulating and complex diagnostics
!  26-aug-13/MR: corrected insertion of imaginary values
!  12-jan-17/MR: undefined diagnostics suppressed in output
!
      use General, only: safe_character_append, compress
      use Cparam, only: max_col_width
      use IO, only: IO_strategy
      use HDF5_IO, only: output_timeseries
      use Sub, only: insert
!
      character (len=1000) :: fform,legend,line
      integer :: iname, nnamel
      real, dimension(2*nname) :: buffer
      integer, parameter :: lun=1
      logical, save :: lfirst_call = .true.
!
!  Add general (not module-specific) quantities for diagnostic output. If the
!  timestep (=dt) is to be written, it is known only after time_step, so the best
!  place to enter it into the save list is here. Use 1.0*(it-1) to have floating
!  point or double precision.
!
      if (lroot) then
        call save_name(tdiagnos,idiag_t)
        call save_name(dt,idiag_dt)
        call save_name(one_real*(it-1),idiag_it)
!
!  Whenever itype_name=ilabel_max_dt, scale result by dt (for printing Courant
!  time).
!
        do iname=1,nname
          if (itype_name(iname)==ilabel_max_dt) fname(iname)=dt*fname(iname)
        enddo

        call gen_form_legend(fform,legend)
!
        if (ldebug) then
          write(0,*) 'PRINTS.prints: format = ', trim(fform)
          write(0,*) 'PRINTS.prints: args   = ', fname(1:nname)
        endif
!
!  This treats all numbers as floating point numbers.  Only those numbers are
!  given (and computed) that are also listed in print.in.
!
        if (lfirst_call) then
          write(*,*)
          write(*,'(" ",A)') trim(legend)
!
!  Write legend to extra file (might want to do only once after each lreset).
!
          open(lun,file=trim(datadir)//'/legend.dat')
          write(lun,'(" ",A)') trim(legend)
          close(lun)
        endif
!
        if (ldebug) write(*,*) 'before writing prints'
        buffer(1:nname) = fname(1:nname)
!
!  Add accumulated values to current ones if existent.
!
        where( fname_keep(1:nname) /= impossible .and. &
               itype_name(1:nname)<ilabel_complex ) &
           buffer(1:nname) = buffer(1:nname)+fname_keep(1:nname)
!
!  Write 'time_series.h5' if output format is HDF5
!
        if (IO_strategy == "HDF5".and.lwrite_ts_hdf5) call output_timeseries (buffer, fname_keep)
        nnamel=nname
        call compress(buffer,cform=='',nnamel)
!
!  Insert imaginary parts behind real ones if quantity is complex.
!
        do iname=nname,1,-1
          if (itype_name(iname)>=ilabel_complex .and. cform(iname)/='') &
            call insert(buffer,(/fname_keep(iname)/),iname+1,nnamel)
        enddo
!
!  Put output line into a string.
!
        write(line,trim(fform)) buffer(1:nnamel)
        call clean_line(line)
!
!  Append to diagnostics file.
!
        open(lun,file=trim(datadir)//'/time_series.dat',position='append')
        if (lfirst_call) then
          write(lun,"('"//comment_char//"',a)") trim(legend)
          lfirst_call = .false.
        endif
        write(lun,'(a)') trim(line)
        !flush(lun)               ! this is a F2003 feature...
        close(lun)
!
!  Write to stdout.
!
        write(*,'(a)') trim(line)
        !flush(6)                 ! this is a F2003 feature....
!
      endif                      ! (lroot)
!
      if (ldebug) write(*,*) 'exit prints'
!
!  reset non-accumulating values (marked with zero in fname_keep)
!
      where( fname_keep(1:nname)==0. .or. &
             itype_name(1:nname)>=ilabel_complex ) fname(1:nname)=0.
      where( itype_name(1:nname)>=ilabel_complex ) fname_keep(1:nname)=0.
!
    endsubroutine prints
!***********************************************************************
    subroutine report_undefined(cname,cform,len,file)
!
!  Reports diagnostics requested in cname, but undefined (their format is empty in cform).
!
!  12-jan-17/MR: coded
!
      use General, only: itoa
      use Syscalls, only: system_cmd

      character(len=*), dimension(:), intent(IN) :: cname,cform
      integer,                        intent(IN) :: len
      character(len=*),               intent(IN) :: file

      integer :: ind,i
      character(LEN=512) :: text
      character(LEN=512) :: sedstring

      sedstring=''
      text='WARNING:'; ind=-1
      do i=1,len
        if (cname(i)/=''.and.cform(i)=='') then
          ind=index(cname(i),'('); if (ind==0) ind=len_trim(cname(i))+1
          text=trim(text)//' '//trim(cname(i)(1:ind-1))//','
          sedstring=trim(sedstring)//" -e'"//trim(itoa(i))//" s/^\(.\)/#\1/'"
        endif
      enddo

      if (ind/=-1) then
!
!  If there are undefined diagnostics, comment them out in *.in.
!
        text=text(1:index(text,',',BACK=.true.)-1)//' diagnostic(s) in '//trim(file)//' undefined or multiply defined!'
        print*, trim(text)
        call system_cmd('sed '//trim(sedstring)//' -isv '//trim(file))
      endif

    endsubroutine report_undefined
!***********************************************************************
    subroutine report_undefined_diagnostics

      if (lroot) then
        call report_undefined(cname,cform,nname,'print.in')
        if (allocated(cname_sound)) &
          call report_undefined(cname_sound,cform_sound,nname_sound,'sound.in')
        if (allocated(cnamex)) &
        call report_undefined(cnamex,cformx,nnamex,'yzaver.in')
        if (allocated(cnamey)) &
        call report_undefined(cnamey,cformy,nnamey,'xzaver.in')
        if (allocated(cnamez)) &
        call report_undefined(cnamez,cformz,nnamez,'xyaver.in')
        if (allocated(cnamer)) &
        call report_undefined(cnamer,cformr,nnamer,'phizaver.in')
        if (allocated(cnamexy)) &
        call report_undefined(cnamexy,cformxy,nnamexy,'zaver.in')
        if (allocated(cnamexz)) &
        call report_undefined(cnamexz,cformxz,nnamexz,'yaver.in')
        if (allocated(cnamerz)) &
        call report_undefined(cnamerz,cformrz,nnamerz,'phiaver.in')
        if (allocated(cnamev)) &
        call report_undefined(cnamev,cformv,nnamev,'video.in')
      endif

    endsubroutine report_undefined_diagnostics
!***********************************************************************
    subroutine gen_form_legend(fform,legend)
!
!  19-aug-13/MR: outsourced from prints
!  10-jan-17/MR: added adjustment of fixed point format to output value
!  12-jan-17/MR: undefined diagnostics suppressed in format
!
      use General, only: safe_character_append, itoa
      use Sub, only: noform
!
      character (len=*)          ,intent(OUT) :: fform
      character (len=*), optional,intent(OUT) :: legend
!
      character, parameter :: comma=','
      character(len=40)    :: tform
      integer              :: iname,index_i,index_d,length
      logical              :: lcompl
      real                 :: rlength
!
!  Produce the format.
!
        fform = '(TL1'
        if (present(legend)) legend=''
!
        do iname=1,nname
!
          if (cform(iname)/='') then

            lcompl=itype_name(iname)>=ilabel_complex

            if (lcompl) then
              tform = comma//'" ("'//comma//trim(cform(iname))//comma &
                      //'","'//comma//trim(cform(iname))//comma//'")"'
            else
              if (fname(iname)/=0.) then
!
!  Check output item for fixed-point formats.
!
                index_i=scan(cform(iname),'fF')
                if (index_i>0) then
                  read(cform(iname)(index_i+1:),*) rlength
                  if (floor(alog10(abs(fname(iname))))+2 > rlength) then
                    index_d=index(cform(iname),'.')
                    length=len(trim(cform(iname)))
                    write(cform(iname)(index_i+1:), &
                          '(f'//trim(itoa(length-index_i))//'.'//trim(itoa(length-index_d))//')') rlength+2
                  endif
                endif
              endif
              if (cform(iname)/='') tform = comma//trim(cform(iname))

            endif
            call safe_character_append(fform, trim(tform))
            if (present(legend)) call safe_character_append(legend, noform(cname(iname),lcompl))

          endif
        enddo
        call safe_character_append(fform, ')')
!
    endsubroutine gen_form_legend
!***********************************************************************
    subroutine clean_line(line)
!
!  removes spurious dots from line, makes it to comment line if unreadable
!
!  14-jan-11/MR: outsourced from prints
!
    character (LEN=*), intent(inout) :: line
    integer :: ind
!
!  If line contains spurious dots, remove them.
!
      ind=index(line,'. ')
      if (ind >= 1) line(ind:ind)=' '
!
!  If the line contains unreadable characters, then comment out line.
!
      ind=index(line,'***') + index(line,'???')
      if (ind > 0) line(1:1)=comment_char
!
    endsubroutine clean_line
!***********************************************************************
    subroutine write_sound_append(legend,sublegend,coorlegend,line,lwrite_legend)
!
!  Append to diagnostics file.
!
!  27-Nov-2014/Bourdin.KIS: cleaned up code from write_sound
!
      character (len=*), intent(in) :: legend, sublegend, coorlegend, line
      logical, intent(in) :: lwrite_legend
!
      integer, parameter :: lun=1
      integer :: len_stroke, len_leg
!
        open(lun,file=trim(directory)//'/sound.dat',position='append')
!
        if (lwrite_legend) then
          if (dimensionality > 0) then
            len_leg=len_trim(legend)
            len_stroke=max(len_leg,len_trim(coorlegend),len_trim(sublegend))
            write(lun,'(a)') trim(legend)//repeat('-',max(0,len_stroke-len_leg))
            write(lun,'(a)') trim(coorlegend)
            if ( ncoords_sound>1 ) write(lun,'(a)') trim(sublegend)
            write(lun,'(a)') comment_char//repeat('-',len_stroke-1)
          else
            write(lun,'(a)') trim(legend)
          endif
        endif
!
        write(lun,'(a)') trim(line)
        close(lun)
!
    endsubroutine write_sound_append
!***********************************************************************
    subroutine write_sound(tout)
!
!  Reads and registers "sound" parameters gathered from the different
!  modules and marked in `sound.in'.
!
!   3-dec-10/dhruba+joern: coded
!  10-jan-11/MR: modified
!
      use General, only: itoa, safe_character_append, safe_character_prepend,compress
      use Sub, only: noform
!
      implicit none
!
      real, intent(in) :: tout
!
      logical :: ldata
      character (len=linelen*3) :: fform,legend,sublegend,coorlegend,line
      character (len=*), parameter :: tform='(f10.4'
      integer, parameter :: ltform=10
      character (len=ltform) :: scoor
      character (len=3*ncoords_sound*max_col_width) :: item
      real    :: coor
      integer :: iname,leng,nc,nch,nleg,nsub,idim,icoor,i,j,len
      logical, save :: lfirst_call = .true.
!
!  Produce the format.
!
      fform = tform//','
      if ( ncoords_sound>1 ) &
        call safe_character_append(fform, trim(itoa(ncoords_sound))//'(')
!
      if (lfirst_call) then
!
        legend = comment_char//noform('t'//tform//')')
!
        if (dimensionality>0) then
!
          coorlegend = comment_char//' Point(s):   '
          nleg = len_trim(coorlegend)+3
!
          do icoor=1,ncoords_sound
!
            if (ncoords_sound>1) coorlegend(nleg+1:) = trim(itoa(icoor))//' = '
            nleg = len_trim(coorlegend)+1
!
            item = '('
            do idim=1,dimensionality
              select case (dim_mask(idim))
              case (1); coor=x(sound_coords_list(icoor,1))
              case (2); coor=y(sound_coords_list(icoor,2))
              case (3); coor=z(sound_coords_list(icoor,3))
              end select
              write(scoor,tform//')') coor
              call safe_character_append(item,trim(adjustl(scoor))//',')
            enddo
            coorlegend(nleg+1:) = item(1:len_trim(item)-1)//'),   '
            nleg = len_trim(coorlegend)+4
!
          enddo
          coorlegend(nleg-4:nleg)=' '
        endif
!
      endif
!
      ldata = .false.
      sublegend = comment_char
!
      do iname=1,nname_sound
!
        if (cform_sound(iname)/=' ') then

          ldata=.true.
          if (lfirst_call) then
!
            item = noform(cname_sound(iname))
            if ( ncoords_sound>1 ) then
!
              nsub = len_trim(legend)
              leng = len_trim(item)
              nc  = leng*(ncoords_sound-1) ; nch = nc/2
              if (nch>0) call safe_character_prepend(item,repeat('-',nch) )
              call safe_character_append(item,repeat('-',nc-nch) )
!
              nch = (leng-2)/2
              do icoor=1,ncoords_sound
                write(sublegend(nsub+nch+1:nsub+nch+2),'(i2)') icoor
                nsub = nsub+leng
              enddo
            endif
            call safe_character_append(legend, item)
!
          endif
          call safe_character_append(fform, trim(cform_sound(iname))//',')
        endif
!
      enddo
!
!  Put output line into a string.
!
      if (ldata) then

        fform = fform(1:len_trim(fform)-1)
        if ( ncoords_sound>1 ) call safe_character_append(fform, ')')
        len=nname_sound
        call compress(fname_sound,cform_sound=='',len)
        write(line,trim(fform)//')') tout, ((fname_sound(i,j), i=1,ncoords_sound), j=1,len)
               !i=1,ncoords_sound), j=1,nname_sound)
               !(1:nname_sound,1:ncoords_sound)
      else
        write(line,tform//')')tout
      endif
!
      call clean_line(line)
      call write_sound_append(legend,sublegend,coorlegend,line,lfirst_call)
!
      if (ldebug) write(*,*) 'exit write_sound'
!
      fname_sound(1:ncoords_sound,1:nname_sound) = 0.0
      lfirst_call = .false.
!
    endsubroutine write_sound
!***********************************************************************
    subroutine get_average_density(mass_per_proc,average_density)
!
!  Calculation of average density.
!
!   1-dec-09/dhruba+piyali: adapted from diagnostic
!   3-may-12/axel+MR: to divide by box_volume, not nw
!
      real :: mass_per_proc,average_density
      intent(out) :: average_density
      intent(in) :: mass_per_proc
!
      real :: mass
!
!  Communicate over all processors.
!
      call mpireduce_sum(mass_per_proc,mass)
!
!  The result is present everywhere
!
      average_density=mass/box_volume
      call mpibcast_real(average_density,comm=MPI_COMM_WORLD)
!
    endsubroutine get_average_density
!**********************************************************************
    subroutine diagnostic(vname,nlname,lcomplex)
!
!  Finalize calculation of diagnostic quantities (0-D).
!
!   2-sep-01/axel: coded
!  14-aug-03/axel: began adding surface integrals
!  26-aug-13/MR:   moved calculation of dVol_rel1 to initialize_diagnostics
!                  added optional parameter lcomplex for use with imaginary part
!
      use General, only: loptest, itoa, safe_character_assign
!
      integer,                 intent(in)   :: nlname
      real, dimension(nlname), intent(inout):: vname
      logical,optional,        intent(in)   :: lcomplex
!
      real, dimension (nlname) :: fmax_tmp, fsum_tmp, fmax, fsum, fweight_tmp
      real :: vol
      integer :: iname, imax_count, isum_count, nmax_count, nsum_count, itype
      logical :: lweight_comm, lalways
      integer, parameter :: lun=1
      character (len=fnlen) :: datadir='data',path=''
      character (len=fnlen) :: fform='(1p30e12.4)'
!
!  Go through all print names, and sort into communicators
!  corresponding to their type.
!
      imax_count=0
      isum_count=0
      lweight_comm=.false.

      lalways=.not.loptest(lcomplex)

      do iname=1,nlname

        itype = itype_name(iname)

        if (lalways.or.itype>=ilabel_complex) then
          if (itype<0) then
            imax_count=imax_count+1
            fmax_tmp(imax_count)=vname(iname)
          elseif (itype>0) then
            itype = mod(itype,ilabel_complex)
            isum_count=isum_count+1
            fsum_tmp(isum_count)=vname(iname)
            if (itype==ilabel_sum_weighted .or. &
                itype==ilabel_sum_weighted_sqrt .or. &
                itype==ilabel_sum_par .or. &
                itype==ilabel_sum_sqrt_par .or. &
                itype==ilabel_sum_log10_par) then
              fweight_tmp(isum_count)=fweight(iname)
              lweight_comm=.true.
            endif
          endif
        endif
      enddo
      nmax_count=imax_count
      nsum_count=isum_count

      if (nmax_count==0 .and. nsum_count==0) return
!
!  Allow for the possibility of writing fsum_tmp/nw for each processor separately.
!
      if (lwrite_fsum) then
        call safe_character_assign(path,trim(datadir)//'/proc'//itoa(iproc))
        open (lun,file=trim(path)//'/fsum.dat',status='unknown',position='append')
        write (lun,fform) t, fsum_tmp/nw
        close(lun)
      endif
!
!  Communicate over all processors.
!
      call mpireduce_max(fmax_tmp,fmax,nmax_count,MPI_COMM_WORLD)
      call mpireduce_sum(fsum_tmp,fsum,nsum_count,comm=MPI_COMM_WORLD)                        ! wrong for Yin-Yang due to overlap
      if (lweight_comm) call mpireduce_sum(fweight_tmp,fweight,nsum_count,comm=MPI_COMM_WORLD)!   ~
!
!  The result is present only on the root processor.
!
      if (lroot) then
!
!  Sort back into original array.
!
        imax_count=0
        isum_count=0
        do iname=1,nlname
!
          itype = itype_name(iname)
          if (lalways.or.itype>=ilabel_complex) then

            if (itype<0) then ! max
              imax_count=imax_count+1
!
              if (itype==ilabel_max)            &
                  vname(iname)=fmax(imax_count)
!
              if (itype==ilabel_max_sqrt)       &
                  vname(iname)=sqrt(fmax(imax_count))
!
              if (itype==ilabel_max_dt)         &
                  vname(iname)=fmax(imax_count)
!
              if (itype==ilabel_max_neg)        &
                  vname(iname)=-fmax(imax_count)
!
              if (itype==ilabel_max_reciprocal) &
                  vname(iname)=1./fmax(imax_count)
!
            elseif (itype>0) then
!
              itype = mod(itype,ilabel_complex)
              isum_count=isum_count+1
!
              if (itype==ilabel_save)            &
                  vname(iname)=fsum(isum_count)
!
              if (itype==ilabel_sum)            &
                  vname(iname)=fsum(isum_count)*dVol_rel1
!
              if (itype==ilabel_sum_sqrt)       &
                  vname(iname)=sqrt(fsum(isum_count)*dVol_rel1)
!
              if (itype==ilabel_sum_log10)       &
                  vname(iname)=log10(fsum(isum_count)*dVol_rel1)
!
              if (itype==ilabel_sum_par)        &
                  vname(iname)=fsum(isum_count)/fweight(isum_count)
!
              if (itype==ilabel_sum_sqrt_par)        &
                  vname(iname)=sqrt(fsum(isum_count))/fweight(isum_count)
!
              if (itype==ilabel_sum_log10_par)        &
                  vname(iname)=log10(fsum(isum_count)/fweight(isum_count))
!
              if (itype==ilabel_integrate.or.itype==ilabel_sum_plain)      &
                  vname(iname)=fsum(isum_count)
!
              if (itype==ilabel_integrate_sqrt)      &
                  vname(iname)=sqrt(fsum(isum_count))
!
              if (itype==ilabel_integrate_log10)      &
                  vname(iname)=log10(fsum(isum_count))
!
              if (itype==ilabel_surf)          &
                  vname(iname)=fsum(isum_count)
!
              if (itype==ilabel_sum_lim) then
                 vol=1.
                 if (lcylinder_in_a_box)  vol=vol*pi*(r_ext**2-r_int**2)
                 if (nzgrid/=1)           vol=vol*Lz
                 if (lsphere_in_a_box)    vol=four_pi_over_three*(r_ext**3-r_int**3)
                 vname(iname)=fsum(isum_count)/vol
              endif
!
              if (itype==ilabel_sum_weighted) then
                if (fweight(isum_count)/=0.0) then
                  vname(iname)=fsum(isum_count)/fweight(isum_count)
                else
                  vname(iname)=0.0
                endif
              endif
!
              if (itype==ilabel_sum_weighted_sqrt) then
                if (fweight(isum_count)/=0.0) then
                  vname(iname)=sqrt(fsum(isum_count)/fweight(isum_count))
                else
                  vname(iname)=0.0
                endif
              endif
!
            endif
          endif
        enddo
!
      endif
!
    endsubroutine diagnostic
!***********************************************************************
    subroutine initialize_time_integrals(f)
!
!  Initialize time_integrals for full chunks.
!
!  28-jun-07/axel+mreinhard: coded
!  29-oct-20/hongzhe: added iuxst, iuyst, iuzst
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(inout) :: f
!
      if (iuut/=0) f(:,:,:,iuxt:iuzt)=0.0
      if (iuust/=0) f(:,:,:,iuxst:iuzst)=0.0
      if (ioot/=0) f(:,:,:,ioxt:iozt)=0.0
      if (ioost/=0) f(:,:,:,ioxst:iozst)=0.0
      if (ibbt/=0) f(:,:,:,ibxt:ibzt)=0.0
      if (ijjt/=0) f(:,:,:,ijxt:ijzt)=0.0
!
    endsubroutine initialize_time_integrals
!***********************************************************************
    subroutine xyaverages_z
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
      real, dimension(nz,nprocz,nnamez) :: fsumz
!
!  Communicate over all processors.
!  The result is only present on the root processor
!
      if (nnamez>0) then
        call mpireduce_sum(fnamez,fsumz,(/nz,nprocz,nnamez/))
        if (lroot) &
            fnamez(:,:,1:nnamez)=fsumz(:,:,1:nnamez)/nxygrid
      endif
!
    endsubroutine xyaverages_z
!***********************************************************************
    subroutine xzaverages_y
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
            fnamey(:,:,1:nnamey)=fsumy(:,:,1:nnamey)/nxzgrid
      endif
!
    endsubroutine xzaverages_y
!***********************************************************************
    subroutine yzaverages_x
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
            fnamex(:,:,1:nnamex)=fsumx(:,:,1:nnamex)/nyzgrid
      endif
!
    endsubroutine yzaverages_x
!***********************************************************************
    subroutine phizaverages_r
!
!  Calculate phiz-averages (still depending on r).
!
!  29-jan-07/wlad: adapted from yzaverages_x and phiaverages_rz
!
      real, dimension (nrcyl,nnamer) :: fsumr
      real, dimension (nrcyl) :: norm
      integer :: in
!
!  Communicate over all processors.
!  The result is only present on the root processor.
!
      if (nnamer>0) then
        !the extra slot is where the normalization is stored
        call mpireduce_sum(fnamer,fsumr,(/nrcyl,nnamer/))
        if (lroot) then
           norm=fsumr(:,nnamer+1)
           do in=1,nnamer
             fnamer(:,in)=fsumr(:,in)/norm
           enddo
        endif
      endif
!
    endsubroutine phizaverages_r
!***********************************************************************
    subroutine yaverages_xz
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
        if (lfirst_proc_y) fnamexz(:,:,1:nnamexz)=fsumxz(:,:,1:nnamexz)/nygrid
      endif
!
    endsubroutine yaverages_xz
!***********************************************************************
    subroutine zaverages_xy
!
!  Calculate z-averages (still depending on x and y).
!
!  19-jun-02/axel: coded
!  25-mar-16/MR: extensions for Yin-Yang grid
!
      use General, only: find_proc
      use Yinyang_mpi, only: reduce_zsum

      real, dimension(:,:,:), allocatable :: fsumxy
      real :: fac
!
      if (nnamexy>0) then

        if (lyang) then
          if (lcaproot) then
            allocate(fsumxy(nnamexy,nx,size(fnamexy_cap,3)))
            fsumxy=0.
          endif
        else
          !if (lfirst_proc_z) allocate(fsumxy(nnamexy,nx,ny))
          allocate(fsumxy(nnamexy,nx,ny))
        endif

        call reduce_zsum(fnamexy,fsumxy)

        fac=1./nzgrid_eff   ! nzgrid_eff approx =4/3*nzgrid for Yin-Yang

        if (.not.lyang) then
          if (lfirst_proc_z) fnamexy=fac*fsumxy
        elseif (lcaproot) then
          fnamexy_cap=fac*fsumxy
        endif

      endif
!
    endsubroutine zaverages_xy
!***********************************************************************
    subroutine phiaverages_rz
!
!  Calculate azimuthal averages (as functions of r_cyl,z).
!  NOTE: these averages depend on (r and) z, so after summation they
!  are still distributed over nprocz CPUs; hence the dimensions of fsumrz
!  (and fnamerz).
!
!  9-dec-02/wolf: coded
!
      integer :: i
      real, dimension (nrcyl,0:nz,nprocz,nnamerz) :: fsumrz
!
!  Communicate over all processors.
!  The result is only present on the root processor
!  normalize by sum of unity which is accumulated in fnamerz(:,0,:,1).
!
      if (nnamerz>0) then
        call mpireduce_sum(fnamerz,fsumrz,(/nrcyl,nz+1,nprocz,nnamerz/))
        if (lroot) then
          do i=1,nnamerz
            fnamerz(:,1:nz,:,i)=fsumrz(:,1:nz,:,i)/spread(fsumrz(:,0,:,1),2,nz)
          enddo
        endif
      endif
!
    endsubroutine phiaverages_rz
!***********************************************************************
    subroutine write_1daverages
!
!  Write 1d averages (z-averages, .., i.e. quantities that are only  1d
!  after averaging). These are written every it1d (default it1) timesteps
!  and appended to their individual files.
!
!   7-aug-03/wolf: coded
!  24-Nov-2018/PABourdin: redesigned
!  09-nov-2022/ccyang: reorganized
!
      use HDF5_IO, only: output_average
      use Mpicomm, only: mpiwtime
!
      logical, save :: lfirst_call = .true.
      logical :: ltimer
      real :: taver
!
      taver = 0.0
      ltimer = ip <= 12 .and. lroot
!
      xya: if (nnamez > 0) then
        if (ltimer) taver = mpiwtime()
        call output_average(datadir, 'xy', nnamez, cnamez, fnamez, nzgrid, t1ddiagnos, lwrite_avg1d_binary, lroot)
        if (ltimer) print *, 'write_1daverages: write xy in ', mpiwtime() - taver, ' seconds'
      endif xya
!
      xza: if (nnamey > 0) then
        if (ltimer) taver = mpiwtime()
        call output_average(datadir, 'xz', nnamey, cnamey, fnamey, nygrid, t1ddiagnos, lwrite_avg1d_binary, lroot)
        if (ltimer) print *, 'write_1daverages: write xz in ', mpiwtime() - taver, ' seconds'
      endif xza
!
      yza: if (nnamex > 0) then
        if (ltimer) taver = mpiwtime()
        call output_average(datadir, 'yz', nnamex, cnamex, fnamex, nxgrid, t1ddiagnos, lwrite_avg1d_binary, lroot)
        if (ltimer) print *, 'write_1daverages: write yz in ', mpiwtime() - taver, ' seconds'
      endif yza
!
      phiz: if (nnamer > 0) then
        if (ltimer) taver = mpiwtime()
        init: if (lfirst_call) then
          call output_average(datadir, 'phi_z', nnamer, cnamer, fnamer, t1ddiagnos, .false., lroot, rcyl)
          lfirst_call = .false.
        else init
          call output_average(datadir, 'phi_z', nnamer, cnamer, fnamer, t1ddiagnos, .false., lroot)
        endif init
        if (ltimer) print *, 'write_1daverages: write phi_z in ', mpiwtime() - taver, ' seconds'
      endif phiz
!
    endsubroutine write_1daverages
!***********************************************************************
    subroutine write_1daverages_prepare(lwrite)
!
!  Prepare l1davg for writing 2D averages.
!  This needs to be done in the beginning of each time step, so
!  the various routines know that they need to calculate averages.
!
!  28-jul-20/joern+axel: adapted from write_2daverages_prepare
!
      use Sub, only: update_snaptime, read_snaptime
!
      logical, intent(IN) :: lwrite
      real, save :: t1davg
      integer, save :: n1davg
      logical, save :: lfirst=.true.
      character (len=fnlen) :: file
      logical :: lwrite_
!
      file=trim(datadir)//'/t1davg.dat'
      if (lfirst) then
        call read_snaptime(trim(file),t1davg,n1davg,d1davg,t)
        lfirst = .false.
      endif
!
!  This routine sets l1davg=T whenever its time to write 2D averages
!    
      lwrite_=lwrite
      call update_snaptime(file,t1davg,n1davg,d1davg,t,lwrite_,ch1davg)
      l1davg=lwrite_
!
    endsubroutine write_1daverages_prepare
!***********************************************************************
    subroutine write_2daverages_prepare(lwrite)
!
!  Prepare l2davg for writing 2D averages.
!  This needs to be done in the beginning of each time step, so
!  the various routines know that they need to calculate averages.
!
!  23-nov-03/axel: adapted from write_2daverages and wvid_prepare
!
      use Sub, only: update_snaptime, read_snaptime
!
      logical, intent(IN) :: lwrite
      real, save :: t2davg
      integer, save :: n2davg
      logical, save :: lfirst=.true.
      character (len=fnlen) :: file
      logical :: lwrite_
!
      file=trim(datadir)//'/t2davg.dat'
      if (lfirst) then
        call read_snaptime(trim(file),t2davg,n2davg,d2davg,t)
        lfirst = .false.
      endif
!
!  This routine sets l2davg=T whenever its time to write 2D averages
!    
      lwrite_=lwrite
      call update_snaptime(file,t2davg,n2davg,d2davg,t,lwrite_,ch2davg)
      l2davg=lwrite_
!
    endsubroutine write_2daverages_prepare
!***********************************************************************
    subroutine write_2daverages
!
!  Write 2d averages (z-averages, phi-averages, .., i.e. quantities that
!  are still 2d after averaging) if it is time.
!  In analogy to 3d output to VARN, the time interval between two writes
!  is determined by a parameter (t2davg) in run.in.
!
!   7-aug-03/wolf: adapted from wsnap
!  24-Nov-2018/PABourdin: redesigned
!  09-nov-2022/ccyang: reorganized
!
      use HDF5_IO, only: output_average
      use IO, only: output_average_2D
      use Mpicomm, only: mpiwtime
!
      logical :: ltimer
      real :: taver
!
      taver = 0.0
      ltimer = ip <= 12 .and. lroot
!
      yavg: if (lwrite_yaverages) then
        if (ltimer) taver = mpiwtime()
        call output_average_2D(directory_dist, 'y', nnamexz, cnamexz, fnamexz, t2davgfirst, lfirst_proc_y)
        if (ltimer) print *, 'write_2daverages: write y averages in ', mpiwtime() - taver, ' seconds'
      endif yavg
!
      zavg: if (lwrite_zaverages .and. (.not. lyang .or. lcaproot)) then
        if (ltimer) taver = mpiwtime()
        if (lcaproot) then
          ! cap root (Yang)
          call output_average_2D(directory_dist, 'z', nnamexy, cnamexy, fnamexy_cap, t2davgfirst, lfirst_proc_z)
        else
          ! z-beam root (Yin)
          call output_average_2D(directory_dist, 'z', nnamexy, cnamexy, fnamexy, t2davgfirst, lfirst_proc_z)
        endif
        if (ltimer) print *, 'write_2daverages: write z averages in ', mpiwtime() - taver, ' seconds'
      endif zavg
!
      phiavg: if (lwrite_phiaverages) then
        if (ltimer) taver = mpiwtime()
        ! normalization is already done in phiaverages_rz
        call output_average(datadir, ch2davg, nrcyl, nnamerz, cnamerz, fnamerz, t2davgfirst, rcyl, drcyl)
        if (ltimer) print *, 'write_2daverages: write phi averages in ', mpiwtime() - taver, ' seconds'
      endif phiavg
!
      if (lroot .and. (ip<=10)) write(*,*) 'write_2daverages: wrote averages (xy,xz,phi#'//trim (ch2davg)//')'
!
    endsubroutine write_2daverages
!***********************************************************************
    subroutine trim_averages
!
!  Trim averages for times past the current time.
!
!  25-apr-16/ccyang: coded
!  23-Nov-2018/PABourdin: redesigned
!
      use IO, only: IO_strategy
      use HDF5_IO, only: trim_average
!
      call trim_average(datadir, 'xy', nzgrid, nnamez)
      call trim_average(datadir, 'xz', nygrid, nnamey)
      call trim_average(datadir, 'yz', nxgrid, nnamex)
!
      if (IO_strategy == "HDF5") then
        call trim_average(datadir, 'y', nxgrid*nzgrid, nnamexz)
        call trim_average(datadir, 'z', nxgrid*nygrid, nnamexy)
        call trim_average(datadir, 'phi', size (rcyl)*nzgrid, nnamerz)
        call trim_average(datadir, 'phi_z', size (rcyl), nnamer)
      endif
!
    endsubroutine trim_averages
!***********************************************************************
    integer function fparse_name(iname,cname,ctest,itest,cform,ncomp)
!
!  Parse name and format of scalar print variable
!  On output, ITEST is set to INAME if CNAME matches CTEST
!  and CFORM is set to the format given as default.
!  E.g. if CTEST='bmax' *i.e. we are testing input line CNAME for 'bmax',
!  CNAME='bmax' will be parsed to ITEST=INAME, CFORM='(1pe10.2)',
!  CNAME='bmax(G5.1)' to ITEST=INAME, CFORM='G5.1',
!  CNAME='brms' to ITEST=<unchanged, normally 0>, CFORM='(1pe10.2)'
!
!   return value is iname if ctest matches
!                   -1    if so and itest/=0 at call (indicates multiple
!                         occurrence of the same diagnostic in print.in)
!                    0    if ctest does not match
!
!   4-may-02/axel: coded
!   6-apr-04/wolf: more liberate format reading
!  26-feb-13/MR  : prepared for ignoring multiple occurrences of
!                  diagnostics in print.in
!  26-aug-13/MR  : removed unneeded 0p setting in format, added handling of D formats
!  27-aug-13/MR  : reinstated 0p
!  10-jan-17/MR  : added correction of floating-point formats if not sufficient to hold sign
!
      use General, only: safe_character_assign, itoa
!
      character (len=*) :: cname, cform
      character (len=*) :: ctest
      integer           :: iname,itest
      integer, optional :: ncomp
!
      intent(in)    :: iname,ctest,ncomp
      intent(inout) :: cname,itest
      intent(out)   :: cform

      integer :: iform0,iform1,iform2,length,index_i,iwidth,idecs,idiff
      character(len=fmtlen) :: tmp
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
        length=iform1-1
      else
        length=iform0-1
      endif
!
!  If the name matches, we keep the name and can strip off the format.
!  The remaining name can then be used for the legend.
!
      if (cname(1:length)==ctest) then

        if (itest==0) then
          if (iform1>0) then
            cform=cname(iform1+1:iform2-1)
          else
            cform='1pE10.2'
          endif
          itest=iname
          fparse_name=iname
        else
          fparse_name=-1
          return          ! diagnostic already defined
        endif

        if (scan(cform(1:1), 'eEdDgG')==1) then
!
!  Increase d in [ED]w.d if insufficient to hold a sign.
!
          if (scan(cform(1:1), "eEdD")==1) then

            index_i=index(cform,'.')
            if (index_i>=3.and.index_i<len(cform)) then
              tmp=cform(2:index_i-1)
              read(tmp,*) iwidth
              tmp=cform(index_i+1:)
              read(tmp,*) idecs
              idiff=iwidth-idecs-7
              if (idiff<0) then
                cform=cform(1:1)//trim(itoa(iwidth-idiff))//cform(index_i:)
!
!  Put changed format back into cname.
!
                if (iform1>0) cname(iform1:)='('//trim(cform)//')'
              endif
            endif
          endif
!
!  Fix annoying Fortran 1p stuff ([EDG]w.d --> 1p[EDG]w.d).
!
          call safe_character_assign(cform, '1p'//trim(cform)//',0p')

        endif
!
!  Integer formats are turned into floating point numbers.
!
        index_i=scan(cform,'iI',.true.)

        if (index_i>0) then
          cform(index_i:index_i)='f'
          cform=trim(cform)//'.0, TL1," "'     ! TL1," " needed, as sometimes, unknown why,
        endif                                  ! the unwanted decimal point does appear.

      else
        fparse_name=0
      endif
!
    endfunction fparse_name
!***********************************************************************
    subroutine parse_name_s(iname,cname,cform,ctest,itest,ncomp)
!
!   subroutine wrapper around fparse_name function: ignores return value
!
!   26-nov-09/MR: coded
!
      character (len=*) :: cname, cform
      character (len=*) :: ctest
      integer :: iname,itest
      integer, optional :: ncomp
!
      intent(in)  :: iname,ctest,ncomp
      intent(out) :: cform
      intent(inout) :: cname,itest
!
      integer :: iret
!
      iret = fparse_name(iname,cname,ctest,itest,cform,ncomp)
!
    endsubroutine parse_name_s
!***********************************************************************
    subroutine parse_name_sa(iname,cname,cform,ctest,itest,ncomp)
!
!   alternative form of parse_name_s: adopts full vectors for cname and cform
!
!   6-mar-13/MR: adapted from parse_name_s
!
      character (len=*), dimension(*) :: cname,cform
      character (len=*) :: ctest
      integer :: iname,itest
      integer, optional :: ncomp
!
      intent(in)  :: iname,ctest,ncomp
      intent(out) :: cform
      intent(inout) :: cname,itest
!
      integer :: iret
!
      iret = fparse_name(iname,cname(iname),ctest,itest,cform(iname),ncomp)
!
    endsubroutine parse_name_sa
!***********************************************************************
    subroutine parse_name_v(cname,cform,cdiag,idiag,ncomp)
!
!   extension to parse_name_s: cname, cform now vectors;
!                              names to look for now in vector cdiag,
!                              indices for found ones in vector idiag;
!                              both do loops incorporated
!
!   29-jun-12/MR: coded
!
      character (len=*), dimension(:) :: cname,cform
      character (len=*), dimension(:) :: cdiag
      integer,           dimension(:) :: idiag
      integer, optional :: ncomp
!
      intent(in)   :: cdiag,ncomp
      intent(inout):: cname
      intent(out)  :: cform,idiag
!
      integer :: i,j
!
      do i=1,size(cname)
        do j=1,size(cdiag)
          if (fparse_name(i,cname(i),cdiag(j),idiag(j),cform(i),ncomp) /= 0) exit
        enddo
      enddo
!
    endsubroutine parse_name_v
!***********************************************************************
    subroutine expand_cname_short(ccname,nname,vlabel,vname,lform)
!
!  Expand string array cname with entries up to index nname such that,
!  if vlabel starts with <vname> or <vname><vname>, it is replaced by the three labels
!  <vname>x..., <vname>y..., <vname>z..., or correspondingly for curvilinear coord systems;
!  updates nname accordingly.
!
!  16-may-12/MR: coded
!
      character (len=*), dimension(:) :: ccname    ! F2003: , allocatable
      character (len=*) :: vlabel
      character (len=*) :: vname
      integer :: nname
      logical, optional :: lform
!
      intent(inout) :: ccname,nname
      intent(in) :: vlabel, vname, lform
!
      character (len=intlen):: tail,form
      integer :: ind, ic,lvn
!
      lvn = len(trim(vname))
      if (vlabel(1:lvn)/=trim(vname)) &
        call fatal_error('expand_cname_short','vlabel does not start with vname')
!
      if (present(lform)) then
        form='T'
      else
        form=''
      endif
!
      ic = name_is_present(ccname,trim(vlabel),form)
      if (ic>0) then
!
        ind=lvn+1
        if ( lvn==1 .and. vlabel(2:2)==vname ) ind=3
!
        tail=vlabel(ind:len(vlabel))
        if ( lspherical_coords ) then
          call expand_cname_full(ccname,nname,ic, &
                                 trim(vname)//'r'//trim(tail), &
                                 trim(vname)//'t'//trim(tail), &
                                 trim(vname)//'p'//trim(tail),form)
        elseif ( lcylindrical_coords ) then
          call expand_cname_full(ccname,nname,ic, &
                                 trim(vname)//'r'//trim(tail), &
                                 trim(vname)//'p'//trim(tail), &
                                 trim(vname)//'z'//trim(tail),form)
        else
          call expand_cname_full(ccname,nname,ic, &
                                 trim(vname)//'x'//trim(tail), &
                                 trim(vname)//'y'//trim(tail), &
                                 trim(vname)//'z'//trim(tail),form)
        endif
      endif
!
    endsubroutine expand_cname_short
!***********************************************************************
    subroutine expand_cname_full(ccname,nname,ic,xlabel,ylabel,zlabel,form)
!
!  Expand string array cname with entries up to index nname such that
!  label at position ic is replaced by the three labels xlabel, ylabel, zlabel;
!  update nname accordingly. Appends format form if present.
!
!   1-apr-04/wolf: coded
!  16-may-12/MR  : new parameters ic = position of label to be expanded in ccname and (optional)
!                  form for a format specification append to the name (for use with
!                  print.in
!
      use General, only : lextend_vector
!
      character (len=*), dimension(:) :: ccname   ! F2003: , allocatable
      character (len=*) :: xlabel,ylabel
      character (len=*), optional :: zlabel
      character (len=intlen), optional :: form
      integer :: nname,ic
!
      intent(inout) :: ccname,nname
      intent(in) :: xlabel,ylabel,zlabel,ic,form
!
      integer :: itot
      logical :: lform
!
      lform=present(form)
      if (lform) then
        if (form=='') lform=.false.
      endif
!
      if (ic<=0) return
!
      if (present(zlabel)) then
        itot=3
      else
        itot=2
      endif
!
      if (.not.lextend_vector(ccname,nname+itot-1)) & ! sanity check
        call fatal_error('expand_cname_full','Not enough memory')
!
      ccname(ic+itot:nname+itot-1) = ccname(ic+1:nname)
      if (lform) then
        ccname(ic) = xlabel//trim(form)
        ccname(ic+1) = ylabel//trim(form)
        if (present(zlabel)) ccname(ic+2) = zlabel//trim(form)
      else
        ccname(ic) = xlabel
        ccname(ic+1) = ylabel
        if (present(zlabel)) ccname(ic+2) = zlabel
      endif
!
      nname = nname+itot-1
!
    endsubroutine expand_cname_full
!***********************************************************************
    subroutine set_type(iname, lsqrt, llog10, lint, lsum, lmax, lmin, lsurf)
!
!  Sets the diagnostic type in itype_name.
!
!  21-sep-17/MR: coded
!
      integer :: iname
      logical, optional :: lsqrt, llog10, lint, lsum, lmax, lmin, lsurf

      if (iname==0) return

      if (present(lsqrt)) &
        itype_name(iname)=ilabel_sum_sqrt
      if (present(lsqrt)) then
        itype_name(iname)=ilabel_sum_sqrt
      elseif (present(llog10)) then
        itype_name(iname)=ilabel_sum_log10
      elseif (present(lint)) then
        itype_name(iname)=ilabel_integrate
      elseif (present(lsurf)) then
        itype_name(iname)=ilabel_surf
      elseif (present(lsum)) then
        itype_name(iname)=ilabel_sum
      elseif (present(lmax)) then
        itype_name(iname)=ilabel_max
      elseif (present(lmin)) then
        itype_name(iname)=ilabel_max_neg
      else
        itype_name(iname)=ilabel_save
      endif

    endsubroutine set_type
!***********************************************************************
    subroutine save_name(a,iname)
!
!  Sets the value of a (must be treated as real) in fname array
!
!  26-may-02/axel: adapted from max_mn_name
!  20-sep-17/MR: removed setting of itype_name as it is set to "save" by default.
!
      real :: a
      integer :: iname
!
      if (iname/=0) fname(iname)=a
!
   endsubroutine save_name
!***********************************************************************
    subroutine save_name_sound(a,iname,iscoord)
!
!  Lists the value of a (must be treated as real) in fname array
!
!  3-Dec-10/dhruba+joern: adapted from max_mn_name
! 25-aug-13/MR: removed unneeded setting of itype, added test of iname.
!
      real :: a
      integer :: iname,iscoord
!
!  This routine is to be called only once per step
!
      if (iname/=0) fname_sound(iscoord,iname)=a
!
   endsubroutine save_name_sound
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
    subroutine max_name_int(a,iname,lneg)
!
!  Successively calculate maximum of a, which is supplied at each call.
!
!  29-aug-05/anders: adapted from save_name
!
      integer :: a, iname
      logical, optional :: lneg
!
      if (iname==0) return

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
    endsubroutine max_name_int
!***********************************************************************
    subroutine max_name_real(a,iname,lneg,l_dt)
!
!  Successively calculate maximum of a, which is supplied at each call.
!
!  13-dec-10/ccyang: adapted from max_name
!
      real, intent(in) :: a
      integer, intent(in) :: iname
      logical, intent(in), optional :: lneg, l_dt
!
      if (iname==0) return

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
      elseif (present(l_dt)) then
        itype_name(iname)=ilabel_max_dt
      else
        itype_name(iname)=ilabel_max
      endif
!
    endsubroutine max_name_real
!***********************************************************************
    subroutine sum_name_real(a,iname)
!
!  Calculate the summation of a, which is supplied at each call.
!
!  19-jun-11/anders: changed to sum single number of all cores
!  17-jun-09/ccyang: adapted from max_name
!  03-sep-09/MR: corrected to real sum
!
      real, intent(in) :: a
      integer, intent(in) :: iname
!
      if (iname==0) return

      fname(iname)=a
!
!  Set corresponding entry in itype_name.
!  Need to set ilabel_surf to avoid multiplication by volume later (to fix).
!
      itype_name(iname)=ilabel_surf
!
    endsubroutine sum_name_real
!***********************************************************************
    subroutine sum_name_int(a,iname)
!
!  Calculate the summation of a, which is supplied at each call.
!
!  19-jun-11/anders: changed to sum single number of all cores
!  17-jun-09/ccyang: adapted from max_name
!  03-sep-09/MR: corrected to real sum
!  12-apr-16/Jørgen+Nils: overloading with int
!
      integer, intent(in) :: a
      integer, intent(in) :: iname
!
      if (iname==0) return

      fname(iname)=a
!
!  Set corresponding entry in itype_name.
!  Need to set ilabel_surf to avoid multiplication by volume later (to fix).
!
      itype_name(iname)=ilabel_surf
!
    endsubroutine sum_name_int
!***********************************************************************
    subroutine max_mn_name(a,iname,lsqrt,l_dt,lneg,lreciprocal)
!
!  Successively calculate maximum of a, which is supplied at each call.
!  Start from zero if lfirstpoint=.true.
!
!   1-apr-01/axel+wolf: coded
!   4-may-02/axel: adapted for fname array
!  23-jun-02/axel: allows for taking square root in the end
!  29-jun-12/MR: incorporated test for iname/=0
!
      real, dimension (nx) :: a
      integer :: iname
      logical, optional :: lsqrt,l_dt,lneg,lreciprocal
!
      if (iname==0) return
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
    subroutine sum_mn_name_arr2(a,iname,lsqrt,llog10,lint,ipart,lplain)
!
!  20-aug-13/MR: derived from sum_mn_name; behaves as before if extent of
!                first dimension of a is 1; if it is 2 considers a to be a complex
!                pencil - stores real(imaginary) part in fname(fname_keep),
!                sets type of diagnostic to complex
!
      real, dimension(:,:), intent(IN) :: a
      integer,              intent(IN) :: iname
      integer, optional,    intent(IN) :: ipart
      logical, optional,    intent(IN) :: lsqrt, llog10, lint, lplain

      if (iname==0) return

      if (size(a,1)==1) then
        call sum_mn_name_std(a(1,:),iname,lsqrt,llog10,lint,ipart,lplain)
      else

        call sum_mn_name_real(a(1,:),iname,fname,lsqrt,llog10,lint,ipart,lplain)
        call sum_mn_name_real(a(2,:),iname,fname_keep)
        if (itype_name(iname) < ilabel_complex ) &
          itype_name(iname)=itype_name(iname)+ilabel_complex

      endif
!
    endsubroutine sum_mn_name_arr2
!***********************************************************************
    subroutine sum_mn_name_std(a,iname,lsqrt,llog10,lint,ipart,lplain)
!
!  20-aug-13/MR: derived from sum_mn_name, behaves as before
!
      use Cdata, only: fname
!
      real, dimension(nx), intent(IN) :: a
      integer,             intent(IN) :: iname
      integer, optional,   intent(IN) :: ipart
      logical, optional,   intent(IN) :: lsqrt, llog10, lint, lplain

      call sum_mn_name_real(a,iname,fname,lsqrt,llog10,lint,ipart,lplain)

    endsubroutine sum_mn_name_std
!***********************************************************************
    subroutine sum_mn_name_real(a,iname,fname,lsqrt,llog10,lint,ipart,lplain)
!
!  Successively calculate sum of a, which is supplied at each call.
!  In subroutine 'diagnostic', the mean is calculated; if 'lint' is
!  explicitly given and true, the integral is calculated, instead.
!  With 'lsqrt=.true.', a square root is applied after building the mean.
!  With 'llog10=.true.', a log10 is applied after building the mean.
!  Start from zero if lfirstpoint=.true.
!  TODO: for nonperiodic arrays we want to multiply boundary data by 1/2.
!  TODO: rename 'ilabel_sum{_sqrt}' as 'ilabel_mean{_sqrt}'.
!
!   1-apr-01/axel+wolf: coded
!   4-may-02/axel: adapted for fname array
!  23-jun-02/axel: allows for taking square root in the end
!  20-jun-07/dhruba: adapted for spherical polar coordinate system
!  30-aug-07/wlad: adapted for cylindrical coordinates
!  22-mar-08/axel: added ladd option, to add to previous values
!  29-aug-11/Bourdin.KIS: added TODO and a comment about building the mean
!  20-aug-13/MR: derived from sum_mn_name, made array of values fname a dummy argument
!  31-mar-15/MR: added switch for proper volume averaging
!
      use Yinyang, only: in_overlap_mask

      real, dimension(nx) :: a,a_scaled
      real, dimension(nname) :: fname
!
      real :: ppart,qpart
      integer :: iname
      integer, optional :: ipart
      logical, optional :: lsqrt, llog10, lint, lplain
!
      intent(in) :: iname

!
!  Only do something if iname is not zero.
!
      if (iname/=0) then
!
!  Only do something if (m,n) is in overlap mask (only relevant for Yin-Yang grid).
!
      if (in_overlap_mask(m,n)) then
!
!  Set corresponding entry in itype_name.
!
        if (present(lsqrt)) then
          itype_name(iname)=ilabel_sum_sqrt
        elseif (present(llog10)) then
          itype_name(iname)=ilabel_sum_log10
        elseif (present(lint)) then
          itype_name(iname)=ilabel_integrate
        elseif (present(lplain)) then
          itype_name(iname)=ilabel_sum_plain
        else
          itype_name(iname)=ilabel_sum
        endif
!
!  Set fraction if old and new stuff.
!
        if (present(ipart)) then
          ppart=1./real(ipart)
          if (lfirstpoint) then
            qpart=1.-ppart
          else
            qpart=1.
          endif
!
!  Add up contributions, taking coordinate system into acount (particles).
!
          if (lspherical_coords) then
            fname(iname)=qpart*fname(iname)+ppart*sum(r2_weight*sinth_weight(m)*a)
          elseif (lcylindrical_coords) then
            fname(iname)=qpart*fname(iname)+ppart*sum(rcyl_weight*a)
          else
            fname(iname)=qpart*fname(iname)+ppart*sum(a)
          endif
!
!  Normal method.
!
        else
          if (lproper_averages) then
            call integrate_mn(a,fname(iname))
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
!  Initialize if one is on the first point, or add up otherwise.
!
            if (lfirstpoint) then
              if (lcartesian_coords.or.lpipe_coords.or.present(lplain)) then
                fname(iname)=sum(a_scaled)
              elseif (lspherical_coords) then
                fname(iname)=sum(r2_weight*a_scaled)*sinth_weight(m)
              elseif (lcylindrical_coords) then
                fname(iname)=sum(rcyl_weight*a_scaled)
              else
                call fatal_error('sum_mn_name_real','not implemented')
              endif
            else
              if (lcartesian_coords.or.lpipe_coords.or.present(lplain)) then
                fname(iname)=fname(iname)+sum(a_scaled)
              elseif (lspherical_coords) then
                fname(iname)=fname(iname)+sum(r2_weight*a_scaled)*sinth_weight(m)
              elseif (lcylindrical_coords) then
                fname(iname)=fname(iname)+sum(rcyl_weight*a_scaled)
              else
                call fatal_error('sum_mn_name_real','not implemented')
              endif
            endif
          endif
        endif
      endif
      endif
!
    endsubroutine sum_mn_name_real
!***********************************************************************
    subroutine sum_mn_name_halfy(a,iname)
!
!  To calculate averages over half the size of box, useful for simulations
!  which includes equator (either cartesian or spherical).
!
!  ??-???-??/dhruba: aped from sum_mn_name
!
      real, dimension (nx) :: a
      real :: sum_name
      integer :: iname
!
      if (iname /= 0) then
        sum_name=0
!
        if (y(m)>=yequator)then
          sum_name=fname_half(iname,2)
        else
          sum_name=fname_half(iname,1)
        endif
!
        if (lfirstpoint) then
          fname_half(iname,1)=0
          fname_half(iname,2)=0
          sum_name=0
          if (lspherical_coords) then
            sum_name=sum(r2_weight*sinth_weight(m)*a)
          elseif (lcylindrical_coords) then
            sum_name=sum(rcyl_weight*a)
          else
            sum_name=sum(a)
          endif
        else
          if (lspherical_coords) then
            sum_name=sum_name+sum(r2_weight*sinth_weight(m)*a)
          elseif (lcylindrical_coords) then
            sum_name=sum_name+sum(rcyl_weight*a)
          else
            sum_name=sum_name+sum(a)
          endif
        endif
!
        if (y(m)>=yequator)then
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
        if (z(n)>=zequator) then
          sum_name=fname_half(iname,2)
        else
          sum_name=fname_half(iname,1)
        endif
!
        if (lfirstpoint) then
          fname_half(iname,1)=0
          fname_half(iname,2)=0
          sum_name=0
          if (lspherical_coords) then
            sum_name=sum(r2_weight*sinth_weight(m)*a)
          elseif (lcylindrical_coords) then
            sum_name=sum(rcyl_weight*a)
          else
            sum_name=sum(a)
          endif
        else
          if (lspherical_coords) then
            sum_name=sum_name+sum(r2_weight*sinth_weight(m)*a)
          elseif (lcylindrical_coords) then
            sum_name=sum_name+sum(rcyl_weight*a)
          else
            sum_name=sum_name+sum(a)
          endif
        endif
!
!  North means 1 and south means 2.
!  However, north corresponds to z < zequator,
!  in analogy to theta < !pi/2 for north.
!
        if (z(n)>=zequator) then
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
!  17-apr-06/anders: coded
!
      real, dimension (:) :: a, weight
      integer :: iname
      logical, optional :: lsqrt
!
      integer, save :: it_save=-1, itsub_save=-1
!
!  Only do something if iname is not zero.
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
        fname(iname)  =fname(iname)  +sum(weight*a)     ! safe_sum(weight*a)
        fweight(iname)=fweight(iname)+sum(weight)       ! safe_sum(weight)
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
    subroutine sum_lim_mn_name(a,iname,p)
!
!  Successively calculate integral of a, which is supplied at each call.
!  Just takes values between r_int < r < r_ext
!  The purpose is to compute quantities just inside a cylinder or sphere
!
!   2-nov-05/wlad: adapted from sum_mn_name
!
      real, dimension (nx) :: a,aux,rlim
      type (pencil_case) :: p
      real :: dv
      integer :: iname,i,isum
!
      if (iname /= 0) then
!
        if (lcylinder_in_a_box) then
          rlim=p%rcyl_mn
        elseif (lsphere_in_a_box) then
          rlim=p%r_mn
        else
          call warning("sum_lim_mn_name","no reason to call it if you are "//&
               "not using a cylinder or a sphere embedded in a "//&
               "Cartesian grid")
        endif
!
         dv=1.
         if (nxgrid/=1) dv=dv*dx
         if (nygrid/=1) dv=dv*dy
         if (nzgrid/=1) dv=dv*dz
!
         do i=1,nx
            if ((rlim(i) <= r_ext).and.(rlim(i) >= r_int)) then
               aux(i) = a(i)
            else
               aux(i) = 0.
            endif
         enddo
!
         if (lfirstpoint) then
            if (lspherical_coords)then
              fname(iname) = 0.
              do isum=l1,l2
                fname(iname)=fname(iname)+ &
                        x(isum)*x(isum)*sinth(m)*aux(isum-nghost)*dv
              enddo
            else
              fname(iname)=sum(aux)*dv
            endif
         else
            if (lspherical_coords)then
              do isum=l1,l2
                fname(iname)=fname(iname)+ &
                      x(isum)*x(isum)*sinth(isum)*aux(isum-nghost)*dv
              enddo
            else
              fname(iname)=fname(iname)+sum(aux)*dv
            endif
         endif
!
         itype_name(iname)=ilabel_sum_lim
!
      endif
!
    endsubroutine sum_lim_mn_name
!*********************************************************
    subroutine surf_mn_name(a,iname)
!
!  Successively calculate surface integral. This routine assumes
!  that "a" contains the partial result for each pencil, so here
!  we just need to add up the contributions from all processors.
!  Start from zero if lfirstpoint=.true.
!
!  14-aug-03/axel: adapted from sum_mn_name
!  15-feb-13/MR: test of iname incorporated
!
      real, intent(in) :: a
      integer, intent(in) :: iname
!
      if (iname>0) then
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
      endif
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
      real, dimension (nx) :: a
      real :: inta
!
      intent(in) :: a
      intent(inout) :: inta
!
      real :: tmp
!
!  initialize by the volume element (which is different for different m and n)
!
      tmp=sum(a*dVol_x(l1:l2))*dVol_y(m)*dVol_z(n)
      !tmp=safe_sum(a*dVol_x(l1:l2))*dVol_y(m)*dVol_z(n)
!
!  initialize if one is on the first point, or add up otherwise
!
      if (lfirstpoint) then
        inta=tmp
      else
        inta=inta+tmp
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
!   31-mar-15/MR: use integrate_mn to avoid doubled code
!
      real, dimension (nx) :: a
      integer :: iname
!
      if (iname==0) return

      call integrate_mn(a,fname(iname))
!
!  set corresponding entry in itype_name
!
      itype_name(iname)=ilabel_integrate
!
    endsubroutine integrate_mn_name
!***********************************************************************
    subroutine xysum_mn_name_z(a,iname)
!
!   3-sep-13/MR: derived from xysum_mn_name_z
!
      use Cdata, only: n
!
      real, dimension(nx), intent(IN) :: a
      integer,             intent(IN) :: iname

      call xysum_mn_name_z_npar(a,n,iname)

    endsubroutine xysum_mn_name_z
!***********************************************************************
    subroutine xysum_mn_name_z_npar(a,n,iname)
!
!  Successively calculate sum over x,y of a, which is supplied at each call.
!  The result fnamez is z-dependent.
!  Start from zero if lfirstpoint=.true.
!
!   5-jun-02/axel: adapted from sum_mn_name
!   3-sep-13/MR: derived from xysum_mn_name_z, index n now parameter
!
      real, dimension(nx), intent(IN) :: a
      integer,             intent(IN) :: n,iname
!
      integer :: isum,lmax,nl
!
!  Only do something if iname is not zero.
!
      if (iname/=0) then
!
!  Initialize to zero, including other parts of the z-array
!  which are later merged with an mpi reduce command.
!
        if (lfirstpoint) fnamez(:,:,iname)=0.0
!
!  n starts with nghost+1=4, so the correct index is n-nghost
!
        lmax=l2; nl=n-nghost
        if (lav_smallx)  lmax=ixav_max
        if (.not.loutside_avg) then
          if (lspherical_coords.or.lcylindrical_coords)then
            do isum=l1,lmax
              fnamez(nl,ipz+1,iname)=fnamez(nl,ipz+1,iname)+ &
                                     x(isum)*a(isum-nghost)
            enddo
          else
            do isum=l1,lmax
              fnamez(nl,ipz+1,iname)=fnamez(nl,ipz+1,iname)+ &
                                     a(isum-nghost)
            enddo
          endif
        endif
      endif
!
    endsubroutine xysum_mn_name_z_npar
!***********************************************************************
    subroutine xzsum_mn_name_y(a,iname)
!
!   3-sep-13/MR: derived from xzsum_mn_name_y
!
      use Cdata, only: m
!
      real, dimension(nx), intent(IN) :: a
      integer,             intent(IN) :: iname

      call xzsum_mn_name_y_mpar(a,m,iname)

    endsubroutine xzsum_mn_name_y
!***********************************************************************
    subroutine xzsum_mn_name_y_mpar(a,m,iname)
!
!  Successively calculate sum over x,z of a, which is supplied at each call.
!  The result fnamey is y-dependent.
!  Start from zero if lfirstpoint=.true.
!
!  12-oct-05/anders: adapted from xysum_mn_name_z
!   3-sep-13/MR: derived from xzsum_mn_name_y, index m now parameter
!
      real, dimension (nx), intent(IN) :: a
      integer,              intent(IN) :: iname,m
!
      integer :: isum,lmax,ml
!
!  Only do something if iname is not zero.
!
      if (iname/=0) then
!
!  Initialize to zero, including other parts of the z-array
!  which are later merged with an mpi reduce command.
!
        if (lfirstpoint) fnamey(:,:,iname)=0.0
!
!  m starts with mghost+1=4, so the correct index is m-nghost.
!
        lmax=l2; ml=m-nghost
!
        if (lav_smallx) lmax=ixav_max
        if (.not.loutside_avg) then
          if (lspherical_coords.and.nxgrid>1)then
            do isum=l1,lmax
              fnamey(ml,ipy+1,iname)=fnamey(ml,ipy+1,iname)+ &
                                     x(isum)*sinth(m)*a(isum-nghost)
            enddo
          else ! also correct for cylindrical
            do isum=l1,lmax
              fnamey(ml,ipy+1,iname)=fnamey(ml,ipy+1,iname)+ &
                                     a(isum-nghost)
            enddo
          endif
        endif
      endif
!
    endsubroutine xzsum_mn_name_y_mpar
!***********************************************************************
    subroutine yzsum_mn_name_x(a,iname)
!
!   3-sep-13/MR: derived from yzsum_mn_name_x
!
      use Cdata, only: m
!
      real, dimension(nx), intent(IN) :: a
      integer,             intent(IN) :: iname

      call yzsum_mn_name_x_mpar(a,m,iname)

    endsubroutine yzsum_mn_name_x
!***********************************************************************
    subroutine yzsum_mn_name_x_mpar(a,m,iname)
!
!  Successively calculate sum over y,z of a, which is supplied at each call.
!  The result fnamex is x-dependent.
!  Start from zero if lfirstpoint=.true.
!
!   2-oct-05/anders: adapted from xysum_mn_name_z
!   3-sep-13/MR: derived from yzsum_mn_name_x, m now parameter
!
      real, dimension (nx), intent(IN) :: a
      integer,              intent(IN) :: iname,m
!
!  Only do something if iname is not zero.
!
      if (iname/=0) then
!
!  Initialize to zero.
!
        if (lfirstpoint) fnamex(:,:,iname)=0.0
!
!  Use different volume differentials for different coordinate systems.
!
        if (lspherical_coords) then
          fnamex(:,ipx+1,iname)=fnamex(:,ipx+1,iname)+x(l1:l2)**2*sinth(m)*a
        elseif (lcylindrical_coords) then
          fnamex(:,ipx+1,iname)=fnamex(:,ipx+1,iname)+x(l1:l2)*a
        else
          fnamex(:,ipx+1,iname)=fnamex(:,ipx+1,iname)+a
        endif
      endif
!
    endsubroutine yzsum_mn_name_x_mpar
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
      integer :: nl
!
      if (iname==0) return
!
!  Initialize to zero, including other parts of the z-array
!  which are later merged with an mpi reduce command.
!
      if (lfirstpoint) fnamez(:,:,iname) = 0.0
!
      fac=1.0
!
      if ((m==m1.and.lfirst_proc_y).or.(m==m2.and.llast_proc_y)) then
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
      nl=n-nghost
      fnamez(nl,ipz+1,iname) = fnamez(nl,ipz+1,iname) + suma
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
      if (iname==0) return
!
!  Initialize to zero, including other parts of the z-array
!  which are later merged with an mpi reduce command.
!
      if (lfirstpoint) fnamey(:,:,iname) = 0.
!
      fac = 1.
!
      if ((n==n1.and.lfirst_proc_z).or.(n==n2.and.llast_proc_z)) then
        if (.not.lperi(3)) fac = .5*fac
      endif
!
      if (lperi(1)) then
        suma = fac*sum(a)
      else
        suma = fac*(sum(a(2:nx-1))+.5*(a(1)+a(nx)))
      endif
!
!  m starts with mghost+1=4, so the correct index is m-nghost.
!
      fnamey(m-nghost,ipy+1,iname) = fnamey(m-nghost,ipy+1,iname) + suma
!
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
      if (iname==0) return
!
!  Initialize to zero.
!
      if (lfirstpoint) fnamex(:,:,iname) = 0.0
!
      fac=1.0
!
      if ((m==m1.and.lfirst_proc_y).or.(m==m2.and.llast_proc_y)) then
        if (.not.lperi(2)) fac = .5*fac
      endif
!
      if ((n==n1.and.lfirst_proc_z).or.(n==n2.and.llast_proc_z)) then
        if (.not.lperi(3)) fac = .5*fac
      endif
!
      fnamex(:,ipx+1,iname) = fnamex(:,ipx+1,iname) + fac*a
!
    endsubroutine yzintegrate_mn_name_x
!***********************************************************************
    subroutine phizsum_mn_name_r(a,iname)
!
!  Successively calculate sum over phi,z of a, which is supplied at each call.
!  Start from zero if lfirstpoint=.true.
!  The fnamer array uses one of its slots in mnamer where we put ones and sum
!  them up in order to get the normalization correct.
!
!  29-jan-07/wlad: adapted from yzsum_mn_name_x and phisum_mn_name
!
      real, dimension (nx) :: a
      integer :: iname,ir,nnghost
!
      if (iname==0) return
!
      if (lfirstpoint) fnamer(:,iname)=0.
      if (lfirstpoint.and.iname==nnamer) fnamer(:,iname+1)=0.
!
      do ir=1,nrcyl
        fnamer(ir,iname) = fnamer(ir,iname) + sum(a*phiavg_profile(ir,:))
      enddo
!
!  Normalization factor, just needs to be done once.
!  As is it a z-average, multiply by nz afterwards.
!
      nnghost=n-nghost
      if ((iname==nnamer).and.(nnghost==1)) then
!  Check if an extra slot is available on fnamer.
        if (nnamer==mnamer) call fatal_error('phizsum_mn_name_r', &
            'no slot for phi-normalization. decrease nnamer')
!
        do ir=1,nrcyl
          fnamer(ir,iname+1)= &
              fnamer(ir,iname+1) + sum(1.*phiavg_profile(ir,:))*nz
        enddo
      endif
!
    endsubroutine phizsum_mn_name_r
!***********************************************************************
    subroutine ysum_mn_name_xz(a,iname)
!
!   3-sep-13/MR: derived from ysum_mn_name_xz
!
      use Cdata, only: n
!
      real, dimension(nx), intent(IN) :: a
      integer,             intent(IN) :: iname

      call ysum_mn_name_xz_npar(a,n,iname)

    endsubroutine ysum_mn_name_xz
!***********************************************************************
    subroutine ysum_mn_name_xz_npar(a,n,iname)
!
!  Successively calculate sum over y of a, which is supplied at each call.
!  The result fnamexz is xz-dependent.
!  Start from zero if lfirstpoint=.true.
!
!   7-jun-05/axel: adapted from zsum_mn_name_xy
!   3-sep-13/MR: derived from ysum_mn_name_xz, n now parameter
!
      real, dimension (nx), intent(IN) :: a
      integer,              intent(IN) :: iname,n
!
      integer :: nl
!
      if (iname==0) return
!
!  Initialize to zero, including other parts of the xz-array
!  which are later merged with an mpi reduce command.
!
      if (lfirstpoint) fnamexz(:,:,iname)=0.
!
!  n starts with nghost+1=4, so the correct index is n-nghost.
!
      nl=n-nghost
!
      if (lspherical_coords.or.lcylindrical_coords)then
        fnamexz(:,nl,iname) = fnamexz(:,nl,iname)+a*x(l1:l2)
      else
        fnamexz(:,nl,iname) = fnamexz(:,nl,iname)+a
      endif
!
    endsubroutine ysum_mn_name_xz_npar
!***********************************************************************
    subroutine zsum_mn_name_xy_scal(a,iname,lint)
!
!  Successively calculate sum over z of a, which is supplied at each call.
!  The result fnamexy is xy-dependent.
!  Start from zero if lfirstpoint=.true.
!     
!  19-jun-02/axel: adapted from xysum_mn_name
!  08-feb-12/ccyang: add option for integration
!   3-sep-13/MR: outsourced zsum_mn_name_xy_mpar
!     
      use General, only: loptest
      use Cdata,   only: n,m,nzgrid_eff 
!
      real, dimension(nx), intent(in) :: a
      integer,             intent(in) :: iname
      logical,   optional, intent(in) :: lint

      if (iname==0) return
!         
!  Scale factor for integration
!       
      if (loptest(lint)) then
        call zsum_mn_name_xy_mpar((real(nzgrid_eff)*zprim(n))*a,m,iname)
      else
        call zsum_mn_name_xy_mpar(a,m,iname)
      endif                                                                                                 
!
    endsubroutine zsum_mn_name_xy_scal
!***********************************************************************
    subroutine zsum_mn_name_xy_arr2(arr,iname)
!
!  Stores multi-component diagnostics. 
!
!  22-sep-20/MR: adapted from zsum_mn_name_xy_mpar
!
      use Cdata, only: m
!
      real,    dimension(:,:,:),      intent(in) :: arr
      integer,                        intent(in) :: iname
!
      integer :: ml,ind,i,j

      if (iname==0) return
!
      if (lfirstpoint) fnamexy(iname,:,:)=0.

      ml=m-nghost

      ind=iname
      do i=1,size(arr,2); do j=1,size(arr,3)
        fnamexy(ind,:,ml)=fnamexy(ind,:,ml)+arr(:,i,j)
        ind=ind+1
      enddo; enddo

    endsubroutine zsum_mn_name_xy_arr2
!***********************************************************************
    subroutine zsum_mn_name_xy_arr(arr,iname)
!
!  Stores multi-component diagnostics. 
!
!  22-sep-20/MR: adapted from zsum_mn_name_xy_mpar
!
      use Cdata, only: m
!
      real,    dimension(:,:),        intent(in) :: arr
      integer,                        intent(in) :: iname
!
      integer :: ml,i

      if (iname==0) return
!
      if (lfirstpoint) fnamexy(iname,:,:)=0.

      ml=m-nghost

      do i=iname,iname+size(arr,2)-1
        fnamexy(i,:,ml)=fnamexy(i,:,ml)+arr(:,i)
      enddo

    endsubroutine zsum_mn_name_xy_arr
!***********************************************************************
    subroutine zsum_mn_name_xy_vec(avec,iname,powers,scal,lint)
!
!  Wrapper for zsum_mn_name_xy_mpar_vec..
!
!  19-jun-02/axel: adapted from xysum_mn_name
!  08-feb-12/ccyang: add option for integration
!   3-sep-13/MR: outsourced zsum_mn_name_xy_mpar
!  31-mar-16/MR: derived from zsum_mn_name_xy
!
      use General, only: loptest
      use Cdata,   only: n,m
!
      real,    dimension(nx,3),        intent(in) :: avec
      integer,                         intent(in) :: iname
      integer, dimension(3),           intent(in) :: powers
      real,    dimension(:), optional, intent(in) :: scal
      logical,               optional, intent(in) :: lint
!
      if (iname==0) return
      if (loptest(lint)) then
!
!  Mulitply with scale factor for integration.
!  Not correct for Yin-Yang with non-equidistant z grid (would this happen?). 
!
        if (present(scal)) then
          call zsum_mn_name_xy_mpar(avec,m,iname,powers,(real(nzgrid_eff)*zprim(n))*scal)
        else
          call zsum_mn_name_xy_mpar(avec,m,iname,powers,(/real(nzgrid_eff)*zprim(n)/))
        endif
      else
        call zsum_mn_name_xy_mpar(avec,m,iname,powers,scal)
      endif
!
    endsubroutine zsum_mn_name_xy_vec
!***********************************************************************
    subroutine zsum_mn_name_xy_mpar_vec(avec,m,iname,powers,scal)
!
!  Calculates a general tensor component from vector avec as 
!  avec(1)**power(1)*avec(2)**power(2)*avec(3)**power(3),
!  optionally the scalar scal (a pencil or one-element vector)
!  is multiplied finally.
!  On Yang procs, avec is first transformed properly.
!
!  31-mar-16/MR: coded
!
      use General, only: transform_thph_yy

      real, dimension(nx,3),        intent(in) :: avec
      integer,                      intent(in) :: iname, m
      integer, dimension(3),        intent(in) :: powers
      real, dimension(:), optional, intent(in) :: scal
!
      real, dimension(nx,3) :: work
      real, dimension(nx) :: quan
      integer :: i
      logical :: lfirst

      if (iname==0) return

      if (lyang) then
!
! On Yang procs: transform theta and phi components if necessary. 
!
        call transform_thph_yy(avec,powers,work)
      else
        work=avec
      endif
!
!  Perform product of powers of vector components..
!
      lfirst=.true.
      do i=1,3
        if (powers(i)/=0) then
          if (lfirst) then
            lfirst=.false.
            if (powers(i)==1) then
              quan=work(:,i)
            else
              quan=work(:,i)**powers(i)
            endif
          else
            if (powers(i)==1) then
              quan=quan*work(:,i)
            else
              quan=quan*work(:,i)**powers(i)
            endif
          endif
        endif
      enddo
!
!  Multiply with scalar if present.
!
      if (present(scal)) then
        if (size(scal)==1) then
          quan=quan*scal(1)
        else
          quan=quan*scal
        endif
      endif
!
!  Sum up result like a scalar.
!
      call zsum_mn_name_xy_mpar(quan,m,iname)
!
    endsubroutine zsum_mn_name_xy_mpar_vec
!***********************************************************************
    subroutine zsum_mn_name_xy_mpar_scal(a,m,iname)
!
!  Accumulates contributions to z-sum within an mn-loop
!  which are later merged with an mpi reduce command.
!  In the Yang part of a Yin-Yang grid, the contributions
!  to the sum along a (extended) phi coordinate line of the Yin 
!  grid are calculated by use of predetermined weights.
!  Likewise for those coordinate lines of the Yin-phi which lie completely
!  within the Yang grid (i.e. the polar caps).
!
!   3-apr-16/MR: derived from zsum_mn_name_xy_mpar; extensions for Yin-Yang grid
!   7-jun-16/MR: outsourced initialize_zaver_yy, reduce_zsum, zsum_y and 
!                corresponding Yin-Yang specific data to Yinyang_mpi
!
      use Cdata, only: n
      use Yinyang_mpi, only: zsum_yy

      real, dimension(nx), intent(in) :: a
      integer,             intent(in) :: iname, m

      integer :: ml
!
      if (iname==0) return
!
      if (lfirstpoint) fnamexy(iname,:,:)=0.

      if (lyang) then
!
!  On Yang grid:
!
        call zsum_yy(fnamexy,iname,m,n,a)
      else
!
! Normal summing-up in Yin procs.
! m starts with nghost+1, so the correct index is m-nghost.
! 
        ml=m-nghost
        fnamexy(iname,:,ml)=fnamexy(iname,:,ml)+a
      endif
!
    endsubroutine zsum_mn_name_xy_mpar_scal
!***********************************************************************
    subroutine calc_phiavg_profile(p)
!
!  Calculate profile for phi-averaging for given pencil.
!
!   2-feb-03/wolf: coded
!
      type (pencil_case) :: p
      real :: r0,width
      integer :: ir
!
      intent(in) :: p
!
!  We use a quartic-Gaussian profile ~ exp(-r^4)
!
!      width = .5*drcyl
      width = .7*drcyl
      do ir=1,nrcyl
        r0 = rcyl(ir)
        phiavg_profile(ir,:) = exp(-0.5*((p%rcyl_mn-r0)/width)**4)
      enddo
!
      if (.not.(lcylinder_in_a_box.or.lsphere_in_a_box)) &
          call warning('calc_phiavg_profile', &
          'no reason to call it if you are '//&
          'not using a cylinder or a sphere embedded in a '//&
          'Cartesian grid')
!
    endsubroutine calc_phiavg_profile
!***********************************************************************
    subroutine phisum_mn_name_rz(a,iname)
!
!  Successively calculate sum over phi of a, which is supplied at each call.
!  Start from zero if lfirstpoint=.true.
!  The fnamerz array has one extra slice in z where we put ones and sum
!  them up in order to get the normalization correct.
!
!   2-feb-03/wolf: adapted from xysum_mn_name_z
!
      real, dimension (nx) :: a
      integer :: iname,n_nghost,ir
!
      if (iname /= 0) then
!
!  Initialize to zero, including other parts of the rz-array
!  which are later merged with an mpi reduce command.
!  At least the root processor needs to reset all ipz slots, as it uses
!  fnamerz(:,:,:,:) for the final averages to write [see
!  phiaverages_rz()]; so we better reset everything:
!      if (lfirstpoint) fnamerz(:,:,ipz+1,iname) = 0.
        if (lfirstpoint) fnamerz(:,:,:,iname) = 0.
!
!  n starts with nghost+1=4, so the correct index is n-nghost
!
        n_nghost=n-nghost
        do ir=1,nrcyl
          fnamerz(ir,n_nghost,ipz+1,iname) &
               = fnamerz(ir,n_nghost,ipz+1,iname) + sum(a*phiavg_profile(ir,:))
        enddo
!
!  sum up ones for normalization; store result in fnamerz(:,0,:,1)
!  Only do this for the first n, or we would sum up nz times too often
!
        if (iname==1 .and. n_nghost==1) then
          do ir=1,nrcyl
            fnamerz(ir,0,ipz+1,iname) &
                 = fnamerz(ir,0,ipz+1,iname) + sum(1.*phiavg_profile(ir,:))
          enddo
        endif
!
      endif
!
    endsubroutine phisum_mn_name_rz
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
    subroutine allocate_sound(nnamel)
!
!  Allocate the variables needed for "sound".
!
!   3-Dec-10/dhruba+joern: coded
!   11-jan-11/MR: parameter nnamel added
!
      use File_io, only : parallel_unit_vec, parallel_open, parallel_close, parallel_count_lines
      use General, only : itoa
      use Sub, only     : location_in_proc
!
      integer, intent(in) :: nnamel
!
      character (LEN=*), parameter :: sound_coord_file = 'sound.coords'
      integer :: isound
      logical :: lval
      integer :: ierr, il
      integer :: mcoords_sound
      integer, allocatable, dimension (:,:) :: sound_inds
      real, dimension(dimensionality) :: coords
      integer :: lsound, msound, nsound, nitems
!
!  Allocate and initialize to zero. Setting it to zero is only
!  necessary because of the pencil test, which doesn't compute these
!  averages, and only evaluates its output for special purposes
!  such as computing mean field energies in calc_bmz, for example,
!
      mcoords_sound = parallel_count_lines(sound_coord_file)

      allocate(sound_inds(mcoords_sound,3),stat=ierr)
      if (ierr>0) call fatal_error('allocate_sound', &
                                   'Could not allocate memory for sound_inds')
!
      ncoords_sound = 0; nitems = 0
!
      call parallel_open(sound_coord_file,nitems=nitems)
!
      do isound=1,mcoords_sound
!
        read(parallel_unit_vec,*) coords
!
        if (location_in_proc(coords,lsound,msound,nsound)) then
!
          lval = .true.
          do il = 1, ncoords_sound
            if ((sound_inds(il,1) == lsound) .and. &
                (sound_inds(il,2) == msound) .and. &
                (sound_inds(il,3) == nsound) ) then
              lval = .false.
              exit
            endif
          enddo
!
          if (lval) then
            ncoords_sound = ncoords_sound + 1
            sound_inds(ncoords_sound,1) = lsound
            sound_inds(ncoords_sound,2) = msound
            sound_inds(ncoords_sound,3) = nsound
          endif
        endif
      enddo
!
      call parallel_close
!
! cname_sound is allocated also for ncoords_sound=0 as needed in read_name_format.
!
      allocate(cname_sound(nnamel),stat=ierr)
      if (ierr>0) call fatal_error('allocate_sound', &
                                   'Could not allocate memory for cname_sound')
!
      lwrite_sound = ncoords_sound>0
      if (lwrite_sound) then
!
        allocate(sound_coords_list(ncoords_sound,3),stat=ierr)
        if (ierr>0) call fatal_error('allocate_sound', &
                                     'Could not allocate memory for sound_coords_list')
        sound_coords_list = sound_inds(1:ncoords_sound,:)
!
        allocate(fname_sound(ncoords_sound,nnamel),stat=ierr)
        if (ierr>0) call fatal_error('allocate_sound', &
                                     'Could not allocate memory for fname_sound')
        if (ldebug) print*, 'allocate_sound: allocated memory for '// &
                            'fname_sound  with nname_sound  =', nnamel
!
        fname_sound = 0.0
!
        allocate(cform_sound(nnamel),stat=ierr)
        if (ierr>0) call fatal_error('allocate_sound', &
                                     'Could not allocate memory for cform_sound')
        cform_sound = ' '
!
      endif
!
!  Now deallocate the temporary memory.
!
      deallocate(sound_inds)
!
    endsubroutine allocate_sound
!***********************************************************************
    subroutine sound_clean_up
!
!  Frees up the memory allocated for sound.
!
!  ??-???-??/dhruba: coded
!
      if (allocated(sound_coords_list)) deallocate(sound_coords_list)
      if (allocated(fname_sound)) deallocate(fname_sound)
      if (allocated(cname_sound)) deallocate(cname_sound)
      if (allocated(cform_sound)) deallocate(cform_sound)
!
    endsubroutine sound_clean_up
!***********************************************************************
    subroutine allocate_fnames(nnamel)
!
!  Allocate arrays needed for diagnostics.
!
!   23-mar-10/Bourdin.KIS: copied from allocate_yaverages
!   11-jan-11/MR: parameter nnamel added
!   18-aug-13/MR: accumulation of diagnostics enabled
!   25-aug-13/MR: added allocation of itype_name
!
      integer, intent(in) :: nnamel
!
      integer :: stat
!
      allocate(cname(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_fnames', &
                                   'Could not allocate memory for cname')
      if (ldebug) print*, 'allocate_fnames    : allocated memory for '// &
                          'cname   with nname   =', nnamel
      cname=''
!
      allocate(fname(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_fnames', &
                                   'Could not allocate memory for fname')
      if (ldebug) print*, 'allocate_fnames    : allocated memory for '// &
                          'fname   with nname   =', nnamel
      allocate(fname_keep(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_fnames', &
                                   'Could not allocate memory for fname_keep')
      fname=0.0
      fname_keep=0.0
!
      allocate(cform(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_fnames', &
                                   'Could not allocate memory for cform')
      if (ldebug) print*, 'allocate_fnames    : allocated memory for '// &
                          'cform   with nname   =', nnamel
      cform=''
!
      allocate(itype_name(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_fnames', &
                                   'Could not allocate memory for itype_name')
      if (ldebug) print*, 'allocate_fnames    : allocated memory for '// &
                          'itype_name with nname   =', nnamel
      itype_name=ilabel_save

    endsubroutine allocate_fnames
!***********************************************************************
    subroutine allocate_vnames(nnamel)
!
!  Allocate arrays needed for video slices.
!
!   23-mar-10/Bourdin.KIS: copied from allocate_yaverages
!   11-jan-11/MR: parameter nnamel added
!
      integer, intent(in) :: nnamel
!
      integer :: stat
!
      allocate(cnamev(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_vnames', &
          'Could not allocate memory for cnamev')
      if (ldebug) print*, 'allocate_vnames    : allocated memory for '// &
          'cnamev  with nnamev  =', nnamel
      cnamev=''
!
      if (allocated(cformv)) deallocate(cformv)
      allocate(cformv(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_vnames', &
          'Could not allocate memory for cformv')
      if (ldebug) print*, 'allocate_vnames    : allocated memory for '// &
          'cformv   with nname   =', nnamel
      cformv=''
!
    endsubroutine allocate_vnames
!***********************************************************************
    subroutine allocate_xyaverages(nnamel)
!
!  Allocate arrays needed for xy-averages.
!
!   24-nov-09/anders: copied from allocate_yaverages
!   11-jan-11/MR: parameter nnamel added
!
      integer, intent(in) :: nnamel
!
      integer :: stat
!
!  Allocate and initialize to zero. Setting it to zero is only
!  necessary because of the pencil test, which doesn't compute these
!  averages, and only evaluates its output for special purposes
!  such as computing mean field energies in calc_bmz, for example,
!
      allocate(cnamez(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_xyaverages', &
                                   'Could not allocate memory for cnamez')
      if (ldebug) print*, 'allocate_xyaverages: allocated memory for '// &
                          'cnamez  with nnamez  =', nnamel
      cnamez=''
!
      allocate(fnamez(nz,nprocz,nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_xyaverages', &
                                   'Could not allocate memory for fnamez')
      if (ldebug) print*, 'allocate_xyaverages: allocated memory for '// &
                          'fnamez  with nnamez  =', nnamel
      fnamez=0.0
!
      allocate(cformz(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_xyaverages', &
                                   'Could not allocate memory for cformz')
      if (ldebug) print*, 'allocate_xyaverages: allocated memory for '// &
                          'cformz  with nnamez  =', nnamel
      cformz=''
!
    endsubroutine allocate_xyaverages
!***********************************************************************
    subroutine allocate_xzaverages(nnamel)
!
!  Allocate arrays needed for xz-averages.
!
!   24-nov-09/anders: copied from allocate_yaverages
!   11-jan-11/MR: parameter nnamel added
!
      integer, intent(in) :: nnamel
!
      integer :: stat
!
      allocate(cnamey(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_xzaverages', &
                                   'Could not allocate memory for cnamey')
      if (ldebug) print*, 'allocate_xzaverages: allocated memory for '// &
                          'cnamey  with nnamey  =', nnamel
      cnamey=''
!
      allocate(fnamey(ny,nprocy,nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_xzaverages', &
                                   'Could not allocate memory for fnamey', .true.)
      if (ldebug) print*, 'allocate_xzaverages: allocated memory for '// &
                          'fnamey  with nnamey  =', nnamel
      fnamey=0.0
!
      allocate(cformy(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_xzaverages', &
                                   'Could not allocate memory for cformy', .true.)
      if (ldebug) print*, 'allocate_xzaverages: allocated memory for '// &
                          'cformy  with nnamey  =', nnamel
      cformy=''
!
    endsubroutine allocate_xzaverages
!***********************************************************************
    subroutine allocate_yzaverages(nnamel)
!
!  Allocate arrays needed for yz-averages.
!
!   24-nov-09/anders: copied from allocate_yaverages
!   11-jan-11/MR: parameter nnamel added
!
      integer, intent(in) :: nnamel
!
      integer :: stat
!
      allocate(cnamex(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_yzaverages', &
                                   'Could not allocate memory for cnamex')
      if (ldebug) print*, 'allocate_yzaverages: allocated memory for '// &
                          'cnamex  with nnamex  =', nnamel
      cnamex=''
!
      allocate(fnamex(nx,nprocx,nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_yzaverages', &
                                   'Could not allocate memory for fnamex')
      if (ldebug) print*, 'allocate_yzaverages: allocated memory for '// &
                          'fnamex  with nnamex  =', nnamel
      fnamex=0.0
!
      allocate(cformx(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_yzaverages', &
                                   'Could not allocate memory for cformx')
      if (ldebug) print*, 'allocate_yzaverages: allocated memory for '// &
                          'cformx  with nnamex  =', nnamel
      cformx=''
!
    endsubroutine allocate_yzaverages
!***********************************************************************
    subroutine allocate_phizaverages(nnamel)
!
!  Allocate arrays needed for phiz-averages.
!
!   24-nov-09/anders: copied from allocate_yaverages
!   11-jan-11/MR: parameter nnamel added
!
      integer, intent(in) :: nnamel
!
      integer :: stat
!
      allocate(cnamer(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_phizaverages', &
                                   'Could not allocate memory for cnamer')
      if (ldebug) print*, 'allocate_phizaverages: allocated memory for '// &
                          'cnamer  with nnamel =', nnamel
      cnamer=''
!
      mnamer=nnamel+1
      allocate(fnamer(nrcyl,mnamer),stat=stat)
      if (stat>0) call fatal_error('allocate_phizaverages', &
                                   'Could not allocate memory for fnamer')
      if (ldebug) print*, 'allocate_phizaverages: allocated memory for '// &
                          'fnamer  with nnamer+1 =', mnamer
      fnamer=0.0
!
      allocate(cformr(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_phizaverages', &
                                   'Could not allocate memory for cformr')
      if (ldebug) print*, 'allocate_phizaverages: allocated memory for '// &
                          'cformr  with nnamel =', nnamel
      cformr=''
!
    endsubroutine allocate_phizaverages
!***********************************************************************
    subroutine allocate_yaverages(nnamel)
!
!  Allocate arrays needed for y-averages.
!
!   12-aug-09/dhruba: coded
!   11-jan-11/MR: parameter nnamel added
!
      integer, intent(in) :: nnamel
!
      integer :: stat
!
      allocate(cnamexz(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_yaverages', &
                                   'Could not allocate memory for cnamexz')
      if (ldebug) print*, 'allocate_yaverages : allocated memory for '// &
                          'cnamexz with nnamexz =', nnamel
      cnamexz=''
!
      allocate(fnamexz(nx,nz,nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_yaverages', &
                                   'Could not allocate memory for fnamexz')
      if (ldebug) print*, 'allocate_yaverages : allocated memory for '// &
                          'fnamexz with nnamexz =', nnamel
      fnamexz=0.0
!
      allocate(cformxz(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_yaverages', &
                                   'Could not allocate memory for cformxz')
      if (ldebug) print*, 'allocate_yaverages : allocated memory for '// &
                          'cformxz with nnamexz =', nnamel
      cformxz=''
!
    endsubroutine allocate_yaverages
!*******************************************************************
    subroutine allocate_zaverages(nnamel)
!
!  Allocate arrays needed for z-averages.
!
!   12-aug-09/dhruba: coded
!   11-jan-11/MR: parameter nnamel added
!   11-mar-16/MR: outsourced allocation of fnamexy to allocate_zaverages_data
!
      integer, intent(in) :: nnamel
!
      integer :: stat
!
      allocate(cnamexy(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_zaverages', &
                                   'Could not allocate memory for cnamexy')
      if (ldebug) print*, 'allocate_zaverages : allocated memory for '// &
                          'cnamexy with nnamexy =', nnamel
      cnamexy=''
!
      allocate(cformxy(nnamel),stat=stat)
      if (stat>0) &
        call fatal_error('allocate_zaverages', &
                         'Could not allocate memory for cformxy')
      if (ldebug) print*, 'allocate_zaverages : allocated memory for '// &
                          'cformxy with nnamel =', nnamel
      cformxy=''

    endsubroutine allocate_zaverages
!*******************************************************************
    subroutine allocate_zaverages_data(nnamel)
!
!  Allocate data array needed for z-averages.
!  Additional fnamexy_cap for collection of z-sums in polar caps of Yang grid.
!
!   11-mar-16/MR: outsourced from allocate_zaverages
!
      use Yinyang_mpi, only: initialize_zaver_yy

      integer, intent(in) :: nnamel
!
      integer :: stat, nyl, nycap
!
      nyl=ny
      if (lyinyang) call initialize_zaver_yy(nyl,nycap)
!
      allocate(fnamexy(nnamel,nx,nyl),stat=stat)
      if (stat>0) call fatal_error('allocate_zaverages_data', &
                                   'Could not allocate memory for fnamexy')
      if (ldebug) print*, 'allocate_zaverages_data: allocated memory for '// &
                          'fnamexy with nnamexy =', nnamel
      fnamexy=0.

      if (lcaproot) then
        allocate(fnamexy_cap(nnamel,nx,nycap),stat=stat)
        if (stat>0) call fatal_error('allocate_zaverages_data', &
                                     'Could not allocate memory for fnamexy_cap')
        if (ldebug) print*, 'allocate_zaverages_data: allocated memory for '// &
                            'fnamexy_cap with nycap, nnamexy =', nycap, nnamel
      endif
!
    endsubroutine allocate_zaverages_data
!*******************************************************************
    subroutine allocate_phiaverages(nnamel)
!
!  Allocate arrays needed for phi-averages.
!
!   24-nov-09/anders: copied from allocate_zaverages
!   11-jan-11/MR: parameter nnamel=iadd+nnamerz instead of nnamerz
!
      integer, intent(in) :: nnamel
!
      integer :: stat
!
      allocate(cnamerz(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_phiaverages', &
                                   'Could not allocate memory for cnamerz')
      if (ldebug) print*, 'allocate_phiaverages : allocated memory for '// &
                          'cnamerz with nnamerz =', nnamel
      cnamerz=''
!
      allocate(fnamerz(nrcyl,0:nz,nprocz,nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_phiaverages', &
                                   'Could not allocate memory for fnamerz')
      if (ldebug) print*, 'allocate_phiaverages : allocated memory for '// &
                          'fnamerz with nnamerz =', nnamel
      fnamerz=0.0
!
      allocate(cformrz(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_phiaverages', &
                                   'Could not allocate memory for cformrz')
      if (ldebug) print*, 'allocate_phiaverages : allocated memory for '// &
                          'cformrz with nnamerz =', nnamel
      cformrz=''
!
    endsubroutine allocate_phiaverages
!***********************************************************************
    subroutine sign_masked_xyaver(quan,idiag,ncount)
!
!  Forms sign-masked averaging over xy-planes (only positive values count).
!  ncount holds the number of points which contribute.
!
!  28-sep-16/MR: coded
!
      use Mpicomm, only: mpiallreduce_sum_int, IXYPLANE

      real   ,dimension(nx),intent(IN)   :: quan
      integer,              intent(IN)   :: idiag
      integer,dimension(nz),intent(INOUT):: ncount

      real, dimension(nx) :: buf
      integer, dimension(:), allocatable :: nsum

      if (idiag>0) then

        where (quan>0.)
          buf=quan
        elsewhere
          buf=0.
        endwhere

        ncount(n-n1+1)=ncount(n-n1+1)+count(buf>0.)
        call xysum_mn_name_z(buf,idiag)

        if (llastpoint) then
!
!  Find total number of contributing points in xy-plane.
!
          allocate(nsum(nz))
          call mpiallreduce_sum_int(ncount,nsum,nz,IXYPLANE)
!
!  Form average by dividing by this number. Multiplication with nxygrid
!  necessary as later the average is divided by that.
!
          where (nsum>0) &
            fnamez(:,ipz,idiag)=fnamez(:,ipz,idiag)/nsum*nxygrid
          ncount=0
        endif
      endif

    endsubroutine sign_masked_xyaver
!***********************************************************************
    subroutine diagnostics_clean_up
!
!   16-jan-17/MR: coded
!
                               call vnames_clean_up
                               call fnames_clean_up
                               call xyaverages_clean_up
                               call xzaverages_clean_up
                               call yzaverages_clean_up
      if (lwrite_phizaverages) call phizaverages_clean_up
      if (lwrite_yaverages)    call yaverages_clean_up
      if (lwrite_zaverages)    call zaverages_clean_up
      if (lwrite_phiaverages)  call phiaverages_clean_up
      if (lwrite_sound)        call sound_clean_up
    
    endsubroutine diagnostics_clean_up
!***********************************************************************
    subroutine fnames_clean_up
!
!  Deallocate space needed for reading the print.in file.
!
!   20-apr-10/Bourdin.KIS: copied from xyaverages_clean_up
!   25-aug-13/MR: added deallocation of itype_name.
!
      if (allocated(fname)) deallocate(fname)
      if (allocated(fname_keep)) deallocate(fname_keep)
      if (allocated(cname)) deallocate(cname)
      if (allocated(cform)) deallocate(cform)
      if (allocated(itype_name)) deallocate(itype_name)
!
    endsubroutine fnames_clean_up
!***********************************************************************
    subroutine vnames_clean_up
!
!  Deallocate space needed for reading the video.in file.
!
!   20-apr-10/Bourdin.KIS: copied from xyaverages_clean_up
!
      if (allocated(cnamev)) deallocate(cnamev,cformv)
!
    endsubroutine vnames_clean_up
!***********************************************************************
    subroutine xyaverages_clean_up
!
!  Deallocate the variables needed for xy-averages.
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
!  Deallocate the variables needed for xz-averages.
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
!  Deallocate the variables needed for yz-averages.
!
!   24-nov-09/anders: copied from yaverages_clean_up
!
      if (allocated(fnamex)) deallocate(fnamex)
      if (allocated(cnamex)) deallocate(cnamex)
      if (allocated(cformx)) deallocate(cformx)
!
    endsubroutine yzaverages_clean_up
!***********************************************************************
    subroutine phizaverages_clean_up
!
!  Deallocate the variables needed for phiz-averages.
!
!   24-nov-09/anders: copied from yaverages_clean_up
!
      if (allocated(fnamer)) deallocate(fnamer)
      if (allocated(cnamer)) deallocate(cnamer)
      if (allocated(cformr)) deallocate(cformr)
!
    endsubroutine phizaverages_clean_up
!***********************************************************************
    subroutine yaverages_clean_up
!
!  Deallocate the variables needed for y-averages.
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
!  Deallocate the variables needed for z-averages.
!
!   12-aug-09/dhruba: coded
!
      if (allocated(fnamexy)) deallocate(fnamexy)
      if (allocated(fnamexy_cap)) deallocate(fnamexy_cap)
      if (allocated(cnamexy)) deallocate(cnamexy)
      if (allocated(cformxy)) deallocate(cformxy)
!
    endsubroutine zaverages_clean_up
!*******************************************************************
    subroutine phiaverages_clean_up
!
!  Dellocate the variables needed for phi-averages.
!
!   24-nov-09/anders: copied from zaverages_clean_up
!
      if (allocated(fnamerz)) deallocate(fnamerz)
      if (allocated(cnamerz)) deallocate(cnamerz)
      if (allocated(cformrz)) deallocate(cformrz)
!
    endsubroutine phiaverages_clean_up
!*******************************************************************
    subroutine init_xaver
!
!  Initialize variables for x-averaging.
!
!  26-oct-09/dhruba: coded
!
     integer :: lx
!
     if (xav_max<x(l1)) loutside_avg=.true.
!
     ixav_max=l1
!
     do lx=l1,l2
       if (x(lx)<xav_max) ixav_max=lx
     enddo
!
   endsubroutine init_xaver
!*******************************************************************
   integer function name_is_present(ccname,vlabel,form)
!
!  Verify if the string vlabel is present or not in the ccname array,
!  return index, put format in form if requested and present
!
!  16-sep-10/dintrans: coded
!  16-may-12/MR: changed function type in integer to return index of found name;
!                new optional parameter form for use with print.in
!                where the name vlabel comes with a format specification
!
      character (len=*), dimension(:) :: ccname
      character (len=*)               :: vlabel
      character (len=*), optional     :: form
!
      intent(inout) :: ccname,form
      intent(in) :: vlabel
!
      integer :: mname, i, ind
      logical :: lform
      character (len=intlen) :: str
!
      if (present(form)) then
        lform = form=='T'
        form  = ''
      else
        lform = .false.
      endif
!
      mname = size(ccname)
      name_is_present=0
!
      do i=1,mname
!
        str=ccname(i)
        if (lform) then
          ind = index(str,'(')
          if (ind>1) then
            form=str(ind:)
            str=str(1:ind-1)
          endif
        endif
        if (str == vlabel) then
          name_is_present=i
          return
        endif
      enddo
!
    endfunction name_is_present
!***********************************************************************
endmodule Diagnostics
