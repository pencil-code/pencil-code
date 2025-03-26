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
  use General, only: safe_sum, loptest
!
  implicit none
!
  public :: initialize_diagnostics, initialize_diagnostic_arrays, prints
  public :: prep_finalize_thread_diagnos
  public :: diagnostic, initialize_time_integrals, get_average_density
  public :: xyaverages_z, xzaverages_y, yzaverages_x
  public :: diagnostics_init_reduc_pointers
  public :: diagnostics_diag_reductions
  !public :: diagnostics_read_diag_accum, diagnostics_write_diag_accum
  !public :: diagnostics_init_private_accumulators
  public :: phizaverages_r, yaverages_xz, zaverages_xy
  public :: phiaverages_rz
  public :: write_1daverages, write_2daverages
  public :: write_sound
  public :: write_1daverages_prepare, write_2daverages_prepare
  public :: name_is_present
  public :: expand_cname, parse_name, fparse_name
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
  public :: yintegrate_mn_name_xz
  public :: zsum_mn_name_xy_mpar_scal, zsum_mn_name_xy_mpar, &
            zsum_mn_name_xy_arr, zsum_mn_name_xy_arr2
  public :: allocate_fnames,allocate_vnames,allocate_sound
  public :: allocate_xyaverages, allocate_xzaverages, allocate_yzaverages
  public :: allocate_phizaverages
  public :: allocate_yaverages, allocate_zaverages, allocate_phiaverages, allocate_zaverages_data
  public :: trim_averages
  public :: fnames_clean_up, vnames_clean_up, diagnostics_clean_up
  public :: get_from_fname
  public :: gen_form_legend
  public :: sign_masked_xyaver
  public :: report_undefined_diagnostics
  public :: allocate_diagnostic_names
  public :: allocate_diagnostic_arrays
  public :: calc_nnames
  public :: save_diagnostic_controls
  public :: restore_diagnostic_controls
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
  real, target, dimension (nrcyl) :: phiavg_norm
  public :: phiavg_norm
  !$omp threadprivate(phiavg_norm)

  private
!
  real, pointer, dimension(:) :: p_phiavg_norm
  real, dimension (nrcyl,nx) :: phiavg_profile=0.0
  real :: dVol_rel1, dA_xy_rel1, dA_yz_rel1, dA_xz_rel1, dL_y_rel1

  character (len=intlen) :: ch1davg, ch2davg
  integer :: ixav_max
!
! Variables for Yin-Yang grid: z-averages.
!
  real, dimension(:,:,:), allocatable :: fnamexy_cap
!
! TP: variables for communication between threads
!
  real ::  dt_save,eps_rkf_save
  integer :: it_save

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
!  22-feb-2024/Kishore: added dA_{xy,yz,xz}_rel1, required for 1D averages in curvilinear coordinates.
!
      real :: dxeff,dyeff,dzeff
      real :: intdr_rel, intdtheta_rel, intdphi_rel, intdz_rel, intrdr_sph
      integer :: i
!
!  Since many of the averaging routines used for diagnostics don't account
!  for nonequidistant coordinates, warn the user.
      if (lroot) then
        if (.not.lproper_averages) then
          if (any(.not.lequidist)) call warning('initialize_diagnostics', &
            'volume averages are calculated wrongly for nonequidistant grids unless lproper_averages=T.')
!
          if (lwrite_xyaverages.and..not.(lequidist(1).and.lequidist(2))) call warning('initialize_diagnostics', &
            '1D average (xyaver) is calculated wrongly for non-equidistant grids unless lproper_averages=T.')
          if (lwrite_xzaverages.and..not.(lequidist(1).and.lequidist(3))) call warning('initialize_diagnostics', &
            '1D average (xzaver) is calculated wrongly for non-equidistant grids unless lproper_averages=T.')
          if (lwrite_yzaverages.and..not.(lequidist(2).and.lequidist(3))) call warning('initialize_diagnostics', &
            '1D average (yzaver) is calculated wrongly for non-equidistant grids unless lproper_averages=T.')
          if (lwrite_phizaverages.and..not.(lequidist(1).and.lequidist(2))) call warning('initialize_diagnostics', &
            '1D average (phizaver) is calculated wrongly for non-equidistant grids unless lproper_averages=T.')
!
          if (lwrite_yaverages.and..not.lequidist(2)) call warning('initialize_diagnostics', &
            '2D average (yaver) is calculated wrongly for non-equidistant grids unless lproper_averages=T.')
        endif
!
        if (lwrite_zaverages.and..not.lequidist(3)) call warning('initialize_diagnostics', &
          '2D average (zaver) is calculated wrongly for non-equidistant grids.')
        if (lwrite_phiaverages.and.any(.not.lequidist)) call warning('initialize_diagnostics', &
          '2D average (phiaver) is calculated wrongly for non-equidistant grids.')
      endif
!
!  Initialize rcyl for the phi-averages grid. Does not need to be
!  done after each reload of run.in, but this is the easiest way
!  of doing it.
!
      if (nrcyl/=0) then
        if ((lcylinder_in_a_box .or. lsphere_in_a_box) .and. &
            .not.(xyz0(1)==-xyz1(1) .and. xyz0(2)==-xyz1(2) .and. Lxyz(1)==Lxyz(2))) then
          call warning("initialize_diagnostics","box not centered at x,y=0 or box not isotropic in x,y,"// &
                       achar(10)//"although this is assumed for [cylinder|star]-in-a-box."// &
                       achar(10)//"We set rcyl range to maximum possible")
          drcyl=minval((/xyz1(1),-xyz0(1),xyz0(2),-xyz1(2)/))/nrcyl
          if (drcyl<=0.) call fatal_error("initialize_diagnostics","drcyl<0 not meaningful")
        else
          drcyl=xyz1(1)/nrcyl
        endif
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
        dA_xy_rel1 = 1./Area_xy
        dA_yz_rel1 = 1./Area_yz
        dA_xz_rel1 = 1./Area_xz
        dL_y_rel1 = 1./Ly
      elseif (lspherical_coords) then
!
!  Prevent zeros from less than 3-dimensional runs
!  (maybe this should be 2pi, but maybe not).
!
        if (nxgrid/=1) then
          intdr_rel = (xyz1(1)**3-xyz0(1)**3)/(3.*dx)
          intrdr_sph = (xyz1(1)**2-xyz0(1)**2)/(2.*dx)
        else
          intdr_rel = 1.
          intrdr_sph = 1.
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
        dA_xy_rel1 = 1./(intrdr_sph*nygrid)
        dA_yz_rel1 = 1./(intdtheta_rel*intdphi_rel)
        dA_xz_rel1 = 1./(intrdr_sph*intdphi_rel)
        dL_y_rel1 = 1./nygrid
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
        dA_xy_rel1 = 1./(intdr_rel*intdphi_rel)
        dA_yz_rel1 = 1./(intdphi_rel*intdz_rel)
        dA_xz_rel1 = 1./(nxgrid*intdz_rel)
        dL_y_rel1 = 1./nygrid
!
      else
        dVol_rel1=1./nwgrid
        dA_xy_rel1 = 1./nxygrid
        dA_xz_rel1 = 1./nxzgrid
        dA_yz_rel1 = 1./nyzgrid
        dL_y_rel1 = 1./nygrid
      endif
!
      if (lroot.and.ip<=10) then
        print*,'dVol_rel1=',dVol_rel1
        print*,'dA_xy_rel1=',dA_xy_rel1
        print*,'dA_xz_rel1=',dA_xz_rel1
        print*,'dA_yz_rel1=',dA_yz_rel1
        print*,'dL_y_rel1=',dL_y_rel1
      endif
!
!  Limits to xaveraging.
!
      call init_xaver
!
    endsubroutine initialize_diagnostics
!***********************************************************************
    subroutine initialize_diagnostic_arrays
!
!  Not needed so far.
!
      if (ldiagnos.and.allocated(fname)) fname=0.
      if (l1davgfirst) then
        if (allocated(fnamex)) fnamex=0.
        if (allocated(fnamey)) fnamey=0.
        if (allocated(fnamez)) fnamez=0.
      endif
      if (l1dphiavg.and.allocated(fnamer)) fnamer=0.
      if (l2davgfirst) then
        if (allocated(fnamexy)) fnamexy=0.
        if (allocated(fnamexz)) fnamexz=0.
        if (allocated(fnamerz)) fnamerz=0.
      endif

    endsubroutine initialize_diagnostic_arrays
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
      use Syscalls, only: system_cmd
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
        call save_name(dtdiagnos,idiag_dt)
        call save_name(eps_rkf_diagnos,idiag_eps_rkf)
        call save_name(one_real*(itdiagnos-1),idiag_it)
!
!  Whenever itype_name=ilabel_max_dt, scale result by dt (for printing Courant
!  time).
!
        do iname=1,nname
          if (itype_name(iname)==ilabel_max_dt) fname(iname)=dtdiagnos*fname(iname)
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
        if (lupdate_cvs) call system_cmd("cvs ci -m 'automatic update' >& /dev/null &")
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

      sedstring=''   !H;1h;$!d;x; '
      text='WARNING:'; ind=-1
      do i=1,len
        if (cname(i)/=''.and.cform(i)=='') then
          ind=index(cname(i),'('); if (ind==0) ind=len_trim(cname(i))+1
          text=trim(text)//' '//trim(cname(i)(1:ind-1))//','
          sedstring=trim(sedstring)//" -e '0,/^ *"//trim(cname(i)(1:ind-1))//"/ s/^ *\("//trim(cname(i)(1:ind-1))//"\)/#\1/'"
        endif
      enddo

      if (ind/=-1) then
!
!  If there are undefined diagnostics, comment them out in *.in.
!
        text=text(1:index(text,',',BACK=.true.)-1)//' diagnostic(s) in '//trim(file)// &
             ' undefined or multiply defined!'
        print*, trim(text)
        call system_cmd('sed '//trim(sedstring)//' -i '//trim(file))
!
!  In the case of duplicate entries, the following would keep the first and comment out the other,
!  but not all unixoid systems have tac.
!
        !call system_cmd('tac '//trim(file)//'|sed '//trim(sedstring)//'| tac > '//trim(file)//'.tmp')
        !call system_cmd('mv '//trim(file)//'.tmp '//trim(file))
      endif

    endsubroutine report_undefined
!***********************************************************************
    subroutine report_undefined_diagnostics

      if (lroot) then
        call report_undefined(cname,cform,nname,'print.in')
        if (allocated(cname_sound)) call report_undefined(cname_sound,cform_sound,nname_sound,'sound.in')
        if (allocated(cnamex))  call report_undefined(cnamex,cformx,nnamex,'yzaver.in')
        if (allocated(cnamey))  call report_undefined(cnamey,cformy,nnamey,'xzaver.in')
        if (allocated(cnamez))  call report_undefined(cnamez,cformz,nnamez,'xyaver.in')
        if (allocated(cnamer))  call report_undefined(cnamer,cformr,nnamer,'phizaver.in')
        if (allocated(cnamexy)) call report_undefined(cnamexy,cformxy,nnamexy,'zaver.in')
        if (allocated(cnamexz)) call report_undefined(cnamexz,cformxz,nnamexz,'yaver.in')
        if (allocated(cnamerz)) call report_undefined(cnamerz,cformrz,nnamerz,'phiaver.in')
        if (allocated(cnamev))  call report_undefined(cnamev,cformv,nnamev,'video.in')
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
      real                 :: rlength,dlength
      integer              :: dlength2
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
                  dlength=floor(alog10(abs(fname(iname))))+2 - rlength
                  if (dlength > 0) then
                    index_d=index(cform(iname),'.')
                    length=len(trim(cform(iname)))
                    dlength2=floor(alog10(rlength+dlength))+1-(length-index_i-2)
                    if (dlength2>0) then
                      length=length+dlength2
                      index_d=min(length,index_d+dlength2)
                    endif
                    write(cform(iname)(index_i+1:), &
                      '(f'//trim(itoa(length-index_i))//'.'//trim(itoa(length-index_d))//')') &
                      rlength+dlength
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
      use General, only: itoa, safe_character_assign
!
      integer,                 intent(in)   :: nlname
      real, dimension(nlname), intent(inout):: vname
      logical,optional,        intent(in)   :: lcomplex
!
      real, dimension (nlname) :: fmax_tmp, fsum_tmp, fmax, fsum, fweight_tmp
      real :: vol
      integer :: iname, imax_count, isum_count, nmax_count, nsum_count, itype, maxreq
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
      fweight_tmp=0.0
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
      call mpireduce_max(fmax_tmp,fmax,nmax_count)
      call mpireduce_sum(fsum_tmp,fsum,nsum_count)  !,nonblock=maxreq)        ! wrong for Yin-Yang due to overlap
      if (lweight_comm) call mpireduce_sum(fweight_tmp,fweight,nsum_count)!   ~
      !call mpiwait(maxreq)
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
    subroutine xyaverages_z(fnamez,ncountsz)
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
!MR: but wastes storage in fnamez
!
!   6-jun-02/axel: coded
!
      real, dimension(nz,nprocz,nnamez) :: fsumz
      integer, dimension(nz) :: nsum, ncount
      integer, dimension(:,:) :: ncountsz
      real, dimension(:,:,:) :: fnamez
      integer :: idiag
!
      if (nnamez>0) then
!
!  Find total number of contributing points in xy-plane for masked averages.
!
        if (any(ncountsz(1,:)>=0)) then
!
!  If there are any sign-masked averages amongst the diagnostics of fnamez.
!
          do idiag=1,nnamez

            ncount=ncountsz(:,idiag)
            if (ncount(1)>=0) then
!
!  If diagnostic with index diag is sign-masked average over xy-planes.
!
              call mpiallreduce_sum_int(ncount,nsum,nz,IXYPLANE)
!
!  Form average by dividing by nsum. Dividing by dA_xy_rel1
!  necessary as below the average is multiplied by that.
!
              where (nsum>0) fnamez(:,ipz+1,idiag)=fnamez(:,ipz+1,idiag)/(dA_xy_rel1*nsum)
              ncountsz(:,idiag)=0

            endif
          enddo
        endif
!
!  Communicate over all processors.
!  The result is only present on the root processor
!
        call mpireduce_sum(fnamez,fsumz,(/nz,nprocz,nnamez/))
        if (lroot) fnamez(:,:,1:nnamez)=fsumz(:,:,1:nnamez)*dA_xy_rel1
      endif
!
    endsubroutine xyaverages_z
!***********************************************************************
    subroutine xzaverages_y(fnamey)
!
!  Calculate xz-averages (still depending on y).
!
!  12-oct-05/anders: adapted from xyaverages_z
!
      real, dimension (ny,nprocy,nnamey) :: fsumy
      real, dimension (:,:,:) :: fnamey
!
!  Communicate over all processors.
!  The result is only present on the root processor.
!
      if (nnamey>0) then
        call mpireduce_sum(fnamey,fsumy,(/ny,nprocy,nnamey/))
        if (lroot) then
          fnamey(:,:,1:nnamey)=fsumy(:,:,1:nnamey)*dA_xz_rel1
        endif
      endif
!
    endsubroutine xzaverages_y
!***********************************************************************
    subroutine yzaverages_x(fnamex)
!
!  Calculate yz-averages (still depending on x).
!
!   2-oct-05/anders: adapted from xyaverages_z
!
      real, dimension (nx,nprocx,nnamex) :: fsumx
      real, dimension (:,:,:) :: fnamex
!
!  Communicate over all processors.
!  The result is only present on the root processor.
!
      if (nnamex>0) then
        call mpireduce_sum(fnamex,fsumx,(/nx,nprocx,nnamex/))
        if (lroot) fnamex(:,:,1:nnamex)=fsumx(:,:,1:nnamex)*dA_yz_rel1
      endif
!
    endsubroutine yzaverages_x
!***********************************************************************
    subroutine phizaverages_r(fnamer)
!
!  Calculate phiz-averages (still depending on r).
!
!  29-jan-07/wlad: adapted from yzaverages_x and phiaverages_rz
!
      real, dimension (nrcyl,nnamer) :: fsumr
      integer :: iname
      real, dimension (nrcyl) :: norm
      real, dimension (:,:) :: fnamer
!
!  Communicate over all processors.
!  The result is only present on the root processor.
!  Finally, normalization is performed with phiavg_norm*nzgrid
!
      if (nnamer>0) then
        call mpireduce_sum(fnamer,fsumr,(/nrcyl,nnamer/))
        if (ipz==0) call mpireduce_sum(phiavg_norm,norm,nrcyl)
        if (lroot) then
          do iname=1,nnamer
            fnamer(:,iname)=fsumr(:,iname)/(norm*nzgrid)
          enddo
        endif
      endif
!
    endsubroutine phizaverages_r
!***********************************************************************
    subroutine yaverages_xz(fnamexz)
!
!  Calculate y-averages (still depending on x and z).
!
!   7-jun-05/axel: adapted from zaverages_xy
!
      real, dimension (nx,nz,nnamexz) :: fsumxz
      real, dimension (:,:,:) :: fnamexz
!
!  Communicate over all processors along y beams.
!  The result is only present on the y-root processors.
!
      if (nnamexz>0) then
        call mpireduce_sum(fnamexz,fsumxz,(/nx,nz,nnamexz/),idir=2)
        if (lfirst_proc_y) fnamexz(:,:,1:nnamexz)=fsumxz(:,:,1:nnamexz)*dL_y_rel1
      endif
!
    endsubroutine yaverages_xz
!***********************************************************************
    subroutine zaverages_xy(fnamexy)
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
      real, dimension(:,:,:) :: fnamexy
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
    subroutine phiaverages_rz(fnamerz)
!
!  Calculate azimuthal averages (as functions of r_cyl,z).
!  NOTE: these averages depend on (r and) z, so after summation they
!  are still distributed over nprocz CPUs; hence the dimensions of fsumrz
!  and fnamerz.
!
!  9-dec-02/wolf: coded
!
      integer :: i
      real, dimension(nrcyl,nz,nprocz,nnamerz) :: fsumrz
      real, dimension(nrcyl) :: norm
      real, dimension(:,:,:,:) :: fnamerz
!
!  Communicate over all processors.
!  The result is only present on the root processor
!  normalize by sum of unity which is accumulated in norm.
!
      if (nnamerz>0) then
        call mpireduce_sum(fnamerz,fsumrz,(/nrcyl,nz,nprocz,nnamerz/))
        if (ipz==0) call mpireduce_sum(phiavg_norm,norm,nrcyl,idir=12)  ! avoid double comm!
        if (lroot) then
          do i=1,nnamerz
            fnamerz(:,:,:,i)=fsumrz(:,:,:,i)/spread(spread(norm,2,nz),3,nprocz)
          enddo
!do i=1,nnamerz
!if (lroot) print*, 'fnamerz(:,:,:,i)=', i,maxval(abs(fnamerz(:,:,:,i)))
!enddo
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
      t1ddiagnos = t
      ltimer = ip <= 12 .and. lroot
!
      if (nnamez > 0) then
        if (ltimer) taver = mpiwtime()
        call output_average(datadir, 'xy', nnamez, cnamez, fnamez, nzgrid, t1ddiagnos, lwrite_avg1d_binary, lroot)
        if (ltimer) print *, 'write_1daverages: write xy in ', mpiwtime() - taver, ' seconds'
      endif
!
      if (nnamey > 0) then
        if (ltimer) taver = mpiwtime()
        call output_average(datadir, 'xz', nnamey, cnamey, fnamey, nygrid, t1ddiagnos, lwrite_avg1d_binary, lroot)
        if (ltimer) print *, 'write_1daverages: write xz in ', mpiwtime() - taver, ' seconds'
      endif
!
      if (nnamex > 0) then
        if (ltimer) taver = mpiwtime()
        call output_average(datadir, 'yz', nnamex, cnamex, fnamex, nxgrid, t1ddiagnos, lwrite_avg1d_binary, lroot)
        if (ltimer) print *, 'write_1daverages: write yz in ', mpiwtime() - taver, ' seconds'
      endif
!
      if (nnamer > 0) then
        if (ltimer) taver = mpiwtime()
        if (lfirst_call) then
          call output_average(datadir, 'phi_z', nnamer, cnamer, fnamer, t1ddiagnos, .false., lroot, rcyl)
          lfirst_call = .false.
        else
          call output_average(datadir, 'phi_z', nnamer, cnamer, fnamer, t1ddiagnos, .false., lroot)
        endif
        if (ltimer) print *, 'write_1daverages: write phi_z in ', mpiwtime() - taver, ' seconds'
      endif
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
      t2davgfirst = t
      ltimer = ip <= 12 .and. lroot
!
      if (lwrite_yaverages) then
        if (ltimer) taver = mpiwtime()
        call output_average_2D('y', nnamexz, cnamexz, fnamexz, t2davgfirst, lfirst_proc_y)
        if (ltimer) print *, 'write_2daverages: write y averages in ', mpiwtime() - taver, ' seconds'
      endif
!
      if (lwrite_zaverages .and. (.not. lyang .or. lcaproot)) then
        if (ltimer) taver = mpiwtime()
        if (lcaproot) then
          ! cap root (Yang)
          call output_average_2D('z', nnamexy, cnamexy, fnamexy_cap, t2davgfirst, lfirst_proc_z)
        else
          ! z-beam root (Yin)
          call output_average_2D('z', nnamexy, cnamexy, fnamexy, t2davgfirst, lfirst_proc_z)
        endif
        if (ltimer) print *, 'write_2daverages: write z averages in ', mpiwtime() - taver, ' seconds'
      endif
!
      if (lwrite_phiaverages) then
        if (ltimer) taver = mpiwtime()
        ! normalization is already done in phiaverages_rz
        call output_average(datadir, ch2davg, nrcyl, nnamerz, cnamerz, fnamerz, t2davgfirst, rcyl, drcyl)

        if (ltimer) print *, 'write_2daverages: write phi averages in ', mpiwtime() - taver, ' seconds'
      endif
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
        call trim_average(datadir, 'phi', nrcyl*nzgrid, nnamerz)
        call trim_average(datadir, 'phi_z', nrcyl, nnamer)
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
    subroutine expand_cname_short(ccname,nname,vlabel,vname,ccform)
!
!  Expand string array cname with entries up to index nname such that,
!  if vlabel starts with <vname> or <vname><vname>, it is replaced by the three labels
!  <vname>x..., <vname>y..., <vname>z..., or correspondingly for curvilinear coord systems;
!  updates nname accordingly.
!
!  16-may-12/MR: coded
!
      character (len=*), dimension(:), allocatable :: ccname,ccform    ! F2003: , allocatable
      character (len=*) :: vlabel
      character (len=*) :: vname
      integer :: nname
!
      intent(inout) :: ccname,ccform,nname
      intent(in) :: vlabel, vname
!
      character (len=intlen):: tail
      integer :: ind,ic,lvn
!
      lvn = len(trim(vname))
      if (vlabel(1:lvn)/=trim(vname)) &
        call fatal_error('expand_cname_short','vlabel does not start with vname')
!
      ic = name_is_present(ccname,trim(vlabel))
      if (ic>0) then
!
        ind=lvn+1
        if ( lvn==1 .and. vlabel(2:2)==vname ) ind=3
!
        tail=vlabel(ind:len(vlabel))
        if ( lspherical_coords ) then
          call expand_cname_full(ccname,ccform,nname, &
                                 trim(vname)//'r'//trim(tail), &
                                 trim(vname)//'t'//trim(tail), &
                                 trim(vname)//'p'//trim(tail),ic_=ic)
        elseif ( lcylindrical_coords ) then
          call expand_cname_full(ccname,ccform,nname, &
                                 trim(vname)//'r'//trim(tail), &
                                 trim(vname)//'p'//trim(tail), &
                                 trim(vname)//'z'//trim(tail),ic_=ic)
        else
          call expand_cname_full(ccname,ccform,nname, &
                                 trim(vname)//'x'//trim(tail), &
                                 trim(vname)//'y'//trim(tail), &
                                 trim(vname)//'z'//trim(tail),ic_=ic)
        endif
      endif
!
    endsubroutine expand_cname_short
!***********************************************************************
    subroutine expand_cname_full(ccname,ccform,nname,xlabel,ylabel,zlabel,name,ic_)
!
!  Expand string array cname with entries up to index nname such that
!  label at position ic is replaced by the three labels xlabel, ylabel, zlabel;
!  update nname accordingly. Takes over format from cform(ic).
!
!   1-apr-04/wolf: coded
!  16-may-12/MR  : new parameter ic = position of label to be expanded in ccname
!
      use General, only : lextend_vector,ioptest
!
      character (len=*), dimension(:), allocatable :: ccname, ccform
      character (len=*) :: xlabel,ylabel
      character (len=*), optional :: name,zlabel
      integer :: nname
      integer, optional :: ic_
!
      intent(inout) :: ccname,ccform,nname
      intent(in) :: name,ic_,xlabel,ylabel,zlabel
!
      integer :: itot,ic,nname_prev

      if (present(name)) then
        ic = name_is_present(ccname,name)
      else
        ic = ioptest(ic_)
      endif
      if (ic<=0) return
!
      if (present(zlabel)) then
        itot=3
      else
        itot=2
      endif

      nname_prev=nname
      nname = nname+itot-1
!
      if (.not.lextend_vector(ccname,nname)) & ! sanity check
        call fatal_error('expand_cname_full','could not allocate ccname')
      if (.not.lextend_vector(ccform,nname)) & ! sanity check
        call fatal_error('expand_cname_full','could not allocate ccform')
!
      ccname(ic+itot:nname) = ccname(ic+1:nname_prev)
      ccname(ic:ic+1) = (/xlabel,ylabel/)
      if (present(zlabel)) ccname(ic+2) = zlabel

      ccform(ic+itot:nname) = ccform(ic+1:nname_prev)
      ccform(ic+1:ic+itot-1) = ccform(ic)
!
    endsubroutine expand_cname_full
!***********************************************************************
    subroutine set_type(iname, lsqrt, llog10, lint, lsum, lmax, lmin, lsurf, ldt)
!
!  Sets the diagnostic type in itype_name.
!
!  21-sep-17/MR: coded
!

      integer :: iname
      logical, optional :: lsqrt, llog10, lint, lsum, lmax, lmin, lsurf, ldt

      if (iname==0) return

      if (loptest(lsqrt)) then
        itype_name(iname)=ilabel_sum_sqrt
      elseif (loptest(llog10)) then
        itype_name(iname)=ilabel_sum_log10
      elseif (loptest(lint)) then
        itype_name(iname)=ilabel_integrate
      elseif (loptest(lsurf)) then
        itype_name(iname)=ilabel_surf
      elseif (loptest(lsum)) then
        itype_name(iname)=ilabel_sum
      elseif (loptest(lmax)) then
        itype_name(iname)=ilabel_max
      elseif (loptest(lmin)) then
        itype_name(iname)=ilabel_max_neg
      elseif (loptest(ldt)) then
        itype_name(iname)=ilabel_max_dt
      else
        itype_name(iname)=ilabel_save
      endif

    endsubroutine set_type
!***********************************************************************
    subroutine save_name(a,iname,l_dt)
!
!  Sets the value of a (must be treated as real) in fname array
!
!  26-may-02/axel: adapted from max_mn_name
!  20-sep-17/MR: removed setting of itype_name as it is set to "save" by default.
!
      real :: a
      integer :: iname
      logical, optional :: l_dt
!
      if (lroot.and.iname/=0) then
        fname(iname)=a
        if (loptest(l_dt)) itype_name(iname)=ilabel_max_dt
      endif
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
      if (iname/=0) fname_half(iname,:)=a(:)
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

      if (loptest(lneg)) then
        if (a<fname(iname)) fname(iname)=a
      else
        if (a>fname(iname)) fname(iname)=a
      endif
!!      fname(iname) = max(fname(iname),real(a))
!
!  Set corresponding entry in itype_name.
!
      if (loptest(lneg)) then
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

!!!      if (loptest(lneg)) then
!!!        if (a<fname(iname)) fname(iname)=a
!!!      else
!!!        if (a>fname(iname)) fname(iname)=a
!!!      endif
      fname(iname) = max(fname(iname),a)
!
!  Set corresponding entry in itype_name.
!
      if (loptest(lneg)) then
        itype_name(iname)=ilabel_max_neg
      elseif (loptest(l_dt)) then
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
      if (loptest(lsqrt)) then
        itype_name(iname)=ilabel_max_sqrt
      elseif (loptest(l_dt)) then
        itype_name(iname)=ilabel_max_dt
      elseif (loptest(lneg)) then
        itype_name(iname)=ilabel_max_neg
      elseif (loptest(lreciprocal)) then
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
        if (itype_name(iname) < ilabel_complex ) itype_name(iname)=itype_name(iname)+ilabel_complex

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
      real, dimension(:), intent(IN) :: a
      integer,            intent(IN) :: iname
      integer, optional,  intent(IN) :: ipart
      logical, optional,  intent(IN) :: lsqrt, llog10, lint, lplain

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

      real, dimension(:) :: a
      real, dimension(nname) :: fname
!
      real :: ppart,qpart
      integer :: iname
      integer, optional :: ipart
      logical, optional :: lsqrt, llog10, lint, lplain
!
      intent(in) :: iname

      real, dimension(size(a)) :: a_scaled
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
        if (loptest(lsqrt)) then
          itype_name(iname)=ilabel_sum_sqrt
        elseif (loptest(llog10)) then
          itype_name(iname)=ilabel_sum_log10
        elseif (loptest(lint)) then
          itype_name(iname)=ilabel_integrate
        elseif (loptest(lplain)) then
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
            if (loptest(lint)) then
              a_scaled=a*xprim(l1:l2)*yprim(m)*zprim(n)
            else
              a_scaled=a
            endif
!
!  Add up contributions, taking coordinate system into acount (fluid).
!  Initialize if one is on the first point, or add up otherwise.
!
            if (lfirstpoint) then
              if (lcartesian_coords.or.lpipe_coords.or.loptest(lplain)) then
                fname(iname)=sum(a_scaled)
              elseif (lspherical_coords) then
                fname(iname)=sum(r2_weight*a_scaled)*sinth_weight(m)
              elseif (lcylindrical_coords) then
                fname(iname)=sum(rcyl_weight*a_scaled)
              else
                call not_implemented('sum_mn_name_real','coordinate system')
              endif
            else
              if (lcartesian_coords.or.lpipe_coords.or.loptest(lplain)) then
                fname(iname)=fname(iname)+sum(a_scaled)
              elseif (lspherical_coords) then
                fname(iname)=fname(iname)+sum(r2_weight*a_scaled)*sinth_weight(m)
              elseif (lcylindrical_coords) then
                fname(iname)=fname(iname)+sum(rcyl_weight*a_scaled)
              else
                call not_implemented('sum_mn_name_real','coordinate system')
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
      integer :: iname
!
      real :: sum_name
!
      if (iname /= 0) then
!
        if (lfirstpoint) fname_half(iname,:)=0
!
        if (lspherical_coords) then
          sum_name=sum(r2_weight*sinth_weight(m)*a)
        elseif (lcylindrical_coords) then
          sum_name=sum(rcyl_weight*a)
        else
          sum_name=sum(a)
        endif
!
        if (y(m)>=yequator)then
          fname_half(iname,2)=fname_half(iname,2)+sum_name
        else
          fname_half(iname,1)=fname_half(iname,1)+sum_name
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
      integer :: iname
!
      real :: sum_name
!
      if (iname /= 0) then
!
!  Note: north means 1 and south means 2.
!  However, north corresponds to z < zequator,
!  in analogy to theta < !pi/2 for north.
!
        if (lfirstpoint) fname_half(iname,:)=0.

        if (lspherical_coords) then
          sum_name=sum(r2_weight*sinth_weight(m)*a)
        elseif (lcylindrical_coords) then
          sum_name=sum(rcyl_weight*a)
        else
          sum_name=sum(a)
        endif
!
!  North means 1 and south means 2.
!  However, north corresponds to z < zequator,
!  in analogy to theta < !pi/2 for north.
!
        if (z(n)>=zequator) then
          fname_half(iname,2)=fname_half(iname,2)+sum_name
        else
          fname_half(iname,1)=fname_half(iname,1)+sum_name
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
        if (loptest(lsqrt)) then
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
      ! logical, save :: lfirsttime=.true.
!
      if (iname /= 0) then
!
        if (lcylinder_in_a_box) then
          rlim=p%rcyl_mn
        elseif (lsphere_in_a_box) then
          rlim=p%r_mn
        ! elseif (lfirsttime) then
          call warning("sum_lim_mn_name","no reason to call it when "// &
               "not using a cylinder or"//achar(10)//"a sphere embedded in a Cartesian grid")
          ! lfirsttime=.false.
        endif
!
        dv=1.
        if (nxgrid/=1) dv=dv*dx
        if (nygrid/=1) dv=dv*dy
        if (nzgrid/=1) dv=dv*dz
!
        where ((rlim <= r_ext).and.(rlim >= r_int))
          aux = a
        elsewhere
          aux = 0.
        endwhere
!
        if (lfirstpoint) then
          if (lspherical_coords)then
            fname(iname)=sinth(m)*sum(x(l1:l2)*x(l1:l2)*aux)*dv
          else
            fname(iname)=sum(aux)*dv
          endif
        else
          if (lspherical_coords)then
            fname(iname)=fname(iname)+sinth(m)*sum(x(l1:l2)*x(l1:l2)*aux)*dv
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
    subroutine surf_mn_name(a,iname,ncontrib)
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
      integer, intent(in) :: iname,ncontrib
!
      if (iname>0) then
!
        if (lfirstpoint) fname(iname)=0.
        if (n==ncontrib) fname(iname)=fname(iname)+a
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
    subroutine xysum_mn_name_z(a,iname,mask)
!
!   3-sep-13/MR: derived from xysum_mn_name_z
!
      use Cdata, only: n
!
      real, dimension(nx), intent(IN) :: a
      integer,             intent(IN) :: iname
      logical, optional, dimension(nx), intent(IN) :: mask
!
      logical, dimension(nx) :: lmask
!
        if (present(mask)) then
          lmask=mask
        else
          lmask=.true.
        endif
!
      if (lproper_averages) then
        if (.not.all(lmask)) call fatal_error('xysum_mn_name_z', &
          'masking not implemented with lproper_averages=T')
        call xyintegrate_mn_name_z(a,iname)
      else
        call xysum_mn_name_z_npar(a,n,iname,MASK=lmask)
      endif
!
    endsubroutine xysum_mn_name_z
!***********************************************************************
    subroutine xysum_mn_name_z_npar(a,n,iname,mask)
!
!  Successively calculate sum over x,y of a, which is supplied at each call.
!  The result fnamez is z-dependent.
!  Start from zero if lfirstpoint=.true.
!
!   5-jun-02/axel: adapted from sum_mn_name
!   3-sep-13/MR: derived from xysum_mn_name_z, index n now parameter
!
      real, dimension(nx), intent(IN) :: a
      logical, optional, dimension(nx), intent(IN) :: mask
      integer,             intent(IN) :: n,iname
!
      integer :: nl
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
        if (.not.loutside_avg) then
!
!  n starts with nghost+1, so the correct index is n-nghost
!
          nl=n-nghost
          if (lspherical_coords.or.lcylindrical_coords)then
            fnamez(nl,ipz+1,iname)=fnamez(nl,ipz+1,iname)+ &
                                   sum(x(l1:ixav_max)*a(:ixav_max-nghost))
          else
            fnamez(nl,ipz+1,iname)=fnamez(nl,ipz+1,iname)+sum(a(:ixav_max-nghost),MASK=mask)
          endif
        endif
!
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

      if (lproper_averages) then
        call xzintegrate_mn_name_y(a,iname)
      else
        call xzsum_mn_name_y_mpar(a,m,iname)
      endif

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
      integer :: ml
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
        if (.not.loutside_avg) then
!
!  m starts with mghost+1, so the correct index is m-nghost.
!
          ml=m-nghost
          if (lspherical_coords.and.nxgrid>1)then
            fnamey(ml,ipy+1,iname)=fnamey(ml,ipy+1,iname)+ &
                                   sinth(m)*sum(x(l1:ixav_max)*a(:ixav_max-nghost))
          else ! also correct for cylindrical
            fnamey(ml,ipy+1,iname)=fnamey(ml,ipy+1,iname)+sum(a(:ixav_max-nghost))
          endif
        endif
!
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

      if (lproper_averages) then
        call yzintegrate_mn_name_x(a,iname)
      else
        call yzsum_mn_name_x_mpar(a,m,iname)
      endif

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
!   22-feb-2024/Kishore: fix for non-Cartesian coordinates
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
!  Note that in spherical/cylindrical coordinates, the line elements for
!  the y and z coordinates depend on x. However, since we are not summing
!  over x, that would simply cancel out when the final normalization is
!  performed in yzaverages_x, and so we simply drop the factors of x in
!  both the places.
!
        if (lspherical_coords) then
          fnamex(:,ipx+1,iname)=fnamex(:,ipx+1,iname)+sinth(m)*a
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
      real, dimension (nx) :: a, tmp
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
      fac=1.
!
      if (.not.lperi(2)) then
        if ((m==m1.and.lfirst_proc_y).or.(m==m2.and.llast_proc_y)) fac = .5
      endif
!
      if (lproper_averages) then
        tmp = a*dAxy_x(l1:l2)*dAxy_y(m)
      else
        tmp  = a
      endif
!
      if (lperi(1)) then
        suma = fac*sum(tmp)
      else
        suma = fac*(sum(tmp(2:nx-1))+.5*(tmp(1)+tmp(nx)))
      endif
!
!  n starts with nghost=4, so the correct index is n-nghost.
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
      real, dimension (nx) :: a, tmp
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
      if (.not.lperi(3)) then
        if ((n==n1.and.lfirst_proc_z).or.(n==n2.and.llast_proc_z)) fac = .5
      endif
!
      if (lproper_averages) then
        tmp = a*dAxz_x(l1:l2)*dAxz_z(n)
      else
        tmp  = a
      endif
!
      if (lperi(1)) then
        suma = fac*sum(tmp)
      else
        suma = fac*(sum(tmp(2:nx-1))+.5*(tmp(1)+tmp(nx)))
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
      real, dimension (nx) :: a, tmp
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
      if (.not.lperi(2)) then
        if ((m==m1.and.lfirst_proc_y).or.(m==m2.and.llast_proc_y)) fac = .5*fac
      endif
!
      if (.not.lperi(3)) then
        if ((n==n1.and.lfirst_proc_z).or.(n==n2.and.llast_proc_z)) fac = .5*fac
      endif
!
!
      if (lproper_averages) then
        tmp = a*dAyz_y(m)*dAyz_z(n)
      else
        tmp  = a
      endif
!
      fnamex(:,ipx+1,iname) = fnamex(:,ipx+1,iname) + fac*tmp
!
    endsubroutine yzintegrate_mn_name_x
!***********************************************************************
    subroutine phizsum_mn_name_r(a,iname)
!
!  Successively calculate sum over phi,z of a, which is supplied at each call.
!  Start from zero if lfirstpoint=.true.
!
!  29-jan-07/wlad: adapted from yzsum_mn_name_x and phisum_mn_name
!
      real, dimension (nx) :: a
      integer :: iname,ir
!
      if (iname==0) return
!
      if (lfirstpoint) fnamer(:,iname)=0.
!
      do ir=1,nrcyl
        fnamer(ir,iname) = fnamer(ir,iname) + sum(a*phiavg_profile(ir,:))
      enddo
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

      if (lproper_averages) then
        call yintegrate_mn_name_xz(a,iname)
      else
        call ysum_mn_name_xz_npar(a,n,iname)
      endif

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
!   22-feb-2024/Kishore: fix for non-Cartesian coordinates
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
!  Note that in spherical/cylindrical coordinates, the line element for
!  the y coordinate is actually x*dy. However, since we are not summing
!  over x, that would simply cancel out when the final normalization is
!  performed in yaverages_xz, and so we simply drop the factor of x in
!  both the places.
      fnamexz(:,nl,iname) = fnamexz(:,nl,iname)+a
!
    endsubroutine ysum_mn_name_xz_npar
!***********************************************************************
    subroutine yintegrate_mn_name_xz(a,iname)
!
!   Integrate over y. Apply trapezoidal rule properly in the case
!   of non-periodic boundaries.
!
!   19-nov-2024/Kishore: adapted from xyintegrate_mn_name_z
!
      real, dimension (nx) :: a, tmp
      integer :: iname, nl
      real :: fac
!
      if (iname==0) return
!
!  Initialize to zero, including other parts of the z-array
!  which are later merged with an mpi reduce command.
!
      if (lfirstpoint) fnamexz(:,:,iname) = 0.0
!
      fac=1.
!
      if (.not.lperi(2)) then
        if ((m==m1.and.lfirst_proc_y).or.(m==m2.and.llast_proc_y)) fac = .5
      endif
!
      if (lproper_averages) then
        tmp = fac*a*yprim(m)
      else
        tmp = fac*a
      endif
!
!  n starts with nghost=4, so the correct index is n-nghost.
!
      nl=n-nghost
      fnamexz(:,nl,iname) = fnamexz(:,nl,iname) + tmp
!
    endsubroutine yintegrate_mn_name_xz
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
      intent(in) :: p
!
      real :: r0,width
      integer :: ir
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
      if (lfirstpoint) phiavg_norm=0.
!
!  Normalization factor, not depending on n, so only done for n=nn(1)
!  As we calculate z-averages, multiply by nzgrid when used.
!
      if (n==nn(1)) phiavg_norm=phiavg_norm+sum(phiavg_profile,2)
!
    endsubroutine calc_phiavg_profile
!***********************************************************************
    subroutine diagnostics_init_reduc_pointers
!
      p_phiavg_norm => phiavg_norm
!
    endsubroutine diagnostics_init_reduc_pointers
!***********************************************************************
    subroutine diagnostics_diag_reductions
!
      p_phiavg_norm = p_phiavg_norm + phiavg_norm
!
    endsubroutine diagnostics_diag_reductions
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
!  n starts with nghost+1, so the correct index is n-nghost
!
        n_nghost=n-nghost
        do ir=1,nrcyl
          fnamerz(ir,n_nghost,ipz+1,iname) &
               = fnamerz(ir,n_nghost,ipz+1,iname) + sum(a*phiavg_profile(ir,:))
        enddo
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
!   22-mar-25/TP: refactored name allocations to separate function
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
      integer, allocatable, dimension (:,:) :: sound_inds
      real, dimension(dimensionality) :: coords
      integer :: lsound, msound, nsound, nitems,mcoords_sound
!
!  Allocate and initialize to zero. Setting it to zero is only
!  necessary because of the pencil test, which doesn't compute these
!  averages, and only evaluates its output for special purposes
!  such as computing mean field energies in calc_bmz, for example,
!
      mcoords_sound = parallel_count_lines(sound_coord_file)

      allocate(sound_inds(mcoords_sound,3),stat=ierr)
      if (ierr>0) call fatal_error('allocate_sound','Could not allocate sound_inds')
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
      if (ierr>0) call fatal_error('allocate_sound','Could not allocate cname_sound')
!
      lwrite_sound = ncoords_sound>0
      if (lwrite_sound) then
!
        allocate(sound_coords_list(ncoords_sound,3),stat=ierr)
        if (ierr>0) call fatal_error('allocate_sound','Could not allocate sound_coords_list')
        sound_coords_list = sound_inds(1:ncoords_sound,:)
!
        if (ierr>0) call fatal_error('allocate_sound','Could not allocate fname_sound')
        if (ldebug) print*, 'allocate_sound: allocated memory for '// &
                            'fname_sound  with nname_sound  =', nnamel
!
!
        allocate(cform_sound(nnamel),stat=ierr)
        if (ierr>0) call fatal_error('allocate_sound','Could not allocate cform_sound')
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
    subroutine allocate_sound_data(nnamel)
!
!   22-mar-25/TP: coded
!
        integer, intent(in) :: nnamel
        integer :: ierr
        allocate(fname_sound(ncoords_sound,nnamel),stat=ierr)
        fname_sound = 0.0

    endsubroutine allocate_sound_data
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
!   25-mar-25/TP: refactored name allocations to their own function
!
      integer, intent(in) :: nnamel
!
      integer :: stat
!
      allocate(fname(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_fnames','Could not allocate fname')
      if (ldebug) print*, 'allocate_fnames    : allocated memory for '// &
                          'fname   with nname   =', nnamel
      allocate(fname_keep(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_fnames','Could not allocate fname_keep')
      fname=0.0
      fname_keep=0.0


    endsubroutine allocate_fnames
!***********************************************************************
    subroutine allocate_cnames(nnamel)
!
!   25-mar-25/TP: separated from allocate_fnames
!
      integer, intent(in) :: nnamel
!
      integer :: stat
      allocate(cname(nnamel),stat=stat)
      if (ldebug) print*, 'allocate_fnames    : allocated memory for '// &
                          'cname   with nname   =', nnamel
      if (ldebug) flush(6)
      if (stat>0) call fatal_error('allocate_fnames','Could not allocate cname')
      cname=''
!
      allocate(cform(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_fnames','Could not allocate cform')
      if (ldebug) print*, 'allocate_fnames    : allocated memory for '// &
                          'cform   with nname   =', nnamel
      cform=''
!
      allocate(itype_name(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_fnames','Could not allocate itype_name')
      if (ldebug) print*, 'allocate_fnames    : allocated memory for '// &
                          'itype_name with nname   =', nnamel
      itype_name=ilabel_save
    endsubroutine
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
      if (stat>0) call fatal_error('allocate_vnames','Could not allocate cnamev')
      if (ldebug) print*, 'allocate_vnames    : allocated memory for '// &
          'cnamev  with nnamev  =', nnamel
      cnamev=''
!
      if (allocated(cformv)) deallocate(cformv)
      allocate(cformv(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_vnames','Could not allocate cformv')
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
!   25-mar-25/TP: refactored name allocations to their own function
!
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
      allocate(fnamez(nz,nprocz,nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_xyaverages','Could not allocate fnamez')
      if (ldebug) print*, 'allocate_xyaverages: allocated memory for '// &
                          'fnamez  with nnamez  =', nnamel
      fnamez=0.0
      allocate(ncountsz(nz,nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_xyaverages','Could not allocate ncountsz')
      if (ldebug) print*, 'allocate_xyaverages: allocated memory for '// &
                          'ncountsz  with nnamez  =', nnamel
      ncountsz=-1
!
    endsubroutine allocate_xyaverages
!***********************************************************************
    subroutine allocate_xyaverages_names(nnamel)
!
!   25-mar-25/TP: carved from allocate_xyaverages
!
      integer, intent(in) :: nnamel
      integer :: stat

      allocate(cnamez(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_xyaverages','Could not allocate cnamez')
      if (ldebug) print*, 'allocate_xyaverages: allocated memory for '// &
                          'cnamez  with nnamez  =', nnamel
      cnamez=''
!
      allocate(cformz(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_xyaverages','Could not allocate cformz')
      if (ldebug) print*, 'allocate_xyaverages: allocated memory for '// &
                          'cformz  with nnamez  =', nnamel
      cformz=''
!
    endsubroutine allocate_xyaverages_names
!***********************************************************************
    subroutine allocate_xzaverages(nnamel)
!
!  Allocate arrays needed for xz-averages.
!
!   24-nov-09/anders: copied from allocate_yaverages
!   11-jan-11/MR: parameter nnamel added
!   21-mar-25/TP: refactored name allocations to their own function
!
      integer, intent(in) :: nnamel
!
      integer :: stat
!
!
      allocate(fnamey(ny,nprocy,nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_xzaverages','Could not allocate fnamey', .true.)
      if (ldebug) print*, 'allocate_xzaverages: allocated memory for '// &
                          'fnamey  with nnamey  =', nnamel
      fnamey=0.0
    endsubroutine allocate_xzaverages
!***********************************************************************
    subroutine allocate_xzaverages_names(nnamel)

!
!   21-mar-25/TP: carved from allocate_xzaverages
!
      integer, intent(in) :: nnamel
      integer :: stat

      allocate(cnamey(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_xzaverages','Could not allocate cnamey')
      if (ldebug) print*, 'allocate_xzaverages: allocated memory for '// &
                          'cnamey  with nnamey  =', nnamel
      cnamey=''
      allocate(cformy(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_xzaverages','Could not allocate cformy', .true.)
      if (ldebug) print*, 'allocate_xzaverages: allocated memory for '// &
                          'cformy  with nnamey  =', nnamel
      cformy=''
    endsubroutine allocate_xzaverages_names
!***********************************************************************
    subroutine allocate_yzaverages(nnamel)
!
!  Allocate arrays needed for yz-averages.
!
!   24-nov-09/anders: copied from allocate_yaverages
!   11-jan-11/MR: parameter nnamel added
!   21-mar-25/TP: refactored name allocations to their own function
!
      integer, intent(in) :: nnamel
!
      integer :: stat
!
      allocate(fnamex(nx,nprocx,nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_yzaverages','Could not allocate fnamex')
      if (ldebug) print*, 'allocate_yzaverages: allocated memory for '// &
                          'fnamex  with nnamex  =', nnamel
      fnamex=0.0
    endsubroutine allocate_yzaverages
!***********************************************************************
    subroutine allocate_yzaverages_names(nnamel)
!
!   21-mar-25/TP: carved from allocate_yzaverages
!

      integer, intent(in) :: nnamel
      integer :: stat

      allocate(cnamex(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_yzaverages','Could not allocate cnamex')
      if (ldebug) print*, 'allocate_yzaverages: allocated memory for '// &
                          'cnamex  with nnamex  =', nnamel
      cnamex=''
!
      allocate(cformx(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_yzaverages','Could not allocate cformx')
      if (ldebug) print*, 'allocate_yzaverages: allocated memory for '// &
                          'cformx  with nnamex  =', nnamel
      cformx=''

    endsubroutine allocate_yzaverages_names
!***********************************************************************
    subroutine allocate_phizaverages(nnamel)
!
!  Allocate arrays needed for phiz-averages.
!
!   24-nov-09/anders: copied from allocate_yaverages
!   11-jan-11/MR: parameter nnamel added
!   21-mar-25/TP: refactored name allocations to their own function
!
      integer, intent(in) :: nnamel
!
      integer :: stat
!
      allocate(fnamer(nrcyl,nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_phizaverages','Could not allocate fnamer')
      if (ldebug) print*, 'allocate_phizaverages: allocated memory for '// &
                          'fnamer  with nnamer'
      fnamer=0.0
    endsubroutine allocate_phizaverages
!***********************************************************************
    subroutine allocate_phizaverages_names(nnamel)
!
!   21-mar-25/TP: carved from allocate_phizaverages
!
      integer, intent(in) :: nnamel
      integer :: stat
      allocate(cnamer(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_phizaverages','Could not allocate cnamer')
      if (ldebug) print*, 'allocate_phizaverages: allocated memory for '// &
                          'cnamer  with nnamer =', nnamel
      cnamer=''
!
!
      allocate(cformr(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_phizaverages','Could not allocate cformr')
      if (ldebug) print*, 'allocate_phizaverages: allocated memory for '// &
                          'cformr  with nnamer =', nnamel
      cformr=''
!
    endsubroutine allocate_phizaverages_names
!***********************************************************************
    subroutine allocate_yaverages(nnamel)
!
!  Allocate arrays needed for y-averages.
!
!   12-aug-09/dhruba: coded
!   11-jan-11/MR: parameter nnamel added
!   21-mar-25/TP: refactored name allocations to their own function
!
      integer, intent(in) :: nnamel
!
      integer :: stat
!
      allocate(fnamexz(nx,nz,nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_yaverages','Could not allocate fnamexz')
      if (ldebug) print*, 'allocate_yaverages : allocated memory for '// &
                          'fnamexz with nnamexz =', nnamel
      fnamexz=0.0
    endsubroutine allocate_yaverages
!*******************************************************************
    subroutine allocate_yaverages_names(nnamel)
!
!   21-mar-25/TP: carved from allocate_phizaverages
!
      integer, intent(in) :: nnamel
      integer :: stat

      allocate(cnamexz(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_yaverages','Could not allocate cnamexz')
      if (ldebug) print*, 'allocate_yaverages : allocated memory for '// &
                          'cnamexz with nnamexz =', nnamel
      cnamexz=''
!
!
      allocate(cformxz(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_yaverages','Could not allocate cformxz')
      if (ldebug) print*, 'allocate_yaverages : allocated memory for '// &
                          'cformxz with nnamexz =', nnamel
      cformxz=''
!
    endsubroutine allocate_yaverages_names
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
      if (stat>0) call fatal_error('allocate_zaverages','Could not allocate cnamexy')
      if (ldebug) print*, 'allocate_zaverages : allocated memory for '// &
                          'cnamexy with nnamexy =', nnamel
      cnamexy=''
!
      allocate(cformxy(nnamel),stat=stat)
      if (stat>0) &
        call fatal_error('allocate_zaverages','Could not allocate cformxy')
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
      if (stat>0) call fatal_error('allocate_zaverages_data','Could not allocate fnamexy')
      if (ldebug) print*, 'allocate_zaverages_data: allocated memory for '// &
                          'fnamexy with nnamexy =', nnamel
      fnamexy=0.

      if (lcaproot) then
        allocate(fnamexy_cap(nnamel,nx,nycap),stat=stat)
        if (stat>0) call fatal_error('allocate_zaverages_data','Could not allocate fnamexy_cap')
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
!   21-mar-25/TP: refactored name allocations to their own function
!
      integer, intent(in) :: nnamel
!
      integer :: stat
!
      allocate(fnamerz(nrcyl,nz,nprocz,nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_phiaverages','Could not allocate fnamerz')
      if (ldebug) print*, 'allocate_phiaverages : allocated memory for '// &
                          'fnamerz with nnamerz =', nnamel
      fnamerz=0.0
    endsubroutine allocate_phiaverages
!***********************************************************************
    subroutine allocate_phiaverages_names(nnamel)
!
!   21-mar-25/TP: carved from allocate_phiaverages
!
      integer, intent(in) :: nnamel
      integer :: stat
      allocate(cnamerz(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_phiaverages','Could not allocate cnamerz')
      if (ldebug) print*, 'allocate_phiaverages : allocated memory for '// &
                          'cnamerz with nnamerz =', nnamel
      cnamerz=''
!
!
      allocate(cformrz(nnamel),stat=stat)
      if (stat>0) call fatal_error('allocate_phiaverages','Could not allocate cformrz')
      if (ldebug) print*, 'allocate_phiaverages : allocated memory for '// &
                          'cformrz with nnamerz =', nnamel
      cformrz=''
!
    endsubroutine allocate_phiaverages_names
!***********************************************************************
    subroutine sign_masked_xyaver(quan,idiag)
!
!  Forms sign-masked averaging over xy-planes (only positive values count).
!
!  28-sep-16/MR: coded
!
      real   ,dimension(nx),intent(IN) :: quan
      integer,              intent(IN) :: idiag

      real, dimension(nx) :: buf
      integer :: nl

      if (idiag>0) then

        where (quan>0.)
          buf=quan
        elsewhere
          buf=0.
        endwhere

        call xysum_mn_name_z(buf,idiag)

        nl=n-nghost
!
!  Number of points which contribute in each z-layer is held ncountsz.
!
        if (lfirstpoint) ncountsz(:,idiag) = 0
        ncountsz(nl,idiag) = ncountsz(nl,idiag)+count(buf>0.)

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
      if (allocated(ncountsz)) deallocate(ncountsz)
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

     if (lav_smallx) then
!
       loutside_avg = xav_max<x(l1)
       if (loutside_avg) return
!
       do lx=l1+1,l2
         if (x(lx)>xav_max) then
           ixav_max=lx-1
           return
         endif
       enddo
       ixav_max=l2

     else
       loutside_avg=.false.
       ixav_max=l2
     endif
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
    subroutine prep_finalize_thread_diagnos
!
!  For accumulated diagnostic variables get which reduction to perform:
!  inds_max_diags/inds_sum_diags hold the indices of the fname entries,
!  which need to be reduced by max/sum. Will be determined only once.
!  Index vectors are threadprivate.
!
!  25-aug-23/TP: modified, bug fix
!
      use General, only: allpos_in_array_int

      integer :: nmax, nsum, nmin, i
      logical :: firstcall=.true., firstcall_from_pencil_check=.false.
      integer, dimension(2) :: max_range, sum_range
      !$omp threadprivate(firstcall)

      max_range(1) = -5
      max_range(2) = -1

      sum_range(1) = 1
      sum_range(2) = 40
      
!  Have to do this ugly workaround since the pencil tests call this function
!  and we want the first non-pencil-test call.
!
      if (firstcall_from_pencil_check .and. .not. lpencil_check_at_work) then
        firstcall = .true.
        deallocate(inds_max_diags)
        deallocate(inds_sum_diags)
      endif

      if (.not.firstcall) return

      nmax = allpos_in_array_int(max_range,itype_name)
      if (nmax>0) then
        allocate(inds_max_diags(nmax))
        nmax = allpos_in_array_int(max_range,itype_name,inds_max_diags)
      endif
      nsum = allpos_in_array_int(sum_range,itype_name)
      if (nsum>0) then
        allocate(inds_sum_diags(nsum))
        nsum = allpos_in_array_int(sum_range,itype_name,inds_sum_diags)
      endif

      firstcall=.false.
      firstcall_from_pencil_check=lpencil_check_at_work

    endsubroutine prep_finalize_thread_diagnos
!***********************************************************************
    subroutine    calc_nnames()
!
!  Calculates the sizes of diagnostic arrays
!  21-mar-25/TP: carved out from rprint_list
!
      use General, only: numeric_precision
      use File_io,         only: parallel_file_exists, parallel_count_lines,read_name_format

      character (LEN=15)           :: print_in_file
      character (LEN=*), parameter :: video_in_file    = 'video.in'
      character (LEN=*), parameter :: sound_in_file    = 'sound.in'
      character (LEN=*), parameter :: xyaver_in_file   = 'xyaver.in'
      character (LEN=*), parameter :: xzaver_in_file   = 'xzaver.in'
      character (LEN=*), parameter :: yzaver_in_file   = 'yzaver.in'
      character (LEN=*), parameter :: phizaver_in_file = 'phizaver.in'
      character (LEN=*), parameter :: yaver_in_file    = 'yaver.in'
      character (LEN=*), parameter :: zaver_in_file    = 'zaver.in'
      character (LEN=*), parameter :: phiaver_in_file  = 'phiaver.in'

      print_in_file = 'print.in'
      if (numeric_precision() == 'D') then
        if (parallel_file_exists(trim(print_in_file)//'.double')) &
            print_in_file = trim(print_in_file)//'.double'
      endif

      nname = max(0,parallel_count_lines(print_in_file,ignore_comments=.true.))

      if ( dvid/=0.0 ) then
        nnamev = max(0,parallel_count_lines(video_in_file))
      endif

      if ( dimensionality>0 .and. dsound/=0.0 ) then
        nname_sound = max(0,parallel_count_lines(sound_in_file))
      endif

      nnamez = parallel_count_lines(xyaver_in_file)
      nnamey = parallel_count_lines(xzaver_in_file)
      nnamex = parallel_count_lines(yzaver_in_file)
      nnamer = max(0,parallel_count_lines(phizaver_in_file))
      nnamexz = parallel_count_lines(yaver_in_file)
      nnamexy = parallel_count_lines(zaver_in_file)
      nnamerz = parallel_count_lines(phiaver_in_file)


    endsubroutine calc_nnames 
!***********************************************************************
    subroutine allocate_diagnostic_names()
!
!  Allocates diagnostic arrays holding the names of the diagnostic outputs
!  Separate from the data allocations because of multithreading concerns
!  21-mar-25/TP: coded
!
      if (nname>0) call allocate_cnames(nname)
      if ( dvid/=0.0 ) then
        if (nnamev>0) call allocate_vnames(nnamev)
      endif
      if ( dimensionality>0 .and. dsound/=0.0 ) then
        if (nname_sound>0) call allocate_sound(nname_sound)
      endif

      if (nnamez>0) call allocate_xyaverages_names(nnamez)

      if (nnamey>0) call allocate_xzaverages_names(nnamey)

      if (nnamex>0) call allocate_yzaverages_names(nnamex)

      if (nnamer>0) then
        if (lcylinder_in_a_box.or.lsphere_in_a_box) then
          call allocate_phizaverages_names(nnamer)
        endif
      endif

      if (nnamexz>0) call allocate_yaverages_names(nnamexz)

!
      if (nnamexy>0) then
        call allocate_zaverages(nnamexy)
      endif

      if (nnamerz>0) then
        if (lcylinder_in_a_box.or.lsphere_in_a_box) then
          call allocate_phiaverages_names(nnamerz)
        endif
      endif
    endsubroutine allocate_diagnostic_names
!***********************************************************************
    subroutine allocate_diagnostic_arrays()
!
!  Allocates diagnostic arrays holding the output data
!  Separate from the name allocations because of multithreading concerns
!  21-mar-25/TP: coded
!

!
!  Read print.in.double if applicable, else print.in.
!  Read in the list of variables to be printed.
!
      if (nname>0) call allocate_fnames(nname)


      if ( dimensionality>0 .and. dsound/=0.0 ) then
        if (nname_sound>0) call allocate_sound_data(nname_sound)
      endif

      if (nnamez>0) call allocate_xyaverages(nnamez)

      if (nnamey>0) call allocate_xzaverages(nnamey)

      if (nnamex>0) call allocate_yzaverages(nnamex)

      if (nnamer>0) then
        if (lcylinder_in_a_box.or.lsphere_in_a_box) then
          call allocate_phizaverages(nnamer)
        endif
      endif

      if (nnamexz>0) call allocate_yaverages(nnamexz)

!
      if (nnamexy>0) then
        if (lwrite_zaverages) call allocate_zaverages_data(nnamexy)
      endif

      if (nnamerz>0) then
        if (lcylinder_in_a_box.or.lsphere_in_a_box) then
          call allocate_phiaverages(nnamerz)
        endif
      endif

    endsubroutine allocate_diagnostic_arrays
!***********************************************************************
   subroutine save_diagnostic_controls
!
!  Saves threadprivate variables to shared ones.
!
!  25-aug-23/TP: Coded
!  19-march-25/TP: moved from Equ to here
!
    l1davgfirst_save = l1davgfirst
    ldiagnos_save = ldiagnos
    l1dphiavg_save = l1dphiavg
    l2davgfirst_save = l2davgfirst

    lout_save = lout
    l1davg_save = l1davg
    l2davg_save = l2davg
    lout_sound_save = lout_sound
    lvideo_save = lvideo
!
!  Record times for diagnostic and 2d average output.
!
    if (l1davgfirst) t1ddiagnos_save=t ! (1-D averages are for THIS time)
    if (l2davgfirst) t2davgfirst_save=t ! (2-D averages are for THIS time)
    if (lvideo     ) tslice_save=t ! (slices are for THIS time)
    if (lout_sound ) tsound_save=t
    if (ldiagnos) then
      t_save  = t ! (diagnostics are for THIS time)
      dt_save = dt
      it_save = it
      eps_rkf_save = eps_rkf
    endif
    lpencil_save = lpencil

    endsubroutine save_diagnostic_controls
!***********************************************************************
    subroutine restore_diagnostic_controls
!
!   Restores the diagnostics flags that were saved when calculating diagnostics started.
!
!   13-nov-23/TP: Written
!   19-march-25/TP: moved from Equ to here
!
    l1davgfirst = l1davgfirst_save
    ldiagnos = ldiagnos_save
    l1dphiavg = l1dphiavg_save
    l2davgfirst = l2davgfirst_save

    lout = lout_save
    l1davg = l1davg_save
    l2davg = l2davg_save
    lout_sound = lout_sound_save
    lvideo = lvideo_save
    t1ddiagnos = t1ddiagnos_save
    t2davgfirst= t2davgfirst_save
    tslice = tslice_save
    tsound = tsound_save

    if (ldiagnos) then
      tdiagnos  = t_save
      dtdiagnos = dt_save
      itdiagnos = it_save
      eps_rkf_diagnos = eps_rkf_save
    endif

    tspec = tspec_save
    lpencil = lpencil_save

    endsubroutine restore_diagnostic_controls
!***********************************************************************
endmodule Diagnostics
