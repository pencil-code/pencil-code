! $Id$
!
!  This module produces slices for animation purposes.
!
module Slices
!
! 16-nov-11/MR: I/O error handling generally introduced
!
  use Cdata
  use Messages
  use Sub, only: xlocation, zlocation, update_snaptime, read_snaptime, position
!
  implicit none
!
  private
!
  public :: wvid, wvid_prepare, setup_slices, wslice
!
  real, target, dimension (nx,ny) :: slice_xy=0.0, slice_xy2=0.0
  real, target, dimension (nx,nz) :: slice_xz=0.0, slice_xz2=0.0
  real, target, dimension (ny,nz) :: slice_yz=0.0
  real, target, dimension (nx,ny) :: slice_xy3=0.0, slice_xy4=0.0
!
  real, public :: tvid
  integer, public :: nvid
  real :: tslice
!
  contains
!***********************************************************************
    subroutine wvid_prepare
!
!  Prepare lvideo for writing slices into video file
!  This is useful for visualization of scalar field (or one component
!  of a vector field) on the perifery of a box.
!  Can be visualized in idl using rvid_box.pro
!
!  20-oct-97/axel: coded
!  08-oct-02/tony: increased size of file to handle datadir//'/tvid.dat'
!  13-nov-02/axel: added more fields, use wslice.
!  18-mar-03/axel: added dust velocity
!
      logical, save :: lfirst_call=.true.
!
      character (len=fnlen) :: file
!
!  Output vid-data in 'tvid' time intervals
!
      file = trim(datadir)//'/tvid.dat'
      if (lfirst_call) then
        call read_snaptime(file,tvid,nvid,dvid,t)
        lfirst_call=.false.
      endif
!
!  This routine sets lvideo=T whenever its time to write a slice
!
      call update_snaptime(file,tvid,nvid,dvid,t,lvideo)
!
!  Save current time so that the time that is written out in wslice() is not
!  from the next time step
!
      if (lvideo) tslice = t
!
    endsubroutine wvid_prepare
!***********************************************************************
    subroutine wvid(f,path)
!
!  Write slices for animation of scalar field
!  (or one component of a vector field) on the perifery of a box.
!  Can be visualized in idl using rvid_box.pro.
!
!  13-nov-02/axel: added more fields, use wslice.
!  22-sep-07/axel: changed Xy to xy2, to be compatible with Mac
!
      use General,         only: itoa
      use Chemistry,       only: get_slices_chemistry
      use Chiral,          only: get_slices_chiral
      use Cosmicray,       only: get_slices_cosmicray
      use Density,         only: get_slices_density,get_slices_pressure
      use Dustdensity,     only: get_slices_dustdensity
      use Dustvelocity,    only: get_slices_dustvelocity
      use EquationOfState, only: get_slices_eos
      use Energy,          only: get_slices_energy
      use Heatflux,        only: get_slices_heatflux
      use Hydro,           only: get_slices_hydro
      use Interstellar,    only: get_slices_interstellar
      use Magnetic,        only: get_slices_magnetic
      use Particles_main,  only: get_slices_particles
      use Pscalar,         only: get_slices_pscalar
      use Radiation,       only: get_slices_radiation
      use Shock,           only: get_slices_shock
      use Special,         only: get_slices_special
      use Testfield,       only: get_slices_testfield
      use Testflow,        only: get_slices_testflow
      use Testscalar,      only: get_slices_testscalar
!
      real, dimension (mx,my,mz,mfarray) :: f
      character(len=*) :: path
!
      type (slice_data) :: slices
      character(len=intlen) :: sindex
      logical, save :: lfirstloop=.true.
      logical :: lnewfile=.true.
      logical :: lslices_legacy=.true.
      integer :: inamev, iostat
!
      slices%index=0
!
!  Loop over slices.
!
      inamev=1
      do while (inamev <= nnamev)
!
        slices%ready=.false.
        slices%xy=>slice_xy
        slices%xz=>slice_xz
        slices%xz2=>slice_xz2
        slices%yz=>slice_yz
        slices%xy2=>slice_xy2
        slices%xy3=>slice_xy3
        slices%xy4=>slice_xy4
!
        slices%name=trim(cnamev(inamev))
!
!  By default assume we're not using module hooks to get the slice contents
!
        lslices_legacy=.true.
!
!  Get slice information from the modules.
!
        lslices_legacy=.false.
        if (lchemistry)    call get_slices_chemistry   (f,slices)
        if (lchiral)       call get_slices_chiral      (f,slices)
        if (lcosmicray)    call get_slices_cosmicray   (f,slices)
        if (ldensity.or.lanelastic)  &
                           call get_slices_density     (f,slices)
        if (lanelastic)    call get_slices_pressure    (f,slices)
        if (ldustdensity)  call get_slices_dustdensity (f,slices)
        if (ldustvelocity) call get_slices_dustvelocity(f,slices)
        if (lenergy)       call get_slices_energy      (f,slices)
        if (leos)          call get_slices_eos         (f,slices)
        if (lheatflux)     call get_slices_heatflux    (f,slices)
        if (lhydro)        call get_slices_hydro       (f,slices)
        if (linterstellar) call get_slices_interstellar(f,slices)
        if (lmagnetic)     call get_slices_magnetic    (f,slices)
        if (lparticles)    call get_slices_particles   (f,slices)
        if (lpscalar)      call get_slices_pscalar     (f,slices)
        if (lradiation)    call get_slices_radiation   (f,slices)
        if (lshock)        call get_slices_shock       (f,slices)
        if (lspecial)      call get_slices_special     (f,slices)
        if (ltestfield)    call get_slices_testfield   (f,slices)
        if (ltestflow)     call get_slices_testflow    (f,slices)
        if (ltestscalar)   call get_slices_testscalar  (f,slices)
!
        if (lslices_legacy) then
          inamev=inamev+1
          cycle
        endif
!
!  Slice, or component of slice, ready for saving.
!
        if (slices%ready) then
          if (slices%index==0) then    ! If this wasn't a multi slice...
            if (associated(slices%yz)) &
              call wslice(path//trim(slices%name)//'.yz',slices%yz, &
                                                     x(ix_loc),ny,nz)
            if (associated(slices%xz)) &
              call wslice(path//trim(slices%name)//'.xz',slices%xz, &
                                                     y(iy_loc),nx,nz)
            if (associated(slices%xz2)) &
              call wslice(path//trim(slices%name)//'.xz2',slices%xz2, &
                                                     y(iy2_loc),nx,nz)
            if (associated(slices%xy)) &
              call wslice(path//trim(slices%name)//'.xy',slices%xy, &
                                                     z(iz_loc),nx,ny)
            if (associated(slices%xy2)) &
              call wslice(path//trim(slices%name)//'.xy2',slices%xy2, &
                                                     z(iz2_loc),nx,ny)
            if (associated(slices%xy3).and.lwrite_slice_xy3) &
              call wslice(path//trim(slices%name)//'.xy3',slices%xy3, &
                                                     z(iz3_loc),nx,ny)
            if (associated(slices%xy4).and.lwrite_slice_xy4) &
              call wslice(path//trim(slices%name)//'.xy4',slices%xy4, &
                                                     z(iz4_loc),nx,ny)
            inamev=inamev+1
          else
            sindex=itoa(slices%index)
            if (associated(slices%yz)) &
              call wslice(path//trim(slices%name)//trim(sindex)//'.yz', &
                                       slices%yz, x(ix_loc),ny,nz)
            if (associated(slices%xz)) &
              call wslice(path//trim(slices%name)//trim(sindex)//'.xz', &
                                       slices%xz, y(iy_loc),nx,nz)
            if (associated(slices%xz2)) &
              call wslice(path//trim(slices%name)//'.xz2',slices%xz2, &
                                                     y(iy2_loc),nx,nz)
            if (associated(slices%xy)) &
              call wslice(path//trim(slices%name)//trim(sindex)//'.xy', &
                                       slices%xy, z(iz_loc),nx,ny)
            if (associated(slices%xy2)) &
              call wslice(path//trim(slices%name)//trim(sindex)//'.xy2', &
                                       slices%xy2, z(iz2_loc),nx,ny)
            if (associated(slices%xy3).and.lwrite_slice_xy3) &
              call wslice(path//trim(slices%name)//trim(sindex)//'.xy3', &
                                       slices%xy3, z(iz3_loc),nx,ny)
            if (associated(slices%xy4).and.lwrite_slice_xy4) &
              call wslice(path//trim(slices%name)//trim(sindex)//'.xy4', &
                                       slices%xy4, z(iz4_loc),nx,ny)
          endif
        else
          if (slices%index==0) then    ! If this wasn't a multi slice...
            if (lfirstloop.and.lroot) then
              if (lnewfile) then
                open(1,file='video.err',iostat=IOSTAT)
                lnewfile=.false.
              else
                open(1,file='video.err',position='append',iostat=IOSTAT)
              endif
              if (.not. outlog(iostat,'open','video.err',dist=-1)) then
                ! video.err is not distributed, backskipping enabled
!
                write(1,*,iostat=IOSTAT) 'unknown slice: ',trim(cnamev(inamev))
                if (.not. outlog(iostat,'write cnamev(inamev)')) then
!
                  close(1,iostat=IOSTAT)
                  if (outlog(iostat,'close')) continue
!
                endif
              endif
            endif
          else
            slices%index=0
          endif
          inamev=inamev+1
        endif
      enddo
!
      lfirstloop=.false.
!
    endsubroutine wvid
!***********************************************************************
    subroutine wslice(filename,a,pos,ndim1,ndim2)
!
!  appending to an existing slice file
!
!  12-nov-02/axel: coded
!  26-jun-06/anders: moved from Slices
!  22-sep-07/axel: changed Xy to xy2, to be compatible with Mac
!
      integer :: ndim1,ndim2
      character (len=*) :: filename
      real, dimension (ndim1,ndim2) :: a
      real, intent(in) :: pos
      integer :: iostat
!
!  check whether we want to write a slice on this processor
!
      if ( (lwrite_slice_xy .and.index(filename,'xy' )>0 &
          .and. index(filename,'xy2')==0 &
          .and. index(filename,'xy3')==0 &
          .and. index(filename,'xy4')==0 ) .or. &
          (lwrite_slice_xy2.and.index(filename,'xy2')>0) .or. &
          (lwrite_slice_xy3.and.index(filename,'xy3')>0) .or. &
          (lwrite_slice_xy4.and.index(filename,'xy4')>0) .or. &
          (lwrite_slice_xz .and.index(filename,'xz' )>0) .or. &
          (lwrite_slice_xz2.and.index(filename,'xz2')>0) .or. &
          (lwrite_slice_yz .and.index(filename,'yz' )>0) ) then
!
        open(1,file=filename,form='unformatted',position='append',IOSTAT=iostat)
!
!  files data/procN/slice*.* are distributed and will be synchronized on I/O error
!
        if (outlog(iostat,'open',filename,dist=1)) return
!
        write(1,IOSTAT=iostat) a,tslice,pos
        if (outlog(iostat,'write a,tslice,pos')) return
!
        close(1,IOSTAT=iostat)
        if (outlog(iostat,'close')) continue
!
      endif
!
    endsubroutine wslice
!***********************************************************************
    subroutine setup_slices_write()
!
!  27-Nov-2014/Bourdin.KIS: cleaned up the actual writing code
!
      integer :: iostat
      integer, parameter :: lun = 1
!
      open(lun,file=trim(directory)//'/slice_position.dat',STATUS='unknown',IOSTAT=iostat)
!
!  file 'data/procN/slice_position.dat' is distributed, but will not be synchronized
!  on I/O error (-> dist=0) as this would make it disfunctional; correct a posteriori if necessary
!
      if (outlog(iostat,'open',trim(directory)//'/slice_position.dat')) return
!
      write(lun,'(l5,i5)',IOSTAT=iostat) lwrite_slice_xy,iz_loc
      if (outlog(iostat,'write lwrite_slice_xy,iz_loc')) return
!
      write(lun,'(l5,i5)',IOSTAT=iostat) lwrite_slice_xy2,iz2_loc
      if (outlog(iostat,'write lwrite_slice_xy2,iz2_loc')) return
!
      write(lun,'(l5,i5)',IOSTAT=iostat) lwrite_slice_xy3,iz3_loc
      if (outlog(iostat,'write lwrite_slice_xy3,iz3_loc')) return
!
      write(lun,'(l5,i5)',IOSTAT=iostat) lwrite_slice_xy4,iz4_loc
      if (outlog(iostat,'write lwrite_slice_xy4,iz4_loc')) return
!
      write(lun,'(l5,i5)',IOSTAT=iostat) lwrite_slice_xz,iy_loc
      if (outlog(iostat,'write lwrite_slice_xz,iy_loc')) return
!
      write(lun,'(l5,i5)',IOSTAT=iostat) lwrite_slice_xz2,iy2_loc
      if (outlog(iostat,'write lwrite_slice_xz2,iy2_loc')) return
!
      write(lun,'(l5,i5)',IOSTAT=iostat) lwrite_slice_yz,ix_loc
      if (outlog(iostat,'write lwrite_slice_yz,ix_loc')) return
!
      close(lun,IOSTAT=iostat)
      if (outlog(iostat,'close')) return
!
    endsubroutine setup_slices_write
!***********************************************************************
    subroutine setup_slices()
!
!  Determine slice positions and whether slices are to be written on this
!  processor.
!
!  29-may-06/tobi: wrapped code from param_io.f90 into this subroutine
!  21-apr-15/MR: corrected i[xyz]_loc determination, see subroutine position
!
!  set slice position. The default for slice_position is 'p' for periphery,
!  although setting ix, iy, iz, iz2 by hand will overwrite this.
!  If slice_position is not 'p', then ix, iy, iz, iz2 are overwritten.
!
      if (slice_position=='p' .or. slice_position=='S') then
        ix_loc=l1; iy_loc=m1; iz_loc=n1; iz2_loc=n2
        lwrite_slice_xy2=llast_proc_z
        lwrite_slice_xy=lfirst_proc_z
        lwrite_slice_xz=lfirst_proc_y
        lwrite_slice_yz=lfirst_proc_x
!
!  slice position in middle of the box in dependence of nprocy,nprocz
!  second horizontal slice is the uppermost layer
!
      elseif (slice_position=='m') then
        if (mod(nprocx,2)==0) then; ix_loc=l1; else; ix_loc=(l1+l2)/2; endif
        if (mod(nprocy,2)==0) then; iy_loc=m1; else; iy_loc=(m1+m2)/2; endif
        if (mod(nprocz,2)==0) then; iz_loc=n1; else; iz_loc=(n1+n2)/2; endif
        iz2_loc=n2
!
!  xy2 is top layer as default.
!  Please set iz2 in run.in to select a different layer
!  where nghost <= iz2 <= mzgrid-nghost
!
        lwrite_slice_xy2=llast_proc_z
        lwrite_slice_xy=(ipz==nprocz/2)
        lwrite_slice_xz=(ipy==nprocy/2)
        lwrite_slice_yz=(ipx==nprocx/2)
!
!  slice positions for spherical coordinates
!  w is for "wedges" since the outputs are
!  the midplane (rphi,theta=y(mpoint)) and four
!  wedges in rtheta (xy)
!
      elseif (slice_position=='w') then
        if (nprocx>1) call warning('setup_slice', &
            'slice_position=w may be wrong for nprocx>1')
        !midplane slices
        !ix_loc=nxgrid/2+nghost
        iy = nygrid/2+nghost
        !meridional wedges, at 4 different
        !equally spaced azimuthal locations
        iz =  0*nzgrid/4+1+nghost
        iz2=  1*nzgrid/4+1+nghost
        iz3=  2*nzgrid/4+1+nghost
        iz4=  3*nzgrid/4+1+nghost
        !yz is not needed
        lwrite_slice_yz =.false.
!
!  Another slice positions for spherical coordinates
!  s is for "surface" meaning theta-phi sections
!  keep iz_loc=n1, corresponding to a meridional slice on n=n1.
!  The lwrite_slice_xz=.true. is needed for the read_video routine.
!
      elseif (slice_position=='s') then
        if (nprocx>1) call warning('setup_slice', &
            'slice_position=s may be wrong for nprocx>1')
        iz_loc=n1; iz2_loc=n2
        call xlocation(xtop_slice,ix_loc,lwrite_slice_yz)
        lwrite_slice_xy2=(ipz==nprocz/4)
        lwrite_slice_xy=lfirst_proc_z
        lwrite_slice_xz=.true.
!
! Another slice position for spherical coordinates, for global disks with
! buffer zones. It will read the midplane (xz2), and three other surfaces:
! the theta-phi wall at constant radius (yz); the meridional plane at constant
! azimuth (xy); and the upper "lid", a radius-azimuth surface at constant
! theta, in the upper disk (xz). Both xy and yz are taken 10 grid cells away from
! the beginning of the grid. This is to avoid the boundary.
!
      elseif (slice_position=='d') then
        ix_loc=l1+10; iy_loc=m1+10; iz_loc=n1
        if (mod(nprocy,2)==0) then; iy2_loc=m1; else; iy2_loc=(m1+m2)/2; endif
        lwrite_slice_xy2=.false.
        lwrite_slice_xy=lfirst_proc_z
        lwrite_slice_xz=lfirst_proc_y
        lwrite_slice_yz=lfirst_proc_x
        lwrite_slice_xz2=(ipy==nprocy/2)
!
!  later we may also want to write other slices
!
        !call xlocation(xtop_slice,ix2_loc,lwrite_slice_yz2)
        !call ylocation(ytop_slice,iy2_loc,lwrite_slice_xz2)
!
!  slice position when the first meshpoint in z is the equator (sphere)
!  For one z-processor, iz remains n1, but iz2 is set to the middle.
!
!  TH: A certain processor layout is implied here
!
      elseif (slice_position=='e') then
        if (nprocx>1) call warning('setup_slice', &
            'slice_position=e may be wrong for nrpocx>1')
        ix_loc=(l1+l2)/2; iy_loc=m1; iz_loc=n1; iz2_loc=n2
        if (nprocy==1) then; iy_loc=(m1+m2)/2; endif
        if (nprocz==1) then; iz2_loc=(iz+n2)/2; endif
        lwrite_slice_xy2=(ipz==nprocz/4)
        lwrite_slice_xy=lfirst_proc_z
        lwrite_slice_xz=(ipy==nprocy/2)
        lwrite_slice_yz=.true.
!
!  slice position similar to periphery (p), but the two z positions
!  can now be given in terms of z (zbot_slice, ztop_slice).
!
      elseif (slice_position=='c') then
        ix_loc=l1; iy_loc=m1
        call zlocation(zbot_slice,iz_loc,lwrite_slice_xy)
        call zlocation(ztop_slice,iz2_loc,lwrite_slice_xy2)
        lwrite_slice_xz=lfirst_proc_y
        lwrite_slice_yz=lfirst_proc_x
!
!  periphery of the box, but the other way around
!
      elseif (slice_position=='q') then
        ix_loc=l2; iy_loc=m2; iz_loc=n2; iz2_loc=n1
        lwrite_slice_xy2=lfirst_proc_z
        lwrite_slice_xy=llast_proc_z
        lwrite_slice_xz=llast_proc_y
        lwrite_slice_yz=llast_proc_x
      else
        if (lroot) then
          call fatal_error('setup_slices', &
                           'No such value for slice_position: '//slice_position)
        endif
      endif
!
!  Spherical admits only position 'w', 's', or 'd'. Break if this is not met.
!  Also, turn extra r-theta slices to false in case of
!  non-spherical coordinates
!
      if (slice_position/='w'.and.slice_position/='s') then
        lwrite_slice_xy3=.false.
        lwrite_slice_xy4=.false.
      endif
      if (slice_position/='d') then
        lwrite_slice_xz2=.false.
      endif
!
!  Overwrite slice positions if any ix,iy,iz,iz2,iz3,iz4 is greater then zero
!
      call position(ix,ipx,nx,ix_loc,lwrite_slice_yz)
      call position(iy,ipy,ny,iy_loc,lwrite_slice_xz)
      call position(iy2,ipy,ny,iy2_loc,lwrite_slice_xz2)
      call position(iz,ipz,nz,iz_loc,lwrite_slice_xy)
      call position(iz2,ipz,nz,iz2_loc,lwrite_slice_xy2)
      call position(iz3,ipz,nz,iz3_loc,lwrite_slice_xy3)
      call position(iz4,ipz,nz,iz4_loc,lwrite_slice_xy4)
!
      call setup_slices_write()
!
!  make sure ix_loc,iy_loc,iy2_loc,iz_loc,iz2_loc,iz3_loc,iz4_loc
!  are not outside the boundaries
!
       ix_loc=min( ix_loc,l2)
       ix_loc=max( ix_loc,l1)
       iy_loc=min( iy_loc,m2) ; iy2_loc=min( iy2_loc,m2)
       iy_loc=max( iy_loc,m1) ; iy2_loc=max( iy2_loc,m1)
       iz_loc=min( iz_loc,n2) ; iz2_loc=min(iz2_loc,n2)
       iz_loc=max( iz_loc,n1) ; iz2_loc=max(iz2_loc,n1)
      iz3_loc=min(iz3_loc,n2) ; iz4_loc=min(iz4_loc,n2)
      iz3_loc=max(iz3_loc,n1) ; iz4_loc=max(iz4_loc,n1)
!
      if (lroot) then
        write (*,*)'read_all_run_pars: slice_position = '//slice_position
        write (*,'(1x,a,4i4)') &
          'read_all_run_pars: ix,iy,iz,iz2 (video files) =',&
          ix,iy,iz,iz2
        if (slice_position=='s'.or.slice_position=='w') write (*,'(1x,a,2i4)') &
          'read_all_run_pars: iz3, iz4 (video files) =',iz3,iz4
        if (slice_position=='d') write (*,'(1x,a,1i4)') &
          'read_all_run_pars: iy2 (video files) =',iy2
      endif
!
    endsubroutine setup_slices
!***********************************************************************
endmodule Slices
