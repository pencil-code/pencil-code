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
  use Mpicomm, only: mpiallreduce_or
  use Sub, only: xlocation, zlocation, update_snaptime, read_snaptime, position
!
  implicit none
!
  public :: wvid, wvid_prepare, setup_slices
!
  real, public :: tvid
  integer, public :: nvid
!
  private
!
  real :: tslice=0.
  real, target, dimension(:,:), allocatable :: slice_xy,slice_xy2,slice_xy3,slice_xy4
  real, target, dimension(:,:), allocatable :: slice_xz,slice_xz2,slice_yz
  logical :: lactive_slice_yz, lactive_slice_xz, lactive_slice_xz2, &
      lactive_slice_xy, lactive_slice_xy2, lactive_slice_xy3, lactive_slice_xy4

!
contains
!***********************************************************************
    subroutine wvid_prepare
!
!  Prepare lvideo for writing slices into video file
!  This is useful for visualization of scalar field (or one component
!  of a vector field) on the periphery of a box.
!  Can be visualized in idl using rvid_box.pro
!
!  20-oct-97/axel: coded
!  08-oct-02/tony: increased size of file to handle datadir//'/tvid.dat'
!  13-nov-02/axel: added more fields
!  18-mar-03/axel: added dust velocity
!
      logical, save :: lfirst_call=.true.
      character (len=fnlen) :: file
!
!  Output vid-data in 'dvid' time intervals
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
!  Save current time so that the time that is written out in
!  output_slice() is not from the next time step
!
      if (lvideo) tslice = t
!
    endsubroutine wvid_prepare
!***********************************************************************
    subroutine wvid(f)
!
!  Write slices for animation of scalar field
!  (or one component of a vector field) on the perifery of a box.
!  Can be visualized in idl using rvid_box.pro.
!
!  13-nov-02/axel: added more fields
!  22-sep-07/axel: changed Xy to xy2, to be compatible with Mac
!  28-oct-18/PABourdin: moved output to IO modules 'output_slice'
!
      use General,         only: itoa
      use IO,              only: output_slice
      use Slices_methods,  only: assign_slices_scal
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
      real, dimension (mx,my,mz,mfarray), intent(IN) :: f
!
      logical :: lslices_legacy=.true.
      integer :: inamev
!
      type (slice_data) :: slices
      character (LEN=labellen) :: sname
!
      slices%index=0
!
!  Loop over slices.
!
      inamev=1
      do while (inamev <= nnamev)
!
        if (trim(cformv(inamev))=='') then
          inamev=inamev+1
          cycle           ! skip undefined slices
        endif

        sname=trim(cnamev(inamev))

        slices%name=sname
        call assign_slices_scal(slices,slice_xy,slice_xz,slice_yz,slice_xy2,slice_xy3,slice_xy4,slice_xz2)
        slices%ready=.false.
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
        if (ldensity .or. lanelastic) &
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
            inamev=inamev+1
          else
            sname=trim(sname)//trim(itoa(slices%index))
          endif
          if (lactive_slice_yz) call output_slice(lwrite_slice_yz, tslice, sname, 'yz', x(ix_loc), ix, slices%yz)
          if (lactive_slice_xz) call output_slice(lwrite_slice_xz, tslice, sname, 'xz', y(iy_loc), iy, slices%xz)
          if (lactive_slice_xz2) call output_slice(lwrite_slice_xz2, tslice, sname, 'xz2', y(iy2_loc), iy2, slices%xz2)
          if (lactive_slice_xy) call output_slice(lwrite_slice_xy, tslice, sname, 'xy', z(iz_loc), iz, slices%xy)
          if (lactive_slice_xy2) call output_slice(lwrite_slice_xy2, tslice, sname, 'xy2', z(iz2_loc), iz2, slices%xy2)
          if (lactive_slice_xy3) call output_slice(lwrite_slice_xy3, tslice, sname, 'xy3', z(iz3_loc), iz3, slices%xy3)
          if (lactive_slice_xy4) call output_slice(lwrite_slice_xy4, tslice, sname, 'xy4', z(iz4_loc), iz4, slices%xy4)
        else
          if (slices%index/=0) slices%index=0
          inamev=inamev+1
        endif
      enddo
!
    endsubroutine wvid
!***********************************************************************
    subroutine setup_slices
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
      !use Slices_methods, only: alloc_slice_buffers
      use General, only: itoa
      use IO, only: output_slice_position
    
      character(LEN=80) :: text, data

      lwrite_slice_xy=.false. 
      lwrite_slice_xz=.false. 
      lwrite_slice_yz=.false. 
      lwrite_slice_xy2=.false.
      lwrite_slice_xy3=.false.
      lwrite_slice_xy4=.false.
      lwrite_slice_xz2=.false.

      if (slice_position=='p' .or. slice_position=='S') then

        lwrite_slice_xy2=llast_proc_z; if (lwrite_slice_xy2)iz2_loc=n2
        lwrite_slice_xy=lfirst_proc_z; if (lwrite_slice_xy) iz_loc=n1
        lwrite_slice_xz=lfirst_proc_y; if (lwrite_slice_xz) iy_loc=m1
        lwrite_slice_yz=lfirst_proc_x; if (lwrite_slice_yz) ix_loc=l1
!
!  slice position in middle of the box in dependence of nprocy,nprocz
!  second horizontal slice is the uppermost layer
!
      elseif (slice_position=='m') then
!
!  xy2 is top layer as default.
!  Please set iz2 in run.in to select a different layer
!  where nghost <= iz2 <= mzgrid-nghost
!
        lwrite_slice_yz=(ipx==nprocx/2)
        if (lwrite_slice_yz) then
          if (mod(nprocx,2)==0) then; ix_loc=l1; else; ix_loc=(l1+l2)/2; endif
        endif
        lwrite_slice_xz=(ipy==nprocy/2)
        if (lwrite_slice_xz) then
          if (mod(nprocy,2)==0) then; iy_loc=m1; else; iy_loc=(m1+m2)/2; endif
        endif
        lwrite_slice_xy=(ipz==nprocz/2)
        if (lwrite_slice_xy) then
          if (mod(nprocz,2)==0) then; iz_loc=n1; else; iz_loc=(n1+n2)/2; endif
        endif
        lwrite_slice_xy2=llast_proc_z; if (lwrite_slice_xy2) iz2_loc=n2
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
        iz =  0*nzgrid/4+1+nghost           !MR: nghost not tb added!
        iz2=  1*nzgrid/4+1+nghost
        iz3=  2*nzgrid/4+1+nghost
        iz4=  3*nzgrid/4+1+nghost
!
!  Another slice positions for spherical coordinates
!  s is for "surface" meaning theta-phi sections
!  keep iz_loc=n1, corresponding to a meridional slice on n=n1.
!
      elseif (slice_position=='s') then

        if (nprocx>1) call warning('setup_slice', &
            'slice_position=s may be wrong for nprocx>1')

        call xlocation(xtop_slice,ix_loc,lwrite_slice_yz)
        lwrite_slice_xy2=(ipz==nprocz/4); if (lwrite_slice_xy2) iz2_loc=n2
        lwrite_slice_xy=lfirst_proc_z; if (lwrite_slice_xy) iz_loc=n1
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

        lwrite_slice_xy=lfirst_proc_z; if (lwrite_slice_xy) iz_loc=n1
        lwrite_slice_xz=lfirst_proc_y; if (lwrite_slice_xz) iy_loc=min(m1+10,m2)
        lwrite_slice_yz=lfirst_proc_x; if (lwrite_slice_yz) ix_loc=min(l1+10,l2)
        lwrite_slice_xz2=(ipy==nprocy/2)
        if (lwrite_slice_xz2) then
          if (mod(nprocy,2)==0) then; iy2_loc=m1; else; iy2_loc=(m1+m2)/2; endif
        endif
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
            'slice_position=e may be wrong for nprocx>1')

        lwrite_slice_xy=lfirst_proc_z; if (lwrite_slice_xy) iz_loc=n1
        lwrite_slice_yz=(ipx==nprocx/2); if (lwrite_slice_yz) ix_loc=(l1+l2)/2

        lwrite_slice_xy2=(ipz==nprocz/4)
        if (lwrite_slice_xy2) then
          if (nprocz==1) then
            iz2_loc=(n1+n2)/2   !MR: not iz2_loc=(iz+n2)/2!  
          else
            iz2_loc=n2
          endif
        endif

        lwrite_slice_xz=(ipy==nprocy/2)
        if (lwrite_slice_xz) then
          if (nprocy==1) then
            iy_loc=(m1+m2)/2
          else
            iy_loc=m1
          endif
        endif
!
!  slice position similar to periphery (p), but the two z positions
!  can now be given in terms of z (zbot_slice, ztop_slice).
!
      elseif (slice_position=='c') then

        call zlocation(zbot_slice,iz_loc,lwrite_slice_xy)
        call zlocation(ztop_slice,iz2_loc,lwrite_slice_xy2)
        lwrite_slice_xz=lfirst_proc_y; if (lwrite_slice_xz) iy_loc=m1
        lwrite_slice_yz=lfirst_proc_x; if (lwrite_slice_yz) ix_loc=l1
!
!  periphery of the box, but the other way around
!
      elseif (slice_position=='q') then

        lwrite_slice_xy2=lfirst_proc_z; if (lwrite_slice_xy2) iz2_loc=n1
        lwrite_slice_xy=llast_proc_z; if (lwrite_slice_xy) iz_loc=n2
        lwrite_slice_xz=llast_proc_y; if (lwrite_slice_xz) iy_loc=m2
        lwrite_slice_yz=llast_proc_x; if (lwrite_slice_yz) ix_loc=l2

      else
        if (lroot) &
          call fatal_error('setup_slices', &
                           'No such value for slice_position: '//slice_position)
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
!  position sets lwrite_slice_* according to whether or not the executing proc 
!  contributes to the slice data
!
      call position(ix,ipx,nx,ix_loc,lwrite_slice_yz)
      call position(iy,ipy,ny,iy_loc,lwrite_slice_xz)
      call position(iy2,ipy,ny,iy2_loc,lwrite_slice_xz2)
      call position(iz,ipz,nz,iz_loc,lwrite_slice_xy)
      call position(iz2,ipz,nz,iz2_loc,lwrite_slice_xy2)
      call position(iz3,ipz,nz,iz3_loc,lwrite_slice_xy3)
      call position(iz4,ipz,nz,iz4_loc,lwrite_slice_xy4)
      call mpiallreduce_or (lwrite_slice_yz, lactive_slice_yz)
      call mpiallreduce_or (lwrite_slice_xz, lactive_slice_xz)
      call mpiallreduce_or (lwrite_slice_xz2, lactive_slice_xz2)
      call mpiallreduce_or (lwrite_slice_xy, lactive_slice_xy)
      call mpiallreduce_or (lwrite_slice_xy2, lactive_slice_xy2)
      call mpiallreduce_or (lwrite_slice_xy3, lactive_slice_xy3)
      call mpiallreduce_or (lwrite_slice_xy4, lactive_slice_xy4)
!
      if (.not.lactive_dimension(3)) then
        lwrite_slice_xy2=.false.; lwrite_slice_xy3=.false.; lwrite_slice_xy4=.false.
        iz2_loc=-1; iz3_loc=-1; iz4_loc=-1
      endif 

      if (.not.lactive_dimension(2)) then
        lwrite_slice_xz2=.false.; iy2_loc=-1
      endif

      call output_slice_position
!
      if (lroot.and.dvid>0.) then  !MR: tbdone: write global position from all procs

        write (*,*)'setup_slices: slice_position = '//slice_position
        text=''; data=''
        if (lwrite_slice_yz) then
          text='ix_loc,'; data=itoa(ix_loc)
        endif
        if (lwrite_slice_xz) then
          text=trim(text)//'iy_loc,'; data=trim(data)//' '//itoa(iy_loc)
        endif
        if (lwrite_slice_xy) then
          text=trim(text)//'iz_loc,'; data=trim(data)//' '//itoa(iz_loc)
        endif
        if (lwrite_slice_xy2) then
          text=trim(text)//'iz2_loc,'; data=trim(data)//' '//itoa(iz2_loc)
        endif
        if (lwrite_slice_xy3) then
          text=trim(text)//'iz3_loc,'; data=trim(data)//' '//itoa(iz3_loc)
        endif
        if (lwrite_slice_xy4) then
          text=trim(text)//'iz4_loc,'; data=trim(data)//' '//itoa(iz4_loc)
        endif
        if (lwrite_slice_xz2) then
          text=trim(text)//'iy2_loc,'; data=trim(data)//' '//itoa(iy2_loc)
        endif
        write (*,*) 'setup_slices: '//trim(text)//' (video files) = '//data
      endif
!
      !call alloc_slice_buffers(slice_xy,slice_xz,slice_yz,slice_xy2,slice_xy3,slice_xy4,slice_xz2)
      if (lwrite_slice_xy .and..not.allocated(slice_xy) ) allocate(slice_xy (nx,ny))
      if (lwrite_slice_xz .and..not.allocated(slice_xz) ) allocate(slice_xz (nx,nz))
      if (lwrite_slice_yz .and..not.allocated(slice_yz) ) allocate(slice_yz (ny,nz))
      if (lwrite_slice_xy2.and..not.allocated(slice_xy2)) allocate(slice_xy2(nx,ny))
      if (lwrite_slice_xy3.and..not.allocated(slice_xy3)) allocate(slice_xy3(nx,ny))
      if (lwrite_slice_xy4.and..not.allocated(slice_xy4)) allocate(slice_xy4(nx,ny))
      if (lwrite_slice_xz2.and..not.allocated(slice_xz2)) allocate(slice_xz2(nx,nz))

    endsubroutine setup_slices
!***********************************************************************
    subroutine prep_xy_slice(izloc)

      use General, only: indgen

      integer, intent(IN) :: izloc

      real, dimension(nygrid/2) :: yloc
      !real, dimension(nygrid/2,1) :: thphprime

      if (ipz<=nprocz/3) then
!
! Line trough Southern cap
!
        yloc=xyz1(2)+indgen(nygrid/2)*dy

        !call yy_transform_strip_other(yloc,(/z(izloc)/),thphprime)
        !nok=prep_interp(thphprime,intcoeffs)

      elseif (ipz>=2*nprocz/3) then
!
! Line trough Nouthern cap
!
        yloc=xyz0(2)-(nygrid/2+1-indgen(nygrid/2))*dy 

        !call yy_transform_strip_other(yloc,(/z(izloc)/),thphprime)
        !nok=prep_interp(thphprime,intcoeffs,iyinyang_intpol_type,thrange_cap)

      endif 

    endsubroutine prep_xy_slice
!***********************************************************************
endmodule Slices
