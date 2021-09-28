! $Id$
!
!  Module for boundary conditions. Extracted from (no)mpicomm, since
!  all non-periodic (external) boundary conditions require the same
!  code for serial and parallel runs.
!
module Boundcond
!
  use Cdata
  use Cparam
  use Messages
  use Mpicomm
  use Deriv, only: set_ghosts_for_onesided_ders, bval_from_neumann, bval_from_3rd
  use General
!
  implicit none
!
  private
!
  public :: update_ghosts, zero_ghosts, finalize_boundcond
  public :: boundconds, boundconds_x, boundconds_y, boundconds_z
  public :: boundconds_x_c, boundconds_y_c, boundconds_z_c
  public :: bc_pencil
  public :: bc_per_x, bc_per_y, bc_per_z
  public :: set_consistent_density_boundary
  public :: set_consistent_vel_boundary
  public :: copy_BCs
  public :: set_periodic_boundcond_on_aux
  public :: jet_x
  public :: initialize_boundcond
!
  interface update_ghosts
    module procedure update_ghosts_all
    module procedure update_ghosts_range
  endinterface
!
  interface zero_ghosts
    module procedure zero_ghosts_all
    module procedure zero_ghosts_range
  endinterface
!
  interface bc_pencil
    module procedure bc_pencil_scalar
    module procedure bc_pencil_vector
  endinterface
!
  integer, parameter :: BOT=1, TOP=2
  logical :: is_vec=.false.
  integer :: jdone=0

  !real, allocatable, dimension(:,:,:) :: slc_dat_xy, slc_dat_xy2, slc_dat_xz, slc_dat_xz2, &
  !                                       slc_dat_yz, slc_dat_yz2
  type(scattered_array), pointer :: slc_dat_xy, slc_dat_xy2, slc_dat_xz, slc_dat_xz2, &
                                    slc_dat_yz, slc_dat_yz2
  integer :: nt_slices=0
  integer, parameter :: sz_slc_chunk=20
!
  contains
!***********************************************************************
    subroutine update_ghosts_all(f)
!
!  Update all ghost zones of f.
!
!  21-sep-02/wolf: extracted from wsnaps
!  28-mar-17/MR: added registration of already communicated variable ranges in f. 
!
      use Grid, only: coarsegrid_interp

      real, dimension (:,:,:,:) :: f
!
      if (ighosts_updated>=0) then
!
!  If registration is activated, register all variables to have been communicated.
!
        ighosts_updated=ighosts_updated+1
        updated_var_ranges(:,ighosts_updated)=(/1,min(mcom,size(f,4))/)
      endif
!
      call boundconds_x(f)
      call initiate_isendrcv_bdry(f)
      call finalize_isendrcv_bdry(f)
      if (lcoarse) call coarsegrid_interp(f)
      call boundconds_y(f)
      call boundconds_z(f)
!
    endsubroutine update_ghosts_all
!***********************************************************************
    subroutine update_ghosts_range(f,ivar1,ivar2_opt)
!
!  Update specific ghost zones of f.
!
!  11-aug-11/wlad: adapted from update_ghosts
!  28-mar-17/MR: added registration of already communicated variable ranges in f. 
!
      use General, only: add_merge_range
      use Grid, only: coarsegrid_interp
!
      real, dimension (:,:,:,:) :: f
      integer  :: ivar1,ivar2
      integer, optional :: ivar2_opt
!
      integer :: nact_ranges,i
!
      ivar2=ivar1
      if (present(ivar2_opt)) ivar2=ivar2_opt
!
      if (ighosts_updated>=0) then
!
!  If registration is activated, figure out which variables out of the range (ivar1,ivar2) have yet to be communicated.
!  These are appended as a set of ranges to the list of ranges in updated_var_ranges after position ighosts_updated.
!  The new total number of variable ranges to be communicated is nact_ranges.
!
        nact_ranges=add_merge_range( updated_var_ranges, ighosts_updated, (/ivar1,ivar2/) )
!
        if (nact_ranges>ighosts_updated) then
          do i=ighosts_updated+1,nact_ranges
            call boundconds_x(f,updated_var_ranges(1,i),updated_var_ranges(2,i))
            call initiate_isendrcv_bdry(f,updated_var_ranges(1,i),updated_var_ranges(2,i))
            call finalize_isendrcv_bdry(f,updated_var_ranges(1,i),updated_var_ranges(2,i))
            if (lcoarse) &
              call coarsegrid_interp(f,updated_var_ranges(1,i),updated_var_ranges(2,i))
            call boundconds_y(f,updated_var_ranges(1,i),updated_var_ranges(2,i))
            call boundconds_z(f,updated_var_ranges(1,i),updated_var_ranges(2,i))
          enddo
          ighosts_updated=nact_ranges
        endif
      else
        call boundconds_x(f,ivar1,ivar2)
        call initiate_isendrcv_bdry(f,ivar1,ivar2)
        call finalize_isendrcv_bdry(f,ivar1,ivar2)
        if (lcoarse) call coarsegrid_interp(f)
        call boundconds_y(f,ivar1,ivar2)
        call boundconds_z(f,ivar1,ivar2)
      endif
!if (lroot) print*, 'update_ghosts_range'
!
    endsubroutine update_ghosts_range
!***********************************************************************
    subroutine zero_ghosts_all(f)
!
!  Zeros the ghost cells for all variables.
!
!  23-oct-13/ccyang: coded.
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      call zero_ghosts_range(f, 1, mfarray)
!
    endsubroutine zero_ghosts_all
!***********************************************************************
    subroutine zero_ghosts_range(f, ivar1, ivar2_opt)
!
!  Zeros the ghost cells for variables ivar1:ivar2.
!
!  23-oct-13/ccyang: coded.
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      integer, intent(in) :: ivar1
      integer, intent(in), optional :: ivar2_opt
!
      integer :: ivar2
!
      ivar2 = ivar1
      if (present(ivar2_opt)) ivar2 = ivar2_opt
!
      xdir: if (nxgrid > 1) then
        f(1:nghost,m1:m2,n1:n2,ivar1:ivar2) = 0.0
        f(mx-nghost+1:mx,m1:m2,n1:n2,ivar1:ivar2) = 0.0
      endif xdir
!
      ydir: if (nygrid > 1) then
        f(:,1:nghost,n1:n2,ivar1:ivar2) = 0.0
        f(:,my-nghost+1:my,n1:n2,ivar1:ivar2) = 0.0
      endif ydir
!
      zdir: if (nzgrid > 1) then
        f(:,:,1:nghost,ivar1:ivar2) = 0.0
        f(:,:,mz-nghost+1:mz,ivar1:ivar2) = 0.0
      endif zdir
!
    endsubroutine zero_ghosts_range
!***********************************************************************
    subroutine initialize_boundcond
!
!  Initialization for reading boundary values from slices:
!  bc_slc_dir - working directory where these reside.
!
      use Sub, only: position
      use File_io, only: file_exists
      use Syscalls, only: directory_exists
      use IO, only: IO_strategy
      use HDF5_IO, only: input_dim

      integer :: ix_bc,ix2_bc,iy_bc,iy2_bc,iz_bc,iz2_bc,idum
      logical :: lread_slice_yz,lread_slice_yz2,lread_slice_xz,lread_slice_xz2, &
                 lread_slice_xy,lread_slice_xy2
      integer ::  mx_in, my_in, mz_in, mvar_in, maux_in, mglobal_in, &
                  nghost_in, nprocx_in, nprocy_in, nprocz_in, nprocz_in_
      logical :: lbcxslc,lbcyslc,lbczslc
      character :: prec_in
      character(LEN=3) :: suff_xy2, suff_xz2, suff_yz2
!
! Set proper BC code for Yin-Yang grid
!
      if (lyinyang) then
        if (lroot) call information('read_all_run_pars', 'all BCs for y and z ignored because of Yin-Yang grid')
        lperi(2:3) = .false.; lpole = .false.
        !bcy='yy'; bcz='yy'    ! not needed when interpolating spherical
        !components of vectors
        bcy='nil'; bcz='nil'
      endif
!
      call check_consistency_of_lperi('initialize_boundcond')
!
      lbcxslc=any(bcx12=='slc'); lbcyslc=any(bcy12=='slc'); lbczslc=any(bcz12=='slc')
      if (lbcxslc.and..not.lactive_dimension(1)) return
      if (lbcyslc.and..not.lactive_dimension(2)) return
      if (lbczslc.and..not.lactive_dimension(3)) return

      lread_slice_xy=.false.; lread_slice_xy2=.false.
      lread_slice_xz=.false.; lread_slice_xz2=.false.
      lread_slice_yz=.false.; lread_slice_yz2=.false.

      if (lbcxslc.or.lbcyslc.or.lbczslc) then

        if (lroot) then

          if (.not.directory_exists(trim(bc_slc_dir)//'/data')) &
            call fatal_error('initialize_boundconds', 'working directory '//trim(bc_slc_dir)//' for slices not found')

          call input_dim(bc_slc_dir, mx_in, my_in, mz_in, mvar_in, maux_in, mglobal_in, &
                         prec_in, nghost_in, nprocx_in, nprocy_in, nprocz_in)
          call input_dim(bc_slc_dir, mx_in, my_in, mz_in, mvar_in, maux_in, mglobal_in, &
                         prec_in, nghost_in, nprocx_in, nprocy_in, nprocz_in_, local=.true.)
          if (mx/=mx_in.or.my/=my_in) &
            call fatal_error('initialize_boundconds', &
                 'data in working directory '//trim(bc_slc_dir)//'have incompatible x or y dimensions')
        endif

        call mpibarrier
        call mpibcast_int(nprocz_in)

        if (lbcxslc) then

          if (any(bcx12(:,1)=='slc')) then
            idum=1
            call position(idum,ipx,nx,ix_bc,lread_slice_yz)
          endif

          if (any(bcx12(:,2)=='slc')) then
            idum=nxgrid
            call position(idum,ipx,nx,ix2_bc,lread_slice_yz2)
            if (lread_slice_yz) then
              suff_yz2='yz2'
            else
              suff_yz2='yz'
            endif
          endif

        endif

        if (lbcyslc) then

          if (any(bcy12(:,1)=='slc')) then
            idum=1
            call position(idum,ipy,ny,iy_bc,lread_slice_xz)
          endif

          if (any(bcy12(:,2)=='slc')) then
            idum=nygrid
            call position(idum,ipy,ny,iy2_bc,lread_slice_xz2)
            if (lread_slice_xz) then
              suff_xz2='xz2'
            else
              suff_xz2='xz'
            endif
          endif

        endif

        if (lbczslc) then

          if (any(bcz12(:,1)=='slc')) then
            idum=1
            call position(idum,ipz,nz,iz_bc,lread_slice_xy)
          endif

          if (any(bcz12(:,2)=='slc')) then
            idum=nzgrid
            call position(idum,ipz,nz,iz2_bc,lread_slice_xy2)
            if (lread_slice_xy) then
              suff_xy2='xy2'
            else 
              suff_xy2='xy'
            endif
          endif

        endif
!
!  Restricted to slice position 'm'!
!
        if (IO_strategy/='HDF5') then
          if (lread_slice_xy ) then
            call init_scattered_array(slc_dat_xy,nx,ny,mvar,sz_slc_chunk,lreloading)
            call get_slice_data(z(n1),find_proc(ipx,ipy,nprocz_in/2-1),'xy',slc_dat_xy,nt_slices)
          endif
          if (lread_slice_xy2) then 
            call init_scattered_array(slc_dat_xy2,nx,ny,mvar,sz_slc_chunk,lreloading)
            call get_slice_data(z(n2),find_proc(ipx,ipy,nprocz_in/2-1),suff_xy2,slc_dat_xy2,nt_slices)
          endif
          if (lread_slice_xz ) then 
            call init_scattered_array(slc_dat_xz ,nx,nz,mvar,sz_slc_chunk,lreloading)
            call get_slice_data(y(m1),find_proc(ipx,nprocy_in/2-1,ipz),'xz',slc_dat_xz,nt_slices)
          endif
          if (lread_slice_xz2) then 
            call init_scattered_array(slc_dat_xz2,nx,nz,mvar,sz_slc_chunk,lreloading)
            call get_slice_data(y(m2),find_proc(ipx,nprocy_in/2-1,ipz),suff_xz2,slc_dat_xz2,nt_slices)
          endif
          if (lread_slice_yz ) then 
            call init_scattered_array(slc_dat_yz ,ny,nz,mvar,sz_slc_chunk,lreloading)
            call get_slice_data(x(l1),find_proc(nprocx_in/2-1,ipy,ipz),'yz',slc_dat_yz,nt_slices)
          endif
          if (lread_slice_yz2) then 
            call init_scattered_array(slc_dat_yz2,ny,nz,mvar,sz_slc_chunk,lreloading)
            call get_slice_data(x(l2),find_proc(nprocx_in/2-1,ipy,ipz),suff_yz2,slc_dat_yz2,nt_slices)
          endif

          call mpibarrier
        else
          call fatal_error('initialize_boundcond','BC set from slice data not implemented for IO_strategy="HDF5"')
        endif      
      endif      

    endsubroutine initialize_boundcond
!***********************************************************************
    subroutine check_consistency_of_lperi(label)
!
!  Check consistency of lperi.
!
!  18-jul-03/axel: coded
!
      character (len=*) :: label
      logical :: lwarning=.true.
      integer :: j
!
!  Identifier.
!
      if (lroot.and.ip<5) print*,'check_consistency_of_lperi: called from',label
!
!  Make the warnings less dramatic looking, if we are only in start
!  and exit this routine altogether if, in addition, ip > 13.
!
      if (label=='check_consistency_of_lperi'.and.ip>13) return
      if (label=='check_consistency_of_lperi') lwarning=.false.
!
      if (nvar > 0) then
!
!  Check x direction.
!
        j=1
        if (any(bcx(1:nvar)=='p'.or. bcx(1:nvar)=='she').and..not.lperi(j).or.&
            any(bcx(1:nvar)/='p'.and.bcx(1:nvar)/='she').and.lperi(j)) &
            call warning_lperi(lwarning,bcx(1:nvar),lperi,j)
!
!  Check y direction.
!
        j=2
        if (any(bcy(1:nvar)=='p').and..not.lperi(j).or.&
            any(bcy(1:nvar)/='p').and.lperi(j)) &
            call warning_lperi(lwarning,bcy(1:nvar),lperi,j)
!
!  Check z direction.
!
        j=3
        if (any(bcz(1:nvar)=='p').and..not.lperi(j).or.&
            any(bcz(1:nvar)/='p').and.lperi(j)) &
            call warning_lperi(lwarning,bcz(1:nvar),lperi,j)
      endif
!
!  Print final warning.
!  Make the warnings less dramatic looking, if we are only in start.
!
      if (lroot .and. (.not. lwarning)) then
        if (label=='check_consistency_of_lperi') then
          print*,'[bad BCs in start.in only affects post-processing' &
               //' of start data, not the run]'
        else
          print*,'check_consistency_of_lperi(run.in): you better stop and check!'
          print*,'------------------------------------------------------'
          print*
        endif
      endif
!
    endsubroutine check_consistency_of_lperi
!***********************************************************************
    subroutine warning_lperi(lwarning,bc,lperi,j)
!
!  Print consistency warning of lperi.
!
!  18-jul-03/axel: coded
!
      character (len=*), dimension(:), intent(in) :: bc
      logical, dimension(3) :: lperi
      logical :: lwarning
      integer :: j
!
      if (lroot) then
        if (lwarning) then
          print*
          print*,'------------------------------------------------------'
          print*,'W A R N I N G'
          lwarning=.false.
        else
          print*
        endif
!
        print*,'warning_lperi: inconsistency, j=', j, ', lperi(j)=',lperi(j)
        print*,'bc=',bc
        print*,"any(bc=='p'.or. bc=='she'), .not.lperi(j) = ", &
          any(bc=='p'.or. bc=='she'), .not.lperi(j)
        print*, "any(bcx/='p'.and.bcx/='she'), lperi(j) = ", &
          any(bc=='p'.or. bc=='she'), .not.lperi(j)
      endif
!
    endsubroutine warning_lperi
!***********************************************************************
    subroutine get_slice_data(pos,iproc_slc,label,slcdat,nt)
!
      use General, only: itoa
      use File_io, only: file_exists
      use IO, only: input_slice

      real, intent(IN) :: pos
      integer, intent(IN) :: iproc_slc
      character(LEN=*) :: label
      type(scattered_array), pointer :: slcdat
      integer :: nt

      character(LEN=fnlen) :: slicedir, file
      real :: pos_slc
     
      slicedir=trim(bc_slc_dir)//'/data/proc'//trim(itoa(iproc_slc))

      if (lhydro) then
        file=trim(slicedir)//'/slice_uu1.'//trim(label)
        call input_slice(file,pos_slc,slcdat,iux,nt)
!print*, 'pos,pos_slc=', pos,pos_slc
!print*, 'iproc, ux:', minval(slcdat(:,:,iux)), maxval(slcdat(:,:,iux))
        if (abs(pos-pos_slc)>dz) &
          call fatal_error_local('get_slice_data', 'slices in '//trim(file)// &
                                 ' at wrong position')
        file=trim(slicedir)//'/slice_uu2.'//trim(label)
        call input_slice(file,pos_slc,slcdat,iuy,nt)
!print*, 'iproc, uy:', minval(slcdat(:,:,iuy)), maxval(slcdat(:,:,iuy))
        file=trim(slicedir)//'/slice_uu3.'//trim(label)
        call input_slice(file,pos_slc,slcdat,iuz,nt)
!print*, 'iproc, uz:', minval(slcdat(:,:,iuz)), maxval(slcdat(:,:,iuz))
      endif

      if (ldensity) then
        file=trim(slicedir)//'/slice_lnrho.'//trim(label)
        call input_slice(file,pos_slc,slcdat,ilnrho,nt)
        if (abs(pos-pos_slc)>dz) &
          call fatal_error_local('get_slice_data', 'slices in '//trim(file)// &
                                 ' at wrong position')
      endif

      if (lentropy) then
        file=trim(slicedir)//'/slice_ss.'//trim(label)
        call input_slice(file,pos_slc,slcdat,iss,nt)
!print*, 'iproc, ss:', minval(slcdat(:,:,iss)), maxval(slcdat(:,:,iss))
        if (abs(pos-pos_slc)>dz) &
          call fatal_error_local('get_slice_data', 'slices in '//trim(file)// &
                                 ' at wrong position')
      endif

      if (lmagnetic) then
        file=trim(slicedir)//'/slice_aa.'//trim(label)
        call input_slice(file,pos_slc,slcdat,iaa,nt)
        if (abs(pos-pos_slc)>dz) &
          call fatal_error_local('get_slice_data', 'slices in '//trim(file)// &
                                 ' at wrong position')
      endif
 
    endsubroutine get_slice_data
!***********************************************************************
    subroutine set_from_slice_x(f,topbot,j)

      use General, only: get_scattered_array

      real, dimension (:,:,:,:) :: f
      character (len=bclen) :: topbot
      integer :: j

      integer, save :: ilayer=-1
      real, save :: last_gettime, timediff
      real, dimension(ny,nz,mvar), save :: ahead_data

      if (lfirst) then
        if (ilayer==-1) then
          last_gettime=t
          ilayer=0
        elseif (t-last_gettime>=timediff) then
          ilayer=mod(ilayer+1,nt_slices)
          last_gettime=t
        endif
      endif

      if (topbot=='bot') then
        call get_scattered_array(j,ilayer,slc_dat_yz,f(l1,m1:m2,n1:n2,j),timediff,ahead_data(:,:,j))
      else
        call get_scattered_array(j,ilayer,slc_dat_yz2,f(l2,m1:m2,n1:n2,j),timediff,ahead_data(:,:,j))
      endif

    endsubroutine set_from_slice_x
!***********************************************************************
    subroutine set_from_slice_y(f,topbot,j)

      use General, only: get_scattered_array

      real, dimension (:,:,:,:) :: f
      character (len=bclen) :: topbot
      integer :: j
    
      integer, save :: ilayer=-1
      real, save :: last_gettime, timediff
      real, dimension(nx,nz,mvar), save :: ahead_data

      if (lfirst) then
        if (ilayer==-1) then
          last_gettime=t
          ilayer=0
        elseif (t-last_gettime>=timediff) then
          ilayer=mod(ilayer+1,nt_slices)
          last_gettime=t
        endif
      endif

      if (topbot=='bot') then
        call get_scattered_array(j,ilayer,slc_dat_xz,f(l1:l2,m1,n1:n2,j),timediff,ahead_data(:,:,j))
      else
        call get_scattered_array(j,ilayer,slc_dat_xz2,f(l1:l2,m2,n1:n2,j),timediff,ahead_data(:,:,j))
      endif

    endsubroutine set_from_slice_y
!***********************************************************************
    subroutine set_from_slice_z(f,topbot,j)

      use General, only: get_scattered_array

      real, dimension (:,:,:,:) :: f
      character (len=bclen) :: topbot
      integer :: j

      logical :: lget,lboth
      integer, dimension(mvar), save :: ilayer=0
      real, dimension(mvar), save :: last_gettime
      real, save :: timediff
      real, dimension(nx,ny,mvar), save :: ahead_data
      real :: w

      lget=.false.
      if (itsub==0.or.lfirst) then   
!
! update only in first substep of integration or before integration has started
!
        if (ilayer(j)==0) then
          last_gettime=t
          ilayer(j)=1
          lget=.true.
        elseif (t-last_gettime(j)>=timediff) then
          if (ilayer(j)==nt_slices) then
            ilayer(j)=1
          else
            ilayer(j)=ilayer(j)+1   !mod(ilayer(j)+1,nt_slices)
          endif
          last_gettime(j)=t
          lget=.true.
        endif
      endif

      lboth = ilayer(j)==1 .or. nt_slices==1
      if (.not.lget) w = (t-last_gettime(j))/timediff

      if (topbot=='bot') then
        if (lget) then
          if (lboth) then
            call get_scattered_array(j,ilayer(j),slc_dat_xy,f(l1:l2,m1:m2,n1,j),timediff,ahead_data(:,:,j))
          else
            f(l1:l2,m1:m2,n1,j)=ahead_data(:,:,j)
            call get_scattered_array(j,ilayer(j)+1,slc_dat_xy,ahead_data(:,:,j),timediff)
          endif
          if (nt_slices>1) f(l1:l2,m1:m2,1,j)=f(l1:l2,m1:m2,n1,j)    ! store obtained slice in unused ghost layer
        else
          if (nt_slices>1) then
            f(l1:l2,m1:m2,n1,j)=(1.-w)*f(l1:l2,m1:m2,1,j)+w*ahead_data(:,:,j)
          else
            f(l1:l2,m1:m2,n1,j)=f(l1:l2,m1:m2,1,j)
          endif
        endif
      else
        if (lget) then
!if (iproc==120) write(103,*) it,itsub,j,ilayer(j)+1,lboth
          if (lboth) then
            call get_scattered_array(j,ilayer(j),slc_dat_xy2,f(l1:l2,m1:m2,n2,j),timediff,ahead_data(:,:,j))
          else
            f(l1:l2,m1:m2,n2,j)=ahead_data(:,:,j)
            call get_scattered_array(j,ilayer(j)+1,slc_dat_xy2,ahead_data(:,:,j),timediff)
          endif
          if (nt_slices>1) f(l1:l2,m1:m2,mz,j)=f(l1:l2,m1:m2,n2,j)    ! store obtained slice in unused ghost layer
        else
!if (iproc==120) write(102,*) it,itsub,j,w
          if (nt_slices>1) then
            f(l1:l2,m1:m2,n2,j)=(1.-w)*f(l1:l2,m1:m2,mz,j)+w*ahead_data(:,:,j)
          else
            f(l1:l2,m1:m2,n2,j)=f(l1:l2,m1:m2,mz,j)
          endif
        endif
      endif

    endsubroutine set_from_slice_z
!***********************************************************************
    subroutine boundconds(f,ivar1_opt,ivar2_opt)
!
!  Apply boundary conditions in all three directions.
!  Note that we _must_ call boundconds_{x,y,z} in this order, or edges and
!  corners will not be OK.
!
!  10-oct-02/wolf: coded
!
      real, dimension (:,:,:,:) :: f
      integer, optional :: ivar1_opt, ivar2_opt
      integer :: ivar1, ivar2
!
      ivar1=1; ivar2=min(mcom,size(f,4))
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
!
      call boundconds_x(f,ivar1,ivar2)
      call boundconds_y(f,ivar1,ivar2)
      call boundconds_z(f,ivar1,ivar2)
!
    endsubroutine boundconds
!***********************************************************************
    subroutine boundconds_x_c(f,ivar1_opt,ivar2_opt)
!
!  Envelope for being called from C code.
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, optional :: ivar1_opt, ivar2_opt

      call boundconds_x(f,ivar1_opt,ivar2_opt)

    endsubroutine boundconds_x_c
!***********************************************************************
    subroutine boundconds_x(f,ivar1_opt,ivar2_opt)
!
!  Boundary conditions in x, except for periodic part handled by communication.
!  Remark: boundconds_x() needs to be called before communicating (because we
!  communicate the x-ghost points), boundconds_[yz] after communication
!  has finished (they need some of the data communicated for the edges
!  (yz-'corners').
!
!   8-jul-02/axel: split up into different routines for x,y and z directions
!  11-nov-02/wolf: unified bot/top, now handled by loop
!  15-dec-06/wolf: Replaced "if (bcx1(1)=='she') then" by "any" command
!  30-sep-16/MR: new BCs 'n1s' = Neumann + 1-sided, "c1s" = heatflux + 1-sided (only for z boundaries)
!
      use EquationOfState
      use Shear
      use Special, only: special_boundconds
!
      real, dimension (:,:,:,:) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      integer :: ivar1, ivar2, j, k
      logical :: ip_ok
      character (len=bclen) :: topbot
      type (boundary_condition) :: bc
!
      if (ldebug) print*, 'boundconds_x: ENTER'
!
      ivar1=1; ivar2=min(mcom,size(f,4))
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
!
!!print*, 'boundconds_x: mcom,mfarray,ivar1,ivar2=', mcom,mfarray,ivar1,ivar2
      select case (nxgrid)
!
      case (1)
        if (ldebug) print*, 'boundconds_x: no x-boundary'
!
!  Boundary conditions in x.
!
      case default
!
!  Use the following construct to keep compiler from complaining if
!  we have no variables (and boundconds) at all (samples/no-modules):
!
        if (all(bcx12(ivar1:ivar2,:)=='she')) then
          call boundcond_shear(f,ivar1,ivar2)
        else
          do k=1,2
            ! loop over 'bot','top'
            if (k==1) then
              topbot='bot'; ip_ok=lfirst_proc_x
            else
              topbot='top'; ip_ok=llast_proc_x
            endif
!
            do j=ivar1,ivar2
!
! Natalia: the next line is for the dustdensity case.
! If ndustspec is large, it is stupid to set bc for all dust species
! in start.in. But if one does not set them, they become 'p' by default
! Since this problem is crucial only for aerosol + chemistry
! the following condition is used. But this place should be modifyed somehow
! Any ideas?
!
!            if ((bcx12(j,k)=='p') .and. lchemistry .and. ldustdensity) bcx12(j,k)=''
!
              if (ldebug) write(*,'(A,I1,A,I2,A,A)') ' bcx',k,'(',j,')=',bcx12(j,k)
              if (bcx12(j,k) == 'she') then
                if (bcx12(j,1) /= bcx12(j,2)) &
                  call fatal_error_local('boundconds_x', 'generalize me to have sheared periodic boundary on only one end.')
                if (k == 1) call boundcond_shear(f, j, j)
              elseif (ip_ok) then
                select case (bcx12(j,k))
                case ('0')
                  ! BCX_DOC: zero value in ghost zones, free value on boundary
                  call bc_zero_x(f,topbot,j)
                case ('p')
                  ! BCX_DOC: periodic
                  call bc_per_x(f,topbot,j)
                case ('s')
                  ! BCX_DOC: symmetry, $f_{N+i}=f_{N-i}$;
                  ! BCX_DOC: implies $f'(x_N)=f'''(x_0)=0$
                  call bc_sym_x(f,+1,topbot,j)
                case ('sf')
                  ! BCX_DOC: symmetry with respect to interface
                  call bc_sf_x(f,+1,topbot,j)
                case ('ss')
                  ! BCX_DOC: symmetry, plus function value given
                  call bc_symset_x(f,+1,topbot,j,val=fbcx(:,k))
                case ('sds')
                  ! BCY_DOC: symmetric-derivative-set
                  call bc_symderset_x(f,topbot,j,val=fbcx(:,k))
                case ('s0d')
                  ! BCX_DOC: symmetry, function value such that df/dx=0
                  call bc_symset0der_x(f,topbot,j)
                case ('a')
                  ! BCX_DOC: antisymmetry, $f_{N+i}=-f_{N-i}$;
                  ! BCX_DOC: implies $f(x_N)=f''(x_0)=0$
                  call bc_sym_x(f,-1,topbot,j)
                case ('af')
                  ! BCX_DOC: antisymmetry with respect to interface
                  call bc_sf_x(f,-1,topbot,j)
                case ('a2')
                  ! BCX_DOC: antisymmetry relative to boundary value,
                  ! BCX_DOC: $f_{N+i}=2 f_{N}-f_{N-i}$;
                  ! BCX_DOC: implies $f''(x_0)=0$
                  call bc_sym_x(f,-1,topbot,j,REL=.true.)
                case ('a2v')
                  ! BCX_DOC: set boundary value and antisymmetry relative to it
                  ! BCX_DOC: $f_{N+i}=2 f_{N}-f_{N-i}$;
                  ! BCX_DOC: implies $f''(x_0)=0$
                  call bc_sym_x(f,-1,topbot,j,REL=.true.,val=fbcx(:,k))
                case ('a2r')
                  ! BCX_DOC: sets $d^2f/dr^2 +2df/dr- 2f/r^2 = 0$
                  ! BCX_DOC: This is the replacement of zero second derivative
                  ! BCX_DOC: in spherical coordinates, in radial direction.
                  call bc_a2r_x(f,topbot,j)
                case ('cpc')
                  ! BCX_DOC: cylindrical perfect conductor
                  ! BCX_DOC: implies $f''+f'/R=0$
                  call bc_cpc_x(f,topbot,j)
                case ('cpp')
                  ! BCX_DOC: cylindrical perfect conductor
                  ! BCX_DOC: implies $f''+f'/R=0$
                  call bc_cpp_x(f,topbot,j)
                case ('cpz')
                  ! BCX_DOC: cylindrical perfect conductor
                  ! BCX_DOC: implies $f''+f'/R=0$
                  call bc_cpz_x(f,topbot,j)
                case ('spr')
                  ! BCX_DOC: spherical perfect conductor
                  ! BCX_DOC: implies $f''+2f'/R=0$ and $f(x_N)=0$
                  call bc_spr_x(f,topbot,j)
                case ('v')
                  ! BCX_DOC: vanishing third derivative
                  call bc_van_x(f,topbot,j)
                case ('cop')
                  ! BCX_DOC: copy value of last physical point to all ghost cells
                  call bc_copy_x(f,topbot,j)
                case ('1s')
                  ! BCX_DOC: onesided
                  call set_ghosts_for_onesided_ders(f,topbot,j,1)
                case ('d1s')
                  ! BCX_DOC: onesided for 1st/2nd derivative in two first inner points, Dirichlet in boundary point
                  if (k==1) then
                    f(l1,:,:,j) = fbcx(j,k)
                  else
                    f(l2,:,:,j) = fbcx(j,k)
                  endif
                  call set_ghosts_for_onesided_ders(f,topbot,j,1,.true.)
                case ('n1s')
                  ! BCX_DOC: onesided for 1st/2nd derivative in two first inner points, Neumann in boundary point
                  call bval_from_neumann(f,topbot,j,1,fbcx(j,k))
                  call set_ghosts_for_onesided_ders(f,topbot,j,1,.true.)
                case ('1so')
                  ! BCX_DOC: onesided
                  call bc_onesided_x_old(f,topbot,j)
                case ('cT')
                  ! BCX_DOC: constant temperature (implemented as
                  ! BCX_DOC: condition for entropy $s$ or temperature $T$)
                  call bc_ss_temp_x(f,topbot)
                case ('c1')
                  ! BCX_DOC: constant conductive flux
                  if (j==iss)   call bc_ss_flux_x(f,topbot)
                  if (j==ilnTT) call bc_lnTT_flux_x(f,topbot)
                case ('Fgs')
                  ! BCX_DOC: Fconv = - chi_t*rho*T*grad(s)
                  if (j==iss) call bc_ss_flux_turb_x(f,topbot)
                case ('Fct')
                  ! BCX_DOC: Fbot = - K*grad(T) - chi_t*rho*T*grad(s)
                  if (j==iss) call bc_ss_flux_condturb_x(f,topbot)
                case ('Fcm')
                  ! BCX_DOC: $Fbot = - K*grad(\overline{T})$
                  ! BCX_DOC: $       - chi_t*\overline{rho}*\overline{T}*grad(\overline{s})$
                  if (j==iss) call bc_ss_flux_condturb_mean_x(f,topbot)
                case ('sT')
                  ! BCX_DOC: symmetric temperature, $T_{N-i}=T_{N+i}$;
                  ! BCX_DOC: implies $T'(x_N)=T'''(x_0)=0$
                  if (j==iss) call bc_ss_stemp_x(f,topbot)
                case ('asT')
                  ! BCX_DOC: select entropy for uniform ghost temperature
                  ! BCX_DOC: matching fluctuating boundary value,
                  ! BCX_DOC: $T_{N-i}=T_{N}=$;
                  ! BCX_DOC: implies $T'(x_N)=T'(x_0)=0$
                  if (j==iss) call bc_ss_a2stemp_x(f,topbot)
                case ('db')
                  ! BCX_DOC:
                  call bc_db_x(f,topbot,j)
                case ('f')
                  ! BCX_DOC: ``freeze'' value, i.e. maintain initial value; antisymm wrt boundary
                  call bc_freeze_var_x(topbot,j)
                  call bc_sym_x(f,-1,topbot,j,REL=.true.)
                case ('fg')
                  ! BCX_DOC: ``freeze'' value, i.e. maintain initial
                  !  value at boundary, also mantaining the
                  !  ghost zones at the initial coded value, i.e.,
                  !  keep the gradient frozen as well
                  call bc_freeze_var_x(topbot,j)
                case ('1')
                  ! BCX_DOC: $f=1$ (for debugging)
                  call bc_one_x(f,topbot,j)
                case ('set')
                  ! BCX_DOC: set boundary value to \var{fbcx}
                  call bc_sym_x(f,-1,topbot,j,REL=.true.,val=fbcx(:,k))
                case ('der')
                  ! BCX_DOC: set derivative on boundary to \var{fbcx}
                  call bc_set_der_x(f,topbot,j,fbcx(j,k))
                case ('slo')
                  ! BCX_DOC: set slope at the boundary = \var{fbcx}
                  call bc_slope_x(f,fbcx(:,k),topbot,j)
                case ('slp')
                  ! BCX_DOC: set slope at the boundary and in ghost cells = \var{fbcx}
                  call bc_ghost_slope_x(f,fbcx(:,k),topbot,j)
                case ('shx')
                  ! BCX_DOC: set shearing boundary proportional to x with slope=\var{fbcx} and abscissa=\var{fbcx2}
                  call bc_shear_x(f,fbcx(:,k),fbcx_2(:,k),topbot,j)
                case ('shy')
                  ! BCX_DOC: set shearing boundary proportional to y with slope=\var{fbcx} and abscissa=\var{fbcx2}
                  call bc_shear_y(f,fbcx(:,k),fbcx_2(:,k),topbot,j)
                case ('shz')
                  ! BCX_DOC: set shearing boundary proportional to z with slope=\var{fbcx} and abscissa=\var{fbcx2}
                  call bc_shear_z(f,fbcx(:,k),fbcx_2(:,k),topbot,j)
                case ('dr0')
                  ! BCX_DOC: set boundary value [really??]
                  call bc_dr0_x(f,fbcx(:,k),topbot,j)
                case ('ovr')
                  ! BCX_DOC: overshoot boundary condition
                  ! BCX_DOC:  ie $(d/dx-1/\mathrm{dist}) f = 0.$
                  call bc_overshoot_x(f,fbcx(:,k),topbot,j)
                case ('out')
                  ! BCX_DOC: allow outflow, but no inflow
                  ! BCX_DOC: forces ghost cells and boundary to not point inwards
                  call bc_outflow_x(f,topbot,j,.true.)
                case ('e1o')
                  ! BCX_DOC: allow outflow, but no inflow
                  ! BCX_DOC: uses the e1 extrapolation scheme
                  call bc_outflow_x_e1(f,topbot,j,.true.)
                case ('ant')
                  ! BCX_DOC: stops and prompts for adding documentation
                  call bc_antis_x(f,fbcx(:,k),topbot,j)
                case ('e1')
                  ! BCX_DOC: extrapolation [describe]
                  call bcx_extrap_2_1(f,topbot,j)
                case ('e2')
                  ! BCX_DOC: extrapolation [describe]
                  call bcx_extrap_2_2(f,topbot,j)
                case ('e3')
                  ! BCX_DOC: extrapolation in log [maintain a power law]
                  call bcx_extrap_2_3(f,topbot,j)
                case ('el')
                  ! BCX_DOC: linear extrapolation from last two active cells
                  call bcx_extrap_linear(f, topbot, j)
                case ('hat')
                  ! BCX_DOC: top hat jet profile in spherical coordinate.
                  !Defined only for the bottom boundary
                  call bc_set_jethat_x(f,j,topbot,fbcx(:,k),fbcx_2(:,k))
                case ('jet')
                  ! BCX_DOC: top hat jet profile in cartezian coordinate.
                  !Defined only for the bottom boundary
                  call bc_set_jet_x(f,j,topbot,fbcx(:,k),fbcx_2(:,k))
                case ('spd')
                  ! BCX_DOC:  sets $d(rA_{\alpha})/dr = \mathtt{fbcx(j)}$
                  call bc_set_spder_x(f,topbot,j,fbcx(j,k))
                case ('sfr')
                  ! BCX_DOC: stress-free boundary condition
                  ! BCX_DOC: for spherical coordinate system.
                  if (j==iux) call fatal_error('boundconds_x', &
                             'stress-free BC at r boundary not allowed for uu_r')
                  call bc_set_sfree_x(f,topbot,j)
                case ('sr1')
                  ! BCX_DOC: Stress-free bc for spherical coordinate system.
                  ! BCX_DOC: Implementation with one-sided derivative.
                  call bc_set_sr1_x(f,topbot,j)
                case ('nfr')
                  ! BCX_DOC: Normal-field bc for spherical coordinate system.
                  ! BCX_DOC: Some people call this the ``(angry) hedgehog bc''.
                  call bc_set_nfr_x(f,topbot,j)
                case ('nr1')
                  ! BCX_DOC: Normal-field bc for spherical coordinate system.
                  ! BCX_DOC: Some people call this the ``(angry) hedgehog bc''.
                  ! BCX_DOC: Implementation with one-sided derivative.
                  call bc_set_nr1_x(f,topbot,j)
                case ('sa2')
                  ! BCX_DOC: $(d/dr)(r B_{\phi}) = 0$ imposes
                  ! BCX_DOC: boundary condition on 2nd derivative of
                  ! BCX_DOC: $r A_{\phi}$. Same applies to $\theta$ component.
                  call bc_set_sa2_x(f,topbot,j)
                case ('pfc')
                  ! BCX_DOC: perfect-conductor in spherical
                  ! BCX_DOC: coordinate: $d/dr( A_r) + 2/r = 0$.
!joern: WARNING, this bc will NOT give a perfect-conductor boundary condition
                  call bc_set_pfc_x(f,topbot,j)
                case ('fix')
                  ! BCX_DOC: set boundary value [really??]
                  call bc_fix_x(f,topbot,j,fbcx(j,k))
                case ('fil')
                  ! BCX_DOC: set boundary value from a file
                  call bc_file_x(f,topbot,j)
                case ('cfb')
                  ! BCX_DOC: radial centrifugal balance
                  if (lcylindrical_coords) then
                    call bc_lnrho_cfb_r_iso(f,topbot)
                  else
                    print*,'not implemented for other than cylindrical'
                    stop
                  endif
                case ('g')
                  ! BCX_DOC: set to given value(s) or function
                  call bc_force_x(f, -1, topbot, j)
                case ('nil')
                  ! BCX_DOC: do nothing; assume that everything is set
                case ('ioc')
                  ! BCX_DOC: inlet/outlet on western/eastern hemisphere
                  ! BCX_DOC: in cylindrical coordinates
                  call bc_inlet_outlet_cyl(f,topbot,j,fbcx(:,k))
                case ('tay')
                  call tayler_expansion(f,topbot,j)
                case ('')
                  ! BCX_DOC: do nothing; assume that everything is set
                case ('slc')
                  call set_from_slice_x(f,topbot,j)
                  call set_ghosts_for_onesided_ders(f,topbot,j,1,.true.)
                case default
                  bc%bcname=bcx12(j,k)
                  bc%ivar=j
                  bc%location=(((k-1)*2)-1)   ! -1/1 for x bot/top
                  bc%value1=fbcx(j,k)
                  bc%value2=fbcx(j,k)
                  bc%done=.false.
!
                  call special_boundconds(f,bc)
!
                  if (.not.bc%done) then
                    write(unit=errormsg,fmt='(A,A4,A,I3)') &
                         "No such boundary condition bcx1/2 = ", &
                         bcx12(j,k), " for j=", j
                    call fatal_error_local("boundconds_x",trim(errormsg))
                  endif
                endselect
              endif
            enddo
          enddo
        endif
      endselect
!
    endsubroutine boundconds_x
!***********************************************************************
    subroutine boundconds_y_c(f,ivar1_opt,ivar2_opt)
!
!  Envelope for being called from C code.
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, optional :: ivar1_opt, ivar2_opt

      call boundconds_y(f,ivar1_opt,ivar2_opt)

    endsubroutine boundconds_y_c
!***********************************************************************
    subroutine boundconds_y(f,ivar1_opt,ivar2_opt)
!
!  Boundary conditions in y, except for periodic part handled by communication.
!  Remark: boundconds_x() needs to be called before communicating (because we
!  communicate the x-ghost points), boundconds_[yz] after communication
!  has finished (they need some of the data communicated for the edges
!  (yz-'corners').
!
!   8-jul-02/axel: split up into different routines for x,y and z directions
!  11-nov-02/wolf: unified bot/top, now handled by loop
!
      use General, only: var_is_vec
      use Special, only: special_boundconds
      use EquationOfState
!
      real, dimension (:,:,:,:) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      integer :: ivar1, ivar2, j, k
      logical :: ip_ok
      character (len=bclen) :: topbot
      type (boundary_condition) :: bc
!
      if (ldebug) print*,'boundconds_y: ENTER'
!
      ivar1=1; ivar2=min(mcom,size(f,4))
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
!
      select case (nygrid)
!
      case (1)
        if (ldebug) print*,'boundconds_y: no y-boundary'
!
!  Boundary conditions in y
!
      case default
        do k=1,2                ! loop over 'bot','top'
          if (k==1) then
            topbot='bot'; 
            ip_ok=lfirst_proc_y
          else
            topbot='top'; 
            ip_ok=llast_proc_y
          endif
!
          jdone=0
          do j=ivar1,ivar2
!
! Natalia: the next line is for the dustdensity case.
! If ndustspec is large, it is stupid to set bc for all dust species
! in start.in. But if one does not set them, they becomes 'p' by default
! Since this problem is crutial  only for aerosol + chemistry
! the following condition is used. But this place should be modifyed somehow
! Any ideas?
!
!            if ((bcy12(j,k)=='p') .and. lchemistry .and. ldustdensity) bcy12(j,k)=''
!
            if (ldebug) write(*,'(A,I1,A,I2,A,A)') ' bcy',k,'(',j,')=',bcy12(j,k)

            if (ip_ok) then
           
              is_vec = var_is_vec(j)
              select case (bcy12(j,k))
              case ('0')
                ! BCY_DOC: zero value in ghost zones, free value on boundary
                call bc_zero_y(f,topbot,j)
              case ('p')
                ! BCY_DOC: periodic
                call bc_per_y(f,topbot,j)
              case ('pp')
                ! BCY_DOC: periodic across the pole
                call bc_pper_y(f,+1,topbot,j)
              case ('yy')
                ! BCY_DOC: Yin-Yang grid
                call bc_yy_y(f,topbot,j)
              case ('ap')
                ! BCY_DOC: anti-periodic across the pole
                call bc_pper_y(f,-1,topbot,j)
              case ('s')
                ! BCY_DOC: symmetry symmetry, $f_{N+i}=f_{N-i}$;
                  ! BCX_DOC: implies $f'(y_N)=f'''(y_0)=0$
                call bc_sym_y(f,+1,topbot,j)
              case ('sf')
                ! BCY_DOC: symmetry with respect to interface
                call bc_sf_y(f,+1,topbot,j)
              case ('ss')
                ! BCY_DOC: symmetry, plus function value given
                call bc_symset_y(f,+1,topbot,j,val=fbcy(:,k))
              case ('sds')
                ! BCY_DOC: symmetric-derivative-set
                call bc_symderset_y(f,topbot,j,val=fbcy(:,k))
              case ('cds')
                ! BCY_DOC: complex symmetric-derivative-set
                call bc_csymderset_y(f,topbot,j,val=fbcy(:,k))
              case ('s0d')
                ! BCY_DOC: symmetry, function value such that df/dy=0
                call bc_symset0der_y(f,topbot,j)
              case ('a')
                ! BCY_DOC: antisymmetry
                call bc_sym_y(f,-1,topbot,j)
              case ('af')
                ! BCY_DOC: antisymmetry with respect to interface
                call bc_sf_y(f,-1,topbot,j)
              case ('a2')
                ! BCY_DOC: antisymmetry relative to boundary value
                call bc_sym_y(f,-1,topbot,j,REL=.true.)
              case ('v')
                ! BCY_DOC: vanishing third derivative
                call bc_van_y(f,topbot,j)
              case ('v3')
                ! BCY_DOC: vanishing third derivative
                call bc_van3rd_y(f,topbot,j)
              case ('out')
                ! BCY_DOC: allow outflow, but no inflow
                ! BCY_DOC: forces ghost cells and boundary to not point inwards
                call bc_outflow_y(f,topbot,j,.true.)
              case ('1s')
                ! BCY_DOC: onesided
                call set_ghosts_for_onesided_ders(f,topbot,j,2)
              case ('d1s')
                ! BCY_DOC: onesided for 1st and 2nd derivative in two first inner points, Dirichlet in boundary point
                if (k==1) then
                  f(:,m1,:,j) = fbcy(j,k)
                else
                  f(:,m2,:,j) = fbcy(j,k)
                endif
                call set_ghosts_for_onesided_ders(f,topbot,j,2,.true.)
              case ('n1s')
                ! BCY_DOC: onesided for 1st and 2nd derivative in two first inner points, Neumann in boundary point
                call bval_from_neumann(f,topbot,j,2,fbcy(j,k))
                call set_ghosts_for_onesided_ders(f,topbot,j,2,.true.)
              case ('cT')
                ! BCY_DOC: constant temp.
                if (j==iss) call bc_ss_temp_y(f,topbot)
              case ('sT')
                ! BCY_DOC: symmetric temp.
                if (j==iss) call bc_ss_stemp_y(f,topbot)
              case ('asT')
                ! BCY_DOC: select entropy for uniform ghost temperature
                ! BCY_DOC: matching fluctuating boundary value,
                ! BCY_DOC: $T_{N-i}=T_{N}=$;
                ! BCY_DOC: implies $T'(x_N)=T'(x_0)=0$
                if (j==iss) call bc_ss_a2stemp_y(f,topbot)
              case ('f')
                ! BCY_DOC: freeze value
                ! tell other modules not to change boundary value
                call bc_freeze_var_y(topbot,j)
                call bc_sym_y(f,-1,topbot,j,REL=.true.) ! antisymm wrt boundary
              case ('s+f')
                ! BCY_DOC: freeze value
                ! tell other modules not to change boundary value
                call bc_freeze_var_y(topbot,j)
                call bc_sym_y(f,+1,topbot,j) ! symm wrt boundary
              case ('fg')
                ! BCY_DOC: ``freeze'' value, i.e. maintain initial
                !  value at boundary, also mantaining the
                !  ghost zones at the initial coded value, i.e.,
                !  keep the gradient frozen as well
                call bc_freeze_var_y(topbot,j)
              case ('fBs')
                ! BCY_DOC: frozen-in B-field (s)
                call bc_frozen_in_bb(topbot,j)
                call bc_sym_y(f,+1,topbot,j) ! symmetry
              case ('fB')
                ! BCY_DOC: frozen-in B-field (a2)
                call bc_frozen_in_bb(topbot,j)
                !call bc_sym_z(f,-1,topbot,j,REL=.true.) ! antisymm wrt boundary
!AB: wasn't this a mistake??
                call bc_sym_y(f,-1,topbot,j,REL=.true.) ! antisymm wrt boundary
              case ('1')
                ! BCY_DOC: f=1 (for debugging)
                call bc_one_y(f,topbot,j)
              case ('set')
                ! BCY_DOC: set boundary value
                call bc_sym_y(f,-1,topbot,j,REL=.true.,val=fbcy(:,k))
              case ('sse')
                ! BCY_DOC:  symmetry, set boundary value
                call bc_sym_y(f,+1,topbot,j,val=fbcy(:,k))
              case ('sep')
                ! BCY_DOC: set boundary value
                call bc_sym_y(f,-1,topbot,j,REL=.true.,val=fbcy(:,k),val2=fbcy_1(:,k),val4=fbcy_2(:,k))
              case ('e1')
                ! BCY_DOC: extrapolation
                call bcy_extrap_2_1(f,topbot,j)
              case ('e2')
                ! BCY_DOC: extrapolation
                call bcy_extrap_2_2(f,topbot,j)
              case ('e3')
                ! BCY_DOC: extrapolation in log [maintain a power law]
                call bcy_extrap_2_3(f,topbot,j)
              case ('der')
                ! BCY_DOC: set derivative on the boundary
                call bc_set_der_y(f,topbot,j,fbcy(j,k))
              case ('cop')
                ! BCY_DOC: outflow: copy value of last physical point to
                ! BCY_DOC: all ghost cells
                call bc_copy_y(f,topbot,j)
              case ('c+k')
                ! BCY_DOC: no-inflow: copy value of last physical point
                ! BCY_DOC: to all ghost cells, but suppressing any inflow
                call bc_copy_y_noinflow(f,topbot,j)
              case ('sfr')
                ! BCY_DOC: stress-free boundary condition for spherical
                ! BCY_DOC: coordinate system.
                  if (j==iux.or.j==iuy) call fatal_error('boundconds_y', &
                             'stress-free BC at theta boundary only allowed for uu_phi')
                call bc_set_sfree_y(f,topbot,j)
              case ('nfr')
                ! BCY_DOC: Normal-field bc for spherical coordinate system.
                ! BCY_DOC: Some people call this the ``(angry) hedgehog bc''.
                call bc_set_nfr_y(f,topbot,j)
              case ('spt')
                ! BCY_DOC: spherical perfect conducting boundary condition
                ! BCY_DOC: along $\theta$ boundary
                ! BCY_DOC: $f''+\cot\theta f'=0$ and $f(x_N)=0$
                call bc_spt_y(f,topbot,j)
              case ('pfc')
                ! BCY_DOC: perfect conducting boundary condition
                ! BCY_DOC: along $\theta$ boundary
!joern: WARNING, this bc will NOT give a perfect-conductor boundary condition
                call bc_set_pfc_y(f,topbot,j)
              case ('str')
                call bc_stratified_y(f,topbot,j)
              case ('nil','')
                ! BCY_DOC: do nothing; assume that everything is set
              case ('slc')
                call set_from_slice_y(f,topbot,j)
                call set_ghosts_for_onesided_ders(f,topbot,j,2,.true.)
              case default
                bc%bcname=bcy12(j,k)
                bc%ivar=j
                bc%value1=fbcy(j,k)
                bc%value2=fbcy(j,k)
                bc%location=(((k-1)*4)-2)   ! -2/2 for y bot/top
                bc%done=.false.
!
                if (lspecial) call special_boundconds(f,bc)
!
                if (.not.bc%done) then
                  write(unit=errormsg,fmt='(A,A4,A,I3)') "No such boundary condition bcy1/2 = ", &
                       bcy12(j,k), " for j=", j
                  call fatal_error_local("boundconds_y",trim(errormsg))
                endif
              endselect
            endif
          enddo
        enddo
      endselect
!
    endsubroutine boundconds_y
!***********************************************************************
    subroutine boundconds_z_c(f,ivar1_opt,ivar2_opt)
!
!  Envelope for being called from C code.
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, optional :: ivar1_opt, ivar2_opt

      call boundconds_z(f,ivar1_opt,ivar2_opt)

    endsubroutine boundconds_z_c
!***********************************************************************
    subroutine boundconds_z(f,ivar1_opt,ivar2_opt)
!
!  Boundary conditions in z, except for periodic part handled by communication.
!  Remark: boundconds_x() needs to be called before communicating (because we
!  communicate the x-ghost points), boundconds_[yz] after communication
!  has finished (they need some of the data communicated for the edges
!  (yz-'corners').
!
!   8-jul-02/axel: split up into different routines for x,y and z directions
!  11-nov-02/wolf: unified bot/top, now handled by loop
!  02-apr-13/MR  : added new boundary condition 'fs' = frozen boundary value
!                  + symmetry about boundary; added 'fa' for alternative reference to
!                  already existing freezing condition (includes antisymmetry)
!  30-dec-16/MR: added BC 'a1s' for constant alpha mean-field model in one dimension
!
      use General, only: var_is_vec
      use Gravity, only: gravz_profile
      use Special, only: special_boundconds
      use EquationOfState
      use Magnetic_meanfield, only: pc_aasb_const_alpha
!
      real, dimension (:,:,:,:) :: f
      integer, optional :: ivar1_opt, ivar2_opt
      real, dimension (size(f,4)) :: fbcz_zero
      integer :: ivar1, ivar2, j, k
      logical :: ip_ok
      character (len=bclen) :: topbot
      type (boundary_condition) :: bc
!
      if (ldebug) print*,'boundconds_z: ENTER'
!
      ivar1=1; ivar2=min(mcom,size(f,4))
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
!
      select case (nzgrid)
!
      case (1)
        if (ldebug) print*,'boundconds_z: no z-boundary'
!
!  Boundary conditions in z
!
      case default
        do k=1,2                ! loop over 'bot','top'
          if (k==1) then
            topbot='bot'
            ip_ok=lfirst_proc_z
          else
            topbot='top'
            ip_ok=llast_proc_z
          endif
!
          jdone=0
          do j=ivar1,ivar2
            if (ldebug) write(*,'(A,I1,A,I2,A,A)') ' bcz',k,'(',j,')=',bcz12(j,k)
            if (ip_ok) then

              is_vec = var_is_vec(j)

              select case (bcz12(j,k))
              case ('0')
                ! BCZ_DOC: zero value in ghost zones, free value on boundary
                call bc_zero_z(f,topbot,j)
              case ('p')
                ! BCZ_DOC: periodic
                call bc_per_z(f,topbot,j)
              case ('yy')
                ! BCZ_DOC: Yin-Yang grid
                call bc_yy_z(f,topbot,j)
              case ('s')
                ! BCZ_DOC: symmetry
                call bc_sym_z(f,+1,topbot,j)
              case ('sf')
                ! BCZ_DOC: symmetry with respect to interface
                call bc_sf_z(f,+1,topbot,j)
              case ('s0d')
                ! BCZ_DOC: symmetry, function value such that df/dz=0
                call bc_symset0der_z(f,topbot,j)
              case ('0ds')
                ! BCZ_DOC: symmetry, function value such that df/dz=0
                call bc_symset0der_z_v2(f,topbot,j)
              case ('a')
                ! BCZ_DOC: antisymmetry
                call bc_sym_z(f,-1,topbot,j)
              case ('a2')
                ! BCZ_DOC: antisymmetry relative to boundary value
                call bc_sym_z(f,-1,topbot,j,REL=.true.)
               case ('a2v')
                  ! BCZ_DOC: set boundary value and antisymmetry relative to it
                  call bc_sym_z(f,-1,topbot,j,REL=.true.,val=fbcz(:,k))
              case ('af')
                ! BCZ_DOC: antisymmetry with respect to interface
                call bc_sf_z(f,-1,topbot,j)
              case ('a0d')
                ! BCZ_DOC: antisymmetry with zero derivative
                fbcz_zero=0.
                call bc_sym_z(f,+1,topbot,j,VAL=fbcz_zero)
              case ('v')
                ! BCZ_DOC: vanishing third derivative
                call bc_van_z(f,topbot,j)
              case ('v3')
                ! BCZ_DOC: vanishing third derivative
                call bc_van3rd_z(f,topbot,j)
              case ('1s')
                ! BCZ_DOC: one-sided
                call set_ghosts_for_onesided_ders(f,topbot,j,3)
              case ('d1s')
                ! BCZ_DOC: onesided for 1st and 2nd derivative in two first inner points, Dirichlet in boundary point
                if (k==1) then
                  f(:,:,n1,j) = fbcz(j,k)
                else
                  f(:,:,n2,j) = fbcz(j,k)
                endif
                call set_ghosts_for_onesided_ders(f,topbot,j,3,.true.)
              case ('n1s')
                ! BCZ_DOC: onesided for 1st and 2nd derivative in two first inner points, Neumann in boundary point
                call bval_from_neumann(f,topbot,j,3,fbcz(j,k))
                call set_ghosts_for_onesided_ders(f,topbot,j,3,.true.)
              case ('a1s')
                ! BCZ_DOC: special for perfect conductor with const alpha and etaT when A considered as B; one-sided for 1st and 2nd derivative in two first inner points
                call pc_aasb_const_alpha(f,topbot,j)
                call set_ghosts_for_onesided_ders(f,topbot,j,3,.true.)
              case ('fg')
                ! BCZ_DOC: ``freeze'' value, i.e. maintain initial
                !  value at boundary, also mantaining the
                !  ghost zones at the initial coded value, i.e.,
                !  keep the gradient frozen as well
                call bc_freeze_var_z(topbot,j)
              case ('c1')
                ! BCZ_DOC: complex
                if (j==iss) call bc_ss_flux(f,topbot)
                if (j==iaa) call bc_aa_pot(f,topbot)
                if (j==ilnTT) call bc_lnTT_flux_z(f,topbot)
              case ('c1s')
                ! BCZ_DOC: complex
                if (j==iss) call bc_ss_flux(f,topbot,.true.)
              case ('Fgs')
                ! BCZ_DOC: Fconv = - chi_t*rho*T*grad(s)
                if (j==iss) call bc_ss_flux_turb(f,topbot)
              case ('Fct')
                ! BCZ_DOC: Fbot = - K*grad(T) - chi_t*rho*T*grad(s)
                if (j==iss) call bc_ss_flux_condturb_z(f,topbot)
              case ('c3')
                ! BCZ_DOC: constant flux at the bottom with a variable hcond
                if (j==ilnTT) call bc_ADI_flux_z(f,topbot)
              case ('pfe')
                ! BCZ_DOC: potential field extrapolation
                if (j==iaa) call bc_aa_pot_field_extrapol(f,topbot)
              case ('p1D')
                ! BCZ_DOC: potential field extrapolation in 1D
                if (j==iay) call bc_aa_pot_1D(f,topbot)
              case ('pot')
                ! BCZ_DOC: potential magnetic field
                if (j==iaa) call bc_aa_pot2(f,topbot)
              case ('pwd')
                ! BCZ_DOC: a variant of 'pot' for nprocx=1
                if (j==iaa) call bc_aa_pot3(f,topbot)
              case ('d2z')
                ! BCZ_DOC:
                call bc_del2zero(f,topbot,j)
              case ('hds')
                ! BCZ_DOC: hydrostatic equilibrium with a high-frequency filter
                call bc_lnrho_hdss_z_iso(f,topbot)
              case ('cT')
                ! BCZ_DOC: constant temp.
                if (j==ilnrho) call bc_lnrho_temp_z(f,topbot)
                call bc_ss_temp_z(f,topbot)
              case ('cT1')
                ! BCZ_DOC: constant temp.
                call bc_ss_temp_z(f,topbot,.true.)
              case ('cT2')
                ! BCZ_DOC: constant temp. (keep lnrho)
                if (j==iss) call bc_ss_temp2_z(f,topbot)
              case ('cT3')
                ! BCZ_DOC: constant temp. (keep lnrho)
                if (j==iss) call bc_ss_temp3_z(f,topbot)
              case ('hs')
                ! BCZ_DOC: hydrostatic equilibrium
                if (.not. lgrav) call fatal_error('boundconds_z', &
                  'hs boundary condition requires gravity')
                if ((.not. ltemperature .or. ltemperature_nolog) .and. (gravz_profile /= 'const')) &
                  call fatal_error('boundconds_z', 'hs boundary condition requires a constant gravity profile')
                if (.not. lequidist(3)) call fatal_error('boundconds_z', &
                  'hs boundary condition requires symmetric grid distances on the z boundary')
                if ((j==ilnrho) .or. (j==irho_b) .or. (j==iss)) then
                  call bc_lnrho_hds_z_iso(f,topbot)
                elseif (j==ipp) then
                  call bc_pp_hds_z_iso(f,topbot)
                else
                  call fatal_error ('boundconds_z', "hs boundary condition requires density or pressure")
                endif
              case ('hse')
                ! BCZ_DOC: hydrostatic extrapolation
                ! BCZ_DOC: rho or lnrho is extrapolated linearily and the
                ! BCZ_DOC: temperature is calculated in hydrostatic equilibrium.
                if (.not. lgrav) &
                    call fatal_error ('boundconds_z', "'hse' requires gravity")
                if (.not. leos) call fatal_error ('boundconds_z', &
                    "'hse' requires an eos module")
                if ((ilnrho == 0) .or. (ilnTT == 0)) &
                    call fatal_error ('boundconds_z', "'hse' requires lnrho and lnTT")
                if (j == ilnTT) then
                  call bcz_hydrostatic_temp(f,topbot)
                elseif (j == ilnrho) then
                  call bcz_hydrostatic_rho(f,topbot)
                else
                  call fatal_error ('boundconds_z', "'hse' works only in lnrho or lnTT")
                endif
              case ('cp')
                ! BCZ_DOC: constant pressure
                ! BCZ_DOC:
                if (j==ilnrho) call bc_lnrho_pressure_z(f,topbot)
              case ('sT')
                ! BCZ_DOC: symmetric temp.
                ! BCZ_DOC:
                if (j==iss) call bc_ss_stemp_z(f,topbot)
              case ('ctz')
                ! BCZ_DOC: for interstellar runs copy T
                if (j==iss) call bc_ctz(f,topbot,iss)
              case ('cdz')
                ! BCZ_DOC: for interstellar runs limit rho
                call bc_cdz(f,topbot,j)
              case ('ism')
                ! BCZ_DOC: for interstellar runs limit rho
                call bc_ism(f,topbot,j)
              case ('asT')
                ! BCZ_DOC: select entropy for uniform ghost temperature
                ! BCZ_DOC: matching fluctuating boundary value,
                ! BCZ_DOC: $T_{N-i}=T_{N}=$;
                ! BCZ_DOC: implies $T'(x_N)=T'(x_0)=0$
                if (j==iss) call bc_ss_a2stemp_z(f,topbot)
              case ('c2')
                ! BCZ_DOC: complex
                ! BCZ_DOC:
                if (j==iss) call bc_ss_temp_old(f,topbot)
              case ('db')
                ! BCZ_DOC: complex
                ! BCZ_DOC:
                call bc_db_z(f,topbot,j)
              case ('ce')
                ! BCZ_DOC: complex
                ! BCZ_DOC:
                if (j==iss) call bc_ss_energy(f,topbot)
              case ('e1')
                ! BCZ_DOC: extrapolation
                call bc_extrap_2_1(f,topbot,j)
              case ('e2')
                ! BCZ_DOC: extrapolation
                call bc_extrap_2_2(f,topbot,j)
              case ('ex')
                ! BCZ_DOC: simple linear extrapolation in first order
                call bcz_extrapol(f,topbot,j)
              case ('exf')
                ! BCZ_DOC: simple linear extrapolation in first order
                !  with a fixed value in the first ghost cell
                call bcz_extrapol_fixed(f,topbot,j)
              case ('exd')
                ! BCZ_DOC: simple linear extrapolation in first order
                !  with an included damping to zero (useful for velocities)
                call bcz_extrapol_damped(f,topbot,j)
              case ('exm')
                ! BCZ_DOC: simple linear extrapolation in first order
                !  with an included local averaging of a 7x7 array
                call bcz_extrapol_mean(f,topbot,j)
              case ('b1')
                ! BCZ_DOC: extrapolation with zero value (improved 'a')
                call bc_extrap0_2_0(f,topbot,j)
              case ('b2')
                ! BCZ_DOC: extrapolation with zero value (improved 'a')
                call bc_extrap0_2_1(f,topbot,j)
              case ('b3')
                ! BCZ_DOC: extrapolation with zero value (improved 'a')
                call bc_extrap0_2_2(f,topbot,j)
              case ('f','fa')
                ! BCZ_DOC: freeze value + antisymmetry
                ! tell other modules not to change boundary value
                call bc_freeze_var_z(topbot,j)
                call bc_sym_z(f,-1,topbot,j,REL=.true.) ! antisymm wrt boundary
              case ('fs')
                ! BCZ_DOC: freeze value + symmetry
                ! tell other modules not to change boundary value
                call bc_freeze_var_z(topbot,j)
                call bc_sym_z(f,+1,topbot,j) ! symmetric wrt boundary
              case ('fBs')
                ! BCZ_DOC: frozen-in B-field (s)
                call bc_frozen_in_bb(topbot,j)
                call bc_sym_z(f,+1,topbot,j) ! symmetry
              case ('fB')
                ! BCZ_DOC: frozen-in B-field (a2)
                call bc_frozen_in_bb(topbot,j)
                call bc_sym_z(f,-1,topbot,j,REL=.true.) ! antisymm wrt boundary
              case ('g')
                ! BCZ_DOC: set to given value(s) or function
                 call bc_force_z(f,-1,topbot,j)
              case ('gs')
                ! BCZ_DOC:
                 call bc_force_z(f,+1,topbot,j)
              case ('1')
                ! BCZ_DOC: f=1 (for debugging)
                call bc_one_z(f,topbot,j)
              case ('StS')
                ! BCZ_DOC: solar surface boundary conditions
                if (j==ilnrho) call bc_stellar_surface(f,topbot)
              case ('set')
                ! BCZ_DOC: set boundary value
                call bc_sym_z(f,-1,topbot,j,REL=.true.,val=fbcz(:,k))
              case ('sep')
                ! BCY_DOC: set boundary value
                call bc_sym_z(f,-1,topbot,j,REL=.true.,val=fbcz(:,k),val2=fbcz_1(:,k),val4=fbcz_2(:,k))
              case ('der')
                ! BCZ_DOC: set derivative on the boundary
                call bc_set_der_z(f,topbot,j,fbcz(j,k))
              case ('div')
                ! BCZ_DOC: set the divergence of $\uv$ to a given value
                ! BCZ_DOC: use bc = 'div' for iuz
                call bc_set_div_z(f,topbot,j,fbcz(j,k))
              case ('ovr')
                ! BCZ_DOC: set boundary value
                call bc_overshoot_z(f,fbcz(:,k),topbot,j)
              case ('inf')
                ! BCZ_DOC: allow inflow, but no outflow
                call bc_inflow_z(f,topbot,j)
              case ('ouf')
                ! BCZ_DOC: allow outflow, but no inflow
                call bc_outflow_z(f,topbot,j)
              case ('in')
                ! BCZ_DOC: allow inflow, but no outflow
                ! BCZ_DOC: forces ghost cells and boundary to not point outwards
                call bc_inflow_z(f,topbot,j,.true.)
              case ('out')
                ! BCZ_DOC: allow outflow, but no inflow
                ! BCZ_DOC: forces ghost cells and boundary to not point inwards
                call bc_outflow_z(f,topbot,j,.true.)
              case ('crk')
                ! BCY_DOC: no-inflow: copy value of last physical point
                ! BCY_DOC: to all ghost cells, but suppressing any inflow
                call bc_copy_z_noinflow(f,topbot,j)
              case ('in0')
                ! BCZ_DOC: allow inflow, but no outflow
                ! BCZ_DOC: forces ghost cells and boundary to not point outwards
                ! BCZ_DOC: relaxes to vanishing 1st derivative at boundary
                call bc_inflow_zero_deriv_z(f,topbot,j)
              case ('ou0')
                ! BCZ_DOC: allow outflow, but no inflow
                ! BCZ_DOC: forces ghost cells and boundary to not point inwards
                ! BCZ_DOC: relaxes to vanishing 1st derivative at boundary
                call bc_outflow_zero_deriv_z(f,topbot,j)
              case ('ind')
                ! BCZ_DOC: allow inflow, but no outflow
                ! BCZ_DOC: forces ghost cells and boundary to not point outwards
                ! BCZ_DOC: creates inwards pointing or zero 1st derivative at boundary
                call bc_inflow_inwards_deriv_z(f,topbot,j)
              case ('oud')
                ! BCZ_DOC: allow outflow, but no inflow
                ! BCZ_DOC: forces ghost cells and boundary to not point inwards
                ! BCZ_DOC: creates outwards pointing or zero 1st derivative at boundary
                call bc_outflow_outwards_deriv_z(f,topbot,j)
              case ('ubs')
                ! BCZ_DOC: copy boundary outflow,
                ! but limit inflow +ve inward gradient (experimental)
                call bc_steady_z(f,topbot,j)
              case ('win')
                ! BCZ_DOC: forces massflux given as
                ! BCZ_DOC: $\Sigma \rho_i ( u_i + u_0)=\textrm{fbcz1/2}(\rho)$
                if (j==ilnrho) then
                   call bc_wind_z(f,topbot,fbcz(j,k))
                   call bc_sym_z(f,+1,topbot,j)           !  's'
                   call bc_sym_z(f,+1,topbot,iuz)         !  's'
                endif
              case ('cop')
                ! BCZ_DOC: copy value of last physical point to all ghost cells
                call bc_copy_z(f,topbot,j)
              case ('nil')
                ! BCZ_DOC: do nothing; assume that everything is set
              case ('slc')
                call set_from_slice_z(f,topbot,j)
                !call set_ghosts_for_onesided_ders(f,topbot,j,3,.true.)
                call bc_sym_z(f,-1,topbot,j,rel=.true.)
              case default
                bc%bcname=bcz12(j,k)
                bc%ivar=j
                bc%location=(((k-1)*6)-3)   ! -3/3 for z bot/top
                bc%value1=fbcz_1(j,k)
                bc%value2=fbcz_2(j,k)
                bc%done=.false.
!
                if (lspecial) call special_boundconds(f,bc)
!
                if (.not.bc%done) then
                  write(unit=errormsg,fmt='(A,A4,A,I3)') "No such boundary condition bcz1/2 = ", &
                       bcz12(j,k), " for j=", j
                  call fatal_error_local("boundconds_z",trim(errormsg))
                endif
              endselect
            endif
          enddo
        enddo
      endselect
!
    endsubroutine boundconds_z
!***********************************************************************
    subroutine bc_pencil_scalar(penc, ncell, nghost, bc, d2_bound, bound)
!
! Apply boundary conditions to a 1D scalar of arbitrary size.
!
! 29-may-12/ccyang: coded
!  2-apr-15/MR: optional parameters d2_bound, bound for use in stress-free 
!               and normal-field BCs added; these BCs, 'a' and 's' implemented
!
! Input/Output Arguments
!   penc - a scalar array to be applied boundary conditions
!
! Input Arguments
!   ncell - number of active cells in penc
!   nghost - number of ghost cells in penc
!   bc1 - boundary condition for the lower boundary
!   bc2 - boundary condition for the upper boundary
!   d2_bound - doubled cumulative cell sizes at boundary
!              (d2_bound(-nghost:-1) - at lower, d2_bound(1:nghost) at upper)
!   bound - boundary coordinates
!
      integer, intent(in) :: ncell, nghost
      real, dimension(1-nghost:ncell+nghost), intent(inout) :: penc
      character(len=*), dimension(2), intent(in) :: bc
      real, dimension(-nghost:nghost), optional :: d2_bound
      real, dimension(2), optional :: bound

      integer :: i
!
! Apply lower boundary condition.
!
      lower: select case (bc(1))
!     Nothing
      case ('') lower
!     Periodic
      case ('p') lower
        penc(1-nghost:0) = penc(ncell-nghost+1:ncell)
!     Zero
      case ('0') lower
        penc(1-nghost:0) = 0.0
!     Zeroth-order extrapolation
      case ('cop') lower
        penc(1-nghost:0) = penc(1)
      case ('s') lower
        penc(1-nghost:0) = penc(2:nghost+1)
      case ('a') lower
        penc(1-nghost:0) = -penc(2:nghost+1)
        penc(1) = 0.
      case ('sfr') lower
        do i=1,nghost
          penc(1-i) = penc(1+i) - penc(1)*(d2_bound(-i)/bound(BOT))
        enddo
      case ('nfr') lower
        do i=1,nghost
          penc(1-i) = penc(1+i) + penc(1)*(d2_bound(-i)/bound(BOT))
        enddo
!     Unknown boundary condition
      case default lower
        call fatal_error('bc_pencil_scalar', 'unknown lower boundary condition')
      endselect lower
!
! Apply upper boundary condition.
!
      upper: select case (bc(2))
!     Nothing
      case ('') upper
!     Periodic
      case ('p') upper
        penc(ncell+1:ncell+nghost) = penc(1:nghost)
!     Zero
      case ('0') upper
        penc(ncell+1:ncell+nghost) = 0.0
!     Zeroth-order extrapolation
      case ('cop') upper
        penc(ncell+1:ncell+nghost) = penc(ncell)
      case ('s') upper
        penc(ncell+1:ncell+nghost) = penc(ncell-nghost:ncell-1) 
      case ('a') upper
        penc(ncell+1:ncell+nghost) = -penc(ncell-nghost:ncell-1)
        penc(ncell) = 0.
      case ('sfr') upper
        do i=1,nghost
          penc(1+i) = penc(1-i) + penc(ncell)*(d2_bound(i)/bound(TOP))
        enddo
      case ('nfr') upper
        do i=1,nghost
          penc(1+i) = penc(1-i) - penc(ncell)*(d2_bound(i)/bound(TOP))
        enddo
!     Unknown boundary condition
      case default upper
        call fatal_error('bc_pencil_scalar', 'unknown upper boundary condition')
      endselect upper
!
    endsubroutine bc_pencil_scalar
!***********************************************************************
    subroutine bc_pencil_vector(penc, ncell, nghost, ncomp, bc, d2_bound, bound)
!
! Apply boundary conditions to a 1D vector of arbitrary size.
!
! 22-may-12/ccyang: coded
!  2-apr-15/MR: optional parameters d2_bound, bound for use in stress-free and normal-field BCs added
!
! Input/Output Arguments
!   penc - a vector array to be applied boundary conditions
!
! Input Arguments
!   ncell - number of active cells in penc
!   nghost - number of ghost cells in penc
!   ncomp - number of components of the vector
!   bc1 - boundary condition for the lower boundary
!   bc2 - boundary condition for the upper boundary
!
      integer, intent(in) :: ncell, nghost, ncomp
      real, dimension(1-nghost:ncell+nghost, ncomp), intent(inout) :: penc
      character(len=bclen), dimension(ncomp,2), intent(in) :: bc
      real, dimension(-nghost:nghost), optional :: d2_bound
      real, dimension(2), optional :: bound
!
      integer :: j
!
      comp: do j = 1, ncomp
        call bc_pencil_scalar(penc(:,j), ncell, nghost, bc(j,:), d2_bound, bound)
      enddo comp
!
    endsubroutine bc_pencil_vector
!***********************************************************************
    subroutine bc_per_x(f,topbot,j)
!
!  Periodic boundary condition
!
!  11-nov-02/wolf: coded
!
      real, dimension (:,:,:,:) :: f
      integer :: j
      character (len=bclen) :: topbot
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        if (nprocx==1) f(1:l1-1,:,:,j) = f(l2i:l2,:,:,j)
!
      case ('top')               ! top boundary
        if (nprocx==1) f(l2+1:,:,:,j) = f(l1:l1i,:,:,j)
!
      case default
        print*, "bc_per_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_per_x
!***********************************************************************
    subroutine bc_per_y(f,topbot,j)
!
!  Periodic boundary condition
!
!  11-nov-02/wolf: coded
!
      real, dimension (:,:,:,:) :: f
      integer :: j
      character (len=bclen) :: topbot
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        if (nprocy==1) f(:,1:m1-1,:,j) = f(:,m2i:m2,:,j)
!
      case ('top')               ! top boundary
        if (nprocy==1) f(:,m2+1:,:,j) = f(:,m1:m1i,:,j)
!
      case default
        print*, "bc_per_y: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_per_y
!***********************************************************************
    subroutine bc_yy_y(f,topbot,j)
!
!  After-communication transformation of vector quantities for Yin-Yang grid.
!
!  30-nov-15/MR: coded
!
      use General, only: transform_cart_spher

      real, dimension (:,:,:,:) :: f
      integer :: j
      character (len=bclen) :: topbot
!
      if (.not.lyinyang) &
        call fatal_error_local('bc_yy_y','BC not legal as no Yin-Yang grid run.')

      if (j<=jdone) return

      if (is_vec) then
!
!  Vector quantities need to be transformed from the Cartesian basis to 
!  the local spherical basis.
!
        jdone=j+2     ! requires adjacent vector components
        if (topbot=='bot') then
          call transform_cart_spher(f,1,nghost,1,mz,j)    ! in-place!
        else
          call transform_cart_spher(f,m2+1,my,1,mz,j)     !  ~
        endif
!
      else
        jdone=0
      endif

    endsubroutine bc_yy_y
!***********************************************************************
    subroutine bc_pper_y(f,sgn,topbot,j)
!
!  Periodic boundary condition across the pole
!
!  15-jun-10/dhruba: aped
!  15-oct-15/fred NB use sgn= 1 for scalars and radial vector components
!                        sgn=-1 for theta and phi vector components
!            In principle similar conditions could apply for R=0
!            for sph/cyl coords, but not yet implemented
!
      real, dimension (:,:,:,:) :: f
      integer :: j,nhalf,sgn
      character (len=bclen) :: topbot
!
      if (.not.lpole(2)) call fatal_error_local('bc_pper_y',&
        "for 'p' lpole=F,T,F , lperi=F,F,T in start.in")
      if (nprocz>1 .and. modulo(nprocz,2)==1) &
        call fatal_error_local('bc_pper_y',&
             "for 'pp' nprocz must be multiple of 2")
!
      nhalf=(n1+n2)/2
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        if (nprocz==1) then
          f(:,:m1-1,n1:nhalf  ,j) = sgn*f(:,m1i:m1:-1,nhalf+1:n2,j)
          f(:,:m1-1,nhalf+1:n2,j) = sgn*f(:,m1i:m1:-1,n1:nhalf  ,j)
        endif
      case ('top')               ! top boundary
        if (nprocz==1) then
          f(:,m2+1:,n1:nhalf  ,j) = sgn*f(:,m2:m2i:-1,nhalf+1:n2,j)
          f(:,m2+1:,nhalf+1:n2,j) = sgn*f(:,m2:m2i:-1,n1:nhalf  ,j)
        endif
      case default
        print*, "bc_pper_y: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_pper_y
!***********************************************************************
    subroutine bc_per_z(f,topbot,j)
!
!  Periodic boundary condition
!
!  11-nov-02/wolf: coded
!
      real, dimension (:,:,:,:) :: f
      integer :: j
      character (len=bclen) :: topbot
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        if (nprocz==1) f(:,:,1:n1-1,j) = f(:,:,n2i:n2,j)
!
      case ('top')               ! top boundary
        if (nprocz==1) f(:,:,n2+1:,j) = f(:,:,n1:n1i,j)
!
      case default
        print*, "bc_per_z: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_per_z
!***********************************************************************
    subroutine bc_yy_z(f,topbot,j)
!
!  After-communication handling of vector quantities for Yin-Yang grid.
!
!  30-nov-15/MR: coded
!  29-feb-16/MR: avoided double transformation in ghost zone corners 
!                which is already done in bc_yy_y
!
      use General, only: transform_cart_spher

      real, dimension (:,:,:,:) :: f
      integer :: j
      character (len=bclen) :: topbot

      integer :: iya, iye
!
      if (.not.lyinyang) &
        call fatal_error_local('bc_yy_z','BC not legal as no Yin-Yang grid run.')

      if (j<=jdone) return

      if (is_vec) then
!
!  Vector quantities need to be transformed from the Cartesian basis to 
!  the local spherical basis.
!
        jdone=j+2     ! requires adjacent vector components

        iya=1; iye=my          
        if (lfirst_proc_y) iya=m1
        if (llast_proc_y) iye=m2
      
        if (topbot=='bot') then
          call transform_cart_spher(f,iya,iye,1,nghost,j)  ! in-place!
        else
          call transform_cart_spher(f,iya,iye,n2+1,mz,j)   ! ~  
        endif
      else
        jdone=0
      endif

    endsubroutine bc_yy_z
!***********************************************************************
    subroutine bc_a2r_x(f,topbot,j)
!
!  Setting d^2f/dr^2 + 2*/r*df/dr - 2*f/r^2 =0,
!  to set del2=0 in spherical coordinates.
!
!  24-nov-12/joern: coded
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      real :: tmp1,tmp2
      integer ::j
!
      select case (topbot)
!
      case ('bot')  ! bottom boundary
!
        tmp1=dx_1(l1)**2
        tmp2=dx_1(l1)*(1./x(l1)+dx_tilde(l1)/2.)
!
        f(l1-1,:,:,j)=(f(l1,:,:,j)*(-2*tmp1-2/x(l1)**2) &
                      +f(l1+1,:,:,j)*(tmp1+tmp2)) &
                      /(-1*tmp1+tmp2)
        f(l1-2,:,:,j)=(f(l1-1,:,:,j)*16*(tmp1-tmp2) &
                      +f(l1,:,:,j)*(-30*tmp1-24/x(l1)**2) &
                      +f(l1+1,:,:,j)*16*(tmp1+tmp2) &
                      +f(l1+2,:,:,j)*(-1*tmp1-2*tmp2)) &
                      /(tmp1-2*tmp2)
        f(l1-3,:,:,j)=(f(l1-2,:,:,j)*27*(-1*tmp1+2*tmp2) &
                      +f(l1-1,:,:,j)*270*(tmp1-tmp2) &
                      +f(l1,:,:,j)*(-490*tmp1-360/x(l1)**2) &
                      +f(l1+1,:,:,j)*270*(tmp1+tmp2) &
                      +f(l1+2,:,:,j)*27*(-1*tmp1-2*tmp2) &
                      +f(l1+3,:,:,j)*(2*tmp1+6*tmp2)) &
                      /(-2*tmp1+6*tmp2)
!
      case ('top')               ! top boundary
!
        tmp1=dx_1(l2)**2
        tmp2=dx_1(l2)*(1/x(l2)+dx_tilde(l2)/2.)
!
        f(l2+1,:,:,j)=(f(l2,:,:,j)*(-2*tmp1-2/x(l2)**2) &
                      +f(l2-1,:,:,j)*(tmp1-tmp2)) &
                      /(-1*tmp1-tmp2)
        f(l2+2,:,:,j)=(f(l2+1,:,:,j)*16*(tmp1+tmp2) &
                      +f(l2,:,:,j)*(-30*tmp1-24/x(l2)**2) &
                      +f(l2-1,:,:,j)*16*(tmp1-tmp2) &
                      +f(l2-2,:,:,j)*(-1*tmp1+2*tmp2)) &
                      /(tmp1+2*tmp2)
        f(l2+3,:,:,j)=(f(l2+2,:,:,j)*27*(-1*tmp1-2*tmp2) &
                      +f(l2+1,:,:,j)*270*(tmp1+tmp2) &
                      +f(l2,:,:,j)*(-490*tmp1-360/x(l2)**2) &
                      +f(l2-1,:,:,j)*270*(tmp1-tmp2) &
                      +f(l2-2,:,:,j)*27*(-1*tmp1+2*tmp2) &
                      +f(l2-3,:,:,j)*(2*tmp1-6*tmp2)) &
                      /(-2*tmp1-6*tmp2)
      case default
        print*, "bc_a2r_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_a2r_x
!***********************************************************************
    subroutine bc_sym_x(f,sgn,topbot,j,rel,val)
!
!  Symmetry boundary conditions.
!  (f,-1,topbot,j)            --> antisymmetry             (f  =0)
!  (f,+1,topbot,j)            --> symmetry                 (f' =0)
!  (f,-1,topbot,j,REL=.true.) --> generalized antisymmetry (f''=0)
!  Don't combine rel=T and sgn=1, that wouldn't make much sense.
!
!  11-nov-02/wolf: coded
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      real, dimension (:), optional :: val
      integer :: sgn,i,j
      logical, optional :: rel
      logical :: relative
!
      if (present(rel)) then; relative=rel; else; relative=.false.; endif
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        if (present(val)) f(l1,:,:,j)=val(j)
        if (relative) then
          do i=1,nghost; f(l1-i,:,:,j)=2*f(l1,:,:,j)+sgn*f(l1+i,:,:,j); enddo
        else
          do i=1,nghost; f(l1-i,:,:,j)=              sgn*f(l1+i,:,:,j); enddo
          if (sgn<0) f(l1,:,:,j) = 0. ! set bdry value=0 (indep of initcond)
        endif
!
      case ('top')               ! top boundary
        if (present(val)) f(l2,:,:,j)=val(j)
        if (relative) then
          do i=1,nghost; f(l2+i,:,:,j)=2*f(l2,:,:,j)+sgn*f(l2-i,:,:,j); enddo
        else
          do i=1,nghost; f(l2+i,:,:,j)=              sgn*f(l2-i,:,:,j); enddo
          if (sgn<0) f(l2,:,:,j) = 0. ! set bdry value=0 (indep of initcond)
        endif
!
      case default
        print*, "bc_sym_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_sym_x
!***********************************************************************
    subroutine bc_cpc_x(f,topbot,j)
!
!  This condition gives A"+A'/R=0.
!  We compute the A1 point using a 2nd-order formula,
!  i.e. A1 = - (1-dx/2R)*A_(-1)/(1+x/2R).
!  Next, we compute A2 using a 4th-order formula,
!  and finally A3 using a 6th-order formula.
!  this can not be used in the setup for -a ..a with cpc on both sides,
!  for both sides A=0 on the boundary does for example not allow a constant Bz
!  removed this restriction in cpp
!  note that for A!=0 boundary conditions for Aphi and Az are not the same,
!  hence cpz
!
!  11-nov-09/axel+koen: coded
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      real, dimension (size(f,2),size(f,3)) :: extra1,extra2
      integer :: i,j
      real :: dxR
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        dxR=-dx/x(l1)
        i=-0; f(l2+i,:,:,j)=0.
        i=-1; f(l2+i,:,:,j)=-(1.-.5*dxR)*f(l2-i,:,:,j)/(1.+.5*dxR)
        extra1=(1.+.5*dxR)*f(l2+i,:,:,j)+(1.-.5*dxR)*f(l2-i,:,:,j)
        i=-2; f(l2+i,:,:,j)=(-(1.-   dxR)*f(l2-i,:,:,j)+16.*extra1)/(1.+dxR)
        extra2=(1.+dxR)*f(l2+i,:,:,j)+(1.-dxR)*f(l2-i,:,:,j)-10.*extra1
        i=-3; f(l2+i,:,:,j)=(-(2.-3.*dxR)*f(l2-i,:,:,j)+27.*extra2)/(2.+3.*dxR)
!
      case ('top')               ! top boundary
        dxR=-dx/x(l2)
        i=0; f(l2+i,:,:,j)=0.
        i=1; f(l2+i,:,:,j)=-(1.-.5*dxR)*f(l2-i,:,:,j)/(1.+.5*dxR)
        extra1=(1.+.5*dxR)*f(l2+i,:,:,j)+(1.-.5*dxR)*f(l2-i,:,:,j)
        i=2; f(l2+i,:,:,j)=(-(1.-   dxR)*f(l2-i,:,:,j)+16.*extra1)/(1.+dxR)
        extra2=(1.+dxR)*f(l2+i,:,:,j)+(1.-dxR)*f(l2-i,:,:,j)-10.*extra1
        i=3; f(l2+i,:,:,j)=(-(2.-3.*dxR)*f(l2-i,:,:,j)+27.*extra2)/(2.+3.*dxR)
!
      case default
        print*, "bc_cpc_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_cpc_x
!***********************************************************************
    subroutine bc_cpz_x(f,topbot,j)
!
!  This condition gives R(RA)"-(RA)'=0, i e perfect conductor condition
!  for Az in cylindrical coordinates.
!  We compute the A1 point using a 2nd-order formula,
!  Next, we compute A2 using a 4th-order formula,
!  and finally A3 using a 6th-order formula.
!
!  28-feb-11/koen: coded
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      real, dimension (size(f,2),size(f,3)) :: f1_co,f2_co
      integer :: i,j
      real :: dxR
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        dxR=dx/x(l1)
        i=-1; f(l1+i,:,:,j)=(f(l1,:,:,j)*2+f(l1-i,:,:,j)*(dxR/2-1))/(dxR/2+1)
        f1_co=(1+dxR/2)*f(l1+i,:,:,j)+(1-dxR/2)*f(l1-i,:,:,j)
        i=-2; f(l1+i,:,:,j)=(-30*f(l1,:,:,j)+16*f1_co+(dxR-1)*f(l1-i,:,:,j))/(dxR+1)
        f2_co=(1+dxR)*f(l1+i,:,:,j)+(1-dxR)*f(l1-i,:,:,j)
        i=-3; f(l1+i,:,:,j)=(490*f(l1,:,:,j)-270*f1_co+27*f2_co+(3*dxR-2)*f(l1-i,:,:,j))/(3*dxR+2)
!
      case ('top')               ! top boundary
        dxR=dx/x(l2)
        i=1; f(l2+i,:,:,j)=(f(l2,:,:,j)*(2+dxR**2)+f(l2-i,:,:,j)*(dxR/2-1))/(dxR/2+1)
        f1_co=(1+dxR/2)*f(l2+i,:,:,j)+(1-dxR/2)*f(l2-i,:,:,j)
        i=2; f(l2+i,:,:,j)=(-30*f(l2,:,:,j)+16*f1_co+(dxR-1)*f(l2-i,:,:,j))/(dxR+1)
        f2_co=(1+dxR)*f(l2+i,:,:,j)+(1-dxR)*f(l2-i,:,:,j)
        i=3; f(l2+i,:,:,j)=(490*f(l2,:,:,j)-270*f1_co+27*f2_co+(3*dxR-2)*f(l2-i,:,:,j))/(3*dxR+2)
!
      case default
        print*, "bc_cpz_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_cpz_x
!***********************************************************************
    subroutine bc_cpp_x(f,topbot,j)
!
!  This condition gives RA"+A'=0, i e perfect conductor condition
!  for Aphi in cylindrical coordinates.
!  We compute the A1 point using a 2nd-order formula,
!  i.e. A1 = - (1-dx/2R)*A_(-1)/(1+x/2R).
!  Next, we compute A2 using a 4th-order formula,
!  and finally A3 using a 6th-order formula.
!
!  28-feb-11/koen: coded
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      real, dimension (size(f,2),size(f,3)) :: f1_co,f2_co
      integer :: i,j
      real :: dxR
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        dxR=dx/x(l1)
        i=-1; f(l1+i,:,:,j)=(f(l1,:,:,j)*(2+dxR**2)+f(l1-i,:,:,j)*(dxR/2-1))/(dxR/2+1)
        f1_co=(1+dxR/2)*f(l1+i,:,:,j)+(1-dxR/2)*f(l1-i,:,:,j)
        i=-2; f(l1+i,:,:,j)=((-30+12*dxR**2)*f(l1,:,:,j)+16*f1_co+(dxR-1)*f(l1-i,:,:,j))/(dxR+1)
        f2_co=(1+dxR)*f(l1+i,:,:,j)+(1-dxR)*f(l1-i,:,:,j)
        i=-3; f(l1+i,:,:,j)=((490+180*dxR**2)*f(l1,:,:,j)-270*f1_co+27*f2_co+(3*dxR-2)*f(l1-i,:,:,j))/(3*dxR+2)
!
      case ('top')               ! top boundary
        dxR=dx/x(l2)
        i=1; f(l2+i,:,:,j)=(f(l2,:,:,j)*(2+dxR**2)+f(l2-i,:,:,j)*(dxR/2-1))/(dxR/2+1)
        f1_co=(1+dxR/2)*f(l2+i,:,:,j)+(1-dxR/2)*f(l2-i,:,:,j)
        i=2; f(l2+i,:,:,j)=((-30+12*dxR**2)*f(l2,:,:,j)+16*f1_co+(dxR-1)*f(l2-i,:,:,j))/(dxR+1)
        f2_co=(1+dxR)*f(l2+i,:,:,j)+(1-dxR)*f(l2-i,:,:,j)
        i=3; f(l2+i,:,:,j)=((490+180*dxR**2)*f(l2,:,:,j)-270*f1_co+27*f2_co+(3*dxR-2)*f(l2-i,:,:,j))/(3*dxR+2)
!
      case default
        print*, "bc_cpp_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_cpp_x

!!***********************************************************************
!    subroutine bc_spr_x(f,topbot,j)
!!
!!  This condition sets values for A_phi and A_theta at the radial boundary.
!!  It solves  A"+2A'/R=0 and A=0 at the boundary.
!!  We compute the A1 point using a 2nd-order formula,
!!  Next, we compute A2 using a 4th-order formula,
!!  and finally A3 using a 6th-order formula.
!!  Has to be used togehter with 's' for A_r.
!!
!!  15-may-13/joern: coded
!!
!      character (len=bclen) :: topbot
!      real, dimension (:,:,:,:) :: f
!      integer :: j
!      real :: tmp
!!
!      select case (topbot)
!!
!      case ('bot')               ! bottom boundary
!        tmp=x(l1)*dx_1(l1)
!!
!        f(l1,:,:,j)  =0
!        f(l1-1,:,:,j)=(f(l1+1,:,:,j)*(-tmp+1))/(tmp+1)
!        f(l1-2,:,:,j)=(f(l1-1,:,:,j)*16*(tmp-1) &
!                      +f(l1+1,:,:,j)*16*(tmp+1) &
!                      +f(l1+2,:,:,j)*(-tmp-2))/(tmp-2)
!        f(l1-3,:,:,j)=(f(l1-2,:,:,j)*27*(0.5*tmp-1) &
!                      +f(l1-1,:,:,j)*135*(-tmp+1) &
!                      +f(l1+1,:,:,j)*135*(-tmp-1) &
!                      +f(l1+2,:,:,j)*27*(0.5*tmp+1) &
!                      +f(l1+3,:,:,j)*(-tmp-3))/(tmp-3)
!!
!      case ('top')               ! top boundary
!        tmp=x(l2)*dx_1(l2)
!!
!        f(l2,:,:,j)  =0
!        f(l2+1,:,:,j)=(f(l2-1,:,:,j)*(tmp+1))/(-tmp+1)
!        f(l2+2,:,:,j)=(f(l2+1,:,:,j)*16*(tmp+1) &
!                      +f(l2-1,:,:,j)*16*(tmp-1) &
!                      +f(l2-2,:,:,j)*(tmp-2))/(-tmp-2)
!        f(l2+3,:,:,j)=(f(l2+2,:,:,j)*27*(0.5*tmp+1) &
!                      +f(l2+1,:,:,j)*135*(-tmp-1) &
!                      +f(l2-1,:,:,j)*135*(-tmp+1) &
!                      +f(l2-2,:,:,j)*27*(0.5*tmp-1) &
!                      +f(l2-3,:,:,j)*(-tmp+3))/(tmp+3)
!!
!      case default
!        print*, "bc_spr_x: ", topbot, " should be 'top' or 'bot'"
!!
!      endselect
!!
!    endsubroutine bc_spr_x
!!***********************************************************************
    subroutine bc_spr_x(f,topbot,j)
!
!  This condition sets values for A_phi and A_theta at the radial boundary.
!  It solves  A"+2A'/R=0 and A=0 at the boundary (A stands for A_phi or A_theta).
!  Has to be used together with 's' for A_r.
!
!  09-may-16/fred: coded
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: ix,j
!
      if (.not.lspherical_coords) &
        call fatal_error('bc_spr_x','only implemented for spherical coordinates')
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
!
        f(l1,:,:,j) = 0.
        do ix=1,nghost
          f(l1-ix,:,:,j) = -f(l1+ix,:,:,j)*x(l1+ix)/x(l1-ix)
        enddo
!
      case ('top')               ! top boundary
!
        f(l2,:,:,j) = 0.
        do ix=1,nghost
          f(l2+ix,:,:,j) = -f(l2-ix,:,:,j)*x(l1-ix)/x(l1+ix)
        enddo
!
      case default
        print*, "bc_spr_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_spr_x
!***********************************************************************
    subroutine bc_symset_x(f,sgn,topbot,j,rel,val)
!
!  This routine works like bc_sym_x, but sets the function value to val
!
!  Symmetry boundary conditions.
!  (f,-1,topbot,j)            --> antisymmetry             (f  =0)
!  (f,+1,topbot,j)            --> symmetry                 (f' =0)
!  (f,-1,topbot,j,REL=.true.) --> generalized antisymmetry (f''=0)
!  Don't combine rel=T and sgn=1, that wouldn't make much sense.
!
!  11-nov-02/wolf: coded
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      real, dimension (:), optional :: val
      integer :: sgn,i,j
      logical, optional :: rel
      logical :: relative
!
      if (present(rel)) then; relative=rel; else; relative=.false.; endif
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        if (present(val)) f(l1,m1:m2,n1:n2,j)=val(j)
        if (relative) then
          do i=1,nghost; f(l1-i,:,:,j)=2*f(l1,:,:,j)+sgn*f(l1+i,:,:,j); enddo
        else
          do i=1,nghost; f(l1-i,:,:,j)=              sgn*f(l1+i,:,:,j); enddo
          !f(l1,:,:,j)=(4.*f(l1+1,:,:,j)-f(l1+2,:,:,j))/3.
        endif
!
      case ('top')               ! top boundary
        if (present(val)) f(l2,m1:m2,n1:n2,j)=val(j)
        if (relative) then
          do i=1,nghost; f(l2+i,:,:,j)=2*f(l2,:,:,j)+sgn*f(l2-i,:,:,j); enddo
        else
          do i=1,nghost; f(l2+i,:,:,j)=              sgn*f(l2-i,:,:,j); enddo
          !f(l2,:,:,j)=(4.*f(l2-1,:,:,j)-f(l2-2,:,:,j))/3.
        endif
!
      case default
        print*, "bc_symset_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_symset_x
!***********************************************************************
    subroutine bc_symderset_x(f,topbot,j,val)
!
!  This routine works like bc_sym_y, but sets the derivative value
!
!  30-may-11/axel: coded
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      real, dimension (:) :: val
      integer :: i,j
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        do i=1,nghost; f(l1-i,:,:,j)=f(l1+i,:,:,j)-dx2_bound(-i)*val(j); enddo
!
      case ('top')               ! top boundary
        do i=1,nghost; f(l2+i,:,:,j)=f(l2-i,:,:,j)+dx2_bound( i)*val(j); enddo
!
      case default
        print*, "bc_symderset_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_symderset_x
!***********************************************************************
    subroutine bc_symset0der_x(f,topbot,j)
!
!  This routine works like bc_sym_x, but sets the function value to what
!  it should be for vanishing one-sided derivative.
!  This is the routine to be used as regularity condition on the axis.
!
!  12-nov-09/axel+koen: coded
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: i,j,i1=1,i2=2,i3=3,i4=4,i5=5,i6=6
!
      select case (topbot)
!
!  bottom (left end of the domain)
!
      case ('bot')               ! bottom boundary
        f(l1,m1:m2,n1:n2,j)=(360.*f(l1+i1,m1:m2,n1:n2,j) &
                            -450.*f(l1+i2,m1:m2,n1:n2,j) &
                            +400.*f(l1+i3,m1:m2,n1:n2,j) &
                            -225.*f(l1+i4,m1:m2,n1:n2,j) &
                             +72.*f(l1+i5,m1:m2,n1:n2,j) &
                             -10.*f(l1+i6,m1:m2,n1:n2,j))/147.
        do i=1,nghost; f(l1-i,:,:,j)=f(l1+i,:,:,j); enddo
!
!  top (right end of the domain)
!
      case ('top')               ! top boundary
        f(l2,m1:m2,n1:n2,j)=(360.*f(l2-i1,m1:m2,n1:n2,j) &
                            -450.*f(l2-i2,m1:m2,n1:n2,j) &
                            +400.*f(l2-i3,m1:m2,n1:n2,j) &
                            -225.*f(l2-i4,m1:m2,n1:n2,j) &
                             +72.*f(l2-i5,m1:m2,n1:n2,j) &
                             -10.*f(l2-i6,m1:m2,n1:n2,j))/147.
        do i=1,nghost; f(l2+i,:,:,j)=f(l2-i,:,:,j); enddo
!
      case default
        print*, "bc_symset0der_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_symset0der_x
!***********************************************************************
    subroutine bc_slope_x(f,slope,topbot,j,rel,val)
!
! FIXME: Documentation is missing => Axel?
! WARNING: the code for "rel=.true." is currently nowhere used.
!
!  25-feb-07/axel: adapted from bc_sym_x
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      real, dimension (:), optional :: val
      real, dimension (:) :: slope
      integer :: i,j
      logical, optional :: rel
      logical :: relative
!
      if (present(rel)) then; relative=rel; else; relative=.false.; endif
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        if (present(val)) f(l1,m1:m2,n1:n2,j)=val(j)
        if (relative) then
          do i=1,nghost
            f(l1-i,:,:,j)=2*f(l1,:,:,j)+slope(j)*f(l1+i,:,:,j)*x(l1+i)/x(l1-i)
          enddo
        else
          do i=1,nghost
            f(l1-i,:,:,j)=f(l1+i,:,:,j)*(x(l1+i)/x(l1-i))**slope(j)
          enddo
!         f(l1,:,:,j)=(2.*x(l1+1)*f(l1+1,:,:,j)-&
!          .5*x(l1+2)*f(l1+2,:,:,j))/(1.5*x(l1))
        endif
!
      case ('top')               ! top boundary
        if (present(val)) f(l2,m1:m2,n1:n2,j)=val(j)
        if (relative) then
          do i=1,nghost
            f(l2+i,:,:,j)=2*f(l2,:,:,j)+slope(j)*f(l2-i,:,:,j)
          enddo
        else
          do i=1,nghost
            f(l2+i,:,:,j)=f(l2-i,:,:,j)*(x(l2-i)/x(l2+i))**slope(j)
          enddo
!         f(l2,:,:,j)=(2.*x(l2-1)*f(l2-1,:,:,j)-&
!           .5*x(l2-2)*f(l2-2,:,:,j))/(1.5*x(l2))
        endif
!
      case default
        print*, "bc_slope_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_slope_x
!***********************************************************************
    subroutine bc_ghost_slope_x(f,slope,topbot,j)
!
!  This maintains a constant slope within the ghost cells.
!
!  02-Sep-2017/PABourdin: coded as a replacement for 'bc_slope_x'
!
      real, dimension(:,:,:,:), intent(inout) :: f
      real, dimension(:), intent(in) :: slope
      character(len=bclen), intent(in) :: topbot
      integer, intent(in) :: j
!
      integer :: i
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        do i = 1, nghost
          f(l1-i,:,:,j) = f(l1,:,:,j) + slope(j) * (x(l1-i) - x(l1))
        enddo
!
      case ('top')               ! top boundary
        do i = 1, nghost
          f(l2+i,:,:,j) = f(l2,:,:,j) + slope(j) * (x(l2+i) - x(l2))
        enddo
!
      case default
        print *, "bc_ghost_slope_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_ghost_slope_x
!***********************************************************************
    subroutine bc_shear_x(f,slope,abscissa,topbot,j)
!
!  This maintains a constant shear proportional to x at the boundary.
!
!  02-Sep-2017/PABourdin: coded
!
      real, dimension(:,:,:,:), intent(inout) :: f
      real, dimension(:), intent(in) :: slope, abscissa
      character(len=bclen), intent(in) :: topbot
      integer, intent(in) :: j
!
      integer :: i
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        do i = 1, nghost
          f(l1-i,:,:,j) = abscissa(j) + slope(j) * x(l1-i)
        enddo
!
      case ('top')               ! top boundary
        do i = 1, nghost
          f(l2+i,:,:,j) = abscissa(j) + slope(j) * x(l2+i)
        enddo
!
      case default
        print *, "bc_shear_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_shear_x
!***********************************************************************
    subroutine bc_shear_y(f,slope,abscissa,topbot,j)
!
!  This maintains a constant shear proportional to y at the boundary.
!
!  04-Sep-2017/PABourdin: coded
!
      real, dimension(:,:,:,:), intent(inout) :: f
      real, dimension(:), intent(in) :: slope, abscissa
      character(len=bclen), intent(in) :: topbot
      integer, intent(in) :: j
!
      integer :: i
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        do i = 1, nghost
          f(:,m1-i,:,j) = abscissa(j) + slope(j) * y(m1-i)
        enddo
!
      case ('top')               ! top boundary
        do i = 1, nghost
          f(:,m2+i,:,j) = abscissa(j) + slope(j) * y(m2+i)
        enddo
!
      case default
        print *, "bc_shear_y: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_shear_y
!***********************************************************************
    subroutine bc_shear_z(f,slope,abscissa,topbot,j)
!
!  This maintains a constant shear proportional to z at the boundary.
!
!  04-Sep-2017/PABourdin: coded
!
      real, dimension(:,:,:,:), intent(inout) :: f
      real, dimension(:), intent(in) :: slope, abscissa
      character(len=bclen), intent(in) :: topbot
      integer, intent(in) :: j
!
      integer :: i
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        do i = 1, nghost
          f(:,:,n1-i,j) = abscissa(j) + slope(j) * z(n1-i)
        enddo
!
      case ('top')               ! top boundary
        do i = 1, nghost
          f(:,:,n2+i,j) = abscissa(j) + slope(j) * z(n2+i)
        enddo
!
      case default
        print *, "bc_shear_z: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_shear_z
!***********************************************************************
    subroutine bc_dr0_x(f,slope,topbot,j,rel,val)
!
! FIXME: This documentation is almost certainly wrong
!
!  Symmetry boundary conditions.
!  (f,-1,topbot,j)            --> antisymmetry             (f  =0)
!  (f,+1,topbot,j)            --> symmetry                 (f' =0)
!  (f,-1,topbot,j,REL=.true.) --> generalized antisymmetry (f''=0)
!  Don't combine rel=T and sgn=1, that wouldn't make much sense.
!
!  25-feb-07/axel: adapted from bc_sym_x
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      real, dimension (:), optional :: val
      real, dimension (:) :: slope
      integer :: i,j
      ! Abbreviations to keep compiler from complaining in 1-d or 2-d:
      integer :: l1_4, l1_5, l1_6
      integer :: l2_4, l2_5, l2_6
      logical, optional :: rel
      logical :: relative
!
      l1_4=l1+4; l1_5=l1+5; l1_6=l1+6
      l2_4=l2-4; l2_5=l2-5; l2_6=l2-6
!
      if (present(rel)) then; relative=rel; else; relative=.false.; endif
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        if (present(val)) f(l1,m1:m2,n1:n2,j)=val(j)
        if (relative) then
          do i=1,nghost
            f(l1-i,:,:,j)=2*f(l1,:,:,j)+slope(j)*f(l1+i,:,:,j)*x(l1+i)/x(l1-i)
          enddo
        else
          f(l1,:,:,j)=(360.*x(l1+1)*f(l1+1,:,:,j)-450.*x(l1+2)*f(l1+2,:,:,j) &
                      +400.*x(l1+3)*f(l1+3,:,:,j)-225.*x(l1_4)*f(l1_4,:,:,j) &
                       +72.*x(l1_5)*f(l1_5,:,:,j)- 10.*x(l1_6)*f(l1_6,:,:,j) &
                      )/(147.*x(l1))
          do i=1,nghost
            f(l1-i,:,:,j)=f(l1+i,:,:,j)+dx2_bound(-i)/x(l1)*f(l1,:,:,j)
          enddo
        endif
!
      case ('top')               ! top boundary
        if (present(val)) f(l2,m1:m2,n1:n2,j)=val(j)
        if (relative) then
          do i=1,nghost
            f(l2+i,:,:,j)=2*f(l2,:,:,j)+slope(j)*f(l2-i,:,:,j)
          enddo
        else
          f(l2,:,:,j)=(360.*x(l2-1)*f(l2-1,:,:,j)-450.*x(l2-2)*f(l2-2,:,:,j) &
                      +400.*x(l2-3)*f(l2-3,:,:,j)-225.*x(l2_4)*f(l2_4,:,:,j) &
                       +72.*x(l2_5)*f(l2_5,:,:,j)- 10.*x(l2_6)*f(l2_6,:,:,j) &
                      )/(147.*x(l2))
          do i=1,nghost
            f(l2+i,:,:,j)=f(l2-i,:,:,j)-dx2_bound(i)/x(l2)*f(l2,:,:,j)
          enddo
        endif
!
      case default
        print*, "bc_dr0_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_dr0_x
!***********************************************************************
    subroutine bc_overshoot_x(f,dist,topbot,j)
!
!  Overshoot boundary conditions, ie (d/dx-1/dist) f = 0.
!  Is implemented as d/dx [ f*exp(-x/dist) ] = 0,
!  so f(l1-i)*exp[-x(l1-i)/dist] = f(l1+i)*exp[-x(l1+i)/dist],
!  or f(l1-i) = f(l1+i)*exp{[x(l1-i)-x(l1+i)]/dist}.
!
!  25-feb-07/axel: adapted from bc_sym_x
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      real, dimension (:) :: dist
      integer :: i,j
!
      select case (topbot)
!
!  bottom
!
      case ('bot')               ! bottom boundary
        do i=1,nghost
          f(l1-i,:,:,j)=f(l1+i,:,:,j)*exp(-dx2_bound(-i)/dist(j))
        enddo
!
!  top
!
      case ('top')               ! top boundary
        do i=1,nghost
          f(l2+i,:,:,j)=f(l2-i,:,:,j)*exp(dx2_bound(i))/dist(j)
        enddo
!
!  default
!
      case default
        print*, "bc_overshoot_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_overshoot_x
!***********************************************************************
    subroutine bc_overshoot_z(f,dist,topbot,j)
!
!  Overshoot boundary conditions, ie (d/dz-1/dist) f = 0.
!  Is implemented as d/dz [ f*exp(-z/dist) ] = 0,
!  so f(n1-i)*exp[-z(n1-i)/dist] = f(n1+i)*exp[-z(n1+i)/dist],
!  or f(n1-i) = f(n1+i)*exp{[z(n1-i)-z(n1+i)]/dist}.
!
!  25-feb-07/axel: adapted from bc_sym_z
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      real, dimension (:) :: dist
      integer :: i,j
!
      select case (topbot)
!
!  bottom
!
      case ('bot')               ! bottom boundary
        do i=1,nghost
          f(:,:,n1-i,j)=f(:,:,n1+i,j)*exp(-dz2_bound(-i)/dist(j))
        enddo
!
!  top
!
      case ('top')               ! top boundary
        do i=1,nghost
          f(:,:,n2+i,j)=f(:,:,n2-i,j)*exp(dz2_bound(i)/dist(j))
        enddo
!
!  default
!
      case default
        print*, "bc_overshoot_z: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_overshoot_z
!***********************************************************************
    subroutine bc_antis_x(f,slope,topbot,j,rel,val)
!
!  Print a warning to prompt potential users to document this.
!  This routine seems an experimental one to me (Axel)
!
!  Symmetry boundary conditions.
!  (f,-1,topbot,j)            --> antisymmetry             (f  =0)
!  (f,+1,topbot,j)            --> symmetry                 (f' =0)
!  (f,-1,topbot,j,REL=.true.) --> generalized antisymmetry (f''=0)
!  Don't combine rel=T and sgn=1, that wouldn't make much sense.
!
!  25-feb-07/axel: adapted from bc_slope_x
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      real, dimension (:), optional :: val
      real, dimension (:) :: slope
      integer :: i,j
      logical, optional :: rel
      logical :: relative
!
!  Print a warning to prompt potential users to document this.
!
      call fatal_error('bc_antis_x','outdated/invalid? Document if needed')
!
      if (present(rel)) then; relative=rel; else; relative=.false.; endif
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        if (present(val)) f(l1,m1:m2,n1:n2,j)=val(j)
        if (relative) then
          do i=1,nghost
            f(l1-i,:,:,j)=2*f(l1,:,:,j)+slope(j)*f(l1+i,:,:,j)*x(l1+i)/x(l1-i)
          enddo
        else
          f(l1,:,:,j)=0.
          do i=1,nghost
            f(l1-i,:,:,j)=-f(l1+i,:,:,j)*(x(l1+i)/x(l1-i))**slope(j)
          enddo
        endif
!
      case ('top')               ! top boundary
        if (present(val)) f(l2,m1:m2,n1:n2,j)=val(j)
        if (relative) then
          do i=1,nghost
            f(l2+i,:,:,j)=2*f(l2,:,:,j)+slope(j)*f(l2-i,:,:,j)
          enddo
        else
          f(l2,:,:,j)=0.
          do i=1,nghost
            f(l2+i,:,:,j)=-f(l2-i,:,:,j)*(x(l2-i)/x(l2+i))**slope(j)
          enddo
        endif
!
      case default
        print*, "bc_antis_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_antis_x
!***********************************************************************
    subroutine bc_sym_y(f,sgn,topbot,j,rel,val,val2,val4)
!
!  Symmetry boundary conditions.
!  (f,-1,topbot,j)            --> antisymmetry             (f  =0)
!  (f,+1,topbot,j)            --> symmetry                 (f' =0)
!  (f,-1,topbot,j,REL=.true.) --> generalized antisymmetry (f''=0)
!  Don't combine rel=T and sgn=1, that wouldn't make much sense.
!
!  11-nov-02/wolf: coded
!  10-apr-05/axel: added val argument
!   9-jun-11/axel: added val2 argument
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      real, dimension (:), optional :: val,val2,val4
      integer :: sgn,i,j
      logical, optional :: rel
      logical :: relative
!
      if (present(rel)) then; relative=rel; else; relative=.false.; endif
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        if (present(val)) f(:,m1,:,j)=val(j)
        if (present(val2)) f(:,m1,:,j)=f(:,m1,:,j)+val2(j)*spread(x**2,2,size(f,3))
        if (present(val4)) f(:,m1,:,j)=f(:,m1,:,j)+val4(j)*spread(x**4,2,size(f,3))
        if (relative) then
          do i=1,nghost; f(:,m1-i,:,j)=2*f(:,m1,:,j)+sgn*f(:,m1+i,:,j); enddo
        else
          do i=1,nghost; f(:,m1-i,:,j)=              sgn*f(:,m1+i,:,j); enddo
          if (sgn<0) f(:,m1,:,j) = 0. ! set bdry value=0 (indep of initcond)
        endif
!
      case ('top')               ! top boundary
        if (present(val)) f(:,m2,:,j)=val(j)
        if (present(val2)) f(:,m2,:,j)=f(:,m2,:,j)+val2(j)*spread(x**2,2,size(f,3))
        if (present(val4)) f(:,m2,:,j)=f(:,m2,:,j)+val4(j)*spread(x**4,2,size(f,3))
        if (relative) then
          do i=1,nghost; f(:,m2+i,:,j)=2*f(:,m2,:,j)+sgn*f(:,m2-i,:,j); enddo
        else
          do i=1,nghost; f(:,m2+i,:,j)=              sgn*f(:,m2-i,:,j); enddo
          if (sgn<0) f(:,m2,:,j) = 0. ! set bdry value=0 (indep of initcond)
        endif
!
      case default
        print*, "bc_sym_y: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_sym_y
!***********************************************************************
    subroutine bc_stratified_y(f,topbot,j)
!
!  Boundary condition that maintains hydrostatic equilibrium in the meriodional direction.
!  This boundary is coded only for spherical coordinates.
!
!  06-oct-13/wlad: coded
!
      use EquationOfState, only: cs0
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      real, dimension(size(f,1)) :: rad,za,zg,H,lnrho
      integer :: i,in,j
!
      if (.not.(j==irho.or.j==ilnrho)) &
           call fatal_error("bc_stratified_y","This boundary condition is specific for density")
      if (.not.lspherical_coords) &
           call fatal_error("bc_stratified_y","This boudary condition is for spherical coordinates only")
!
      rad=x
!
      select case (topbot)
      case ('bot')
        za=rad*costh(m1)
        H=cs0*rad
        do i=1,nghost
          zg=rad*costh(m1-i)
          do in=1,size(f,3)
            if (ldensity_nolog) then 
              lnrho = alog(f(:,m1,in,j)) - (zg**2-za**2)/(2*H**2)
              f(:,m1-i,in,j) = exp(lnrho)
            else
              lnrho = f(:,m1,in,j) - (zg**2-za**2)/(2*H**2)
              f(:,m1-i,in,j) = lnrho
            endif
          enddo
        enddo
!
      case ('top')
        za=rad*costh(m2)
        H=cs0*rad
        do i=1,nghost
          zg=rad*costh(m2+i)
          do in=1,size(f,3)
            if (ldensity_nolog) then
              lnrho = alog(f(:,m2,in,j)) - (zg**2-za**2)/(2*H**2)
              f(:,m2+i,in,j) = exp(lnrho)
            else
              lnrho = f(:,m2,in,j) - (zg**2-za**2)/(2*H**2)
              f(:,m2+i,in,j) = lnrho
            endif
          enddo
        enddo
!
      case default
        print*, "bc_sym_y: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_stratified_y
!***********************************************************************
    subroutine bc_symset_y(f,sgn,topbot,j,rel,val)
!
!  This routine works like bc_sym_y, but sets the function value to what
!  it should be for vanishing one-sided derivative.
!  At the moment the derivative is only 2nd order accurate.
!
!  Symmetry boundary conditions.
!  (f,-1,topbot,j)            --> antisymmetry             (f  =0)
!  (f,+1,topbot,j)            --> symmetry                 (f' =0)
!  (f,-1,topbot,j,REL=.true.) --> generalized antisymmetry (f''=0)
!  Don't combine rel=T and sgn=1, that wouldn't make much sense.
!
!  11-nov-02/wolf: coded
!  10-apr-05/axel: added val argument
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      real, dimension (:), optional :: val
      integer :: sgn,i,j
      logical, optional :: rel
      logical :: relative
!
      if (present(rel)) then; relative=rel; else; relative=.false.; endif
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        if (present(val)) f(l1:l2,m1,n1:n2,j)=val(j)
        if (relative) then
          do i=1,nghost; f(:,m1-i,:,j)=2*f(:,m1,:,j)+sgn*f(:,m1+i,:,j); enddo
        else
          do i=1,nghost; f(:,m1-i,:,j)=              sgn*f(:,m1+i,:,j); enddo
          !f(:,m1,:,j)=(4.*f(:,m1+1,:,j)-f(:,m1+2,:,j))/3.
        endif
!
      case ('top')               ! top boundary
        if (present(val)) f(l1:l2,m2,n1:n2,j)=val(j)
        if (relative) then
          do i=1,nghost; f(:,m2+i,:,j)=2*f(:,m2,:,j)+sgn*f(:,m2-i,:,j); enddo
        else
          do i=1,nghost; f(:,m2+i,:,j)=              sgn*f(:,m2-i,:,j); enddo
          !f(:,m2,:,j)=(4.*f(:,m2-1,:,j)-f(:,m2-2,:,j))/3.
        endif
!
      case default
        print*, "bc_symset_y: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_symset_y
!***********************************************************************
    subroutine bc_symderset_y(f,topbot,j,val)
!
!  This routine works like bc_sym_y, but sets the derivative value
!
!  30-may-11/axel: coded
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      real, dimension (:) :: val
      integer :: i,j
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        do i=1,nghost; f(:,m1-i,:,j)=f(:,m1+i,:,j)-dy2_bound(-i)*val(j); enddo
!
      case ('top')               ! top boundary
        do i=1,nghost; f(:,m2+i,:,j)=f(:,m2-i,:,j)+dy2_bound(i)*val(j); enddo
!
      case default
        print*, "bc_symderset_y: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_symderset_y
!***********************************************************************
    subroutine bc_csymderset_y(f,topbot,j,val)
!
!  This routine works like bc_sym_y, but sets the derivative value
!
!  30-may-11/axel: coded
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      real, dimension (size(f,1),size(f,3)) :: derval
      real, dimension (:) :: val
      integer :: i,j
!
      derval=spread((xyz1(1)-x)*val(j),2,size(f,3))
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        do i=1,nghost; f(:,m1-i,:,j)=f(:,m1+i,:,j)-dy2_bound(-i)*derval; enddo
!
      case ('top')               ! top boundary
        do i=1,nghost; f(:,m2+i,:,j)=f(:,m2-i,:,j)+dy2_bound(i)*derval; enddo
!
      case default
        print*, "bc_csymderset_y: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_csymderset_y
!***********************************************************************
    subroutine bc_symset0der_y(f,topbot,j)
!
!  This routine works like bc_sym_y, but sets the function value to what
!  it should be for vanishing one-sided derivative.
!  This is the routine to be used as regularity condition on the axis.
!
!  19-nov-09/axel: adapted from bc_symset0der_x
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: i,j,i1=1,i2=2,i3=3,i4=4,i5=5,i6=6
!
      select case (topbot)
!
!  bottom (left end of the domain)
!
      case ('bot')               ! bottom boundary
        f(:,m1,:,j)=(360.*f(:,m1+i1,:,j) &
                    -450.*f(:,m1+i2,:,j) &
                    +400.*f(:,m1+i3,:,j) &
                    -225.*f(:,m1+i4,:,j) &
                     +72.*f(:,m1+i5,:,j) &
                     -10.*f(:,m1+i6,:,j))/147.
        do i=1,nghost; f(:,m1-i,:,j)=f(:,m1+i,:,j); enddo
!
!  top (right end of the domain)
!
      case ('top')               ! top boundary
        f(:,m2,:,j)=(360.*f(:,m2-i1,:,j) &
                    -450.*f(:,m2-i2,:,j) &
                    +400.*f(:,m2-i3,:,j) &
                    -225.*f(:,m2-i4,:,j) &
                     +72.*f(:,m2-i5,:,j) &
                     -10.*f(:,m2-i6,:,j))/147.
        do i=1,nghost; f(:,m2+i,:,j)=f(:,m2-i,:,j); enddo
!
      case default
        print*, "bc_symset0der_y: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_symset0der_y
!***********************************************************************
    subroutine bc_spt_y(f,topbot,j)
!
!  This condition sets values for A_r or/and A_phi at the theta boundary.
!  It solves A"+\cot(theta)A'=0 and A=0 at the boundary.
!  We compute the A1 point using a 2nd-order formula,
!  Next, we compute A2 using a 4th-order formula,
!  and finally A3 using a 6th-order formula.
!  is has to be used togehter with 'sse' with 'fbcy_top' or 'fbcy_bot'
!  where A_theta=0
!
!  23-may-13/joern: coded
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      real :: tmp
      integer :: j
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        tmp=cotth(m1)/dy_1(m1)
!
        f(:,m1,:,j)  =0
        f(:,m1-1,:,j)=(f(:,m1+1,:,j)*(-1-0.5*tmp))/(1-0.5*tmp)
        f(:,m1-2,:,j)=(f(:,m1-1,:,j)*8*(2-tmp) &
                      +f(:,m1+1,:,j)*8*(2+tmp) &
                      +f(:,m1+2,:,j)*(-1-tmp))/(1-tmp)
        f(:,m1-3,:,j)=(f(:,m1-2,:,j)*13.5*(1-tmp) &
                      +f(:,m1-1,:,j)*135*(-1+0.5*tmp) &
                      +f(:,m1+1,:,j)*135*(-1-0.5*tmp) &
                      +f(:,m1+2,:,j)*13.5*(1+tmp) &
                      +f(:,m1+3,:,j)*(-1-1.5*tmp))/(1-1.5*tmp)
!
      case ('top')               ! top boundary
        tmp=cotth(m2)/dy_1(m2)
!
        f(:,m2,:,j)  =0
        f(:,m2+1,:,j)=(f(:,m2-1,:,j)*(-1-0.5*tmp))/(1-0.5*tmp)
        f(:,m2+2,:,j)=(f(:,m2+1,:,j)*8*(2-tmp) &
                      +f(:,m2+1,:,j)*8*(2+tmp) &
                      +f(:,m2+2,:,j)*(-1-tmp))/(1-tmp)
        f(:,m2+3,:,j)=(f(:,m2+2,:,j)*13.5*(1-tmp) &
                      +f(:,m2+1,:,j)*135*(-1+0.5*tmp) &
                      +f(:,m2-1,:,j)*135*(-1-0.5*tmp) &
                      +f(:,m2-2,:,j)*13.5*(1+tmp) &
                      +f(:,m2-3,:,j)*(-1-1.5*tmp))/(1-1.5*tmp)
!
      case default
        print*, "bc_spt_y: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_spt_y
!***********************************************************************
    subroutine bc_sym_z(f,sgn,topbot,j,rel,val,val2,val4)
!
!  Symmetry boundary conditions.
!  (f,-1,topbot,j)            --> antisymmetry             (f  =0)
!  (f,+1,topbot,j)            --> symmetry                 (f' =0)
!  (f,-1,topbot,j,REL=.true.) --> generalized antisymmetry (f''=0)
!  Don't combine rel=T and sgn=1, that wouldn't make much sense.
!
!  11-nov-02/wolf: coded
!  10-apr-05/axel: added val argument
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      real, dimension (:), optional :: val,val2,val4
      integer :: sgn,i,j
      logical, optional :: rel
      logical :: relative
!
      if (present(rel)) then; relative=rel; else; relative=.false.; endif
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        if (present(val)) f(:,:,n1,j)=val(j)
        if (present(val2)) f(:,:,n1,j)=f(:,:,n1,j)+val2(j)*spread(x**2,2,size(f,2))
        if (present(val4)) f(:,:,n1,j)=f(:,:,n1,j)+val4(j)*spread(x**4,2,size(f,2))
        if (relative) then
          do i=1,nghost; f(:,:,n1-i,j)=2*f(:,:,n1,j)+sgn*f(:,:,n1+i,j);
            if (.false..and.j==3) then
!if (i==1) print*, f(4,4:9,n1,j)
!if (i==1) print*, f(4,4:23,n1-1,j)
              if (any(f(4:131,4:131,n1-i,j)/=f(4:131,4:131,n1+i,j))) &
                print'(a,i2,1x,e20.12)','boundcond ghost:i=', i, maxval(abs(f(4:131,4:131,n1-i,j)-f(4:131,4:131,n1+i,j)))
              if (any(f(4:131,4:131,n1-i,j)/=f(4:131,4:131,n1,j))) &
                print'(a,i2,1x,e20.12)','boundcond bound-:i=', i, maxval(abs(f(4:131,4:131,n1-i,j)-f(4:131,4:131,n1,j)))
              if (any(f(4:131,4:131,n1+i,j)/=f(4:131,4:131,n1,j))) &
                print'(a,i2,1x,e20.12)','boundcond bound+:i=', i, maxval(abs(f(4:131,4:131,n1+i,j)-f(4:131,4:131,n1,j)))
            endif
          enddo
        else
!if (ldownsampling) print*, 'size,n1,j=', size(f,1), size(f,2), size(f,3), size(f,4),n1,j 
          do i=1,nghost; f(:,:,n1-i,j)=              sgn*f(:,:,n1+i,j); enddo
          if (sgn<0) f(:,:,n1,j) = 0. ! set bdry value=0 (indep of initcond)
        endif
!
      case ('top')               ! top boundary
        if (present(val)) f(:,:,n2,j)=val(j)
        if (present(val2)) f(:,:,n2,j)=f(:,:,n2,j)+val2(j)*spread(x**2,2,size(f,2))
        if (present(val4)) f(:,:,n2,j)=f(:,:,n2,j)+val4(j)*spread(x**4,2,size(f,2))
        if (relative) then
          do i=1,nghost; f(:,:,n2+i,j)=f(:,:,n2,j)+(f(:,:,n2,j)+sgn*f(:,:,n2-i,j));
          enddo
        else
          do i=1,nghost; f(:,:,n2+i,j)=              sgn*f(:,:,n2-i,j); enddo
          if (sgn<0) f(:,:,n2,j) = 0. ! set bdry value=0 (indep of initcond)
        endif
!
      case default
        print*, "bc_sym_z: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_sym_z
!***********************************************************************
    subroutine bc_sf_x(f,sgn,topbot,j)
!
!  Symmetric/antisymmetric boundary conditions with respect to the interface.
!  i.e. where the reflection plane is between the last mesh point and first
!  ghost point
!
!    sgn = +1  -->  symmetric
!    sgn = -1  -->  antisymmetric
!
!  12-nov-16/ccyang: coded
!
      real, dimension(:,:,:,:), intent(inout) :: f
      integer, intent(in) :: sgn, j
      character(3), intent(in) :: topbot
!
      integer :: i
!
      select case(topbot)
      case('bot')               ! bottom boundary
        forall (i=1:nghost) f(l1-i,:,:,j) = real(sgn) * f(l1+i-1,:,:,j)
      case('top')               ! top boundary
        forall (i=1:nghost) f(l2+i,:,:,j) = real(sgn) * f(l2-i+1,:,:,j)
      case default
        print *, 'bc_sf_x: unknown input; topbot = ', topbot
      endselect
!
    endsubroutine bc_sf_x
!***********************************************************************
    subroutine bc_sf_y(f,sgn,topbot,j)
!
!  Symmetric/antisymmetric boundary conditions with respect to the interface.
!  i.e. where the reflection plane is between the last mesh point and first
!  ghost point
!
!    sgn = +1  -->  symmetric
!    sgn = -1  -->  antisymmetric
!
!  12-nov-16/ccyang: coded
!
      real, dimension(:,:,:,:), intent(inout) :: f
      integer, intent(in) :: sgn, j
      character(3), intent(in) :: topbot
!
      integer :: i
!
      select case(topbot)
      case('bot')               ! bottom boundary
        forall (i=1:nghost) f(:,m1-i,:,j) = real(sgn) * f(:,m1+i-1,:,j)
      case('top')               ! top boundary
        forall (i=1:nghost) f(:,m2+i,:,j) = real(sgn) * f(:,m2-i+1,:,j)
      case default
        print *, 'bc_sf_y: unknown input; topbot = ', topbot
      endselect
!
    endsubroutine bc_sf_y
!***********************************************************************
    subroutine bc_sf_z(f,sgn,topbot,j)
!
!  Symmetric/antisymmetric boundary conditions with respect to the interface.
!  i.e. where the reflection plane is between the last mesh point and first
!  ghost point
!
!    sgn = +1  -->  symmetric
!    sgn = -1  -->  antisymmetric
!
!  14-feb-09/ccyang: coded
!
      real, dimension(:,:,:,:), intent(inout) :: f
      integer, intent(in) :: sgn, j
      character(3), intent(in) :: topbot
!
      integer :: i
!
      select case(topbot)
      case('bot')               ! bottom boundary
        forall (i=1:nghost) f(:,:,n1-i,j) = real(sgn) * f(:,:,n1+i-1,j)
      case('top')               ! top boundary
        forall (i=1:nghost) f(:,:,n2+i,j) = real(sgn) * f(:,:,n2-i+1,j)
      case default
        print *, 'bc_sf_z: unknown input; topbot = ', topbot
      endselect
!
    endsubroutine bc_sf_z
!***********************************************************************
    subroutine bc_symset0der_z(f,topbot,j)
!
!  This routine works like bc_sym_z, but sets the function value to what
!  it should be for vanishing one-sided derivative.
!  This is the routine to be used as regularity condition on the axis.
!
!  22-nov-09/axel: adapted from bc_symset0der_y
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: i,j,i1=1,i2=2,i3=3,i4=4,i5=5,i6=6
!
      select case (topbot)
!
!  bottom (left end of the domain)
!
      case ('bot')               ! bottom boundary
        f(:,:,n1,j)=(360.*f(:,:,n1+i1,j) &
                    -450.*f(:,:,n1+i2,j) &
                    +400.*f(:,:,n1+i3,j) &
                    -225.*f(:,:,n1+i4,j) &
                     +72.*f(:,:,n1+i5,j) &
                     -10.*f(:,:,n1+i6,j))/147.
        do i=1,nghost; f(:,:,n1-i,j)=f(:,:,n1+i,j); enddo
!
!  top (right end of the domain)
!
      case ('top')               ! top boundary
        f(:,:,n2,j)=(360.*f(:,:,n2-i1,j) &
                    -450.*f(:,:,n2-i2,j) &
                    +400.*f(:,:,n2-i3,j) &
                    -225.*f(:,:,n2-i4,j) &
                     +72.*f(:,:,n2-i5,j) &
                     -10.*f(:,:,n2-i6,j))/147.
        do i=1,nghost; f(:,:,n2+i,j)=f(:,:,n2-i,j); enddo
!
      case default
        print*, "bc_symset0der_z: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_symset0der_z
!***********************************************************************
    subroutine bc_set_der_x(f,topbot,j,val)
!
!  Sets the derivative on the boundary to a given value.
!
!  14-may-2006/tobi: coded
!
      character (len=bclen), intent (in) :: topbot
      real, dimension (:,:,:,:), intent (inout) :: f
      integer, intent (in) :: j
      real, intent (in) :: val
!
      integer :: i
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        do i=1,nghost; f(l1-i,:,:,j) = f(l1+i,:,:,j) - dx2_bound(-i)*val; enddo
!
      case ('top')               ! top boundary
        do i=1,nghost; f(l2+i,:,:,j) = f(l2-i,:,:,j) + dx2_bound(i)*val; enddo
!
      case default
        call warning('bc_set_der_x',topbot//" should be 'top' or 'bot'")
!
      endselect
!
    endsubroutine bc_set_der_x
!***********************************************************************
    subroutine bc_fix_x(f,topbot,j,val)
!
!  Sets the value of f, particularly:
!    A_{\alpha}= <val>
!  on the boundary to a given value
!
!  27-apr-2007/dhruba: coded
!
      character (len=bclen), intent (in) :: topbot
      real, dimension (:,:,:,:), intent (inout) :: f
      integer, intent (in) :: j
!
      real, intent (in) :: val
      integer :: i
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        do i=1,nghost; f(l1-i,:,:,j)=val; enddo
      case ('top')               ! top boundary
        do i=1,nghost; f(l2+i,:,:,j)=val; enddo
      case default
        call warning('bc_fix_x',topbot//" should be 'top' or 'bot'")
!
      endselect
!
    endsubroutine bc_fix_x
!***********************************************************************
    subroutine bc_file_x(f,topbot,j)
!
!  Sets the value of f from a file
!
!   9-jan-2008/axel+nils+natalia: coded
!
      character (len=bclen), intent (in) :: topbot
      real, dimension (:,:,:,:), intent (inout) :: f
      integer, intent (in) :: j
!
      real, dimension (:,:,:,:), allocatable :: bc_file_x_array
      integer :: i,lbc0,lbc1,lbc2,stat,iszx,io_code
      real :: lbc,frac
      logical, save :: lbc_file_x=.true.
!
      if (ldownsampling) then
        call warning('bc_file_x','Not available for downsampling')
        return
      endif
!
!  Allocate memory for large array.
!
      allocate(bc_file_x_array(mx,my,mz,mvar),stat=stat)
      if (stat>0) call fatal_error('bc_file_x', &
          'Could not allocate memory for bc_file_x_array')
!
      if (lbc_file_x) then
        if (lroot) then
          print*,'opening bc_file_x.dat'
          open(9,file=trim(directory_dist)//'/bc_file_x.dat',form='unformatted')
          read(9,iostat=io_code) bc_file_x_array
          if (io_code < 0) then
            ! end of file
            if (lroot) print*,'need file with dimension: ',mx,my,mz,mvar
            deallocate(bc_file_x_array)
            call stop_it("boundary file bc_file_x.dat has incorrect size")
          endif
          close(9)
        endif
        lbc_file_x=.false.
      endif
      
      iszx=size(f,1)
!
      select case (topbot)
!
!  x - Udrift_bc*t = dx * (ix - Udrift_bc*t/dx)
!
      case ('bot')               ! bottom boundary
        lbc=Udrift_bc*t*dx_1(1)+1.
        lbc0=int(lbc)
        frac=mod(lbc,real(lbc0))
        lbc1=iszx+mod(-lbc0,iszx)
        lbc2=iszx+mod(-lbc0-1,iszx)
        do i=1,nghost
          f(l1-i,:,:,j)=(1-frac)*bc_file_x_array(lbc1,:,:,j) &
                           +frac*bc_file_x_array(lbc2,:,:,j)
        enddo
      case ('top')               ! top boundary
!
!  note: this "top" thing hasn't been adapted or tested yet.
!  The -lbc0-1 has been changed to +lbc0+1, but has not been tested yet.
!
        lbc=Udrift_bc*t*dx_1(1)+1.
        lbc0=int(lbc)
        frac=mod(lbc,real(lbc0))
        lbc1=iszx+mod(+lbc0,iszx)
        lbc2=iszx+mod(+lbc0+1,iszx)
        do i=1,nghost
          f(l2+i,:,:,j)=(1-frac)*bc_file_x_array(lbc1,:,:,j) &
                           +frac*bc_file_x_array(lbc2,:,:,j)
        enddo
      case default
        call warning('bc_fix_x',topbot//" should be 'top' or 'bot'")
!
      endselect
!
      deallocate(bc_file_x_array)
!
    endsubroutine bc_file_x
!***********************************************************************
    subroutine bc_set_spder_x(f,topbot,j,val)
!
!  Sets the derivative, particularly:
!    d(rA_{\alpha})/dr = <val>
!  on the boundary to a given value
!
!  27-apr-2007/dhruba: coded
!
      character (len=bclen), intent (in) :: topbot
      real, dimension (:,:,:,:), intent (inout) :: f
      integer, intent (in) :: j
!
      real, intent (in) :: val
      integer :: i
!
      if (lspherical_coords)then
        select case (topbot)
        case ('bot')               ! bottom boundary
          do i=1,nghost
            f(l1-i,:,:,j)=f(l1+i,:,:,j)-dx2_bound(-i)*(val-f(l1,:,:,j)*r1_mn(1))
          enddo
        case ('top')               ! top boundary
          do i=1,nghost
            f(l2+i,:,:,j)=f(l2-i,:,:,j)+dx2_bound(i)*(val-f(l2,:,:,j)*r1_mn(nx))
          enddo
!
        case default
          call warning('bc_set_spder_x',topbot//" should be 'top' or 'bot'")
!
        endselect
      else
        call stop_it('bc_set_spder_x valid only in spherical coordinate system')
      endif
!
    endsubroutine bc_set_spder_x
! **********************************************************************
    subroutine bc_set_pfc_x(f,topbot,j)
!
!joern: WARNING, this bc will NOT give a perfect-conductor boundary condition
!
!  In spherical polar coordinate system,
!  at a radial boundary set : $A_{\theta} = 0$ and $A_{phi} = 0$,
!  and demand $div A = 0$ gives the condition on $A_r$ to be
!  $d/dr( A_r) + 2 A_r/r = 0$ . This subroutine sets this condition of
!  $j$ the component of f. As this is related to setting the
!  perfect conducting boundary condition we call this "pfc".
!
!  25-Aug-2007/dhruba: coded
!
      character (len=bclen), intent (in) :: topbot
      real, dimension (:,:,:,:), intent (inout) :: f
      integer, intent (in) :: j
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
! The coding assumes we are using 6-th order centered finite difference for our
! derivatives.
        f(l1-1,:,:,j)= f(l1+1,:,:,j) +  2.*60.*f(l1,:,:,j)*dx/(45.*x(l1))
        f(l1-2,:,:,j)= f(l1+2,:,:,j) +  2.*60.*f(l1,:,:,j)*dx/(9.*x(l1))
        f(l1-3,:,:,j)= f(l1+3,:,:,j) +  2.*60.*f(l1,:,:,j)*dx/x(l1)
      case ('top')               ! top boundary
        f(l2+1,:,:,j)= f(l2-1,:,:,j) -  2.*60.*f(l2,:,:,j)*dx/(45.*x(l2))
        f(l2+2,:,:,j)= f(l2-2,:,:,j) -  2.*60.*f(l2,:,:,j)*dx/(9.*x(l2))
        f(l2+3,:,:,j)= f(l2-3,:,:,j) -  2.*60.*f(l2,:,:,j)*dx/(x(l2))
!
      case default
        call warning('bc_set_pfc_x',topbot//" should be 'top' or 'bot'")
!
      endselect
!
    endsubroutine bc_set_pfc_x
!***********************************************************************
    subroutine bc_set_nfr_x(f,topbot,j)
!
!  Normal-field (or angry-hedgehog) boundary condition for spherical
!  coordinate system.
!  d_r(A_{\theta}) = -A_{\theta}/r  with A_r = 0 sets B_{r} to zero
!  in spherical coordinate system.
!  (compare with next subroutine sfree )
!
!  25-Aug-2007/dhruba: coded
!
      character (len=bclen), intent (in) :: topbot
      real, dimension (:,:,:,:), intent (inout) :: f
      integer, intent (in) :: j
      integer :: k
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        do k=1,nghost
          f(l1-k,:,:,j)= f(l1+k,:,:,j)*(x(l1+k)/(x(l1+k)-dx2_bound(-k)))
        enddo
!
     case ('top')               ! top boundary
       do k=1,nghost
         f(l2+k,:,:,j)= f(l2-k,:,:,j)*(x(l2-k)/(x(l2-k)+dx2_bound(k)))
       enddo
!
      case default
        call warning('bc_set_nfr_x',topbot//" should be 'top' or 'bot'")
!
      endselect
!
    endsubroutine bc_set_nfr_x
!***********************************************************************
    subroutine bc_set_nr1_x(f,topbot,j)
!
!  Normal-field (or angry-hedgehog) boundary condition for spherical
!  coordinate system.
!  d_r(A_{\theta}) = -A_{\theta}/r  with A_r = 0 sets B_{r} to zero
!  in spherical coordinate system.
!  Implementation with one-sided derivative.
!
!  13-Dec-2016/MR: coded
!
      character (len=bclen), intent (in) :: topbot
      real, dimension (:,:,:,:), intent (inout) :: f
      integer, intent (in) :: j
!
      if (topbot=='bot') then
        call bval_from_3rd(f,topbot,j,1,-1./x(l1))
      else
        call bval_from_3rd(f,topbot,j,1,-1./x(l2)) 
      endif
      call set_ghosts_for_onesided_ders(f,topbot,j,1,.true.)
!
    endsubroutine bc_set_nr1_x
!***********************************************************************
    subroutine bc_set_sr1_x(f,topbot,j)
!
!  Stress-free boundary condition for spherical
!  coordinate system: \partial_r u_(\theta,\phi)|_r_(i,a) = u_(r_(i,a),\theta,\phi)/r_(i,a)
!  Implementation with one-sided derivative.
!
!  4-Sep-2017/MR: coded
!
      character (len=bclen), intent (in) :: topbot
      real, dimension (:,:,:,:), intent (inout) :: f
      integer, intent (in) :: j
!
      if (topbot=='bot') then
        call bval_from_3rd(f,topbot,j,1,1./x(l1))
      else
        call bval_from_3rd(f,topbot,j,1,1./x(l2))
      endif
      call set_ghosts_for_onesided_ders(f,topbot,j,1,.true.)
!
    endsubroutine bc_set_sr1_x
! **********************************************************************
    subroutine bc_set_sa2_x(f,topbot,j)
!
!  To set the boundary condition:
!  d_r(r B_{\phi} = 0 we need to se
!  (d_r)^2(r A_{\theta}) = 0 which sets the condition 'a2'
!  on r A_{\theta} and vice-versa for A_{\phi}
!
!  03-Dec-2009/dhruba: coded
!
      character (len=bclen), intent (in) :: topbot
      real, dimension (:,:,:,:), intent (inout) :: f
      integer, intent (in) :: j
      integer :: k
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        do k=1,nghost
          f(l1-k,:,:,j)= f(l1,:,:,j)*2.*(x(l1)/x(l1-k))&
                         -f(l1+k,:,:,j)*(x(l1+k)/x(l1-k))
        enddo
!
     case ('top')               ! top boundary
       do k=1,nghost
         f(l2+k,:,:,j)= f(l2,:,:,j)*2.*(x(l2)/x(l2+k))&
                        -f(l2-k,:,:,j)*(x(l2-k)/x(l2+k))
       enddo
!
      case default
        call warning('bc_set_sa2_x',topbot//" should be 'top' or 'bot'")
!
      endselect
!
    endsubroutine bc_set_sa2_x
! **********************************************************************
    subroutine bc_set_sfree_x(f,topbot,j)
!
!  Details are given in an appendix in the manual.
!  Lambda effect : stresses due to Lambda effect are added to the stress-tensor.
!  For rotation along the z direction and also for not very strong rotation such
!  that the breaking of rotational symmetry is only due to gravity, the only
!  new term appears in the r-phi component. This implies that this term
!  affects only the boundary condition of u_{\phi} for the radial boundary.
!
!  25-Aug-2007/dhruba: coded
!  21-Mar-2009/axel: get llambda_effect using get_shared_variable
!
      use SharedVariables, only : get_shared_variable
!
      character (len=bclen), intent (in) :: topbot
      real, dimension (:,:,:,:), intent (inout) :: f
      integer, intent (in) :: j
!
      real, pointer :: nu,Lambda_V0t,Lambda_V0b,Lambda_V1t,Lambda_V1b
      logical, pointer :: llambda_effect
      integer :: iy, k
      real :: fac,sth,lambda_exp
!
! -------- Either case get the lambda variables first -----------
!
      call get_shared_variable('nu',nu,caller='bc_set_sfree_x')
      call get_shared_variable('llambda_effect',llambda_effect)
      if (llambda_effect) then
         call get_shared_variable('Lambda_V0t',Lambda_V0t)
         call get_shared_variable('Lambda_V1t',Lambda_V1t)
         call get_shared_variable('Lambda_V0b',Lambda_V0b)
         call get_shared_variable('Lambda_V1b',Lambda_V1b)
      endif
!
      select case (topbot)
!
! Bottom boundary
!
      case ('bot')
!
        if ((llambda_effect).and.(j==iuz)) then
          do iy=1,size(f,2)
            sth=sinth(iy)
            lambda_exp=1.+(Lambda_V0b+Lambda_V1b*sth*sth)/nu
            do k=1,nghost
               fac=(1.-dx2_bound(-k)/x(l1+k))**lambda_exp
               if (Omega==0) then
                 f(l1-k,iy,:,j) = f(l1+k,iy,:,j)*fac
               else
                 f(l1-k,iy,:,j) = (f(l1+k,iy,:,j)+Omega*x(l1+k)*sth)*fac &
                                 -Omega*(x(l1+k)-dx2_bound(-k))*sth
               endif
            enddo
          enddo
        else
          do k=1,nghost
            f(l1-k,:,:,j) = f(l1+k,:,:,j)*(1.-dx2_bound(-k)/x(l1+k))
!
!  Alternative formulation
!
            !f(l1-k,:,:,j)= f(l1+k,:,:,j) - f(l1,:,:,j)*(dx2_bound(-k)/x(l1))
          enddo
        endif
!
! Top boundary
!
      case ('top')
        if ((llambda_effect).and.(j==iuz)) then
          do iy=1,size(f,2)
            sth=sinth(iy)
            lambda_exp=1.+(Lambda_V0t+Lambda_V1t*sth*sth)/nu
            do k=1,nghost
              fac=(1.+dx2_bound(k)/x(l2-k))**lambda_exp
              if (Omega==0) then
                f(l2+k,iy,:,j) = f(l2-k,iy,:,j)*fac
              else
                f(l2+k,iy,:,j) = (f(l2-k,iy,:,j)+Omega*x(l2-k)*sth)*fac &
                                -Omega*(x(l2-k)+dx2_bound(k))*sth
              endif
            enddo
          enddo
        else
          do k=1,nghost
            f(l2+k,:,:,j)= f(l2-k,:,:,j)*(1.+dx2_bound(k)/x(l2-k))
!
!  Alternative formulation
!
            !f(l2+k,:,:,j)= f(l2-k,:,:,j) + f(l2,:,:,j)*(dx2_bound(k)/x(l2))
          enddo
        endif
!
      case default
        call warning('bc_set_sfree_x',topbot//" should be 'top' or 'bot'")
!
      endselect
!
    endsubroutine bc_set_sfree_x
! **********************************************************************
    subroutine bc_set_jethat_x(f,jj,topbot,fracall,uzeroall)
!
!  Sets tophat velocity profile at the inner (bot) boundary
!
!  03-jan-2008/dhruba: coded
!
      use Sub, only: step
!
      character (len=bclen), intent (in) :: topbot
      real, dimension (:,:,:,:), intent (inout) :: f
      integer, intent(in) :: jj
      integer :: i,j,k
      real, dimension(:),intent(in) :: fracall,uzeroall
      real :: frac,uzero,ylim,ymid,y1,zlim,zmid,z1
      real :: yhat_min,yhat_max,zhat_min,zhat_max
      real, parameter :: width_hat=0.01
      real, dimension (m2-m1+1) :: hatprofy
      real, dimension (n2-n1+1) :: hatprofz
!
      y1 = xyz1(2)
      z1 = xyz1(3)
      frac = fracall(jj)
      uzero = uzeroall(jj)
!
      if (lspherical_coords)then
!
        select case (topbot)
        case ('bot')               ! bottom boundary
          ylim = (y1-y0)*frac
          ymid = y0+(y1-y0)/2.
          yhat_min=ymid-ylim/2.
          yhat_max=ymid+ylim/2
          hatprofy=step(y(m1:m2),yhat_min,width_hat)*(1.-step(y(m1:m2),yhat_max,width_hat))
          zlim = (z1-z0)*frac
          zmid = z0+(z1-z0)/2.
          zhat_min=zmid-zlim/2.
          zhat_max=zmid+zlim/2
          hatprofz=step(z(n1:n2),zhat_min,width_hat)*(1.-step(z(n1:n2),zhat_max,width_hat))
          do j=m1,m2
            do k=n1,n2
                f(l1,j,k,iux)= uzero*hatprofy(j)*hatprofz(k)
                do i=1,nghost
                  f(l1-i,j,k,iux)= uzero*hatprofy(j)*hatprofz(k)
                enddo
            enddo
          enddo
        case ('top')               ! top boundary
          call warning('bc_set_jethat_x','Jet flowing out of the exit boundary ?')
          do i=1,nghost
            f(l2+i,:,:,jj)=0.
          enddo
!
        case default
          call warning('bc_set_jethat_x',topbot//" should be 'top' or 'bot'")
        endselect

     else
        call stop_it('Boundary condition jethat is valid only in spherical coordinate system')
     endif
!
    endsubroutine bc_set_jethat_x
! **********************************************************************
    subroutine bc_set_jet_x(f,jj,topbot,velocity,radius)
!
!  Sets tophat velocity profile at the inner (bot) boundary
!
!  06-nov-2013/nils: adapted from bc_set_jethat_x. Made this new routine
!                    because there some awckward choices made in the
!                    other one, and the other one is for spherical geometries.
!
      use Sub, only: step
!
      character (len=bclen), intent (in) :: topbot
      real, dimension (m2-m1+1,n2-n1+1) :: prof
      real, dimension (:,:,:,:), intent (inout) :: f
      integer, intent(in) :: jj
      integer :: i,j,k
      real, dimension(:),intent(in) :: velocity,radius
      real :: vel,rad
!
      vel = velocity(jj)
      rad = radius(jj)
!
      if (lcartesian_coords) then
         select case (topbot)
         case ('bot')               ! bottom boundary
            call jet_x(prof,vel,rad)
            do j=m1,m2
               do k=n1,n2
                  f(l1,j,k,iux)= prof(j-nghost,k-nghost)
                  do i=1,nghost
                     f(l1-i,j,k,iux)= prof(j-nghost,k-nghost)
                  enddo
               enddo
            enddo
!
         case ('top')               ! top boundary
            call warning('bc_set_jet_x',&
                 'Jet flowing out of the exit boundary ?')
            do i=1,nghost
               f(l2+i,:,:,jj)=0.
            enddo
!
         case default
            call warning('bc_set_jethat_x',topbot//" should be 'top' or 'bot'")
         endselect
!
      else
         call stop_it('Boundary condition jethat is valid only in spherical coordinate system')
      endif
!
    endsubroutine bc_set_jet_x
! **********************************************************************
    subroutine jet_x(prof,vel,rad)
!
!  06-nov-2013/nils: Set jet profile
!
      use Sub, only: step
!
      real, dimension (:,:), intent (out) :: prof
      integer :: j,k
      real :: vel,rad,ymid,y1,zlim,zmid,z1
      real :: yhat_min,yhat_max,zhat_min,zhat_max
      real :: width_hat
      real, dimension (size(prof,1)) :: hatprofy
      real, dimension (size(prof,2)) :: hatprofz
!
      y1 = xyz1(2)
      z1 = xyz1(3)
      width_hat=dy*2
!
      ymid = y0+(y1-y0)/2.
      yhat_min=ymid-rad/2.
      yhat_max=ymid+rad/2
      hatprofy= step(y(m1:m2),yhat_min,width_hat) &
               -step(y(m1:m2),yhat_max,width_hat)
!
      if (nzgrid>1) then
         zlim = (z1-z0)*rad
         zmid = z0+(z1-z0)/2.
         zhat_min=zmid-zlim/2.
         zhat_max=zmid+zlim/2
         hatprofz= step(z(n1:n2),zhat_min,width_hat) &
                  -step(z(n1:n2),zhat_max,width_hat)
      endif
!
      do j=1,size(prof,1)
         do k=1,size(prof,2)
            if (nzgrid>1) then
               prof(j,k)= vel*hatprofy(j)*hatprofz(k)
            else
               prof(j,k)= vel*hatprofy(j)
            endif
         enddo
      enddo
!
    endsubroutine jet_x
! **********************************************************************
    subroutine bc_set_nfr_y(f,topbot,j)
!
!  Normal field boundary condition for spherical coordinate system.
!  d_{\theta}(A_{\phi}) = -A_{\phi}cot(\theta)/r with A_{\theta} = 0 sets
!  B_r = B_{\phi} = 0 in spherical polar coordinate system. This subroutine
!  sets only the first part of this boundary condition for 'j'-th component
!  of f.
!
!  25-Aug-2007/dhruba: coded
!
      character (len=bclen), intent (in) :: topbot
      real, dimension (:,:,:,:), intent (inout) :: f
      integer, intent (in) :: j
      integer :: k
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        do k=1,nghost
          f(:,m1-k,:,j)= f(:,m1+k,:,j)*(sinth(m1+k)/sin(y(m1+k)-dy2_bound(-k)))
        enddo
      case ('top')               ! top boundary
        do k=1,nghost
          f(:,m2+k,:,j)= f(:,m2-k,:,j)*(sinth(m2-k)/sin(y(m2-k)+dy2_bound(k)))
        enddo
!
      case default
        call warning('bc_set_nfr_y',topbot//" should be 'top' or 'bot'")
!
      endselect
!
    endsubroutine bc_set_nfr_y
! **********************************************************************
    subroutine bc_set_sfree_y(f,topbot,j)
!
!  Stress-free boundary condition for spherical coordinate system.
!  d_{\theta}(u_{\phi}) = u_{\phi}cot(\theta)  with u_{\theta} = 0 sets
!  S_{\theta \phi} component of the strain matrix to be zero in spherical
!  coordinate system. This subroutine sets only the first part of this
!  boundary condition for 'j'-th component of f.
!
!  25-Aug-2007/dhruba: coded
!
      use SharedVariables, only : get_shared_variable
!
      character (len=bclen), intent (in) :: topbot
      real, dimension (:,:,:,:), intent (inout) :: f
      integer, intent (in) :: j
      real, pointer :: Lambda_H1,nu
      real, pointer :: LH1_rprof(:)
      logical, pointer :: llambda_effect
      integer :: k,ix
      real :: cos2thm_k,cos2thmpk,somega
      real,dimension(size(f,1)):: LH1
!
! -------- Either case get the lambda variables first -----------
!
      call get_shared_variable('nu',nu,caller='bc_set_sfree_y')
      call get_shared_variable('llambda_effect',llambda_effect)
      if (llambda_effect) then
         call get_shared_variable('Lambda_H1',Lambda_H1)
         call get_shared_variable('LH1_rprof',LH1_rprof)
         LH1=Lambda_H1*LH1_rprof
      endif
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        if ((llambda_effect).and.(j==iuz).and.(Lambda_H1/=0.)) then
          do k=1,nghost
              cos2thm_k= costh(m1-k)**2-sinth(m1-k)**2
              cos2thmpk= costh(m1+k)**2-sinth(m1+k)**2
            if (Omega==0) then
               do ix=1,size(f,1)
                  f(ix,m1-k,:,j)= f(ix,m1+k,:,j)* &
                       (exp(LH1(ix)*cos2thmpk/(4.*nu))*sin1th(m1+k)) &
                       *(exp(-LH1(ix)*cos2thm_k/(4.*nu))*sinth(m1-k))
               enddo
            else
              do ix=1,size(f,1)
! DM+GG: temporally commented out
!                somega=x(ix)*Omega*sinth(m1-k)*( &
!                   exp(2*cos2thm_k*LH1(ix)/(4.*nu))&
!                        -exp((cos2thmpk+cos2thm_k)*LH1(ix)/(4.*nu)) )
                somega=x(ix)*Omega*sinth(m1-k)*( &
                   exp(cos2thmpk*LH1(ix)/(4.*nu))&
                        /exp((cos2thm_k)*LH1(ix)/(4.*nu)) -1.)
                f(ix,m1-k,:,j)= f(ix,m1+k,:,j)* &
                   (exp(LH1(ix)*cos2thmpk/(4.*nu))*sin1th(m1+k)) &
                   *(exp(-LH1(ix)*cos2thm_k/(4.*nu))*sinth(m1-k)) &
                      +somega
              enddo
            endif
          enddo
        else
          do k=1,nghost
            f(:,m1-k,:,j)= f(:,m1+k,:,j)*(sin(y(m1+k)-dy2_bound(-k))*sin1th(m1+k))
          enddo
        endif
      case ('top')               ! top boundary
        if ((llambda_effect).and.(j==iuz).and.(Lambda_H1/=0)) then
          do k=1,nghost
            cos2thm_k= costh(m2-k)**2-sinth(m2-k)**2
            cos2thmpk= costh(m2+k)**2-sinth(m2+k)**2
            if (Omega==0)then
               do ix=1,size(f,1)
                  f(ix,m2+k,:,j)= f(ix,m2-k,:,j)* &
                   (exp(LH1(ix)*cos2thm_k/(4.*nu))*sin1th(m2-k)) &
                  *(exp(-LH1(ix)*cos2thmpk/(4.*nu))*sinth(m2+k))
               enddo
             else
              do ix=1,size(f,1)
! DM+GG: Temporally comented out
!                somega=x(ix)*Omega*sinth(m2+k)*( &
!                   exp(2*cos2thmpk*LH1(ix)/(4.*nu))&
!                        -exp((cos2thmpk+cos2thm_k)*LH1(ix)/(4.*nu)) )
                 somega=x(ix)*Omega*sinth(m2+k)*( &
                      exp(cos2thm_k*LH1(ix)/(4.*nu))    &
                      / exp(cos2thmpk*LH1(ix)/(4.*nu))-1.)
                f(ix,m2+k,:,j)= f(ix,m2-k,:,j)* &
                     (exp(LH1(ix)*cos2thm_k/(4.*nu))*sin1th(m2-k)) &
                    *(exp(-LH1(ix)*cos2thmpk/(4.*nu))*sinth(m2+k)) &
                      +somega
              enddo
             endif
          enddo
        else
          do k=1,nghost
            f(:,m2+k,:,j)= f(:,m2-k,:,j)*(sin(y(m2-k)+dy2_bound(k))*sin1th(m2-k))
          enddo
        endif
!
     case default
        call warning('bc_set_sfree_y',topbot//" should be 'top' or 'bot'")
!
      endselect
!
    endsubroutine bc_set_sfree_y
! **********************************************************************
    subroutine bc_set_pfc_y(f,topbot,j)
!
!joern: WARNING, this bc will NOT give a perfect-conductor boundary condition
!
!  In spherical polar coordinate system,
!  at a theta boundary set : $A_{r} = 0$ and $A_{\phi} = 0$,
!  and demand $div A = 0$ gives the condition on $A_{\theta}$ to be
!  $d/d{\theta}( A_{\theta}) + \cot(\theta)A_{\theta} = 0$ .
!  This subroutine sets this condition on
!  $j$ the component of f. As this is related to setting the
!  perfect conducting boundary condition we call this "pfc".
!
!  25-Aug-2007/dhruba: coded
!
      character (len=bclen), intent (in) :: topbot
      real, dimension (:,:,:,:), intent (inout) :: f
      integer, intent (in) :: j
      real :: cottheta
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
!
!  The coding assumes we are using 6-th order centered finite difference for our
!  derivatives.
!
        cottheta= cotth(m1)
        f(:,m1-1,:,j)= f(:,m1+1,:,j) +  60.*dy*cottheta*f(:,m1,:,j)/45.
        f(:,m1-2,:,j)= f(:,m1+2,:,j) -  60.*dy*cottheta*f(:,m1,:,j)/9.
        f(:,m1-3,:,j)= f(:,m1+3,:,j) +  60.*dy*cottheta*f(:,m1,:,j)
      case ('top')               ! top boundary
        cottheta= cotth(m2)
        f(:,m2+1,:,j)= f(:,m2-1,:,j) -  60.*dy*cottheta*f(:,m2,:,j)/45.
        f(:,m2+2,:,j)= f(:,m2-2,:,j) +  60.*dy*cottheta*f(:,m2,:,j)/9.
        f(:,m2+3,:,j)= f(:,m2-3,:,j) -  60.*dy*cottheta*f(:,m2,:,j)
!
      case default
        call warning('bc_set_pfc_y',topbot//" should be 'top' or 'bot'")
!
      endselect
!
    endsubroutine bc_set_pfc_y
!***********************************************************************
    subroutine bc_set_der_y(f,topbot,j,val)
!
!  Sets the derivative on the boundary to a given value.
!
!  14-may-2006/tobi: coded
!
      character (len=bclen), intent (in) :: topbot
      real, dimension (:,:,:,:), intent (inout) :: f
      integer, intent (in) :: j
      real, intent (in) :: val
!
      integer :: i
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        do i=1,nghost; f(:,m1-i,:,j) = f(:,m1+i,:,j) - dy2_bound(-i)*val; enddo
!
      case ('top')               ! top boundary
        do i=1,nghost; f(:,m2+i,:,j) = f(:,m2-i,:,j) + dy2_bound(i)*val; enddo
!
      case default
        call warning('bc_set_der_y',topbot//" should be 'top' or 'bot'")
!
      endselect
!
    endsubroutine bc_set_der_y
!***********************************************************************
    subroutine bc_set_der_z(f,topbot,j,val)
!
!  Sets the derivative on the boundary to a given value.
!
!  14-may-2006/tobi: coded
!
      character (len=bclen), intent (in) :: topbot
      real, dimension (:,:,:,:), intent (inout) :: f
      integer, intent (in) :: j
      real, intent (in) :: val
!
      integer :: i
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        do i=1,nghost; f(:,:,n1-i,j) = f(:,:,n1+i,j) - dz2_bound(-i)*val; enddo
!
      case ('top')               ! top boundary
        do i=1,nghost; f(:,:,n2+i,j) = f(:,:,n2-i,j) + dz2_bound(i)*val; enddo
!
      case default
        call warning('bc_set_der_z',topbot//" should be 'top' or 'bot'")
!
      endselect
!
    endsubroutine bc_set_der_z
!***********************************************************************
    subroutine bc_set_div_z(f,topbot,j,val)
!
!  Sets the derivative on the boundary to a given value
!
!  17-may-2010/bing: coded
!
      character (len=bclen), intent (in) :: topbot
      real, dimension (:,:,:,:), intent (inout) :: f
      real, dimension (l2-l1+1,m2-m1+1) :: fac,duz_dz
      real, intent(in) :: val
!
      integer, intent (in) :: j
!
      integer :: iref=-1,pos,nxl,nyl
!
      if (j/=iuz) call fatal_error_local('bc_set_div_z','please set div for uz only')
!
      nxl=l2-l1+1; nyl=m2-m1+1
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        iref = n1
!
      case ('top')               ! top boundary
        iref = n2
!
      case default
        call warning('bc_set_der_x',topbot//" should be 'top' or 'bot'")
!
      endselect
!
! take the x derivative of ux
      if (nxgrid/=1) then
        fac=(1./60)*spread(dx_1(l1:l2),2,nyl)
        duz_dz= fac*(+45.0*(f(l1+1:l2+1,m1:m2,iref,iux)-f(l1-1:l2-1,m1:m2,iref,iux)) &
                     - 9.0*(f(l1+2:l2+2,m1:m2,iref,iux)-f(l1-2:l2-2,m1:m2,iref,iux)) &
                     +     (f(l1+3:l2+3,m1:m2,iref,iux)-f(l1-3:l2-3,m1:m2,iref,iux)))
      else
        if (ip<=5) print*, 'bc_set_div_z: Degenerate case in x-direction'
      endif
!
! take the y derivative of uy and add to dux/dx
      if (nygrid/=1) then
        fac=(1./60)*spread(dy_1(m1:m2),1,nxl)
        duz_dz=duz_dz + fac*(+45.0*(f(l1:l2,m1+1:m2+1,iref,iuy)-f(l1:l2,m1-1:m2-1,iref,iuy)) &
                             - 9.0*(f(l1:l2,m1+2:m2+2,iref,iuy)-f(l1:l2,m1-2:m2-2,iref,iuy)) &
                             +     (f(l1:l2,m1+3:m2+3,iref,iuy)-f(l1:l2,m1-3:m2-3,iref,iuy)))
      else
        if (ip<=5) print*, 'bc_set_div_z: Degenerate case in y-direction'
      endif
!
! add given number to set div(u)=val; default val=0
! duz/dz = val - dux/dx - duy/dy
!
      duz_dz = val - duz_dz
!
! set the derivative of uz at the boundary
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        do pos=1,nghost
          f(l1:l2,m1:m2,n1-pos,j) = f(l1:l2,m1:m2,n1+pos,j) - dz2_bound(-pos)*duz_dz
        enddo
!
      case ('top')               ! top boundary
        do pos=1,nghost
          f(l1:l2,m1:m2,n2+pos,j) = f(l1:l2,m1:m2,n2-pos,j) + dz2_bound(pos)*duz_dz
        enddo
!
      case default
        call warning('bc_set_div_z',topbot//" should be 'top' or 'bot'")
!
      endselect
!
    endsubroutine bc_set_div_z
!***********************************************************************
    subroutine bc_van_x(f,topbot,j)
!
!  Vanishing boundary conditions.
!  (TODO: clarify what this means)
!
!  26-apr-06/tobi: coded
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: i,j
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
          do i=1,nghost
            f(l1-i,:,:,j)=((nghost+1-i)*f(l1,:,:,j))/(nghost+1)
          enddo
!
      case ('top')               ! top boundary
          do i=1,nghost
            f(l2+i,:,:,j)=((nghost+1-i)*f(l2,:,:,j))/(nghost+1)
          enddo
!
      case default
        print*, "bc_van_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_van_x
!***********************************************************************
    subroutine bc_van_y(f,topbot,j)
!
!  Vanishing boundary conditions.
!  (TODO: clarify what this means)
!
!  26-apr-06/tobi: coded
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: i,j
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
          do i=1,nghost
            f(:,m1-i,:,j)=((nghost+1-i)*f(:,m1,:,j))/(nghost+1)
          enddo
!
      case ('top')               ! top boundary
          do i=1,nghost
            f(:,m2+i,:,j)=((nghost+1-i)*f(:,m2,:,j))/(nghost+1)
          enddo
!
      case default
        print*, "bc_van_y: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_van_y
!***********************************************************************
    subroutine bc_van_z(f,topbot,j)
!
!  Vanishing boundary conditions.
!  (TODO: clarify what this means)
!
!  26-apr-06/tobi: coded
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: i,j
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
          do i=1,nghost
            f(:,:,n1-i,j)=((nghost+1-i)*f(:,:,n1,j))/(nghost+1)
          enddo
!
      case ('top')               ! top boundary
          do i=1,nghost
            f(:,:,n2+i,j)=((nghost+1-i)*f(:,:,n2,j))/(nghost+1)
          enddo
!
      case default
        print*, "bc_van_z: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_van_z
!***********************************************************************
    subroutine bc_van3rd_y(f,topbot,j)
!
!  Boundary condition with vanishing 3rd derivative
!  (useful for vertical hydrostatic equilibrium in discs)
!
!  30-jul-13/wlad: copied from z
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j
!
      real, dimension (size(f,1),size(f,3)) :: cpoly0,cpoly1,cpoly2
      integer :: i
!
      select case (topbot)
!
      case ('bot')
        cpoly0(:,:)=alog(f(:,m1,:,j))
        cpoly1(:,:)=-(3*alog(f(:,m1,:,j))-4*alog(f(:,m1+1,:,j))+alog(f(:,m1+2,:,j)))/(2*dy)
        cpoly2(:,:)=-(-alog(f(:,m1,:,j))+2*alog(f(:,m1+1,:,j))-alog(f(:,m1+2,:,j))) /(2*dy**2)
        do i=1,nghost
          f(:,m1-i,:,j) = exp(cpoly0(:,:) - cpoly1(:,:)*i*dy + cpoly2(:,:)*(i*dy)**2)
        enddo
!
      case ('top')
        cpoly0(:,:)=alog(f(:,m2,:,j))
        cpoly1(:,:)=-(-3*alog(f(:,m2,:,j))+4*alog(f(:,m2-1,:,j))-alog(f(:,m2-2,:,j)))/(2*dy)
        cpoly2(:,:)=-(-alog(f(:,m2,:,j))+2*alog(f(:,m2-1,:,j))-alog(f(:,m2-2,:,j)))/(2*dy**2)
        do i=1,nghost
          f(:,m2+i,:,j) = exp(cpoly0(:,:) + cpoly1(:,:)*i*dy + cpoly2(:,:)*(i*dy)**2)
        enddo
!
      endselect
!
    endsubroutine bc_van3rd_y
!***********************************************************************
    subroutine bc_van3rd_z(f,topbot,j)
!
!  Boundary condition with vanishing 3rd derivative
!  (useful for vertical hydrostatic equilibrium in discs)
!
!  19-aug-03/anders: coded
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j
!
      real, dimension (size(f,1),size(f,2)) :: cpoly0,cpoly1,cpoly2
      integer :: i
!
      select case (topbot)
!
      case ('bot')
        cpoly0=f(:,:,n1,j)
        cpoly1=-(3*f(:,:,n1,j)-4*f(:,:,n1+1,j)+f(:,:,n1+2,j))/(2*dz)
        cpoly2=-(-f(:,:,n1,j)+2*f(:,:,n1+1,j)-f(:,:,n1+2,j)) /(2*dz**2)
        do i=1,nghost
          f(:,:,n1-i,j) = cpoly0 - cpoly1*i*dz + cpoly2*(i*dz)**2
        enddo
!
      case ('top')
        cpoly0=f(:,:,n2,j)
        cpoly1=-(-3*f(:,:,n2,j)+4*f(:,:,n2-1,j)-f(:,:,n2-2,j))/(2*dz)
        cpoly2=-(-f(:,:,n2,j)+2*f(:,:,n2-1,j)-f(:,:,n2-2,j))/(2*dz**2)
        do i=1,nghost
          f(:,:,n2+i,j) = cpoly0 + cpoly1*i*dz + cpoly2*(i*dz)**2
        enddo
!
      endselect
!
    endsubroutine bc_van3rd_z
!***********************************************************************
    subroutine bc_onesided_x_old(f,topbot,j)
!
!  One-sided conditions.
!  These expressions result from combining Eqs(207)-(210), astro-ph/0109497,
!  corresponding to (9.207)-(9.210) in Ferriz-Mas proceedings.
!
!  05-apr-03/axel: coded
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: i,j,k
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        do i=1,nghost
          k=l1-i
          f(k,:,:,j)=7*f(k+1,:,:,j) &
                   -21*f(k+2,:,:,j) &
                   +35*f(k+3,:,:,j) &
                   -35*f(k+4,:,:,j) &
                   +21*f(k+5,:,:,j) &
                    -7*f(k+6,:,:,j) &
                      +f(k+7,:,:,j)
        enddo
!
      case ('top')               ! top boundary
        do i=1,nghost
          k=l2+i
          f(k,:,:,j)=7*f(k-1,:,:,j) &
                   -21*f(k-2,:,:,j) &
                   +35*f(k-3,:,:,j) &
                   -35*f(k-4,:,:,j) &
                   +21*f(k-5,:,:,j) &
                    -7*f(k-6,:,:,j) &
                      +f(k-7,:,:,j)
        enddo
!
      case default
        print*, "bc_onesided_x_old ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_onesided_x_old
!***********************************************************************
    subroutine bc_onesided_z_orig(f,topbot,j)
!
!  One-sided conditions.
!  These expressions result from combining Eqs(207)-(210), astro-ph/0109497,
!  corresponding to (9.207)-(9.210) in Ferriz-Mas proceedings.
!
!  05-apr-03/axel: coded
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: i,j,k
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        do i=1,nghost
          k=n1-i
          f(:,:,k,j)=7*f(:,:,k+1,j) &
                   -21*f(:,:,k+2,j) &
                   +35*f(:,:,k+3,j) &
                   -35*f(:,:,k+4,j) &
                   +21*f(:,:,k+5,j) &
                    -7*f(:,:,k+6,j) &
                      +f(:,:,k+7,j)
        enddo
!
      case ('top')               ! top boundary
        do i=1,nghost
          k=n2+i
          f(:,:,k,j)=7*f(:,:,k-1,j) &
                   -21*f(:,:,k-2,j) &
                   +35*f(:,:,k-3,j) &
                   -35*f(:,:,k-4,j) &
                   +21*f(:,:,k-5,j) &
                    -7*f(:,:,k-6,j) &
                      +f(:,:,k-7,j)
        enddo
!
      case default
        print*, "bc_onesided_z ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_onesided_z_orig
!***********************************************************************
    subroutine bc_extrap_2_1(f,topbot,j)
!
!  Extrapolation boundary condition.
!  Correct for polynomials up to 2nd order, determined 1 further degree
!  of freedom by minimizing L2 norm of coefficient vector.
!
!  19-jun-03/wolf: coded
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        f(:,:,n1-1,j)=0.25*(  9*f(:,:,n1,j)- 3*f(:,:,n1+1,j)- 5*f(:,:,n1+2,j)+ 3*f(:,:,n1+3,j))
        f(:,:,n1-2,j)=0.05*( 81*f(:,:,n1,j)-43*f(:,:,n1+1,j)-57*f(:,:,n1+2,j)+39*f(:,:,n1+3,j))
        f(:,:,n1-3,j)=0.05*(127*f(:,:,n1,j)-81*f(:,:,n1+1,j)-99*f(:,:,n1+2,j)+73*f(:,:,n1+3,j))
!
      case ('top')               ! top boundary
        f(:,:,n2+1,j)=0.25*(  9*f(:,:,n2,j)- 3*f(:,:,n2-1,j)- 5*f(:,:,n2-2,j)+ 3*f(:,:,n2-3,j))
        f(:,:,n2+2,j)=0.05*( 81*f(:,:,n2,j)-43*f(:,:,n2-1,j)-57*f(:,:,n2-2,j)+39*f(:,:,n2-3,j))
        f(:,:,n2+3,j)=0.05*(127*f(:,:,n2,j)-81*f(:,:,n2-1,j)-99*f(:,:,n2-2,j)+73*f(:,:,n2-3,j))
!
      case default
        print*, "bc_extrap_2_1: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_extrap_2_1
!***********************************************************************
    subroutine bcx_extrap_2_1(f,topbot,j)
!
!  Extrapolation boundary condition for x.
!  Correct for polynomials up to 2nd order, determined 1 further degree
!  of freedom by minimizing L2 norm of coefficient vector.
!
!  19-jun-03/wolf: coded
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        f(l1-1,:,:,j)=0.25*(  9*f(l1,:,:,j)- 3*f(l1+1,:,:,j)- 5*f(l1+2,:,:,j)+ 3*f(l1+3,:,:,j))
        f(l1-2,:,:,j)=0.05*( 81*f(l1,:,:,j)-43*f(l1+1,:,:,j)-57*f(l1+2,:,:,j)+39*f(l1+3,:,:,j))
        f(l1-3,:,:,j)=0.05*(127*f(l1,:,:,j)-81*f(l1+1,:,:,j)-99*f(l1+2,:,:,j)+73*f(l1+3,:,:,j))
!
      case ('top')               ! top boundary
        f(l2+1,:,:,j)=0.25*(  9*f(l2,:,:,j)- 3*f(l2-1,:,:,j)- 5*f(l2-2,:,:,j)+ 3*f(l2-3,:,:,j))
        f(l2+2,:,:,j)=0.05*( 81*f(l2,:,:,j)-43*f(l2-1,:,:,j)-57*f(l2-2,:,:,j)+39*f(l2-3,:,:,j))
        f(l2+3,:,:,j)=0.05*(127*f(l2,:,:,j)-81*f(l2-1,:,:,j)-99*f(l2-2,:,:,j)+73*f(l2-3,:,:,j))
!
      case default
        print*, "bcx_extrap_2_1: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bcx_extrap_2_1
!***********************************************************************
    subroutine bcy_extrap_2_1(f,topbot,j)
!
!  Extrapolation boundary condition for y.
!  Correct for polynomials up to 2nd order, determined 1 further degree
!  of freedom by minimizing L2 norm of coefficient vector.
!
!  19-jun-03/wolf: coded
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        f(:,m1-1,:,j)=0.25*(  9*f(:,m1,:,j)- 3*f(:,m1+1,:,j)- 5*f(:,m1+2,:,j)+ 3*f(:,m1+3,:,j))
        f(:,m1-2,:,j)=0.05*( 81*f(:,m1,:,j)-43*f(:,m1+1,:,j)-57*f(:,m1+2,:,j)+39*f(:,m1+3,:,j))
        f(:,m1-3,:,j)=0.05*(127*f(:,m1,:,j)-81*f(:,m1+1,:,j)-99*f(:,m1+2,:,j)+73*f(:,m1+3,:,j))
!
      case ('top')               ! top boundary
        f(:,m2+1,:,j)=0.25*(  9*f(:,m2,:,j)- 3*f(:,m2-1,:,j)- 5*f(:,m2-2,:,j)+ 3*f(:,m2-3,:,j))
        f(:,m2+2,:,j)=0.05*( 81*f(:,m2,:,j)-43*f(:,m2-1,:,j)-57*f(:,m2-2,:,j)+39*f(:,m2-3,:,j))
        f(:,m2+3,:,j)=0.05*(127*f(:,m2,:,j)-81*f(:,m2-1,:,j)-99*f(:,m2-2,:,j)+73*f(:,m2-3,:,j))
!
      case default
        print*, "bcy_extrap_2_1: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bcy_extrap_2_1
!***********************************************************************
    subroutine bc_extrap_2_2(f,topbot,j)
!
!  Extrapolation boundary condition.
!  Correct for polynomials up to 2nd order, determined 2 further degrees
!  of freedom by minimizing L2 norm of coefficient vector.
!
!  19-jun-03/wolf: coded
!  01-jul-03/axel: introduced abbreviations n1p4,n2m4
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j,n1p4,n2m4
!
!  abbreviations, because otherwise the ifc compiler complains
!  for 1-D runs without vertical extent
!
      n1p4=n1+4
      n2m4=n2-4
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        f(:,:,n1-1,j)=0.2   *(  9*f(:,:,n1,j)                 -  4*f(:,:,n1+2,j)- 3*f(:,:,n1+3,j)+ 3*f(:,:,n1p4,j))
        f(:,:,n1-2,j)=0.2   *( 15*f(:,:,n1,j)- 2*f(:,:,n1+1,j)-  9*f(:,:,n1+2,j)- 6*f(:,:,n1+3,j)+ 7*f(:,:,n1p4,j))
        f(:,:,n1-3,j)=1./35.*(157*f(:,:,n1,j)-33*f(:,:,n1+1,j)-108*f(:,:,n1+2,j)-68*f(:,:,n1+3,j)+87*f(:,:,n1p4,j))
!
      case ('top')               ! top boundary
        f(:,:,n2+1,j)=0.2   *(  9*f(:,:,n2,j)                 -  4*f(:,:,n2-2,j)- 3*f(:,:,n2-3,j)+ 3*f(:,:,n2m4,j))
        f(:,:,n2+2,j)=0.2   *( 15*f(:,:,n2,j)- 2*f(:,:,n2-1,j)-  9*f(:,:,n2-2,j)- 6*f(:,:,n2-3,j)+ 7*f(:,:,n2m4,j))
        f(:,:,n2+3,j)=1./35.*(157*f(:,:,n2,j)-33*f(:,:,n2-1,j)-108*f(:,:,n2-2,j)-68*f(:,:,n2-3,j)+87*f(:,:,n2m4,j))
!
      case default
        print*, "bc_extrap_2_2: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_extrap_2_2
!***********************************************************************
    subroutine bcx_extrap_2_2(f,topbot,j)
!
!  Extrapolation boundary condition.
!  Correct for polynomials up to 2nd order, determined 2 further degrees
!  of freedom by minimizing L2 norm of coefficient vector.
!
!  19-jun-03/wolf: coded
!  01-jul-03/axel: introduced abbreviations n1p4,n2m4
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j,l1p4,l2m4
!
!  abbreviations, because otherwise the ifc compiler complains
!  for 1-D runs without vertical extent
!
      l1p4=l1+4
      l2m4=l2-4
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        f(l1-1,:,:,j)=0.2   *(  9*f(l1,:,:,j)                 -  4*f(l1+2,:,:,j)- 3*f(l1+3,:,:,j)+ 3*f(l1p4,:,:,j))
        f(l1-2,:,:,j)=0.2   *( 15*f(l1,:,:,j)- 2*f(l1+1,:,:,j)-  9*f(l1+2,:,:,j)- 6*f(l1+3,:,:,j)+ 7*f(l1p4,:,:,j))
        f(l1-3,:,:,j)=1./35.*(157*f(l1,:,:,j)-33*f(l1+1,:,:,j)-108*f(l1+2,:,:,j)-68*f(l1+3,:,:,j)+87*f(l1p4,:,:,j))
!
      case ('top')               ! top boundary
        f(l2+1,:,:,j)=0.2   *(  9*f(l2,:,:,j)                 -  4*f(l2-2,:,:,j)- 3*f(l2-3,:,:,j)+ 3*f(l2m4,:,:,j))
        f(l2+2,:,:,j)=0.2   *( 15*f(l2,:,:,j)- 2*f(l2-1,:,:,j)-  9*f(l2-2,:,:,j)- 6*f(l2-3,:,:,j)+ 7*f(l2m4,:,:,j))
        f(l2+3,:,:,j)=1./35.*(157*f(l2,:,:,j)-33*f(l2-1,:,:,j)-108*f(l2-2,:,:,j)-68*f(l2-3,:,:,j)+87*f(l2m4,:,:,j))
!
      case default
        print*, "bcx_extrap_2_2: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bcx_extrap_2_2
!***********************************************************************
    subroutine bcy_extrap_2_2(f,topbot,j)
!
!  Extrapolation boundary condition.
!  Correct for polynomials up to 2nd order, determined 2 further degrees
!  of freedom by minimizing L2 norm of coefficient vector.
!
!   19-jun-03/wolf: coded
!   01-jul-03/axel: introduced abbreviations n1p4,n2m4
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j,m1p4,m2m4
!
!  abbreviations, because otherwise the ifc compiler complains
!  for 1-D runs without vertical extent
!
      m1p4=m1+4
      m2m4=m2-4
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        f(:,m1-1,:,j)=0.2   *(  9*f(:,m1,:,j)                 -  4*f(:,m1+2,:,j)- 3*f(:,m1+3,:,j)+ 3*f(:,m1p4,:,j))
        f(:,m1-2,:,j)=0.2   *( 15*f(:,m1,:,j)- 2*f(:,m1+1,:,j)-  9*f(:,m1+2,:,j)- 6*f(:,m1+3,:,j)+ 7*f(:,m1p4,:,j))
        f(:,m1-3,:,j)=1./35.*(157*f(:,m1,:,j)-33*f(:,m1+1,:,j)-108*f(:,m1+2,:,j)-68*f(:,m1+3,:,j)+87*f(:,m1p4,:,j))
!
      case ('top')               ! top boundary
        f(:,m2+1,:,j)=0.2   *(  9*f(:,m2,:,j)                 -  4*f(:,m2-2,:,j)- 3*f(:,m2-3,:,j)+ 3*f(:,m2m4,:,j))
        f(:,m2+2,:,j)=0.2   *( 15*f(:,m2,:,j)- 2*f(:,m2-1,:,j)-  9*f(:,m2-2,:,j)- 6*f(:,m2-3,:,j)+ 7*f(:,m2m4,:,j))
        f(:,m2+3,:,j)=1./35.*(157*f(:,m2,:,j)-33*f(:,m2-1,:,j)-108*f(:,m2-2,:,j)-68*f(:,m2-3,:,j)+87*f(:,m2m4,:,j))
!
      case default
        print*, "bcy_extrap_2_2: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bcy_extrap_2_2
!***********************************************************************
    subroutine bcy_extrap_2_3(f,topbot,j)
!
!  Extrapolation boundary condition in logarithm:
!  It maintains a power law
!
!  18-dec-08/wlad: coded
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j,l,i
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        do i=1,nghost
          do n=1,size(f,3)
            do l=1,size(f,1)
              if (f(l,m1+i,n,j)/=0.) then
                f(l,m1-i,n,j)=f(l,m1,n,j)**2/f(l,m1+i,n,j)
              else
                f(l,m1-i,n,j)=0.
              endif
            enddo
          enddo
        enddo
!
      case ('top')               ! top boundary
        do i=1,nghost
          do n=1,size(f,3)
            do l=1,size(f,1)
              if (f(l,m2-i,n,j)/=0.) then
                f(l,m2+i,n,j)=f(l,m2,n,j)**2/f(l,m2-i,n,j)
              else
                f(l,m2+i,n,j)=0.
              endif
            enddo
          enddo
        enddo
!
      case default
        print*, "bcy_extrap_2_3: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bcy_extrap_2_3
!***********************************************************************
    subroutine bc_extrap0_2_0(f,topbot,j)
!
!  Extrapolation boundary condition for f(bdry)=0.
!  Correct for polynomials up to 2nd order, determined no further degree
!  of freedom by minimizing L2 norm of coefficient vector.
!
!  09-oct-03/wolf: coded
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j
!
      select case (topbot)
!
! Nyquist-filtering
!
      case ('bot')
! bottom boundary
        f(:,:,n1  ,j)=0.        ! set bdry value=0 (indep of initcond)
        f(:,:,n1-1,j)=(1/11.)*&
             (-17*f(:,:,n1+1,j)- 9*f(:,:,n1+2,j)+ 8*f(:,:,n1+3,j))
        f(:,:,n1-2,j)=      2*&
             (- 2*f(:,:,n1+1,j)-   f(:,:,n1+2,j)+   f(:,:,n1+3,j))
        f(:,:,n1-3,j)=(3/11.)*&
             (-27*f(:,:,n1+1,j)-13*f(:,:,n1+2,j)+14*f(:,:,n1+3,j))
!
      case ('top')
! top boundary
        f(:,:,n2  ,j)=0.        ! set bdry value=0 (indep of initcond)
        f(:,:,n2+1,j)=(1/11.)*&
             (-17*f(:,:,n2-1,j)- 9*f(:,:,n2-2,j)+ 8*f(:,:,n2-3,j))
        f(:,:,n2+2,j)=      2*&
             (- 2*f(:,:,n2-1,j)-   f(:,:,n2-2,j)+   f(:,:,n2-3,j))
        f(:,:,n2+3,j)=(3/11.)*&
             (-27*f(:,:,n2-1,j)-13*f(:,:,n2-2,j)+14*f(:,:,n2-3,j))
!
      case default
        print*, "bc_extrap0_2_0: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_extrap0_2_0
!***********************************************************************
    subroutine bc_extrap0_2_1(f,topbot,j)
!
!  Extrapolation boundary condition for f(bdry)=0.
!  Correct for polynomials up to 2nd order, determined 1 further degree
!  of freedom by minimizing L2 norm of coefficient vector.
!
!  NOTE: This is not the final formula, but just bc_extrap_2_1() with f(bdry)=0
!
!  09-oct-03/wolf: coded
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        f(:,:,n1  ,j)=0.        ! set bdry value=0 (indep of initcond)
        f(:,:,n1-1,j)=0.25*(- 3*f(:,:,n1+1,j)- 5*f(:,:,n1+2,j)+ 3*f(:,:,n1+3,j))
        f(:,:,n1-2,j)=0.05*(-43*f(:,:,n1+1,j)-57*f(:,:,n1+2,j)+39*f(:,:,n1+3,j))
        f(:,:,n1-3,j)=0.05*(-81*f(:,:,n1+1,j)-99*f(:,:,n1+2,j)+73*f(:,:,n1+3,j))
!
      case ('top')               ! top boundary
        f(:,:,n2  ,j)=0.        ! set bdry value=0 (indep of initcond)
        f(:,:,n2+1,j)=0.25*(- 3*f(:,:,n2-1,j)- 5*f(:,:,n2-2,j)+ 3*f(:,:,n2-3,j))
        f(:,:,n2+2,j)=0.05*(-43*f(:,:,n2-1,j)-57*f(:,:,n2-2,j)+39*f(:,:,n2-3,j))
        f(:,:,n2+3,j)=0.05*(-81*f(:,:,n2-1,j)-99*f(:,:,n2-2,j)+73*f(:,:,n2-3,j))
!
      case default
        print*, "bc_extrap0_2_1: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_extrap0_2_1
!***********************************************************************
    subroutine bc_extrap0_2_2(f,topbot,j)
!
!  Extrapolation boundary condition for f(bdry)=0.
!  Correct for polynomials up to 2nd order, determined 1 further degree
!  of freedom by minimizing L2 norm of coefficient vector.
!
!  NOTE: This is not the final formula, but just bc_extrap_2_2() with f(bdry)=0
!
!  09-oct-03/wolf: coded
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j,n1p4,n2m4
!
!  abbreviations, because otherwise the ifc compiler complains
!  for 1-D runs without vertical extent
!
      n1p4=n1+4
      n2m4=n2-4
!
      select case (topbot)
!
      case ('bot')
! bottom boundary
        f(:,:,n1  ,j)= 0.       ! set bdry value=0 (indep of initcond)
        f(:,:,n1-1,j)=0.2   *(                 -  4*f(:,:,n1+2,j)- 3*f(:,:,n1+3,j)+ 3*f(:,:,n1p4,j))
        f(:,:,n1-2,j)=0.2   *(- 2*f(:,:,n1+1,j)-  9*f(:,:,n1+2,j)- 6*f(:,:,n1+3,j)+ 7*f(:,:,n1p4,j))
        f(:,:,n1-3,j)=1./35.*(-33*f(:,:,n1+1,j)-108*f(:,:,n1+2,j)-68*f(:,:,n1+3,j)+87*f(:,:,n1p4,j))
!
      case ('top')
! top boundary
        f(:,:,n2  ,j)= 0.       ! set bdry value=0 (indep of initcond)
        f(:,:,n2+1,j)=0.2   *(                 -  4*f(:,:,n2-2,j)- 3*f(:,:,n2-3,j)+ 3*f(:,:,n2m4,j))
        f(:,:,n2+2,j)=0.2   *(- 2*f(:,:,n2-1,j)-  9*f(:,:,n2-2,j)- 6*f(:,:,n2-3,j)+ 7*f(:,:,n2m4,j))
        f(:,:,n2+3,j)=1./35.*(-33*f(:,:,n2-1,j)-108*f(:,:,n2-2,j)-68*f(:,:,n2-3,j)+87*f(:,:,n2m4,j))
!
      case default
        print*, "bc_extrap0_2_2: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_extrap0_2_2
!***********************************************************************
    subroutine bcx_extrap_2_3(f,topbot,j)
!
!  Extrapolation boundary condition in logarithm:
!  It maintains a power law.
!
!  y_{b+i} = y_b + a * (x_{b+1} - x_b)
!
!  where a = (y_b - y_{b-1})/(x_b-x_{b-1})
!
!  18-dec-08/wlad: coded
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j,i
      real :: yl1,ypi,ymi,xl1,xmi,xpi,yyi,xl2,yl2
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        do i=1,nghost
          do n=1,size(f,3);do m=1,size(f,2)
            yl1=alog(f(l1,m,n,j)) ; ypi=alog(f(l1+i,m,n,j))
            xl1=alog(x(l1)) ; xmi=alog(x(l1-i)) ; xpi=alog(x(l1+i))
!
            yyi = yl1 - (ypi-yl1)*(xl1-xmi)/(xpi-xl1)
            f(l1-i,m,n,j) = exp(yyi)
          enddo;enddo
        enddo
!
      case ('top')               ! top boundary
        do i=1,nghost
          do n=1,size(f,3);do m=1,size(f,2)
            yl2=alog(f(l2,m,n,j)) ; ymi=alog(f(l2-i,m,n,j))
            xpi=alog(x(l2+i)) ; xl2=alog(x(l2)) ; xmi=alog(x(l2-i))
!
            yyi = yl2 + (yl2-ymi)*(xpi-xl2)/(xl2-xmi)
            f(l2+i,m,n,j) = exp(yyi)
          enddo;enddo
        enddo
!
      case default
        print*, "bcx_extrap_2_3: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bcx_extrap_2_3
!***********************************************************************
    subroutine bcx_extrap_linear(f, topbot, j)
!
!  Applies linear extrapolation to the ghost cells.
!
!  05-jun-18/ccyang: coded.
!
      real, dimension(:,:,:,:), intent(inout) :: f
      character(len=bclen), intent(in) :: topbot
      integer, intent(in) :: j
!
      integer :: i
      real :: dx1
!
      select case (topbot)
!
      case ('bot')
        ! bottom (left end of the domain)
        dx1 = 1.0 / (x(l1+1) - x(l1))
        do i = 1, nghost
          f(l1-i,:,:,j) = (dx1 * (x(l1+1) - x(l1-i))) * f(l1,:,:,j) + (dx1 * (x(l1-i) - x(l1))) * f(l1+1,:,:,j)
        enddo
!
      case ('top')
        ! top (right end of the domain)
        dx1 = 1.0 / (x(l2) - x(l2-1))
        do i = 1, nghost
          f(l2+i,:,:,j) = (dx1 * (x(l2) - x(l2+i))) * f(l2-1,:,:,j) + (dx1 * (x(l2+i) - x(l2-1))) * f(l2,:,:,j)
        enddo
!
      case default
        call fatal_error('bcx_extrap_linear', 'invalid argument', lfirst_proc_xy)
!
      endselect
!
    endsubroutine bcx_extrap_linear
!***********************************************************************
    subroutine bcz_extrapol(f,topbot,j)
!
!  Simple linear extrapolation in first order.
!  The last two grid points are used to determine the slope.
!
!  23-nov-10/Bourdin.KIS: coded
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: i, j
!
      real, dimension (size(f,1),size(f,2)) :: slope
!
!
      select case (topbot)
      case ('bot')
        ! bottom (left end of the domain)
        slope = (f(:,:,n1+1,j) - f(:,:,n1,j)) / dz2_bound(-1) 
        do i = 1, nghost
          f(:,:,n1-i,j) = f(:,:,n1,j) - slope * dz2_bound(-i) 
        enddo
      case ('top')
        ! top (right end of the domain)
        slope = (f(:,:,n2,j) - f(:,:,n2-1,j)) / dz2_bound(1)
        do i = 1, nghost
          f(:,:,n2+i,j) = f(:,:,n2,j) + slope * dz2_bound(i) 
        enddo
      case default
        call fatal_error ('bcz_extrapol', 'invalid argument', lfirst_proc_xy)
      endselect
!
    endsubroutine bcz_extrapol
!***********************************************************************
    subroutine bcz_extrapol_fixed(f,topbot,j)
!
!  Simple linear extrapolation in first order
!  with a fixed value in the first ghost cell.
!  The last two grid points are used to determine the slope.
!
!  23-nov-10/Bourdin.KIS: coded
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: i, j
!
      real, dimension (size(f,1),size(f,2)) :: m
!
!
      select case (topbot)
      case ('bot')
        ! bottom (left end of the domain)
        m = (f(:,:,n1+1,j) - f(:,:,n1,j)) / (z(n1+1) - z(n1))
        do i = 2, nghost
          f(:,:,n1-i,j) = f(:,:,n1-1,j) + m * (z(n1-i) - z(n1-1))
        enddo
      case ('top')
        ! top (right end of the domain)
        m = (f(:,:,n2,j) - f(:,:,n2-1,j)) / (z(n2) - z(n2-1))
        do i = 2, nghost
          f(:,:,n2+i,j) = f(:,:,n2+1,j) + m * (z(n2+i) - z(n2+1))
        enddo
      case default
        call fatal_error ('bcz_extrapol_fixed', 'invalid argument', lfirst_proc_xy)
      endselect
!
    endsubroutine bcz_extrapol_fixed
!***********************************************************************
    subroutine bcz_extrapol_damped(f,topbot,j)
!
!  Simple linear extrapolation in first order
!  with an included damping to zero (useful for velocities).
!  The last two grid points are used to determine the slope.
!  The parameters 'fbcz_bot' and 'fbcz_top' are used as damping factors,
!  which should have values between 0.0 (no damping) and 1.0 (full damping).
!  A typical value of 0.001 corresponds to a half-value time of ~1000 timesteps.
!  Negative values let the damping be persistent, even if 'lfade_damp' is true.
!
!  23-nov-10/Bourdin.KIS: coded
!
      use SharedVariables, only: get_shared_variable
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: i, j
!
      real, dimension (size(f,1),size(f,2)) :: slope
      real :: gamma_bot, gamma_top, tau, fade_fact
      real, pointer :: tdamp, tfade_start
      logical, pointer :: ldamp_fade
!
!
      ! bottom and top damping factors:
      gamma_bot = 1.0 - abs (fbcz_bot(j))
      gamma_top = 1.0 - abs (fbcz_top(j))
!
      call get_shared_variable ('ldamp_fade', ldamp_fade, caller='bcz_extrapol_damped')
      if (ldamp_fade) then
        ! fading of damping is active
        call get_shared_variable ('tdamp', tdamp)
        call get_shared_variable ('tfade_start', tfade_start)
        if (t > tfade_start) then
          if (t < tdamp) then
            ! tau is a normalized t, the transition interval is [-0.5, 0.5]:
            tau = (t-tfade_start) / (tdamp-tfade_start) - 0.5
            fade_fact = 0.5 * (1 - tau * (3 - 4*tau**2))
            ! apply damping with fading:
            gamma_bot = 1.0 - abs (fbcz_bot(j)) * fade_fact
            gamma_top = 1.0 - abs (fbcz_top(j)) * fade_fact
          else
            ! damping has already ended (t >= tdamp)
            gamma_bot = 1.0
            gamma_top = 1.0
          endif
          ! apply fading-persistent damping:
          if (fbcz_bot(j) < 0.0) gamma_bot = 1.0 - abs (fbcz_bot(j))
          if (fbcz_top(j) < 0.0) gamma_top = 1.0 - abs (fbcz_top(j))
        endif
      endif
!
      select case (topbot)
      case ('bot')
        ! bottom (left end of the domain)
        slope = (f(:,:,n1+1,j) - f(:,:,n1,j)) / dz2_bound(-1)
        do i = 1, nghost
          f(:,:,n1-i,j) = (f(:,:,n1,j) - slope * dz2_bound(-i)) * gamma_bot**i
        enddo
        f(:,:,n1,j) = 0.5*(f(:,:,n1-1,j) + f(:,:,n1+1,j))
      case ('top')
        ! top (right end of the domain)
        slope = (f(:,:,n2,j) - f(:,:,n2-1,j)) / dz2_bound(1)
        do i = 1, nghost
          f(:,:,n2+i,j) = (f(:,:,n2,j) + slope * dz2_bound(i)) * gamma_top**i
        enddo
        f(:,:,n2,j) = 0.5*(f(:,:,n2-1,j) + f(:,:,n2+1,j))
      case default
        call fatal_error ('bcz_extrapol_damped', 'invalid argument', lfirst_proc_xy)
      endselect
!
    endsubroutine bcz_extrapol_damped
!***********************************************************************
    subroutine bcz_extrapol_mean (f, topbot, j)
!
!  Simple linear extrapolation in first order
!  with an included local averaging of a 7x7 array.
!  The last two grid points in z are used to determine the slope.
!  The 3 neighbouring grid points in x and y contribute to the local average.
!
!  11-apr-11/Bourdin.KIS: coded
!
      character (len=bclen), intent(in) :: topbot
      real, dimension (:,:,:,:), intent(inout) :: f
      integer, intent(in) :: j
!
      integer :: i
      real, dimension (size(f,1),size(f,2)) :: slope, rho_ref
!
      select case (topbot)
      case ('bot')
        ! bottom (left end of the domain)
        rho_ref = f(:,:,n1,j)
        call average_xy (rho_ref, 3)
        slope = (f(:,:,n1+1,j) - rho_ref) / dz2_bound(-1)
        do i = 1, nghost
          f(:,:,n1-i,j) = rho_ref - slope * dz2_bound(-i)
        enddo
        f(:,:,n1,j) = 0.5*(f(:,:,n1-1,j) + f(:,:,n1+1,j))
      case ('top')
        ! top (right end of the domain)
        rho_ref = f(:,:,n2,j)
        call average_xy (rho_ref, 3)
        slope = (rho_ref - f(:,:,n2-1,j)) / dz2_bound(1)
        do i = 1, nghost
          f(:,:,n2+i,j) = rho_ref + slope * dz2_bound(i)
        enddo
        f(:,:,n2,j) = 0.5*(f(:,:,n2-1,j) + f(:,:,n2+1,j))
      case default
        call fatal_error ('bcz_extrapol_mean', 'invalid argument', lfirst_proc_xy)
      endselect
!
    endsubroutine bcz_extrapol_mean
!***********************************************************************
    subroutine average_xy (data, num)
!
!  Simple averaging over a num*num array in x and y direction.
!
!  11-apr-11/Bourdin.KIS: coded
!
      use Mpicomm, only: communicate_xy_ghosts
!
      real, dimension (:,:), intent(inout) :: data
      integer, intent(in) :: num
!
      real, dimension (size(data,1),size(data,2)) :: out
      integer :: px, py
!
!
      out = 0.0
      do px = l1, l2
        do py = m1, m2
          out(px,py) = sum (data(px-num:px+num,py-num:py+num))
        enddo
      enddo
      data = out / (num*2+1)**2
!
      call communicate_xy_ghosts (data)
!
    endsubroutine average_xy
!***********************************************************************
    subroutine bc_db_z(f,topbot,j)
!
!  "One-sided" boundary condition for density.
!  Set ghost zone to reproduce one-sided boundary condition
!  (2nd order):
!  Finding the derivatives on the boundary using a one
!  sided final difference method. This derivative is being
!  used to calculate the boundary points. This will probably
!  only be used for ln(rho)
!
!  may-2002/nils: coded
!  11-jul-2002/nils: moved into the density module
!  13-aug-2002/nils: moved into boundcond
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j
!
      real, dimension (size(f,1),size(f,2)) :: fder
      integer :: i
!
      select case (topbot)
!
! Bottom boundary
!
      case ('bot')
        do i=1,nghost
          fder=(-3*f(:,:,n1-i+1,j)+4*f(:,:,n1-i+2,j)&
               -f(:,:,n1-i+3,j))/(2*dz)
          f(:,:,n1-i,j)=f(:,:,n1-i+2,j)-2*dz*fder
        enddo
      case ('top')
        do i=1,nghost
          fder=(3*f(:,:,n2+i-1,j)-4*f(:,:,n2+i-2,j)&
               +f(:,:,n2+i-3,j))/(2*dz)
          f(:,:,n2+i,j)=f(:,:,n2+i-2,j)+2*dz*fder
        enddo
      case default
        print*,"bc_db_z: invalid argument for 'bc_db_z'"
      endselect
!
    endsubroutine bc_db_z
!***********************************************************************
    subroutine bc_db_x(f,topbot,j)
!
!  "One-sided" boundary condition for density.
!  Set ghost zone to reproduce one-sided boundary condition
!  (2nd order):
!  Finding the derivatives on the boundary using a one
!  sided final difference method. This derivative is being
!  used to calculate the boundary points. This will probably
!  only be used for ln(rho)
!
!  may-2002/nils: coded
!  11-jul-2002/nils: moved into the density module
!  13-aug-2002/nils: moved into boundcond
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j
!
      real, dimension (size(f,2),size(f,3)) :: fder
      integer :: i
!
      select case (topbot)
!
! Bottom boundary
!
      case ('bot')
        do i=1,nghost
          fder=(-3*f(l1-i+1,:,:,j)+4*f(l1-i+2,:,:,j)&
               -f(l1-i+3,:,:,j))/(2*dx)
          f(l1-i,:,:,j)=f(l1-i+2,:,:,j)-2*dx*fder
        enddo
      case ('top')
        do i=1,nghost
          fder=(3*f(l2+i-1,:,:,j)-4*f(l2+i-2,:,:,j)&
               +f(l2+i-3,:,:,j))/(2*dx)
          f(l2+i,:,:,j)=f(l2+i-2,:,:,j)+2*dx*fder
        enddo
      case default
        print*,"bc_db_x: invalid argument for 'bc_db_x'"
      endselect
!
    endsubroutine bc_db_x
!***********************************************************************
    subroutine bc_force_z(f,sgn,topbot,j)
!
!  Force values of j-th variable on vertical boundary topbot.
!  This can either be used for freezing variables at the boundary, or for
!  enforcing a certain time-dependent function of (x,y).
!
!  Currently this is hard-coded for velocity components (ux,uy) and quite
!  useless. Plan is to read time-dependent velocity field from disc and
!  apply it as boundary condition here.
!
!  26-apr-2004/wolf: coded
!
      use EquationOfState, only: gamma_m1, cs2top, cs2bot
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: sgn,i,j
!
      select case (topbot)
!
!  lower boundary
!
      case ('bot')
         select case (force_lower_bound)
         case ('uxy_sin-cos')
            call bc_force_uxy_sin_cos(f,n1,j)
         case ('axy_sin-cos')
            call bc_force_axy_sin_cos(f,n1,j)
         case ('uxy_convection')
            call uu_driver(f)
         !case ('kepler')
         !   call bc_force_kepler(f,n1,j)
         case ('mag_time')
            call bc_force_aa_time(f)
         case ('mag_convection')
            call bc_force_aa_time(f)
            call uu_driver(f)
         case ('cT')
            f(:,:,n1,j) = log(cs2bot/gamma_m1)
         case ('vel_time')
            call bc_force_ux_time(f,n1,j)
         case default
            if (lroot) print*, "No such value for force_lower_bound: <", &
                 trim(force_lower_bound),">"
            call stop_it("")
         endselect
         !
         !  Now fill ghost zones imposing antisymmetry w.r.t. the values just set:
         !
         do i=1,nghost; f(:,:,n1-i,j)=2*f(:,:,n1,j)+sgn*f(:,:,n1+i,j); enddo
!
!  upper boundary
!
      case ('top')
         select case (force_upper_bound)
         case ('uxy_sin-cos')
            call bc_force_uxy_sin_cos(f,n2,j)
         case ('axy_sin-cos')
            call bc_force_axy_sin_cos(f,n2,j)
         case ('uxy_convection')
            call uu_driver(f)
         !case ('kepler')
         !   call bc_force_kepler(f,n2,j)
         case ('cT')
            f(:,:,n2,j) = log(cs2top/gamma_m1)
         case ('vel_time')
            call bc_force_ux_time(f,n2,j)
         case default
            if (lroot) print*, "No such value for force_upper_bound: <", &
                 trim(force_upper_bound),">"
            call stop_it("")
         endselect
         !
         !  Now fill ghost zones imposing antisymmetry w.r.t. the values just set:
         !
         do i=1,nghost; f(:,:,n2+i,j)=2*f(:,:,n2,j)+sgn*f(:,:,n2-i,j); enddo
      case default
        print*,"bc_force_z: invalid argument topbot=",topbot
      endselect
!
    endsubroutine bc_force_z
!***********************************************************************
    subroutine bc_force_x(f, sgn, topbot, j)
!
!  Force values of j-th variable on x-boundaries topbot.
!
!  09-mar-2007/dintrans: coded
!
      use SharedVariables, only : get_shared_variable
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      real, pointer :: ampl_forc, k_forc, w_forc
      integer :: sgn, i, j
!
      select case (topbot)
!
!  lower boundary
!
      case ('bot')
         select case (force_lower_bound)
         case ('vel_time')
           if (j /= iuy) call stop_it("BC_FORCE_X: only valid for uy")
           call get_shared_variable('ampl_forc', ampl_forc, caller='bc_force_x')
           call get_shared_variable('k_forc', k_forc)
           call get_shared_variable('w_forc', w_forc)
           if (headtt) print*, 'BC_FORCE_X: ampl_forc, k_forc, w_forc=',&
               ampl_forc, k_forc, w_forc
           f(l1,:,:,iuy) = spread(ampl_forc*sin(k_forc*y)*cos(w_forc*t), 2, size(f,3))
         case default
            if (lroot) print*, "No such value for force_lower_bound: <", &
                 trim(force_lower_bound),">"
            call stop_it("")
         endselect
         !
         !  Now fill ghost zones imposing antisymmetry w.r.t. the values just set:
         !
         do i=1,nghost; f(l1-i,:,:,j)=2*f(l1,:,:,j)+sgn*f(l1+i,:,:,j); enddo
!
!  upper boundary
!
      case ('top')
         select case (force_upper_bound)
         case ('vel_time')
            if (j /= iuy) call stop_it("BC_FORCE_X: only valid for uy")
            call get_shared_variable('ampl_forc', ampl_forc, caller='bc_force_x')
            call get_shared_variable('k_forc', k_forc)
            call get_shared_variable('w_forc', w_forc)
            if (headtt) print*, 'BC_FORCE_X: ampl_forc, k_forc, w_forc=',&
                   ampl_forc, k_forc, w_forc
            f(l2,:,:,iuy) = spread(ampl_forc*sin(k_forc*y)*cos(w_forc*t), 2, size(f,3))
         case default
            if (lroot) print*, "No such value for force_upper_bound: <", &
                 trim(force_upper_bound),">"
            call stop_it("")
         endselect
         !
         !  Now fill ghost zones imposing antisymmetry w.r.t. the values just set:
         !
         do i=1,nghost; f(l2+i,:,:,j)=2*f(l2,:,:,j)+sgn*f(l2-i,:,:,j); enddo
      case default
        print*,"bc_force_x: invalid argument topbot=",topbot
      endselect
!
    endsubroutine bc_force_x
!***********************************************************************
    subroutine bc_force_uxy_sin_cos(f,idz,j)
!
!  Set (ux, uy) = (cos y, sin x) in vertical layer
!
!  26-apr-2004/wolf: coded
!
      real, dimension (:,:,:,:) :: f
      integer :: idz,j
      real :: kx,ky
!
      if (iuz == 0) call stop_it("BC_FORCE_UXY_SIN_COS: Bad idea...")
!
      if (j==iux) then
        if (Ly>0) then; ky=2*pi/Ly; else; ky=0.; endif
        f(:,:,idz,j) = spread(cos(ky*y),1,size(f,1))
      elseif (j==iuy) then
        if (Lx>0) then; kx=2*pi/Lx; else; kx=0.; endif
        f(:,:,idz,j) = spread(sin(kx*x),2,size(f,2))
      elseif (j==iuz) then
        f(:,:,idz,j) = 0.
      endif
!
    endsubroutine bc_force_uxy_sin_cos
!***********************************************************************
    subroutine bc_force_axy_sin_cos(f,idz,j)
!
!  Set (ax, ay) = (cos y, sin x) in vertical layer
!
!  26-apr-2004/wolf: coded
!  10-apr-2005/axel: adapted for A
!
      real, dimension (:,:,:,:) :: f
      integer :: idz,j
      real :: kx,ky
!
      if (iaz == 0) call stop_it("BC_FORCE_AXY_SIN_COS: Bad idea...")
!
      if (j==iax) then
        if (Ly>0) then; ky=2*pi/Ly; else; ky=0.; endif
        f(:,:,idz,j) = spread(cos(ky*y),1,size(f,1))
      elseif (j==iay) then
        if (Lx>0) then; kx=2*pi/Lx; else; kx=0.; endif
        f(:,:,idz,j) = spread(sin(kx*x),2,size(f,2))
      elseif (j==iaz) then
        f(:,:,idz,j) = 0.
      endif
!
    endsubroutine bc_force_axy_sin_cos
!***********************************************************************
    subroutine bc_one_x(f,topbot,j)
!
!  Set bdry values to 1 for debugging purposes
!
!  11-jul-02/wolf: coded
!
      real, dimension (:,:,:,:) :: f
      integer :: j
      character (len=bclen) :: topbot
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
          f(1:l1-1,:,:,j)=1.
!
      case ('top')               ! top boundary
          f(l2+1:,:,:,j)=1.
!
      case default
        print*, "bc_one_x: ",topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_one_x
!***********************************************************************
    subroutine bc_one_y(f,topbot,j)
!
!  Set bdry values to 1 for debugging purposes
!
!  11-jul-02/wolf: coded
!
      real, dimension (:,:,:,:) :: f
      integer :: j
      character (len=bclen) :: topbot
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
          f(:,1:m1-1,:,j)=1.
!
      case ('top')               ! top boundary
          f(:,m2+1:,:,j)=1.
!
      case default
        print*, "bc_one_y: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_one_y
!***********************************************************************
    subroutine bc_one_z(f,topbot,j)
!
!  Set bdry values to 1 for debugging purposes
!
!  11-jul-02/wolf: coded
!
      real, dimension (:,:,:,:) :: f
      integer :: j
      character (len=bclen) :: topbot
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
          f(:,:,1:n1-1,j)=1.
!
      case ('top')               ! top boundary
          f(:,:,n2+1:,j)=1.
!
      case default
        print*, "bc_one_z: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_one_z
!***********************************************************************
    subroutine bc_freeze_var_x(topbot,j)
!
!  Tell other modules that variable with slot j is to be frozen in on
!  given boundary
!
      integer :: j
      character (len=bclen) :: topbot
!
      lfrozen_bcs_x = .true.    ! set flag
!
      select case (topbot)
      case ('bot')               ! bottom boundary
        lfrozen_bot_var_x(j) = .true.
      case ('top')               ! top boundary
        lfrozen_top_var_x(j) = .true.
      case default
        print*, "bc_freeze_var_x: ", topbot, " should be 'top' or 'bot'"
      endselect
!
    endsubroutine bc_freeze_var_x
!***********************************************************************
    subroutine bc_freeze_var_y(topbot,j)
!
!  Tell other modules that variable with slot j is to be frozen in on
!  given boundary
!
      integer :: j
      character (len=bclen) :: topbot
!
      lfrozen_bcs_y = .true.    ! set flag
!
      select case (topbot)
      case ('bot')               ! bottom boundary
        lfrozen_bot_var_y(j) = .true.
      case ('top')               ! top boundary
        lfrozen_top_var_y(j) = .true.
      case default
        print*, "bc_freeze_var_y: ", topbot, " should be 'top' or 'bot'"
      endselect
!
    endsubroutine bc_freeze_var_y
!***********************************************************************
    subroutine bc_freeze_var_z(topbot,j)
!
!  Tell other modules that variable with slot j is to be frozen in on
!  given boundary
!
      integer :: j
      character (len=bclen) :: topbot
!
      lfrozen_bcs_z = .true.    ! set flag
!
      select case (topbot)
      case ('bot')               ! bottom boundary
        lfrozen_bot_var_z(j) = .true.
      case ('top')               ! top boundary
        lfrozen_top_var_z(j) = .true.
      case default
        print*, "bc_freeze_var_z: ", topbot, " should be 'top' or 'bot'"
      endselect
!
    endsubroutine bc_freeze_var_z
!***********************************************************************
     subroutine uu_driver(f,quenching)
!
!  Simulated velocity field used as photospherec motions
!  Use of velocity field produced by Boris Gudiksen
!
!  27-mai-04/bing: coded
!  11-aug-06/axel: make it compile with nprocx>0, renamed quenching -> quen
!  18-jun-08/bing: quenching depends on B^2, not only Bz^2
!
       use EquationOfState, only : gamma,gamma_m1,gamma1,cs20,lnrho0
       use File_io, only : file_exists
       use Mpicomm, only : mpisend_real, mpirecv_real
!
       real, dimension (:,:,:,:), intent (inout) :: f
       logical, optional :: quenching
!
       real, dimension (nx,ny), save :: uxl,uxr,uyl,uyr
       real, dimension (:,:), allocatable :: tmp
       real, dimension (nx,ny) :: uxd,uyd,quen,pp,betaq,fac,bbx,bby,bbz,bb2

       integer :: tag_xl=321,tag_yl=322,tag_xr=323,tag_yr=324
       integer :: tag_tl=345,tag_tr=346,tag_dt=347
       integer :: lend=0,ierr,frame=0,pos,iref,px,py
       real, save :: tl=0.,tr=0.,delta_t=0.
       real  :: zmin
       logical :: quench
!
       character (len=*), parameter :: vel_times_dat = 'driver/vel_times.dat'
       character (len=*), parameter :: vel_field_dat = 'driver/vel_field.dat'
       integer :: unit=1
!
       if (ldownsampling) then
         call warning('uu_driver','Not available for downsampling')
         return
       endif
!
       if (lroot .and. .not. file_exists(vel_times_dat)) &
           call fatal_error_local('uu_driver','Could not find file "'//trim(vel_times_dat)//'"')
       if (lroot .and. .not. file_exists(vel_field_dat)) &
           call fatal_error_local('uu_driver', 'Could not find file "'//trim(vel_field_dat)//'"')
!
       if (present(quenching)) then
         quench = quenching
       else
         ! Right now quenching is per default active
         quench=.true.
       endif
!
!  Read the time table
!
       if ((t*unit_time<tl+delta_t) .or. (t*unit_time>=tr+delta_t)) then
!
         if (lroot) then
           inquire(IOLENGTH=lend) tl
           open (unit,file=vel_times_dat,form='unformatted',status='unknown',recl=lend,access='direct')
!
           ierr = 0
           frame = 0
           do while (ierr == 0)
             frame=frame+1
             read (unit,rec=frame,iostat=ierr) tl
             read (unit,rec=frame+1,iostat=ierr) tr
             if (ierr /= 0) then
               frame=1
               delta_t = t*unit_time                  ! EOF is reached => read again
               read (unit,rec=frame,iostat=ierr) tl
               read (unit,rec=frame+1,iostat=ierr) tr
               ierr=-1
             else
               if (t*unit_time>=tl+delta_t .and. t*unit_time<tr+delta_t) ierr=-1
               ! correct time step is reached
             endif
           enddo
           close (unit)
!
           do px=0, nprocx-1
             do py=0, nprocy-1
               if ((px /= 0) .or. (py /= 0)) then
                 call mpisend_real (tl, px+py*nprocx, tag_tl)
                 call mpisend_real (tr, px+py*nprocx, tag_tr)
                 call mpisend_real (delta_t, px+py*nprocx, tag_dt)
               endif
             enddo
           enddo
         else
           call mpirecv_real (tl, 0, tag_tl)
           call mpirecv_real (tr, 0, tag_tr)
           call mpirecv_real (delta_t, 0, tag_dt)
         endif
!
! Read velocity field
!
         if (lroot) then

           allocate(tmp(nxgrid,nygrid),stat=ierr)
!
           if (ierr>0) call fatal_error('uu_driver', &
              'Could not allocate memory for array please check', .true.)

           open (unit,file=vel_field_dat,form='unformatted',status='unknown',recl=lend*nxgrid*nygrid,access='direct')
!
           read (unit,rec=2*frame-1) tmp
           do px=0, nprocx-1
             do py=0, nprocy-1
               if ((px /= 0) .or. (py /= 0)) then
                 uxl = tmp(px*nx+1:(px+1)*nx,py*ny+1:(py+1)*ny)
                 call mpisend_real (uxl, (/ nx, ny /), px+py*nprocx, tag_xl)
               endif
             enddo
           enddo
           uxl = tmp(1:nx,1:ny)
!
           read (unit,rec=2*frame) tmp
           do px=0, nprocx-1
             do py=0, nprocy-1
               if ((px /= 0) .or. (py /= 0)) then
                 uyl = tmp(px*nx+1:(px+1)*nx,py*ny+1:(py+1)*ny)
                 call mpisend_real (uyl, (/ nx, ny /), px+py*nprocx, tag_yl)
               endif
             enddo
           enddo
           uyl = tmp(1:nx,1:ny)
!
           read (unit,rec=2*frame+1) tmp
           do px=0, nprocx-1
             do py=0, nprocy-1
               if ((px /= 0) .or. (py /= 0)) then
                 uxr = tmp(px*nx+1:(px+1)*nx,py*ny+1:(py+1)*ny)
                 call mpisend_real (tmp(px*nx+1:(px+1)*nx,py*ny+1:(py+1)*ny), (/ nx, ny /), px+py*nprocx, tag_xr)
               endif
             enddo
           enddo
           uxr = tmp(1:nx,1:ny)
!
           read (unit,rec=2*frame+2) tmp
           uyr = tmp(1:nx,1:ny)
           do px=0, nprocx-1
             do py=0, nprocy-1
               if ((px /= 0) .or. (py /= 0)) then
                 uyr = tmp(px*nx+1:(px+1)*nx,py*ny+1:(py+1)*ny)
                 call mpisend_real (tmp(px*nx+1:(px+1)*nx,py*ny+1:(py+1)*ny), (/ nx, ny /), px+py*nprocx, tag_yr)
               endif
             enddo
           enddo
           uyr = tmp(1:nx,1:ny)
!
           close (unit)
           deallocate(tmp)
         else
           call mpirecv_real (uxl, (/ nx, ny /), 0, tag_xl)
           call mpirecv_real (uyl, (/ nx, ny /), 0, tag_yl)
           call mpirecv_real (uxr, (/ nx, ny /), 0, tag_xr)
           call mpirecv_real (uyr, (/ nx, ny /), 0, tag_yr)
         endif
!
         uxl = uxl / 10. / unit_velocity
         uxr = uxr / 10. / unit_velocity
         uyl = uyl / 10. / unit_velocity
         uyr = uyr / 10. / unit_velocity
!
       endif
!
!   simple linear interploation between timesteps
!
       if (tr /= tl) then
         uxd  = (t*unit_time - (tl+delta_t)) * (uxr - uxl) / (tr - tl) + uxl
         uyd  = (t*unit_time - (tl+delta_t)) * (uyr - uyl) / (tr - tl) + uyl
       else
         uxd = uxl
         uyd = uyl
       endif
!
!   suppress footpoint motion at low plasma beta
!
       zmin = minval(abs(z(n1:n2)))
       iref = n1
       do pos=n1,n2
         if (abs(z(pos))==zmin) iref=pos; exit
       enddo
!
!   Calculate B^2 for plasma beta
!
       if (quench) then
!-----------------------------------------------------------------------
         if (nygrid/=1) then
           fac=(1./60)*spread(dy_1(m1:m2),1,nx)
           bbx= fac*(+ 45.0*(f(l1:l2,m1+1:m2+1,iref,iaz)-f(l1:l2,m1-1:m2-1,iref,iaz)) &
               -  9.0*(f(l1:l2,m1+2:m2+2,iref,iaz)-f(l1:l2,m1-2:m2-2,iref,iaz)) &
               +      (f(l1:l2,m1+3:m2+3,iref,iaz)-f(l1:l2,m1-3:m2-3,iref,iaz)))
         else
           if (ip<=5) print*, 'uu_driver: Degenerate case in y-direction'
         endif
         if (nzgrid/=1) then
           fac=(1./60)*spread(spread(dz_1(iref),1,nx),2,ny)
           bbx= bbx -fac*(+ 45.0*(f(l1:l2,m1:m2,iref+1,iay)-f(l1:l2,m1:m2,iref-1,iay)) &
               -  9.0*(f(l1:l2,m1:m2,iref+2,iay)-f(l1:l2,m1:m2,iref-2,iay)) &
               +      (f(l1:l2,m1:m2,iref+3,iay)-f(l1:l2,m1:m2,iref-2,iay)))
         else
           if (ip<=5) print*, 'uu_driver: Degenerate case in z-direction'
         endif
!-----------------------------------------------------------------------
         if (nzgrid/=1) then
           fac=(1./60)*spread(spread(dz_1(iref),1,nx),2,ny)
           bby= fac*(+ 45.0*(f(l1:l2,m1:m2,iref+1,iax)-f(l1:l2,m1:m2,iref-1,iax)) &
               -  9.0*(f(l1:l2,m1:m2,iref+2,iax)-f(l1:l2,m1:m2,iref-2,iax)) &
               +      (f(l1:l2,m1:m2,iref+3,iax)-f(l1:l2,m1:m2,iref-3,iax)))
         else
           if (ip<=5) print*, 'uu_driver: Degenerate case in z-direction'
         endif
         if (nxgrid/=1) then
           fac=(1./60)*spread(dx_1(l1:l2),2,ny)
           bby=bby-fac*(+45.0*(f(l1+1:l2+1,m1:m2,iref,iaz)-f(l1-1:l2-1,m1:m2,iref,iaz)) &
               -  9.0*(f(l1+2:l2+2,m1:m2,iref,iaz)-f(l1-2:l2-2,m1:m2,iref,iaz)) &
               +      (f(l1+3:l2+3,m1:m2,iref,iaz)-f(l1-3:l2-3,m1:m2,iref,iaz)))
         else
           if (ip<=5) print*, 'uu_driver: Degenerate case in x-direction'
         endif
!-----------------------------------------------------------------------
         if (nxgrid/=1) then
           fac=(1./60)*spread(dx_1(l1:l2),2,ny)
           bbz= fac*(+ 45.0*(f(l1+1:l2+1,m1:m2,iref,iay)-f(l1-1:l2-1,m1:m2,iref,iay)) &
               -  9.0*(f(l1+2:l2+2,m1:m2,iref,iay)-f(l1-2:l2-2,m1:m2,iref,iay)) &
               +      (f(l1+3:l2+3,m1:m2,iref,iay)-f(l1-3:l2-3,m1:m2,iref,iay)))
         else
           if (ip<=5) print*, 'uu_driver: Degenerate case in x-direction'
         endif
         if (nygrid/=1) then
           fac=(1./60)*spread(dy_1(m1:m2),1,nx)
           bbz=bbz-fac*(+45.0*(f(l1:l2,m1+1:m2+1,iref,iax)-f(l1:l2,m1-1:m2-1,iref,iax)) &
               -  9.0*(f(l1:l2,m1+2:m2+2,iref,iax)-f(l1:l2,m1-2:m2-2,iref,iax)) &
               +      (f(l1:l2,m1+3:m2+3,iref,iax)-f(l1:l2,m1-3:m2-3,iref,iax)))
         else
           if (ip<=5) print*, 'uu_driver: Degenerate case in y-direction'
         endif
!-----------------------------------------------------------------------
!
         bb2 = bbx*bbx + bby*bby + bbz*bbz
         bb2 = bb2/(2.*mu0)
!
         if (ltemperature) then
           pp=gamma_m1*gamma1*exp(f(l1:l2,m1:m2,iref,ilnrho)+f(l1:l2,m1:m2,iref,ilnTT))
         else if (lentropy) then
           if (pretend_lnTT) then
             pp=gamma_m1*gamma1*exp(f(l1:l2,m1:m2,iref,ilnrho)+f(l1:l2,m1:m2,iref,iss))
           else
             pp=gamma*(f(l1:l2,m1:m2,iref,iss)+ &
                 f(l1:l2,m1:m2,iref,ilnrho))-gamma_m1*lnrho0
             pp=exp(pp) * cs20*gamma1
           endif
         else
           pp=gamma1*cs20*exp(lnrho0)
         endif
!
!   limit plasma beta
!
         betaq = pp / max(tini,bb2)*1e-3
!
         quen=(1.+betaq**2)/(10.+betaq**2)
       else
         quen(:,:)=1.
       endif
!
!   Fill z=0 layer with velocity field
!
       f(l1:l2,m1:m2,iref,iux)=uxd*quen
       f(l1:l2,m1:m2,iref,iuy)=uyd*quen
       if (iref/=n1) f(l1:l2,m1:m2,n1,iux:iuz)=0.
!
     endsubroutine uu_driver
!***********************************************************************
    subroutine bc_force_aa_time(f)
!
!  Reads in time series of magnetograms
!
!  17-feb-10/bing: coded
!  25-jul-10/Bourdin.KIS: parallelized
!
      use Fourier, only : setup_extrapol_fact, field_extrapol_z_parallel
      use Mpicomm, only : mpisend_real, mpirecv_real, &
                          mpisend_logical, mpirecv_logical
!
      real, dimension (:,:,:,:) :: f
      real, save :: t_l=0., t_r=0., delta_t=0.
      integer :: ierr, lend, frame, stat, rec_l, rec_r
      integer :: rec_vxl, rec_vxr, rec_vyl, rec_vyr ! l- and r-record position if file
      integer, parameter :: bnx=nxgrid, bny=ny/nprocx ! data in pencil shape
      integer, parameter :: enx=nygrid, eny=nx/nprocy ! transposed data in pencil shape
      integer :: px, py, partner
      integer, parameter :: tag_l=208, tag_r=209, tag_dt=210
      logical, save :: luse_vel_field = .false., first_run = .true.
      logical :: ex
!
      ! temporal storage for frames before (l) and after (r) current time step:
      real, dimension (:,:), allocatable, save :: Bz0_l, Bz0_r
      real, dimension (:,:), allocatable, save :: vx_l, vx_r, vy_l, vy_r
      real, dimension (:,:), allocatable :: vx_tmp, vy_tmp
      ! current magnetic field z-component (interpolated between l and r):
      real, dimension (:,:), allocatable :: Bz0
      ! current velocity x- and y-component (interpolated between l and r):
      real, dimension (:,:), allocatable, save :: vx, vy
!
      real, dimension (:,:,:), allocatable, save :: exp_fact ! exponential factor
      integer :: i
      real, parameter :: reduce_factor=0.25
!
      real :: time_SI
!
      character (len=*), parameter :: mag_field_dat = 'driver/mag_field.dat'
      character (len=*), parameter :: mag_times_dat = 'driver/mag_times.dat'
      character (len=*), parameter :: mag_vel_field_dat = 'driver/mag_vel_field.dat'
!
      if (ldownsampling) then
        call warning('bc_force_aa_time','Not available for downsampling')
        return
      endif
!
      if (first_run) then
!
!  Check for consistency:
        if ((.not. lequidist(1)) .or. (.not. lequidist(2))) &
            call fatal_error ('bc_force_aa_time', 'not yet implemented for non-equidistant grids', lfirst_proc_xy)
        if (mod (nygrid, nprocxy) /= 0) &
            call fatal_error ('bc_force_aa_time', 'nygrid needs to be an integer multiple of nprocx*nprocy', lfirst_proc_xy)
        if (mod (nxgrid, nprocxy) /= 0) &
            call fatal_error ('bc_force_aa_time', 'nxgrid needs to be an integer multiple of nprocx*nprocy', lfirst_proc_xy)
!
!  Check for existence of necessary driver files:
        if (lfirst_proc_xy) then
          inquire (file=mag_field_dat, exist=ex)
          if (.not. ex) call fatal_error ('bc_force_aa_time', 'File does not exists: '//trim(mag_field_dat), .true.)
          inquire (file=mag_times_dat, exist=ex)
          if (.not. ex) call fatal_error ('bc_force_aa_time', 'File does not exists: '//trim(mag_times_dat), .true.)
          inquire (file=mag_vel_field_dat, exist=ex)
          if (ex) then
            luse_vel_field = .true.
            print *, 'bc_force_aa_time: using time dependant magnetogram _with_ corresponding horizontal velocities.'
          else
            print *, 'bc_force_aa_time: using time dependant magnetogram _without_ corresponding horizontal velocities.'
          endif
        endif
!
        do i=1,nprocxy-1
          if (lroot) then
            call mpisend_logical(luse_vel_field,i,i)
          elseif (iproc==i) then
            call mpirecv_logical(luse_vel_field,0,iproc)
          endif
        enddo
!
!  Initialization of magnetograms and velocity fields.
        allocate(Bz0_l(bnx,bny),stat=stat)
        if (stat>0) call fatal_error('bc_force_aa_time','Could not allocate memory for Bz0_l',.true.)
        allocate(Bz0_r(bnx,bny),stat=stat)
        if (stat>0) call fatal_error('bc_force_aa_time','Could not allocate memory for Bz0_r',.true.)
        if (luse_vel_field) then
          allocate(vx_l(nx,ny),stat=stat)
          if (stat>0) call fatal_error('bc_force_aa_time','Could not allocate memory for vx_l',.true.)
          allocate(vx_r(nx,ny),stat=stat)
          if (stat>0) call fatal_error('bc_force_aa_time','Could not allocate memory for vx_r',.true.)
          allocate(vy_l(nx,ny),stat=stat)
          if (stat>0) call fatal_error('bc_force_aa_time','Could not allocate memory for vy_l',.true.)
          allocate(vy_r(nx,ny),stat=stat)
          if (stat>0) call fatal_error('bc_force_aa_time','Could not allocate memory for vy_r',.true.)
          allocate(vx(nx,ny),stat=stat)
          if (stat>0) call fatal_error('bc_force_aa_time','Could not allocate memory for vx',.true.)
          allocate(vy(nx,ny),stat=stat)
          if (stat>0) call fatal_error('bc_force_aa_time','Could not allocate memory for vy',.true.)
        endif
!
        first_run = .false.
!
      endif
!
      allocate(Bz0(bnx,bny),stat=stat)
      if (stat>0) call fatal_error('bc_force_aa_time','Could not allocate memory for Bz0',.true.)
!
      time_SI = t*unit_time
!
      if (t_r+delta_t <= time_SI) then
!
        if (lfirst_proc_xy) then
!  Read and distribute Bz data (in pencil shape)
!
          inquire (IOLENGTH=lend) t_l
          open (10,file=mag_times_dat,form='unformatted',status='unknown', &
              recl=lend,access='direct')
!
          ierr = 0
          t_l = 0.
          frame = 0
          do while (ierr == 0)
            frame=frame+1
            read (10,rec=frame,iostat=ierr) t_l
            read (10,rec=frame+1,iostat=ierr) t_r
            if (ierr /= 0) then
              ! EOF is reached => read again
              delta_t = time_SI
              frame=1
              read (10,rec=frame,iostat=ierr) t_l
              read (10,rec=frame+1,iostat=ierr) t_r
              ierr=-1
            else
              ! test, if correct time step is reached
              if ((t_l+delta_t < time_SI) .and. (t_r+delta_t > time_SI)) ierr = -1
            endif
          enddo
          close (10)
!
          if (luse_vel_field) then
            allocate (vx_tmp(nxgrid,nygrid), stat=stat)
            if (stat>0) call fatal_error ('bc_force_aa_time', 'Could not allocate memory for vx_tmp', .true.)
            allocate (vy_tmp(nxgrid,nygrid), stat=stat)
            if (stat>0) call fatal_error ('bc_force_aa_time', 'Could not allocate memory for vy_tmp', .true.)
            open (10, file=mag_vel_field_dat, form='unformatted', status='unknown', &
                recl=lend*nxgrid*nygrid, access='direct')
!
            rec_vxl = 1 + (frame-1)*nprocxy*2
            rec_vxr = 1 + frame*nprocxy*2
            rec_vyl = rec_vxl + frame*nprocxy
            rec_vyr = rec_vxr + frame*nprocxy
!
            ! read _l data in the order of occurence in file
            read (10,rec=rec_vxl) vx_l
            read (10,rec=rec_vyl) vy_l
!
            ! send _l data to remote
            do py = 1, nprocy
              do px = 1, nprocx
                partner = px + py*nprocx + ipz*nprocxy
                if (partner == iproc) cycle
                vx_l = vx_tmp(1+(px-1)*nprocx:px*nprocx,1+(py-1)*nprocy:py*nprocy)
                vy_l = vy_tmp(1+(px-1)*nprocx:px*nprocx,1+(py-1)*nprocy:py*nprocy)
                call mpisend_real (vx_l, (/ nx, ny /), partner, tag_l)
                call mpisend_real (vy_l, (/ nx, ny /), partner, tag_r)
              enddo
            enddo
            ! read local _l data
            vx_l = vx_tmp(1:nprocx,1:nprocy)
            vy_l = vy_tmp(1:nprocx,1:nprocy)
!
            ! read _r data in the order of occurence in file
            read (10,rec=rec_vxr) vx_r
            read (10,rec=rec_vyr) vy_r
!
            ! send _r data to remote
            do py = 1, nprocy
              do px = 1, nprocx
                partner = px + py*nprocx + ipz*nprocxy
                if (partner == iproc) cycle
                vx_r = vx_tmp(1+(px-1)*nprocx:px*nprocx,1+(py-1)*nprocy:py*nprocy)
                vy_r = vy_tmp(1+(px-1)*nprocx:px*nprocx,1+(py-1)*nprocy:py*nprocy)
                call mpisend_real (vx_r, (/ nx, ny /), partner, tag_l)
                call mpisend_real (vy_r, (/ nx, ny /), partner, tag_r)
              enddo
            enddo
            ! read local _r data
            vx_r = vx_tmp(1:nprocx,1:nprocy)
            vy_r = vy_tmp(1:nprocx,1:nprocy)
!
            close (10)
            if (allocated (vx_tmp)) deallocate (vx_tmp)
            if (allocated (vy_tmp)) deallocate (vy_tmp)
          endif
!
          open (10,file=mag_field_dat,form='unformatted',status='unknown', &
              recl=lend*bnx*bny,access='direct')
          rec_l = 1 + (frame-1)*nprocxy
          rec_r = 1 + frame*nprocxy
          do py=1, nprocxy-1
            partner = py + ipz*nprocxy
            ! read Bz data for remote processors
            read (10,rec=rec_l+py) Bz0_l
            read (10,rec=rec_r+py) Bz0_r
            ! send Bz data to remote
            call mpisend_real (Bz0_l, (/ bnx, bny /), partner, tag_l)
            call mpisend_real (Bz0_r, (/ bnx, bny /), partner, tag_r)
            call mpisend_real (t_l, partner, tag_l)
            call mpisend_real (t_r, partner, tag_r)
            call mpisend_real (delta_t, partner, tag_dt)
          enddo
          ! read local Bz data
          read (10,rec=rec_l) Bz0_l
          read (10,rec=rec_r) Bz0_r
          close (10)
!
        else
!
          if (luse_vel_field) then
            ! wait for vx and vy data from root processor
            call mpirecv_real (vx_l, (/ nx, ny /), ipz*nprocxy, tag_l)
            call mpirecv_real (vy_l, (/ nx, ny /), ipz*nprocxy, tag_r)
            call mpirecv_real (vx_r, (/ nx, ny /), ipz*nprocxy, tag_l)
            call mpirecv_real (vy_r, (/ nx, ny /), ipz*nprocxy, tag_r)
          endif
!
          ! wait for Bz data from root processor
          call mpirecv_real (Bz0_l, (/ bnx, bny /), ipz*nprocxy, tag_l)
          call mpirecv_real (Bz0_r, (/ bnx, bny /), ipz*nprocxy, tag_r)
          call mpirecv_real (t_l, ipz*nprocxy, tag_l)
          call mpirecv_real (t_r, ipz*nprocxy, tag_r)
          call mpirecv_real (delta_t, ipz*nprocxy, tag_dt)
!
        endif
!
        ! Gauss to Tesla and SI to PENCIL units
        Bz0_l = Bz0_l * 1e-4 / unit_magnetic
        Bz0_r = Bz0_r * 1e-4 / unit_magnetic
!
        if (luse_vel_field) then
          vx_l = vx_l / unit_velocity
          vy_l = vy_l / unit_velocity
          vx_r = vx_r / unit_velocity
          vy_r = vy_r / unit_velocity
        endif
!
      endif
!
      Bz0 = (time_SI - (t_l+delta_t)) * (Bz0_r - Bz0_l) / (t_r - t_l) + Bz0_l
!
      if (luse_vel_field) then
        vx = (time_SI - (t_l+delta_t)) * (vx_r - vx_l) / (t_r - t_l) + vx_l
        vy = (time_SI - (t_l+delta_t)) * (vy_r - vy_l) / (t_r - t_l) + vy_l
      endif
!
!  Fourier Transform of Bz0:
!
      if (.not. allocated (exp_fact)) then
        ! Setup exponential factor for bottom boundary
        allocate (exp_fact(enx,eny,nghost+1), stat=stat)
        if (stat > 0) call fatal_error ('bc_force_aa_time', 'Could not allocate memory for exp_fact', .true.)
        call setup_extrapol_fact (z(1:n1), z(n1), exp_fact, reduce_factor)
      endif
!
      call field_extrapol_z_parallel (Bz0, f(l1:l2,m1:m2,n1-nghost:n1,iax:iay), exp_fact)
      call communicate_vect_field_ghosts (f, 'bot')
      if (luse_vel_field) call communicate_vect_field_ghosts (f, 'bot', iux)
!
      if (allocated(Bz0)) deallocate(Bz0)
!
    endsubroutine bc_force_aa_time
!***********************************************************************
    subroutine bc_lnTT_flux_x(f,topbot)
!
!  Constant flux boundary condition for temperature (called when bcx='c1')
!
!  12-Mar-2007/dintrans: coded
!
      use SharedVariables, only: get_shared_variable
!
      real, dimension (:,:,:,:) :: f
      character (len=bclen) :: topbot
!
      real, pointer :: hcond0, hcond1, Fbot
      real, dimension (size(f,2),size(f,3)) :: tmp_yz
      integer :: i
!
!  Do the 'c1' boundary condition (constant heat flux) for lnTT.
!  check whether we want to do top or bottom (this is processor dependent)
!
      call get_shared_variable('hcond0',hcond0,caller='bc_lnTT_flux_x')
      call get_shared_variable('hcond1',hcond1)
      call get_shared_variable('Fbot',Fbot)
!
      if (headtt) print*,'bc_lnTT_flux_x: Fbot,hcond,dx=',Fbot,hcond0*hcond1,dx
!
      select case (topbot)
!
!  bottom boundary
!  ===============
!
      case ('bot')
        tmp_yz=-Fbot/(hcond0*hcond1)/exp(f(l1,:,:,ilnTT))
!
!  enforce dlnT/dx = - Fbot/(K*T)
!
        do i=1,nghost
          f(l1-i,:,:,ilnTT)=f(l1+i,:,:,ilnTT)-dx2_bound(-i)*tmp_yz
        enddo
!
      case default
        call fatal_error('bc_lnTT_flux_x','invalid argument')
!
      endselect
!
    endsubroutine bc_lnTT_flux_x
!***********************************************************************
    subroutine bc_lnTT_flux_z(f,topbot)
!
!  Constant flux boundary condition for temperature
!  (called when bcz='c1')
!
!  12-May-07/dintrans: coded
!
      use SharedVariables, only: get_shared_variable
!
      real, dimension (:,:,:,:) :: f
      character (len=bclen) :: topbot
!
      real, dimension (size(f,1),size(f,2)) :: tmp_xy
      real, pointer :: hcond0, Fbot
      integer :: i
!
!  Do the 'c1' boundary condition (constant heat flux) for lnTT or TT (if
!  ltemperature_nolog=.true.) at the bottom _only_.
!  lnTT version: enforce dlnT/dz = - Fbot/(K*T)
!    TT version: enforce   dT/dz = - Fbot/K
!
      call get_shared_variable('hcond0',hcond0,caller='bc_lnTT_flux_z')
      call get_shared_variable('Fbot',Fbot)
!
      if (headtt) print*,'bc_lnTT_flux_z: Fbot,hcond,dz=',Fbot,hcond0,dz
!
      select case (topbot)
      case ('bot')
        if (ltemperature_nolog) then
          tmp_xy=-Fbot/hcond0
        else
          tmp_xy=-Fbot/hcond0/exp(f(:,:,n1,ilnTT))
        endif
        do i=1,nghost
          f(:,:,n1-i,ilnTT)=f(:,:,n1+i,ilnTT)-dz2_bound(-i)*tmp_xy
        enddo
!
      case default
        call fatal_error('bc_lnTT_flux_z','invalid argument')
!
      endselect
!
    endsubroutine bc_lnTT_flux_z
!***********************************************************************
    subroutine bc_ss_flux_x(f,topbot)
!
!  Constant flux boundary condition for entropy (called when bcx='c1')
!
!  17-mar-07/dintrans: coded
!  16-apr-12/MR: eliminated cs2_yz; allocation of rho_yz -> work_yz only if necessary;
!                introduced heatflux_boundcond_x (necessary for nonequidistant grid)
!   5-feb-15/MR: added reference state
!  11-feb-15/MR: corrected use of reference state
!
      use EquationOfState, only: gamma, gamma_m1, lnrho0, cs20
      use SharedVariables, only: get_shared_variable
!
      real, dimension (:,:,:,:) :: f
      character (len=bclen) :: topbot
!
      real, dimension (:,:), allocatable :: tmp_yz,work_yz
      real, pointer :: FbotKbot, FtopKtop, Fbot, Ftop, cp
      real, pointer :: hcond0_kramers, nkramers
      logical, pointer :: lheatc_kramers
      integer :: i,stat
      real, dimension (:,:), pointer :: reference_state
      real :: fac
!
!  Do the 'c1' boundary condition (constant heat flux) for entropy.
!
      call get_shared_variable('lheatc_kramers',lheatc_kramers, caller='bc_ss_flux_x')
!
!  Allocate memory for large arrays.
!
      allocate(tmp_yz(size(f,2),size(f,3)),stat=stat)
      if (stat>0) call fatal_error('bc_ss_flux_x', &
          'Could not allocate memory for tmp_yz')
!
      if (lheatc_kramers) then
!
        call get_shared_variable('hcond0_kramers',hcond0_kramers)
        call get_shared_variable('nkramers',nkramers)
        call get_shared_variable('cp',cp)
!
      endif
!
      if (lheatc_kramers.or.lreference_state) then
        allocate(work_yz(size(f,2),size(f,3)),stat=stat)
        if (stat>0) call fatal_error('bc_ss_flux_x', &
                                     'Could not allocate memory for work_yz')
      endif
!
      if (lreference_state) &
        call get_shared_variable('reference_state',reference_state)
! 
      fac=gamma_m1/gamma
!
! Check whether we want to do top or bottom (this is processor dependent)
!
      select case (topbot)
!
!  bottom boundary
!  ===============
!
      case ('bot')
!
        call get_shared_variable('FbotKbot',FbotKbot)
        if ((headtt) .and. (lroot)) print*,'bc_ss_flux_x: FbotKbot=',FbotKbot
!
!  Deal with the simpler pretend_lnTT=T case first. Now ss is actually
!  lnTT and the boundary condition reads glnTT=FbotKbot/T
!
        if (pretend_lnTT) then      ! TODO: non-equidistant grid
          do i=1,nghost
            f(l1-i,:,:,iss)=f(l1+i,:,:,iss)+dx2_bound(-i)*FbotKbot/exp(f(l1,:,:,iss))
          enddo
        else
!
!  calculate Fbot/(K*cs2)
!
!       cs2_yz=cs20*exp(gamma_m1*(f(l1,:,:,ilnrho)-lnrho0)+cv1*f(l1,:,:,iss))
!
!  Both, bottom and top boundary conditions are corrected for linear density
!
          if (ldensity_nolog) then
            if (lheatc_kramers) work_yz=f(l1,:,:,irho)
            if (lreference_state) then
              tmp_yz= cs20*exp(gamma_m1*(log(f(l1,:,:,irho)+reference_state(1,iref_rho))-lnrho0)  &
                     +gamma*(f(l1,:,:,iss)+reference_state(1,iref_s)))
            else
              tmp_yz=cs20*exp(gamma_m1*(log(f(l1,:,:,irho))-lnrho0)+gamma*f(l1,:,:,iss))
            endif
          else
            if (lheatc_kramers) work_yz=exp(f(l1,:,:,ilnrho))
!print*, 'bc_ss_flux_x: iproc, lnrho, ss=', iproc, maxval(f(l1,:,:,ilnrho)), &
!minval(f(l1,:,:,ilnrho)), maxval(f(l1,:,:,iss)), minval(f(l1,:,:,iss))
            tmp_yz=cs20*exp(gamma_m1*(f(l1,:,:,ilnrho)-lnrho0)+gamma*f(l1,:,:,iss))
          endif
          if (lheatc_kramers) then
!
            call get_shared_variable('Fbot',Fbot)
            if ((headtt) .and. (lroot)) print*,'bc_ss_flux_x: Fbot=',Fbot
!
            tmp_yz = Fbot*work_yz**(2*nkramers)*(cp*gamma_m1)**(6.5*nkramers)/ &
                     (hcond0_kramers*tmp_yz**(6.5*nkramers+1.))
!
          else
            tmp_yz=FbotKbot/tmp_yz
          endif
!
!  enforce ds/dx           + gamma_m1/gamma*dlnrho/dx = - gamma_m1/gamma*Fbot/(K*cs2)
!  or with reference state:
!                + ds_0/dx + gamma_m1/gamma*d/dx(ln(rho'+rho_0)) = - gamma_m1/gamma*Fbot/(K*cs2)
!
          if (lreference_state) then
            work_yz= 1./(f(l1,:,:,irho)+reference_state(1,iref_rho))
            tmp_yz = tmp_yz + reference_state(1,iref_gs)/fac + reference_state(1,iref_grho)*work_yz
            call heatflux_boundcond_x( f, tmp_yz, fac, BOT, work_yz )
          else
            call heatflux_boundcond_x( f, tmp_yz, fac, BOT )
          endif
!
        endif
!
!  top boundary
!  ============
!
      case ('top')
!
        call get_shared_variable('FtopKtop',FtopKtop)
         if ((headtt) .and. (lroot)) print*,'bc_ss_flux_x: FtopKtop=',FtopKtop
!
!  Deal with the simpler pretend_lnTT=T case first. Now ss is actually
!  lnTT and the boundary condition reads glnTT=FtopKtop/T
!
        if (pretend_lnTT) then
          do i=1,nghost
            f(l2+i,:,:,iss)=f(l2-i,:,:,iss)-dx2_bound(i)*FtopKtop/exp(f(l2,:,:,iss))
          enddo
        else
!
!  calculate Ftop/(K*cs2)
!
          if (ldensity_nolog) then
            if (lheatc_kramers) work_yz=f(l2,:,:,irho)
            if (lreference_state) then
              tmp_yz=cs20*exp(gamma_m1*(log(f(l2,:,:,irho)+reference_state(nx,iref_rho))-lnrho0) &
                     +gamma*(f(l2,:,:,iss)+reference_state(nx,iref_s)))
            else
              tmp_yz=cs20*exp(gamma_m1*(log(f(l2,:,:,irho))-lnrho0)+gamma*f(l2,:,:,iss))
            endif
          else
            if (lheatc_kramers) work_yz=exp(f(l2,:,:,ilnrho))
            tmp_yz=cs20*exp(gamma_m1*(f(l2,:,:,ilnrho)-lnrho0)+gamma*f(l2,:,:,iss))
          endif
          if (lheatc_kramers) then
!
            call get_shared_variable('Ftop',Ftop)
            if ((headtt) .and. (lroot)) print*,'bc_ss_flux_x: Ftop=',Ftop
!
            tmp_yz = Ftop*work_yz**(2*nkramers)*(cp*gamma_m1)**(6.5*nkramers)/ &
                     (hcond0_kramers*tmp_yz**(6.5*nkramers+1.))
          else
            tmp_yz=FtopKtop/tmp_yz
          endif
!
          if (lreference_state) &
            tmp_yz = tmp_yz + reference_state(nx,iref_gs)
!
!  enforce ds/dx + gamma_m1/gamma*dlnrho/dx = gamma_m1/gamma*Ftop/(K*cs2)   ! check sign
!
          if (lreference_state) then
            work_yz= 1./(f(l2,:,:,irho)+reference_state(nx,iref_rho))
            tmp_yz = tmp_yz + reference_state(nx,iref_gs)/fac + reference_state(nx,iref_grho)*work_yz
            call heatflux_boundcond_x( f, -tmp_yz, fac, TOP, work_yz )     ! check sign of tmp_yz
          else
            call heatflux_boundcond_x( f, -tmp_yz, fac, TOP )
          endif
!
        endif
!
        case default
          call fatal_error('bc_ss_flux_x','invalid argument')
!
        endselect
!
!  Deallocate large arrays.
!
      if (allocated(tmp_yz)) deallocate(tmp_yz)
      if (allocated(work_yz)) deallocate(work_yz)
!
    endsubroutine bc_ss_flux_x
!************************************************************************
    subroutine heatflux_boundcond_x( f, inh, fac, topbot, coef )
!
!  encapsules BC 'prescribed heat flux at x boundary'
!
!  17-apr-12/MR: outsourced from bc_ss_flux_x
!  11-feb-15/MR: optional parameter coef for correct use of reference state
!
      use Deriv, only: heatflux_deriv_x
!
      real, dimension(:,:,:,:),       intent(INOUT):: f
      real, dimension(:,:),           intent(IN)   :: inh
      real, dimension(:,:), optional, intent(IN)   :: coef
      real                          , intent(IN)   :: fac
      integer                       , intent(IN)   :: topbot
!
      integer :: i,ll,ia,ie
!
      if ( .not.lequidist(1) ) then
        if ( heatflux_deriv_x( f, inh, fac, topbot ) ) return
      endif
!
      if (topbot==BOT) then
        ll=l1; ia=1; ie=nghost
      else
        ll=l2; ia=-nghost; ie=-1
      endif

      do i=ia,ie
        if (ldensity_nolog) then
          if (present(coef)) then
            f(ll-i,:,:,iss)=f(ll+i,:,:,iss)+fac* &
                ( (f(ll+i,:,:,irho)-f(ll-i,:,:,irho))*coef + dx2_bound(-i)*inh )
          else
            f(ll-i,:,:,iss)=f(ll+i,:,:,iss)+fac* &
                (log(f(ll+i,:,:,irho)/f(ll-i,:,:,irho)) + dx2_bound(-i)*inh)
          endif
        else
          f(ll-i,:,:,iss)=f(ll+i,:,:,iss)+fac* &
              (f(ll+i,:,:,ilnrho)-f(ll-i,:,:,ilnrho) + dx2_bound(-i)*inh)
        endif
      enddo
!
  endsubroutine heatflux_boundcond_x
!***********************************************************************
    subroutine bc_del2zero(f,topbot,j)
!
!  Potential field boundary condition
!
!  11-oct-06/wolf: Adapted from Tobi's bc_aa_pot2
!
      use Fourier, only: fourier_transform_xy_xy, kx_fft, ky_fft
!
      real, dimension (:,:,:,:), intent (inout) :: f
      character (len=bclen), intent (in) :: topbot
      integer, intent (in) :: j
!
      real, dimension (l2-l1+1,m2-m1+1) :: kx,ky,kappa,exp_fact,tmp_re,tmp_im
      integer :: i,nxl,nyl
!
!  Allocate memory for large arrays.
!
      nxl=l2-l1+1; nyl=m2-m1+1
!
!  Get local wave numbers
!
      kx = spread(kx_fft(ipx*nxl+1:ipx*nxl+nxl),2,nyl)
      ky = spread(ky_fft(ipy*nyl+1:ipy*nyl+nyl),1,nxl)
!
!  Calculate 1/k^2, zero mean
!
      if (lshear) then
        kappa = sqrt((kx+deltay*ky/Lx)**2+ky**2)
      else
        kappa = sqrt(kx**2 + ky**2)
      endif
!
!  Check whether we want to do top or bottom (this is processor dependent)
!
      select case (topbot)
!
!  Potential field condition at the bottom
!
      case ('bot')
!
        do i=1,nghost
!
! Calculate delta_z based on z(), not on dz to improve behavior for
! non-equidistant grid (still not really correct, but could be OK)
!
          exp_fact = exp(-kappa*dz2_bound(-i))
!
!  Determine potential field in ghost zones
!
          !  Fourier transforms of x- and y-components on the boundary
          tmp_re = f(l1:l2,m1:m2,n1+i,j)
          tmp_im = 0.0
          call fourier_transform_xy_xy(tmp_re,tmp_im)
          tmp_re = tmp_re*exp_fact
          tmp_im = tmp_im*exp_fact
          ! Transform back
          call fourier_transform_xy_xy(tmp_re,tmp_im,linv=.true.)
          f(l1:l2,m1:m2,n1-i,j) = tmp_re
!
        enddo
!
!  Potential field condition at the top
!
      case ('top')
!
        do i=1,nghost
!
! Calculate delta_z based on z(), not on dz to improve behavior for
! non-equidistant grid (still not really correct, but could be OK)
!
          exp_fact = exp(-kappa*dz2_bound(i))
!
!  Determine potential field in ghost zones
!
          !  Fourier transforms of x- and y-components on the boundary
          tmp_re = f(l1:l2,m1:m2,n2-i,j)
          tmp_im = 0.0
          call fourier_transform_xy_xy(tmp_re,tmp_im)
          tmp_re = tmp_re*exp_fact
          tmp_im = tmp_im*exp_fact
          ! Transform back
          call fourier_transform_xy_xy(tmp_re,tmp_im,linv=.true.)
          f(l1:l2,m1:m2,n2+i,j) = tmp_re
!
        enddo
!
      case default
!
        if (lroot) print*,"bc_del2zero: invalid argument"
!
      endselect
!
    endsubroutine bc_del2zero
!***********************************************************************
    subroutine bc_zero_x(f,topbot,j)
!
!  Zero value in the ghost zones.
!
!  11-aug-2009/anders: implemented
!
      real, dimension (:,:,:,:) :: f
      character (len=bclen) :: topbot
      integer :: j
!
      select case (topbot)
!
!  Bottom boundary.
!
      case ('bot')
        f(1:l1-1,:,:,j)=0.0
!
!  Top boundary.
!
      case ('top')
        f(l2+1:,:,:,j)=0.0
!
!  Default.
!
      case default
        print*, "bc_zero_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_zero_x
!***********************************************************************
    subroutine bc_zero_y(f,topbot,j)
!
!  Zero value in the ghost zones.
!
!  13-jul-2011/Tijmen: adapted from bc_zero_x
!
      real, dimension (:,:,:,:) :: f
      character (len=bclen) :: topbot
      integer :: j
!
      select case (topbot)
!
!  Bottom boundary.
!
      case ('bot')
        f(:,1:m1-1,:,j)=0.0
!
!  Top boundary.
!
      case ('top')
        f(:,m2+1:,:,j)=0.0
!
!  Default.
!
      case default
        print*, "bc_zero_y: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_zero_y
!***********************************************************************
    subroutine bc_zero_z(f,topbot,j)
!
!  Zero value in the ghost zones.
!
!  13-aug-2007/anders: implemented
!
      real, dimension (:,:,:,:) :: f
      character (len=bclen) :: topbot
      integer :: j
!
      select case (topbot)
!
!  Bottom boundary.
!
      case ('bot')
        f(:,:,1:n1-1,j)=0.0
!
!  Top boundary.
!
      case ('top')
        f(:,:,n2+1:,j)=0.0
!
!  Default.
!
      case default
        print*, "bc_zero_z: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_zero_z
!***********************************************************************
    subroutine bc_inflow_z(f,topbot,j,lforce_ghost)
!
!  Inflow boundary conditions.
!
!  If the velocity vector points out of the box, the velocity boundary
!  condition is set to 's', otherwise it is set to 'a'.
!  If 'lforce_ghost' is true, the boundary and ghost cell values are forced
!  to not point outwards. Otherwise the boundary value is forced to be 0.
!
!  25-dec-2010/Bourdin.KIS: adapted from 'bc_outflow_z'
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j
      logical, optional :: lforce_ghost
!
      integer :: i, ix, iy
      logical :: lforce
!
      lforce = .false.
      if (present (lforce_ghost)) lforce = lforce_ghost
!
      select case (topbot)
!
!  Bottom boundary.
!
      case ('bot')
        do iy=1,size(f,2); do ix=1,size(f,1)
          if (f(ix,iy,n1,j)>0.0) then  ! 's'
            do i=1,nghost; f(ix,iy,n1-i,j)=+f(ix,iy,n1+i,j); enddo
          else                         ! 'a'
            do i=1,nghost; f(ix,iy,n1-i,j)=-f(ix,iy,n1+i,j); enddo
            f(ix,iy,n1,j)=0.0
          endif
          if (lforce) then
            do i = 0, nghost
              if (f(ix,iy,n1-i,j) < 0.0) f(ix,iy,n1-i,j) = 0.0
            enddo
          endif
        enddo; enddo
!
!  Top boundary.
!
      case ('top')
        do iy=1,size(f,2); do ix=1,size(f,1)
          if (f(ix,iy,n2,j)<0.0) then  ! 's'
            do i=1,nghost; f(ix,iy,n2+i,j)=+f(ix,iy,n2-i,j); enddo
          else                         ! 'a'
            do i=1,nghost; f(ix,iy,n2+i,j)=-f(ix,iy,n2-i,j); enddo
            f(ix,iy,n2,j)=0.0
          endif
          if (lforce) then
            do i = 0, nghost
              if (f(ix,iy,n2+i,j) > 0.0) f(ix,iy,n2+i,j) = 0.0
            enddo
          endif
        enddo; enddo
!
!  Default.
!
      case default
        print*, "bc_inflow_z: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_inflow_z
!***********************************************************************
    subroutine bc_outflow_x(f,topbot,j,lforce_ghost)
!
!  Outflow boundary conditions.
!
!  If the velocity vector points out of the box, the velocity boundary
!  condition is set to 's', otherwise it is set to 'a'.
!  If 'lforce_ghost' is true, the boundary and ghost cell values are forced
!  to not point inwards. Otherwise the boundary value is forced to be 0.
!
!  14-jun-2011/axel: adapted from bc_outflow_z
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j
      logical, optional :: lforce_ghost
!
      integer :: i, iy, iz
      logical :: lforce
!
      lforce = .false.
      if (present (lforce_ghost)) lforce = lforce_ghost
!
      select case (topbot)
!
!  Bottom boundary.
!
      case ('bot')
        do iy=1,size(f,2); do iz=1,size(f,3)
          if (f(l1,iy,iz,j)<0.0) then  ! 's'
            do i=1,nghost; f(l1-i,iy,iz,j)=+f(l1+i,iy,iz,j); enddo
          else                         ! 'a'
            do i=1,nghost; f(l1-i,iy,iz,j)=-f(l1+i,iy,iz,j); enddo
            f(l1,iy,iz,j)=0.0
          endif
          if (lforce) then
            do i = 0, nghost
              if (f(l1-i,iy,iz,j) > 0.0) f(l1-i,iy,iz,j) = 0.0
            enddo
          endif
        enddo; enddo
!
!  Top boundary.
!
      case ('top')
        do iy=1,size(f,2); do iz=1,size(f,3)
          if (f(l2,iy,iz,j)>0.0) then  ! 's'
            do i=1,nghost; f(l2+i,iy,iz,j)=+f(l2-i,iy,iz,j); enddo
          else                         ! 'a'
            do i=1,nghost; f(l2+i,iy,iz,j)=-f(l2-i,iy,iz,j); enddo
            f(l2,iy,iz,j)=0.0
          endif
          if (lforce) then
            do i = 0, nghost
              if (f(l2+i,iy,iz,j) < 0.0) f(l2+i,iy,iz,j) = 0.0
            enddo
          endif
        enddo; enddo
!
!  Default.
!
      case default
        print*, "bc_outflow_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_outflow_x
!***********************************************************************
    subroutine bc_outflow_x_e1(f,topbot,j,lforce_ghost)
!
!  Outflow boundary conditions.
!
!  If the velocity vector points out of the box, the velocity boundary
!  condition is set to 's', otherwise it is set to 'a'.
!  If 'lforce_ghost' is true, the boundary and ghost cell values are forced
!  to not point inwards. Otherwise the boundary value is forced to be 0.
!
!  14-jun-2011/axel: adapted from bc_outflow_x
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j
      logical, optional :: lforce_ghost
!
      integer :: i, iy, iz
      logical :: lforce
!
      lforce = .false.
      if (present (lforce_ghost)) lforce = lforce_ghost
!
      select case (topbot)
!
!  Bottom boundary.
!
      case ('bot')
        do iy=1,size(f,2); do iz=1,size(f,3)
          if (f(l1,iy,iz,j)<0.0) then  ! 's'
            do i=1,nghost; f(l1-i,iy,iz,j)=+f(l1+i,iy,iz,j); enddo
            f(l1-1,iy,iz,j)=0.25*(  9*f(l1,iy,iz,j)- 3*f(l1+1,iy,iz,j)- 5*f(l1+2,iy,iz,j)+ 3*f(l1+3,iy,iz,j))
            f(l1-2,iy,iz,j)=0.05*( 81*f(l1,iy,iz,j)-43*f(l1+1,iy,iz,j)-57*f(l1+2,iy,iz,j)+39*f(l1+3,iy,iz,j))
            f(l1-3,iy,iz,j)=0.05*(127*f(l1,iy,iz,j)-81*f(l1+1,iy,iz,j)-99*f(l1+2,iy,iz,j)+73*f(l1+3,iy,iz,j))
          else                         ! 'a'
            do i=1,nghost; f(l1-i,iy,iz,j)=-f(l1+i,iy,iz,j); enddo
            f(l1,iy,iz,j)=0.0
          endif
          if (lforce) then
            do i = 0, nghost
              if (f(l1-i,iy,iz,j) > 0.0) f(l1-i,iy,iz,j) = 0.0
            enddo
          endif
        enddo; enddo
!
!  Top boundary.
!
      case ('top')
        do iy=1,size(f,2); do iz=1,size(f,3)
          if (f(l2,iy,iz,j)>0.0) then  ! 's'
            f(l2+1,iy,iz,j)=0.25*(  9*f(l2,iy,iz,j)- 3*f(l2-1,iy,iz,j)- 5*f(l2-2,iy,iz,j)+ 3*f(l2-3,iy,iz,j))
            f(l2+2,iy,iz,j)=0.05*( 81*f(l2,iy,iz,j)-43*f(l2-1,iy,iz,j)-57*f(l2-2,iy,iz,j)+39*f(l2-3,iy,iz,j))
            f(l2+3,iy,iz,j)=0.05*(127*f(l2,iy,iz,j)-81*f(l2-1,iy,iz,j)-99*f(l2-2,iy,iz,j)+73*f(l2-3,iy,iz,j))
          else                         ! 'a'
            do i=1,nghost; f(l2+i,iy,iz,j)=-f(l2-i,iy,iz,j); enddo
            f(l2,iy,iz,j)=0.0
          endif
          if (lforce) then
            do i = 0, nghost
              if (f(l2+i,iy,iz,j) < 0.0) f(l2+i,iy,iz,j) = 0.0
            enddo
          endif
        enddo; enddo
!
!  Default.
!
      case default
        print*, "bc_outflow_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_outflow_x_e1
!***********************************************************************
    subroutine bc_outflow_y(f,topbot,j,lforce_ghost)
!
!  Outflow boundary conditions.
!
!  If the velocity vector points out of the box, the velocity boundary
!  condition is set to 's', otherwise it is set to 'a'.
!  If 'lforce_ghost' is true, the boundary and ghost cell values are forced
!  to not point inwards. Otherwise the boundary value is forced to be 0.
!
!  08-oct-2013/wlad: copied from z
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j
      logical, optional :: lforce_ghost
!
      integer :: i, ix, iz
      logical :: lforce
!
      lforce = .false.
      if (present (lforce_ghost)) lforce = lforce_ghost
!
      select case (topbot)
!
!  Bottom boundary.
!
      case ('bot')
        do iz=1,size(f,3); do ix=1,size(f,1)
          if (f(ix,m1,iz,j)<0.0) then  ! 's'
            do i=1,nghost; f(ix,m1-i,iz,j)=+f(ix,m1+i,iz,j); enddo
          else                         ! 'a'
            do i=1,nghost; f(ix,m1-i,iz,j)=-f(ix,m1+i,iz,j); enddo
            f(ix,m1,iz,j)=0.0
          endif
          if (lforce) then
            do i = 0, nghost
              if (f(ix,m1-i,iz,j) > 0.0) f(ix,m1-i,iz,j) = 0.0
            enddo
          endif
        enddo; enddo
!
!  Top boundary.
!
      case ('top')
        do iz=1,size(f,3); do ix=1,size(f,1)
          if (f(ix,m2,iz,j)>0.0) then  ! 's'
            do i=1,nghost; f(ix,m2+i,iz,j)=+f(ix,m2-i,iz,j); enddo
          else                         ! 'a'
            do i=1,nghost; f(ix,m2+i,iz,j)=-f(ix,m2-i,iz,j); enddo
            f(ix,m2,iz,j)=0.0
          endif
          if (lforce) then
            do i = 0, nghost
              if (f(ix,m2+i,iz,j) < 0.0) f(ix,m2+i,iz,j) = 0.0
            enddo
          endif
        enddo; enddo
!
!  Default.
!
      case default
        print*, "bc_outflow_y: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_outflow_y
!***********************************************************************
    subroutine bc_outflow_z(f,topbot,j,lforce_ghost)
!
!  Outflow boundary conditions.
!
!  If the velocity vector points out of the box, the velocity boundary
!  condition is set to 's', otherwise it is set to 'a'.
!  If 'lforce_ghost' is true, the boundary and ghost cell values are forced
!  to not point inwards. Otherwise the boundary value is forced to be 0.
!
!  12-aug-2007/anders: implemented
!  25-dec-2010/Bourdin.KIS: added forcing of boundary and ghost cell values
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j
      logical, optional :: lforce_ghost
!
      integer :: i, ix, iy
      logical :: lforce
!
      lforce = .false.
      if (present (lforce_ghost)) lforce = lforce_ghost
!
      select case (topbot)
!
!  Bottom boundary.
!
      case ('bot')
        do iy=1,size(f,2); do ix=1,size(f,1)
          if (f(ix,iy,n1,j)<0.0) then  ! 's'
            do i=1,nghost; f(ix,iy,n1-i,j)=+f(ix,iy,n1+i,j); enddo
          else                         ! 'a'
            do i=1,nghost; f(ix,iy,n1-i,j)=-f(ix,iy,n1+i,j); enddo
            f(ix,iy,n1,j)=0.0
          endif
          if (lforce) then
            do i = 0, nghost
              if (f(ix,iy,n1-i,j) > 0.0) f(ix,iy,n1-i,j) = 0.0
            enddo
          endif
        enddo; enddo
!
!  Top boundary.
!
      case ('top')
        do iy=1,size(f,2); do ix=1,size(f,1)
          if (f(ix,iy,n2,j)>0.0) then  ! 's'
            do i=1,nghost; f(ix,iy,n2+i,j)=+f(ix,iy,n2-i,j); enddo
          else                         ! 'a'
            do i=1,nghost; f(ix,iy,n2+i,j)=-f(ix,iy,n2-i,j); enddo
            f(ix,iy,n2,j)=0.0
          endif
          if (lforce) then
            do i = 0, nghost
              if (f(ix,iy,n2+i,j) < 0.0) f(ix,iy,n2+i,j) = 0.0
            enddo
          endif
        enddo; enddo
!
!  Default.
!
      case default
        print*, "bc_outflow_z: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_outflow_z
!***********************************************************************
    subroutine bc_inflow_zero_deriv_z(f,topbot,j)
!
!  Inflow boundary condition, tries to create zero 1st derivative at boundary.
!
!  If the velocity vector points out of the box, the velocity boundary
!  condition is set to 's', otherwise it is set to relax to zero derivative.
!  The boundary and ghost cell values are forced to not point outwards.
!
!  27-dec-2010/Bourdin.KIS: adapted from 'bc_outflow_const_deriv_z'
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j
!
      integer :: i, ix, iy
!
      select case (topbot)
!
      ! bottom boundary
      case ('bot')
        do iy = 1, size(f,2)
          do ix = 1, size(f,1)
            ! 's' boundary condition
            do i = 1, nghost
              f(ix,iy,n1-i,j) = f(ix,iy,n1+i,j)
            enddo
            ! force inflow
            do i = 0, nghost
              if (f(ix,iy,n1-i,j) < 0.0) f(ix,iy,n1-i,j) = 0.0
            enddo
          enddo
        enddo
!
      ! top boundary
      case ('top')
        do iy = 1, size(f,2)
          do ix = 1, size(f,1)
            ! 's' boundary condition
            do i = 1, nghost
              f(ix,iy,n2+i,j) = f(ix,iy,n2-i,j)
            enddo
            ! force inflow
            do i = 0, nghost
              if (f(ix,iy,n2+i,j) > 0.0) f(ix,iy,n2+i,j) = 0.0
            enddo
          enddo
        enddo
!
      case default
        print*, "bc_inflow_zero_deriv_z: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_inflow_zero_deriv_z
!***********************************************************************
    subroutine bc_outflow_zero_deriv_z(f,topbot,j)
!
!  Outflow boundary condition, tries to create zero 1st derivative at boundary.
!
!  If the velocity vector points out of the box, the velocity boundary
!  condition is set to 's', otherwise it is set to relax to zero derivative.
!  The boundary and ghost cell values are forced to not point inwards.
!
!  27-dec-2010/Bourdin.KIS: adapted from 'bc_outflow_z'
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j
!
      integer :: i, ix, iy
!
      select case (topbot)
!
      ! bottom boundary
      case ('bot')
        do iy = 1, size(f,2)
          do ix = 1, size(f,1)
            ! 's' boundary condition
            do i = 1, nghost
              f(ix,iy,n1-i,j) = f(ix,iy,n1+i,j)
            enddo
            ! force outflow
            do i = 0, nghost
              if (f(ix,iy,n1-i,j) > 0.0) f(ix,iy,n1-i,j) = 0.0
            enddo
          enddo
        enddo
!
      ! top boundary
      case ('top')
        do iy = 1, size(f,2)
          do ix = 1, size(f,1)
            ! 's' boundary condition
            do i = 1, nghost
              f(ix,iy,n2+i,j) = f(ix,iy,n2-i,j)
            enddo
            ! force outflow
            do i = 0, nghost
              if (f(ix,iy,n2+i,j) < 0.0) f(ix,iy,n2+i,j) = 0.0
            enddo
          enddo
        enddo
!
      case default
        print*, "bc_outflow_zero_deriv_z: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_outflow_zero_deriv_z
!***********************************************************************
    subroutine bc_inflow_inwards_deriv_z(f,topbot,j)
!
!  Inflow boundary condition with inwards 1st derivative at boundary.
!
!  The velocity boundary condition is set to 's' or 'a' for steady flows.
!  The boundary and ghost cell values are forced to not point outwards.
!
!  10-jul-2012/Bourdin.KIS: adapted from 'bc_inflow_zero_deriv_z'
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j
!
      integer :: i, ix, iy
!
      select case (topbot)
!
      ! bottom boundary
      case ('bot')
        do iy = 1, size(f,2)
          do ix = 1, size(f,1)
            ! 's' or 'a' boundary condition for forced inflow
            do i = 1, nghost
              f(ix,iy,n1-i,j) = abs (f(ix,iy,n1+i,j))
            enddo
            if (f(ix,iy,n1,j) < 0.0) f(ix,iy,n1,j) = 0.0
          enddo
        enddo
!
      ! top boundary
      case ('top')
        do iy = 1, size(f,2)
          do ix = 1, size(f,1)
            ! 's' or 'a' boundary condition for forced inflow
            do i = 1, nghost
              f(ix,iy,n2+i,j) = -abs (f(ix,iy,n2-i,j))
            enddo
            if (f(ix,iy,n2,j) > 0.0) f(ix,iy,n2,j) = 0.0
          enddo
        enddo
!
      case default
        print*, "bc_inflow_inwards_deriv_z: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_inflow_inwards_deriv_z
!***********************************************************************
    subroutine bc_outflow_outwards_deriv_z(f,topbot,j)
!
!  Outflow boundary condition with outwards 1st derivative at boundary.
!
!  The velocity boundary condition is set to 's' or 'a' for steady flows.
!  The boundary and ghost cell values are forced to not point inwards.
!
!  10-jul-2012/Bourdin.KIS: adapted from 'bc_outflow_zero_deriv_z'
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j
!
      integer :: i, ix, iy
!
      select case (topbot)
!
      ! bottom boundary
      case ('bot')
        do iy = 1, size(f,2)
          do ix = 1, size(f,1)
            ! 's' or 'a' boundary condition for forced outflow
            do i = 1, nghost
              f(ix,iy,n1-i,j) = -abs (f(ix,iy,n1+i,j))
            enddo
            if (f(ix,iy,n1,j) > 0.0) f(ix,iy,n1,j) = 0.0
          enddo
        enddo
!
      ! top boundary
      case ('top')
        do iy = 1, size(f,2)
          do ix = 1, size(f,1)
            ! 's' or 'a' boundary condition for forced outflow
            do i = 1, nghost
              f(ix,iy,n2+i,j) = abs (f(ix,iy,n2-i,j))
            enddo
            if (f(ix,iy,n2,j) < 0.0) f(ix,iy,n2,j) = 0.0
          enddo
        enddo
!
      case default
        print*, "bc_outflow_outwards_deriv_z: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_outflow_outwards_deriv_z
!***********************************************************************
    subroutine bc_steady_z(f,topbot,j)
!
!  Steady in/outflow boundary conditions.
!
!  Match ghost to outward velocity on boundary. Impose positive inward
!  gradient in ghost zones for inflow on boundary.
!
!  06-nov-2010/fred: implemented
!  14-mar-2011/fred: amended
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j
!
      integer :: i, ix, iy
!
      select case (topbot)
!
!  Bottom boundary.
!
      case ('bot')
        do iy=1,size(f,2); do ix=1,size(f,1)
          if (f(ix,iy,n1,j) <= 0.0) then
            do i=1,nghost; f(ix,iy,n1-i,j)=f(ix,iy,n1,j); enddo
          else
            if (f(ix,iy,n1,j) > f(ix,iy,n1+1,j)) then
              f(ix,iy,n1-1,j)=0.5*(f(ix,iy,n1,j)    +f(ix,iy,n1+1,j))
            else
              f(ix,iy,n1-1,j)=2.0* f(ix,iy,n1,j)    -f(ix,iy,n1+1,j)
            endif
            do i=2,nghost
              f(ix,iy,n1-i,j)=2.0* f(ix,iy,n1-i+1,j)-f(ix,iy,n1-i+2,j)
            enddo
          endif
        enddo; enddo
!
!  Top boundary.
!
      case ('top')
        do iy=1,size(f,2); do ix=1,size(f,1)
          if (f(ix,iy,n2,j) >= 0.0) then
            do i=1,nghost; f(ix,iy,n2+i,j)=f(ix,iy,n2,j); enddo
          else
            if (f(ix,iy,n2,j) < f(ix,iy,n2-1,j)) then
              f(ix,iy,n2+1,j)=0.5*(f(ix,iy,n2,j)    +f(ix,iy,n2-1,j))
            else
              f(ix,iy,n2+1,j)=2.0* f(ix,iy,n2,j)    -f(ix,iy,n2-1,j)
            endif
            do i=2,nghost
              f(ix,iy,n2+i,j)=2.0* f(ix,iy,n2+i-1,j)-f(ix,iy,n2+i-2,j)
            enddo
          endif
        enddo; enddo
!
!  Default.
!
      case default
        print*, "bc_steady_z: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_steady_z
!***********************************************************************
    subroutine bc_copy_x(f,topbot,j)
!
!  Copy value in last grid point to all ghost cells.
!
!  11-aug-2009/anders: implemented
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j
!
      integer :: i
!
      select case (topbot)
!
!  Bottom boundary.
!
      case ('bot')
        do i=1,nghost; f(l1-i,:,:,j)=f(l1,:,:,j); enddo
!
!  Top boundary.
!
      case ('top')
        do i=1,nghost; f(l2+i,:,:,j)=f(l2,:,:,j); enddo
!
!  Default.
!
      case default
        print*, "bc_copy_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_copy_x
!***********************************************************************
    subroutine bc_copy_y(f,topbot,j)
!
!  Copy value in last grid point to all ghost cells.
!
!  08-june-2010/wlyra: implemented
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j
!
      integer :: i
!
      select case (topbot)
!
!  Bottom boundary.
!
      case ('bot')
        do i=1,nghost; f(:,m1-i,:,j)=f(:,m1,:,j); enddo
!
!  Top boundary.
!
      case ('top')
        do i=1,nghost; f(:,m2+i,:,j)=f(:,m2,:,j); enddo
!
!  Default.
!
      case default
        print*, "bc_copy_y: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_copy_y
!***********************************************************************
    subroutine bc_copy_y_noinflow(f,topbot,j)
!
!  Copy value in last grid point to all ghost cells. Set to zero if
!  the sign is wrong. This bc is different from outflow (cop). Outflow
!  just copies the last point to the ghost cells, thus permitting both
!  outflow (uy pointing out of the box) and inflow (uy pointing back to
!  the box). 'c+k' is a no-inflow, purely outflow boundary. It sets the
!  velocity to zero if that was pointing back to the box. The 'k' means
!  "kill". "copy if outflow, kill if inflow".
!
!  08-june-2010/wlyra: implemented
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      real :: value
      integer :: j,l,n
!
      integer :: i
!
      select case (topbot)
!
!  Bottom boundary.
!
      case ('bot')
        do l=1,size(f,1); do n=1,size(f,3)
          value=0.
          if (f(l,m1,n,j)<0) value=f(l,m1,n,j)
          do i=1,nghost
            f(l,m1-i,n,j)=value
          enddo
        enddo;enddo
!
!  Top boundary.
!
      case ('top')
        do l=1,size(f,1); do n=1,size(f,3)
          value=0.
          if (f(l,m2,n,j) > 0) value=f(l,m2,n,j)
          do i=1,nghost
            f(l,m2+i,n,j)=value
          enddo
        enddo; enddo
!
!  Default.
!
      case default
        print*, "bc_copy_y_noinflow: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_copy_y_noinflow
!***********************************************************************
    subroutine bc_copy_z(f,topbot,j)
!
!  Copy value in last grid point to all ghost cells.
!
!  15-aug-2007/anders: implemented
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j
!
      integer :: i
!
      select case (topbot)
!
!  Bottom boundary.
!
      case ('bot')
        forall(i=1:nghost) f(:,:,n1-i,j) = f(:,:,n1,j)
!
!  Top boundary.
!
      case ('top')
        forall(i=1:nghost) f(:,:,n2+i,j) = f(:,:,n2,j)
!
!  Default.
!
      case default
        print*, "bc_copy_z: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_copy_z
!***********************************************************************
    subroutine bc_copy_z_noinflow(f,topbot,j)
!
!  Copy value in last grid point to all ghost cells. Set to zero if
!  the sign is wrong. This bc is different from outflow (cop). Outflow
!  just copies the last point to the ghost cells, thus permitting both
!  outflow (uy pointing out of the box) and inflow (uy pointing back to
!  the box). 'crk' is a no-inflow, purely outflow boundary. It sets the
!  velocity to zero if that was pointing back to the box. The 'k' means
!  "kill". "copy amd reduce if outflow, kill if inflow". Additionally the velocity
!  in the ghost zones are reduced by a factor 
!  2i, where i is the i-th ghost zone
!
!  22-mar-2018/piyali: copied from bc_copy_z_noinflow
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j
!
      real :: value
      integer :: i,l
!
      select case (topbot)
!
!  Bottom boundary.
!
      case ('bot')
        do l=1,size(f,1); do m=1,size(f,2)
          value=0.
          if (f(l,m,n1,j)<0) value=f(l,m,n1,j)
          do i=1,nghost
            f(l,m,n1-i,j)=value/(1.0*i)
          enddo
        enddo;enddo
!
!  Top boundary.
!
      case ('top')
        do l=1,size(f,1); do m=1,size(f,2)
          do i=1,nghost
            value=0.
            if (f(l,m,n2,j) > 0) value=f(l,m,n2-i,j)
            f(l,m,n2+i,j)=value/(1.0*i)
          enddo
        enddo; enddo
!
!  Default.
!
      case default
        print*, "bc_copy_z_noinflow: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_copy_z_noinflow
!***********************************************************************
    subroutine bc_frozen_in_bb(topbot,j)
!
!  Set flags to indicate that magnetic flux is frozen-in at the
!  boundary. The implementation occurs in daa_dt where magnetic
!  diffusion is switched off in the corresponding layer.
!
      use SharedVariables, only: get_shared_variable
!
      character (len=bclen) :: topbot
      integer :: j
!
      logical, save :: lfirstcall=.true.
      logical, pointer, save, dimension (:) :: lfrozen_bb_bot, lfrozen_bb_top
!
      if (lfirstcall) then
        call get_shared_variable('lfrozen_bb_bot',lfrozen_bb_bot,caller='bc_frozen_in_bb')
        call get_shared_variable('lfrozen_bb_top',lfrozen_bb_top)
      endif
!
      select case (topbot)
      case ('bot')               ! bottom boundary
        lfrozen_bb_bot(j-iax+1) = .true.    ! set flag
      case ('top')               ! top boundary
        lfrozen_bb_top(j-iax+1) = .true.    ! set flag
      case default
        print*, "bc_frozen_in_bb: ", topbot, " should be 'top' or 'bot'"
      endselect
!
      lfirstcall=.false.
!
    endsubroutine bc_frozen_in_bb
!***********************************************************************
    subroutine bcz_hydrostatic_temp(f,topbot)
!
!  The logarithmic density in the ghost cells is used to calculate the
!  logarithmic temperature under the assumption of a hydrostatic equilibrium.
!
!  19-nov-2010/Bourdin.KIS: coded
!
      use EquationOfState, only: gamma, gamma_m1, get_cp1
      use SharedVariables, only: get_shared_variable
!
      real, dimension (:,:,:,:), intent (inout) :: f
      character (len=bclen), intent (in) :: topbot
!
      integer :: i
      real, dimension (size(f,1),size(f,2)) :: T_inv, grad_rho
      real :: g_ref, delta_z, inv_cp_cv, cp_inv
      real, dimension (:), pointer :: gravz_zpencil
!
!
      call get_shared_variable ('gravz_zpencil', gravz_zpencil, caller='bcz_hydrostatic_temp')
      call get_cp1 (cp_inv)
      inv_cp_cv = gamma / gamma_m1 * cp_inv
!
      select case (topbot)
      case ('bot')
        ! bottom (left end of the domain)
        do i = 1, nghost
          delta_z = z(n1-i) - z(n1-i+1)
          g_ref = gravz_zpencil(n1-i+1)
          T_inv = exp (-f(:,:,n1-i+1,ilnTT))
          grad_rho = f(:,:,n1-i+1,ilnrho) - f(:,:,n1-i,ilnrho)
          f(:,:,n1-i,ilnTT) = f(:,:,n1-i+1,ilnTT) + grad_rho + g_ref*delta_z*inv_cp_cv*T_inv
        enddo
      case ('top')
        ! top (right end of the domain)
        do i = 1, nghost
          delta_z = z(n2+i) - z(n2+i-1)
          g_ref = gravz_zpencil(n2+i-1)
          T_inv = exp (-f(:,:,n2+i-1,ilnTT))
          grad_rho = f(:,:,n2+i-1,ilnrho) - f(:,:,n2+i,ilnrho)
          f(:,:,n2+i,ilnTT) = f(:,:,n2+i-1,ilnTT) + grad_rho + g_ref*delta_z*inv_cp_cv*T_inv
        enddo
      case default
        call fatal_error ('bcz_hydrostatic_temp', 'invalid argument', lfirst_proc_xy)
      endselect
!
    endsubroutine bcz_hydrostatic_temp
!***********************************************************************
    subroutine bcz_hydrostatic_rho(f,topbot)
!
!  The logarithmic temperature in the ghost cells is used to calculate the
!  logarithmic density under the assumption of a hydrostatic equilibrium.
!
!  24-May-2019/PABourdin: adapted from bcz_hydrostatic_temp
!
      use EquationOfState, only: gamma, gamma_m1, get_cp1
      use SharedVariables, only: get_shared_variable
!
      real, dimension (:,:,:,:), intent (inout) :: f
      character (len=bclen), intent (in) :: topbot
!
      integer :: i
      real, dimension (size(f,1),size(f,2)) :: T_inv
      real :: g_ref, delta_z, inv_cp_cv, cp_inv
      real, dimension (:), pointer :: gravz_zpencil
!
!
      call get_shared_variable ('gravz_zpencil', gravz_zpencil, caller='bcz_hydrostatic_rho')
      call get_cp1 (cp_inv)
      inv_cp_cv = gamma / gamma_m1 * cp_inv
!
      select case (topbot)
      case ('bot')
        ! bottom (left end of the domain)
        do i = 1, nghost
          delta_z = z(n1-i) - z(n1-i+1)
          g_ref = 0.5 * (gravz_zpencil(n1-i) + gravz_zpencil(n1-i+1))
          T_inv = exp (-0.5 * (f(:,:,n1-i,ilnTT) + f(:,:,n1-i+1,ilnTT)))
          f(:,:,n1-i,ilnrho) = f(:,:,n1-i+1,ilnrho) + g_ref*delta_z*inv_cp_cv*T_inv
        enddo
      case ('top')
        ! top (right end of the domain)
        do i = 1, nghost
          delta_z = z(n2+i) - z(n2+i-1)
          g_ref = 0.5 * (gravz_zpencil(n2+i) + gravz_zpencil(n2+i-1))
          T_inv = exp (-0.5 * (f(:,:,n2+i,ilnTT) + f(:,:,n2+i-1,ilnTT)))
          f(:,:,n2+i,ilnrho) = f(:,:,n2+i-1,ilnrho) + g_ref*delta_z*inv_cp_cv*T_inv
        enddo
      case default
        call fatal_error ('bcz_hydrostatic_rho', 'invalid argument', lfirst_proc_xy)
      endselect
!
    endsubroutine bcz_hydrostatic_rho
!***********************************************************************
    subroutine finalize_boundcond(f)
!
!  Call finalization routines, i.e. freeing allocated memory.
!
! 15-aug-2011/Bourdin.KIS: adapted from finalize_modules
!
      use General, only: keep_compiler_quiet
!
      real, dimension(:,:,:,:) :: f
!
      call bc_aa_pot_field_extrapol(f,"all",.true.)
!
    endsubroutine finalize_boundcond
!***********************************************************************
    subroutine bc_aa_pot_field_extrapol(f,topbot,lfinalize)
!
!  Potential field extrapolation in z-direction for the ghost cells.
!  To have a smooth transition at the boundary from non-force-free to a
!  force-free field, the Az component is also extrapolated in the same way.
!  This reduces strong currents at the top boundary.
!  At the bottom boundary the extrapolation increases contrasts in A.
!  This mimicks flux tubes that become narrower below the photosphere.
!
!  9-jul-2010/Bourdin.KIS: coded
!
      use Fourier, only: vect_pot_extrapol_z_parallel, kx_fft, ky_fft
!
      real, dimension (:,:,:,:), intent (inout) :: f
      character (len=bclen), intent (in) :: topbot
      logical, optional :: lfinalize
!
      real, dimension (:,:,:), allocatable, save :: exp_fact_top, exp_fact_bot
      integer, parameter :: bnx=nygrid, bny=nx/nprocy
      integer :: kx_start, stat, pos_z
      real :: delta_z, reduce_factor=1.
!
      if (ldownsampling) then
        call warning('bc_force_aa_time','Not available for downsampling')       
        return
      endif

      if (present (lfinalize)) then
        if (lfinalize) then
          if (allocated (exp_fact_bot)) deallocate (exp_fact_bot)
          if (allocated (exp_fact_top)) deallocate (exp_fact_top)
          return
        endif
      endif
!
      ! reduce_factor reduces the structure increase at the bottom boundary
      ! to help numerically resolving the strong gradients in the ghost cells.
      ! Set reduce_factor to [0,1] by using fbcz_bot(iaa) in run_pars.
      ! A value of 1 (default) switches this reducing-mechanism off.
      ! A value of 0 just clones the vector field A from the n1-layer.
      if (fbcz_bot(iaa) /= 0.) reduce_factor = fbcz_bot(iaa)
!
      if (.not. ((lfirst_proc_z .and. (topbot == 'bot')) .or. (llast_proc_z .and. (topbot == 'top')))) &
          call fatal_error ('bc_aa_pot_field_extrapol', 'Only implemented for topmost or downmost z-layer.', lfirst_proc_xy)
!
      if (mod (nx, nprocy) /= 0) &
          call fatal_error ('bc_aa_pot_field_extrapol', 'nx needs to be an integer multiple of nprocy.', lfirst_proc_xy)
!
!  Check whether we want to do top or bottom z boundary
!
      select case (topbot)
      case ('bot')
        if (.not. allocated (exp_fact_bot)) then
          ! Setup exponential factor for bottom boundary
          allocate (exp_fact_bot(bnx,bny,nghost), stat=stat)
          if (stat > 0) call fatal_error ('bc_aa_pot_field_extrapol', 'Could not allocate memory for exp_fact_bot', .true.)
          ! Get wave numbers already in transposed pencil shape and calculate exp(|k|)
          kx_start = (ipx+ipy*nprocx)*bny
          exp_fact_bot = spread (exp (sqrt (spread (ky_fft(1:bnx), 2, bny) ** 2 + &
                                            spread (kx_fft(kx_start+1:kx_start+bny), 1, bnx) ** 2)), 3, nghost)
          do pos_z = 1, nghost
            ! dz is positive => enhance structures or contrast
            delta_z = reduce_factor * (z(n1) - z(n1-nghost+pos_z-1))
            ! Include normalization factor for fourier transform: 1/(nxgrid*nygrid)
            exp_fact_bot(:,:,pos_z) = exp_fact_bot(:,:,pos_z) ** delta_z / (nxgrid*nygrid)
          enddo
        endif
        call vect_pot_extrapol_z_parallel &
             (f(l1:l2,m1:m2,n1,iax:iaz), f(l1:l2,m1:m2,n1-nghost:n1-1,iax:iaz), exp_fact_bot)
      case ('top')
        if (.not. allocated (exp_fact_top)) then
          ! Setup exponential factor for top boundary
          allocate (exp_fact_top(bnx,bny,nghost), stat=stat)
          if (stat > 0) call fatal_error ('bc_aa_pot_field_extrapol', 'Could not allocate memory for exp_fact_top', .true.)
          ! Get wave numbers already in transposed pencil shape and calculate exp(|k|)
          kx_start = (ipx+ipy*nprocx)*bny
          exp_fact_top = spread (exp (sqrt (spread (ky_fft(1:bnx), 2, bny) ** 2 + &
                                            spread (kx_fft(kx_start+1:kx_start+bny), 1, bnx) ** 2)), 3, nghost)
          do pos_z = 1, nghost
            ! dz is negative => decay of structures
            delta_z = z(n2) - z(n2+pos_z)
            ! Include normalization factor for fourier transform: 1/(nxgrid*nygrid)
            exp_fact_top(:,:,pos_z) = exp_fact_top(:,:,pos_z) ** delta_z / (nxgrid*nygrid)
          enddo
        endif
        call vect_pot_extrapol_z_parallel &
             (f(l1:l2,m1:m2,n2,iax:iaz), f(l1:l2,m1:m2,n2+1:n2+nghost,iax:iaz), exp_fact_top)
      case default
        call fatal_error ('bc_aa_pot_field_extrapol', 'invalid argument', lfirst_proc_xy)
      endselect
!
!  The vector potential needs to be known outside of (l1:l2,m1:m2) as well
!
      call communicate_vect_field_ghosts(f,topbot)
!
    endsubroutine bc_aa_pot_field_extrapol
!***********************************************************************
    subroutine bc_aa_pot3(f,topbot)
!
!  Potential field boundary condition
!
!  11-oct-06/wolf: Adapted from Tobi's bc_aa_pot2
!
      use Fourier, only: fourier_transform_xy_xy, kx_fft, ky_fft
!
      real, dimension (:,:,:,:), intent (inout) :: f
      character (len=bclen), intent (in) :: topbot
!
      real, dimension (l2-l1+1,m2-m1+1,iax:iaz) :: aa_re,aa_im
      real, dimension (l2-l1+1,m2-m1+1) :: kx,ky,kappa,exp_fact,tmp_re,tmp_im
      real    :: delta_z
      integer :: i,j,nxl,nyl
!
!  Get local wave numbers
!
      kx = spread(kx_fft(ipx*nxl+1:ipx*nxl+nxl),2,nyl)
      ky = spread(ky_fft(ipy*nyl+1:ipy*nyl+nyl),1,nxl)
!
!  Calculate 1/k^2, zero mean
!
      kappa = sqrt(kx**2 + ky**2)
!
!  Fourier transforms of x- and y-components on the boundary
!  Check whether we want to do top or bottom (this is processor dependent)
!
      select case (topbot)
!
      case ('bot')
        ! Potential field condition at the bottom
        do j=1,nghost
!
! Calculate delta_z based on z(), not on dz to improve behavior for
! non-equidistant grid (still not really correct, but could be OK)
!
          delta_z  = z(n1+j) - z(n1-j)
          exp_fact = exp(-kappa*delta_z)
          ! Determine potential field in ghost zones
          do i=iax,iaz
            tmp_re = f(l1:l2,m1:m2,n1+j,i)
            tmp_im = 0.0
            call fourier_transform_xy_xy(tmp_re,tmp_im)
            aa_re(:,:,i) = tmp_re*exp_fact
            aa_im(:,:,i) = tmp_im*exp_fact
          enddo
         ! Transform back
          do i=iax,iaz
            tmp_re = aa_re(:,:,i)
            tmp_im = aa_im(:,:,i)
            call fourier_transform_xy_xy(tmp_re,tmp_im,linv=.true.)
            f(l1:l2,m1:m2,n1-j,i) = tmp_re
          enddo
        enddo
!
      case ('top')
        ! Potential field condition at the top
        do j=1,nghost
!
! Calculate delta_z based on z(), not on dz to improve behavior for
! non-equidistant grid (still not really correct, but could be OK)
!
          delta_z  = z(n2+j) - z(n2-j)
          exp_fact = exp(-kappa*delta_z)
          ! Determine potential field in ghost zones
          do i=iax,iaz
            tmp_re = f(l1:l2,m1:m2,n2-j,i)
            tmp_im = 0.0
            call fourier_transform_xy_xy(tmp_re,tmp_im)
            aa_re(:,:,i) = tmp_re*exp_fact
            aa_im(:,:,i) = tmp_im*exp_fact
          enddo
          ! Transform back
          do i=iax,iaz
            tmp_re = aa_re(:,:,i)
            tmp_im = aa_im(:,:,i)
            call fourier_transform_xy_xy(tmp_re,tmp_im,linv=.true.)
            f(l1:l2,m1:m2,n2+j,i) = tmp_re
          enddo
        enddo
!
      case default
        call fatal_error('bc_aa_pot3', 'invalid argument', lfirst_proc_xy)
!
      endselect
!
!  The vector potential needs to be known outside of (l1:l2,m1:m2) as well
!
      call communicate_vect_field_ghosts(f,topbot)
!
    endsubroutine bc_aa_pot3
!***********************************************************************
    subroutine bc_aa_pot2(f,topbot)
!
!  Potential field boundary condition
!
!  10-oct-06/tobi: Coded
!
      use Fourier, only: fourier_transform_xy_xy, fourier_transform_y_y, kx_fft, ky_fft
!
      real, dimension (:,:,:,:), intent (inout) :: f
      character (len=bclen), intent (in) :: topbot
!
      real, dimension (l2-l1+1,m2-m1+1,iax:iaz) :: aa_re,aa_im
      real, dimension (l2-l1+1,m2-m1+1) :: kx,ky,kappa
      real, dimension (l2-l1+1,m2-m1+1) :: tmp_re,tmp_im,fac
      integer :: i,j,nxl,nyl
!
      nxl=l2-l1+1; nyl=m2-m1+1
!
!  Get local wave numbers
!
      if (nxgrid>1) then
        kx = spread(kx_fft(ipx*nxl+1:ipx*nxl+nxl),2,nyl)
        ky = spread(ky_fft(ipy*nyl+1:ipy*nyl+nyl),1,nxl)
      else
        kx(1,:) = 0.0
        ky(1,:) = ky_fft(ipy*nyl+1:ipy*nyl+nyl)
      endif
!
!  Calculate 1/k^2, zero mean
!
      kappa = sqrt(kx**2 + ky**2)
!
!  Fourier transforms of x- and y-components on the boundary
!  Check whether we want to do top or bottom (this is processor dependent)
!
      select case (topbot)
!
      case ('bot')
        ! Potential field condition at the bottom
        do i=iax,iaz
          tmp_re = f(l1:l2,m1:m2,n1,i)
          tmp_im = 0.0
          if (nxgrid>1) then
            call fourier_transform_xy_xy(tmp_re,tmp_im)
          else
            call fourier_transform_y_y(tmp_re,tmp_im)
          endif
          aa_re(:,:,i) = tmp_re
          aa_im(:,:,i) = tmp_im
        enddo
        ! Determine potential field in ghost zones
        do j=1,nghost
          fac = exp(-j*kappa*dz)
          do i=iax,iaz
            tmp_re = fac*aa_re(:,:,i)
            tmp_im = fac*aa_im(:,:,i)
            if (nxgrid>1) then
              call fourier_transform_xy_xy(tmp_re,tmp_im,linv=.true.)
            else
              call fourier_transform_y_y(tmp_re,tmp_im,linv=.true.)
            endif
            f(l1:l2,m1:m2,n1-j,i) = tmp_re
          enddo
        enddo
!
      case ('top')
        ! Potential field condition at the top
        do i=iax,iaz
          tmp_re = f(l1:l2,m1:m2,n2,i)
          tmp_im = 0.0
          if (nxgrid>1) then
            call fourier_transform_xy_xy(tmp_re,tmp_im)
          else
            call fourier_transform_y_y(tmp_re,tmp_im)
          endif
          aa_re(:,:,i) = tmp_re
          aa_im(:,:,i) = tmp_im
        enddo
        ! Determine potential field in ghost zones
        do j=1,nghost
          fac = exp(-j*kappa*dz)
          do i=iax,iaz
            tmp_re = fac*aa_re(:,:,i)
            tmp_im = fac*aa_im(:,:,i)
            if (nxgrid>1) then
              call fourier_transform_xy_xy(tmp_re,tmp_im,linv=.true.)
            else
              call fourier_transform_y_y(tmp_re,tmp_im,linv=.true.)
            endif
            f(l1:l2,m1:m2,n2+j,i) = tmp_re
          enddo
        enddo
!
      case default
        call fatal_error('bc_aa_pot2', 'invalid argument', lfirst_proc_xy)
!
      endselect
!
!  The vector potential needs to be known outside of (l1:l2,m1:m2) as well
!
        call communicate_vect_field_ghosts(f,topbot)
!
    endsubroutine bc_aa_pot2
!***********************************************************************
      subroutine bc_aa_pot(f,topbot)
!
!  Potential field boundary condition for magnetic vector potential at
!  bottom or top boundary (in z).
!
!  14-jun-2002/axel: adapted from similar
!   8-jul-2002/axel: introduced topbot argument
!
      real, dimension (:,:,:,:) :: f
      character (len=bclen) :: topbot
!
      real, dimension (l2-l1+1,m2-m1+1) :: f2,f3
      real, dimension (l2-l1+1,m2-m1+1,nghost+1) :: fz
      integer :: j
!
!  potential field condition
!  check whether we want to do top or bottom (this is processor dependent)
!
      select case (topbot)
!
!  potential field condition at the bottom
!
      case ('bot')
        if (headtt) print*,'bc_aa_pot: pot-field bdry cond at bottom'
        if (mod(nxgrid,nygrid)/=0) &
             call fatal_error("bc_aa_pot", "pot-field doesn't work "//&
                 "with mod(nxgrid,nygrid)/=0", lfirst_proc_xy)
        do j=0,1
          f2=f(l1:l2,m1:m2,n1+1,iax+j)
          f3=f(l1:l2,m1:m2,n1+2,iax+j)
          call potential_field(fz,f2,f3,-1)
          f(l1:l2,m1:m2,1:n1,iax+j)=fz
        enddo
!
        f2=f(l1:l2,m1:m2,n1,iax)
        f3=f(l1:l2,m1:m2,n1,iay)
        call potentdiv(fz,f2,f3,-1)
        f(l1:l2,m1:m2,1:n1,iaz)=-fz
!
!  potential field condition at the top
!
      case ('top')
        if (headtt) print*,'bc_aa_pot: pot-field bdry cond at top'
        if (mod(nxgrid,nygrid)/=0) &
             call fatal_error("bc_aa_pot", "pot-field doesn't work "//&
                 "with mod(nxgrid,nygrid)/=0", lfirst_proc_xy)
        do j=0,1
          f2=f(l1:l2,m1:m2,n2-1,iax+j)
          f3=f(l1:l2,m1:m2,n2-2,iax+j)
          call potential_field(fz,f2,f3,+1)
          f(l1:l2,m1:m2,n2:,iax+j)=fz
        enddo
!
        f2=f(l1:l2,m1:m2,n2,iax)
        f3=f(l1:l2,m1:m2,n2,iay)
        call potentdiv(fz,f2,f3,+1)
        f(l1:l2,m1:m2,n2:,iaz)=-fz
      case default
        call fatal_error('bc_aa_pot', 'invalid argument', lfirst_proc_xy)
      endselect
!
      call communicate_vect_field_ghosts(f,topbot)
!
      endsubroutine bc_aa_pot
!***********************************************************************
      subroutine potential_field(fz,f2,f3,irev)
!
!  solves the potential field boundary condition;
!  fz is the boundary layer, and f2 and f3 are the next layers inwards.
!  The condition is the same on the two sides.
!
!  20-jan-00/axel+wolf: coded
!  22-mar-00/axel: corrected sign (it is the same on both sides)
!  29-sep-06/axel: removed multiple calls, removed normalization, non-para
!
      use Fourier, only: fourier_transform_xy_xy, kx_fft, ky_fft
!
      real, dimension (:,:,:) :: fz
      real, dimension (:,:) :: f2,f3
      integer :: irev
!
      real, dimension(l2-l1+1,m2-m1+1) :: fac,kk,f1r,f1i,g1r,g1i,f2r,f2i,f3r,f3i
      real :: delz
      integer :: i,nxl,nyl
!
      nxl=l2-l1+1; nyl=m2-m1+1
!
!  initialize workspace
!
      f2r=f2; f2i=0
      f3r=f3; f3i=0
!
!  Transform; real and imaginary parts
!
      call fourier_transform_xy_xy(f2r,f2i)
      call fourier_transform_xy_xy(f3r,f3i)
!
!  define wave vector
!  calculate sqrt(k^2)
!
      kk=sqrt(spread(kx_fft(ipx*nxl+1:ipx*nxl+nxl)**2,2,nyl)+spread(ky_fft(ipy*nyl+1:ipy*nyl+nyl)**2,1,nxl))
!
!  one-sided derivative
!
      fac=1./(3.+2.*dz*kk)
      f1r=fac*(4.*f2r-f3r)
      f1i=fac*(4.*f2i-f3i)
!
!  set ghost zones
!
      do i=0,nghost
        delz=i*dz
        fac=exp(-kk*delz)
        g1r=fac*f1r
        g1i=fac*f1i
!
!  Transform back
!
        call fourier_transform_xy_xy(g1r,g1i,linv=.true.)
!
!  reverse order if irev=-1 (if we are at the bottom)
!
        if (irev==+1) fz(:,:,       i+1) = g1r
        if (irev==-1) fz(:,:,nghost-i+1) = g1r
      enddo
!
    endsubroutine potential_field
!***********************************************************************
    subroutine potentdiv(fz,f2,f3,irev)
!
!  solves the divA=0 for potential field boundary condition;
!  f2 and f3 correspond to Ax and Ay (input) and fz corresponds to Ax (out)
!  In principle we could save some ffts, by combining with the potential
!  subroutine above, but this is now easier
!
!  22-mar-02/axel: coded
!  29-sep-06/axel: removed multiple calls, removed normalization, non-para
!   7-oct-06/axel: corrected sign for irev==+1.
!
      use Fourier, only: fourier_transform_xy_xy
!
      real, dimension (:,:,:) :: fz
      real, dimension (:,:) :: f2,f3
      integer :: irev
!
      real, dimension (l2-l1+1,m2-m1+1) :: fac,kk,kkkx,kkky, &
                                           f1r,f1i,g1r,g1i,f2r,f2i,f3r,f3i
      real, dimension (nygrid) :: ky
      real, dimension (nx) :: kx
      real :: delz
      integer :: i, nxl, nyl
!
      if (ldownsampling) then
        call warning('bc_force_aa_time','Not available for downsampling')      
        return
      endif
!
      nxl=l2-l1+1; nyl=m2-m1+1
      f2r=f2; f2i=0
      f3r=f3; f3i=0
!
!  Transform
!
      call fourier_transform_xy_xy(f2r,f2i)
      call fourier_transform_xy_xy(f3r,f3i)
!
!  define wave vector
!
      kx=cshift((/(i-nxl/2,i=0,nxl-1)/),+nxl/2)*2*pi/Lx
      ky=cshift((/(i-nygrid/2,i=0,nygrid-1)/),+nygrid/2)*2*pi/Ly
!
!  calculate 1/k^2, zero mean
!
      kk=sqrt(spread(kx**2,2,nyl)+spread(ky(ipy*nyl+1:(ipy+1)*nyl)**2,1,nxl))
      kkkx=spread(kx,2,nyl)
      kkky=spread(ky(ipy*nyl+1:(ipy+1)*nyl),1,nxl)
!
!  calculate 1/kk
!
      kk(1,1)=1.
      fac=1./kk
      fac(1,1)=0.
!
      f1r=fac*(-kkkx*f2i-kkky*f3i)
      f1i=fac*(+kkkx*f2r+kkky*f3r)
!
!  set ghost zones
!
      do i=0,nghost
        delz=i*dz
        fac=exp(-kk*delz)
        g1r=fac*f1r
        g1i=fac*f1i
!
!  Transform back
!
        call fourier_transform_xy_xy(g1r,g1i,linv=.true.)
!
!  reverse order if irev=-1 (if we are at the bottom)
!  but reverse sign if irev=+1 (if we are at the top)
!
        if (irev==+1) fz(:,:,       i+1) = -g1r
        if (irev==-1) fz(:,:,nghost-i+1) = +g1r
      enddo
!
    endsubroutine potentdiv
!***********************************************************************
    subroutine bc_wind_z(f,topbot,massflux)
!
!  Calculates u_0 so that rho*(u+u_0)=massflux.
!  Set 'win' for rho and
!  massflux can be set as fbcz1/2(rho) in run.in.
!
!  18-06-2008/bing: coded
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: i,j,ipt,ntb=-1
      real :: massflux,u_add
      real :: local_flux,local_mass
      real :: total_flux,total_mass
      real :: get_lf,get_lm
      integer :: nroot
!
      if (ldownsampling) then
        call warning('bc_force_aa_time','Not available for downsampling')    
        return
      endif
!
      if (headtt) then
        print*,'bc_wind: Massflux',massflux
!
!   check wether routine can be implied
!
        if (.not.(lequidist(1) .and. lequidist(2))) &
            call fatal_error('bc_wind_z', &
            'non equidistant grid in x and y not implemented')
!
!   check for warnings
!
        if (.not. ldensity)  &
            call warning('bc_wind',"no defined density, using rho=1 ?")
      endif
!
      select case (topbot)
!
!  Bottom boundary.
!
      case ('bot')
        ntb = n1
        nroot = 0
!
!  Top boundary.
!
      case ('top')
        ntb = n2
        nroot = ipz*nprocx*nprocy
!
!  Default.
!
      case default
        print*, "bc_wind: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
      local_flux=sum(exp(f(l1:l2,m1:m2,ntb,ilnrho))*f(l1:l2,m1:m2,ntb,iuz))
      local_mass=sum(exp(f(l1:l2,m1:m2,ntb,ilnrho)))
!
!  One  processor has to collect the data
!
      if (iproc/=nroot) then
        ! send to first processor at given height
        !
        call mpisend_real(local_flux,nroot,111+iproc)
        call mpisend_real(local_mass,nroot,211+iproc)
      else
        total_flux=local_flux
        total_mass=local_mass
        do i=0,nprocx-1
          do j=0,nprocy-1
            ipt = i+nprocx*j+ipz*nprocx*nprocy
            if (ipt/=nroot) then
              call mpirecv_real(get_lf,ipt,111+ipt)
              call mpirecv_real(get_lm,ipt,211+ipt)
              total_flux=total_flux+get_lf
              total_mass=total_mass+get_lm
            endif
          enddo
        enddo
!
!  Get u0 addition rho*(u+u0) = wind
!  rho*u + u0 *rho =wind
!  u0 = (wind-rho*u)/rho
!
        u_add = (massflux-total_flux) / total_mass
      endif
!
!  now distribute u_add
!
      if (iproc/=nroot) then
        call mpirecv_real(u_add,nroot,311+iproc)
      else
        do i=0,nprocx-1
          do j=0,nprocy-1
            ipt = i+nprocx*j+ipz*nprocx*nprocy
            if (ipt/=nroot) then
              call mpisend_real(u_add,ipt,311+ipt)
            endif
          enddo
        enddo
      endif
!
!  Set boundary
!
      f(l1:l2,m1:m2,ntb,iuz) =  f(l1:l2,m1:m2,ntb,iuz)+u_add
!
    endsubroutine bc_wind_z
!***********************************************************************
    subroutine bc_ADI_flux_z(f,topbot)
!
!  Constant flux boundary condition for temperature (called when bcz='c3')
!  at the bottom _only_ in the ADI case where hcond(n1)=hcond(x)
!  TT version: enforce dT/dz = - Fbot/K
!  30-jan-2009/dintrans: coded
!
      use SharedVariables, only: get_shared_variable
!
      real, pointer :: Fbot
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      real, dimension (size(f,1)) :: tmp_x
      integer :: i
!
      call get_shared_variable('Fbot', Fbot, caller='bc_ADI_flux_z')
!
      if (headtt) print*,'bc_ADI_flux_z: Fbot, hcondADI, dz=', &
           Fbot, hcondADI, dz
!
      if (topbot=='bot') then
        tmp_x=-Fbot/hcondADI
        do i=1,nghost
          f(:,4,n1-i,ilnTT)=f(:,4,n1+i,ilnTT)-dz2_bound(-i)*tmp_x
        enddo
      else
        call fatal_error('bc_ADI_flux_z', 'invalid argument')
      endif
!
    endsubroutine bc_ADI_flux_z
!***********************************************************************
    subroutine bc_force_ux_time(f, idz, j)
!
!  Set ux = ampl_forc*sin(k_forc*x)*cos(w_forc*t)
!
!  05-jun-2009/dintrans: coded from bc_force_uxy_sin_cos
!  Note: the ampl_forc, k_forc & w_forc run parameters are set in
!  'hydro' and shared using the 'shared_variables' module
!
      use SharedVariables, only : get_shared_variable
!
      real, dimension (:,:,:,:) :: f
      integer :: idz, j
      real    :: kx
      real, pointer, save :: ampl_forc, k_forc, w_forc, x_forc, dx_forc
      logical, save :: l1st=.true.
!
      if (headtt) then
        if (iuz == 0) call stop_it("BC_FORCE_UX_TIME: Bad idea...")
        if (Lx  == 0) call stop_it("BC_FORCE_UX_TIME: Lx cannot be 0")
        if (j /= iux) call stop_it("BC_FORCE_UX_TIME: only valid for ux")
      endif
!
      if (l1st) then
        call get_shared_variable('ampl_forc', ampl_forc, caller='bc_force_ux_time')
        call get_shared_variable('k_forc', k_forc)
        call get_shared_variable('w_forc', w_forc)
        call get_shared_variable('x_forc', x_forc)
        call get_shared_variable('dx_forc', dx_forc)
        if (headtt) print*, 'bc_force_ux_time: ampl_forc, k_forc, '//&
             'w_forc, x_forc, dx_forc=', ampl_forc, k_forc, w_forc, &
             x_forc, dx_forc
        l1st=.false.
      endif
!
      if (k_forc /= impossible) then
        kx=2*pi/Lx*k_forc
        f(:,:,idz,j) = spread(ampl_forc*sin(kx*x)*cos(w_forc*t), 2, size(f,2))
      else
        f(:,:,idz,j) = spread(ampl_forc*exp(-((x-x_forc)/dx_forc)**2)*cos(w_forc*t), 2, size(f,2))
      endif
!
    endsubroutine bc_force_ux_time
!***********************************************************************
    subroutine bc_inlet_outlet_cyl(f,topbot,j,val)
!
! For pi/2 < y < 3pi/4,
! set r and theta velocity corresponding to a constant x-velocity
! and symmetric for lnrho/rho.
!
! Otherwise, set symmetric for velocities, and constant
! for lnrho/rho.
!
! NB! Assumes y to have the range 0 < y < 2pi
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j,i
      real, dimension(:) :: val
!
      select case (topbot)
      case ('bot')
        call fatal_error('bc_inlet_outlet_cyl', &
          'this boundary condition is not allowed for bottom boundary')
      case ('top')
        do m=m1,m2
          if (      (y(m)>=xyz0(2) +   Lxyz(2)/4)&
              .and. (y(m)<=xyz0(2) + 3*Lxyz(2)/4)) then
            if (j==iux) then
              f(l2,m,:,j) = cos(y(m))*val(j)
              do i=1,nghost; f(l2+i,m,:,j) = 2*f(l2,m,:,j) - f(l2-i,m,:,j); enddo
            elseif (j==iuy) then
              f(l2,m,:,j) = -sin(y(m))*val(j)
              do i=1,nghost; f(l2+i,m,:,j) = 2*f(l2,m,:,j) - f(l2-i,m,:,j); enddo
            elseif ((j==ilnrho) .or. (j==irho)) then
              do i=1,nghost; f(l2+i,m,:,j) = f(l2-i,m,:,j); enddo
            endif
!
          else
            if (j==iux) then
              do i=1,nghost; f(l2+i,m,:,j) = f(l2-i,m,:,j); enddo
            elseif (j==iuy) then
              do i=1,nghost; f(l2+i,m,:,j) = f(l2-i,m,:,j); enddo
            elseif ((j==ilnrho) .or. (j==irho)) then
              f(l2,m,:,j) = val(j)
              do i=1,nghost; f(l2+i,m,:,j) = 2*f(l2,m,:,j) - f(l2-i,m,:,j); enddo
            endif
          endif
        enddo
      endselect
!
    endsubroutine bc_inlet_outlet_cyl
!***********************************************************************
    subroutine bc_pp_hds_z_iso(f,topbot)
!
!  Boundary condition for pressure
!
!  This sets \partial_{z} p = \rho g_{z},
!  i.e. it enforces hydrostatic equlibrium at the boundary for the
!  pressure with an isothermal EOS.
!
!  16-dec-2009/dintrans: coded
!
      use Gravity, only: gravz
      use EquationOfState, only : cs20
!
      real, dimension (:,:,:,:), intent (inout) :: f
      character (len=bclen), intent (in) :: topbot
      real    :: haut
      integer :: i
!
      haut=cs20/gravz
      if (topbot=='bot') then
        do i=1,nghost
          f(:,:,n1-i,ipp) = f(:,:,n1+i,ipp)-dz2_bound(-i)*f(:,:,n1,ipp)/haut
        enddo
      else
        do i=1,nghost
          f(:,:,n2+i,ipp) = f(:,:,n2-i,ipp)+dz2_bound(i)*f(:,:,n2,ipp)/haut
        enddo
      endif
!
    endsubroutine bc_pp_hds_z_iso
!***********************************************************************
    subroutine bc_symset0der_z_v2(f,topbot,j)
!
!  This routine modified from bc_sym_z, but to a lower order.
!  Only available for z axis, activate with "0ds"
!  This is the routine to be used as regularity condition on the axis.
!
!  25-Oct-10/tijmen & bing: coded
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: i,j
!
      select case (topbot)
!
!  bottom (left end of the domain)
      case ('bot')
        f(:,:,n1,j)=(-18.*f(:,:,n1+1,j) &
                     +9.*f(:,:,n1+2,j) &
                     -2.*f(:,:,n1+3,j))/11.
!
        do i=1,nghost; f(:,:,n1-i,j)=f(:,:,n1+i,j); enddo
!
!  top (right end of the domain)
      case ('top')
        f(:,:,n2,j)=(+18.*f(:,:,n2-1,j) &
                     -9.*f(:,:,n2-2,j) &
                     +2.*f(:,:,n2-3,j))/11.
!
        do i=1,nghost; f(:,:,n2+i,j)=f(:,:,n2-i,j); enddo
!
      case default
        print*, "bc_symset0der_z_v2: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_symset0der_z_v2
!***********************************************************************
    subroutine bc_aa_pot_1D(f,topbot)
!
!  Computes a potential field extrapolation for a
!  1D magnetic field boundary with nprocx >= 1
!
!  27-Oct-10/bing: coded
!
      use Fourier, only: fourier_transform_other, kx_fft
!
      real, dimension (:,:,:,:), intent (inout) :: f
      character (len=bclen), intent (in) :: topbot
      real, dimension (nxgrid) :: fft_az_r,fft_az_i,A_r,A_i,exp_fact
      real, dimension (nxgrid) :: iay_global
      integer :: i,j,ipos,dir
!
      if (ldownsampling) then
        call warning('bc_force_aa_time','Not available downsampling')
        return
      endif
!
      select case (topbot)
!
!  bottom (left end of the domain)
      case ('bot')
        ipos = n1
        dir = -1
!
!  top (right end of the domain)
      case ('top')
        ipos = n2
        dir = 1
!
      case default
        print*, "bc_aa_pot_1D: ", topbot, " should be 'top' or 'bot'"
        ipos=1
        dir=0
!
      endselect
!
      if (nygrid>1) call fatal_error('bc_aa_pot_1D','only for nygrid=1')
!
      if (iproc==0) then
        iay_global(1:nx) = f(l1:l2,m1,ipos,iay)
        if (nprocx>1) then
          do j=1,nprocx-1
            call mpirecv_real(iay_global(j*nx+1:(j+1)*nx),nx,j,j*100)
          enddo
        endif
        fft_az_r=iay_global
        call fourier_transform_other(fft_az_r,fft_az_i)
      else
        if (nprocx>1) call mpisend_real(f(l1:l2,m1,ipos,iay),nx,0,iproc*100)
      endif
!
      do i=1,nghost
        if (iproc==0) then
!
          exp_fact = exp(-abs(kx_fft)*(z(ipos+dir*i)-z(ipos)))
!
          A_r = exp_fact*fft_az_r
          A_i = exp_fact*fft_az_i
!
          call fourier_transform_other(A_r,A_i,linv=.true.)
!
          f(l1:l2,m1,ipos+dir*i,iay) = A_r(1:nx)
!
          if (nprocx>1) then
            do j=1,nprocx-1
              call mpisend_real(A_r(j*nx+1:(j+1)*nx),nx,j,j*100)
            enddo
          endif
        else
          if (nprocx>1) call mpirecv_real(f(l1:l2,m1,ipos+dir*i,iay),nx,0,iproc*100)
        endif
!
      enddo
!
    endsubroutine bc_aa_pot_1D
!***********************************************************************
    subroutine bc_ctz(f,topbot,j)
!
!  Set entropy to match temperature in the ghost zones to boundary value
!  value with small increment. Density ghost zones need to be calculated
!  again here and corners must be included to avoid NAN's.
!
!  14-mar-11/fred: check that 'cdz' is also set for bcz density.
!
      use EquationOfState, only: get_cv1,get_cp1
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j,k
      real :: cv1,cp1,cv,cp
      real, dimension (size(f,1),size(f,2),size(f,3)) :: lnrho_
!
      call get_cv1(cv1); cv=1./cv1
      call get_cp1(cp1); cp=1./cp1
!
      call bc_cdz(f,topbot,j-1)
!
      lnrho_=f(:,:,:,j-1)
      if (ldensity_nolog) then
        where (lnrho_<=0) lnrho_=tini
        lnrho_=log(lnrho_)
      endif
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        do k=1,3
          f(:,:,n1-k,j)=f(:,:,n1-k+1,j)+(cp-cv)*&
              (lnrho_(:,:,n1-k+1)-lnrho_(:,:,n1-k))
        enddo
!
      case ('top')               ! top boundary
        do k=1,3
          f(:,:,n2+k,j)=f(:,:,n2+k-1,j)+(cp-cv)*&
              (lnrho_(:,:,n2+k-1)-lnrho_(:,:,n2+k))
        enddo
!
      case default
        print*, "bc_ctz ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_ctz
!***********************************************************************
    subroutine bc_cdz(f,topbot,j)
!
!  Set ghost values to diminishing amplitude of boundary value.
!  Motivation density spikes in 'ism' runs leading to temp spikes that
!  crash the code on outflows, but halo density much lower so 'cop'
!  induces mass inflows which are too high to be physically sustainable.
!
!  13-feb-11/fred: adapted from bc_ctz
!  16-apr-14/fred: revised constant from 10. to 1.11 to preserve exponential
!                  reduction in density with height ref Ferriere Review 2001
!                  Eq.(5)
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j,k
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
          do k=1,3
            f(:,:,n1-k,j)=f(:,:,n1-k+1,j)*(1.0-1.11*dz)
          enddo
!
      case ('top')               ! top boundary
          do k=1,3
            f(:,:,n2+k,j)=f(:,:,n2+k-1,j)*(1.0-1.11*dz)
          enddo
!
      case default
        print*, "bc_cdz ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_cdz
!***********************************************************************
    subroutine bc_ism(f,topbot,j)
!
!  30-nov-15/fred: Replaced bc_ctz and bc_cdz.
!  Apply observed scale height locally from Reynolds 1991, Manchester & Taylor
!  1981 for warm ionized gas - dominant scale height above 500 parsecs.
!  Apply constant local temperature across boundary for entropy.
!  Motivation to prevent numerical spikes in shock fronts, which cannot be 
!  absorbed in only three ghost cells, but boundary thermodynamics still 
!  responsive to interior dynamics.
!
      use EquationOfState, only: get_cv1,get_cp1
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j,k
      real, parameter :: density_scale_cgs=2.7774e21 !900pc Reynolds 91, etc
      real :: density_scale1, density_scale
      real :: cv1,cp1,cv,cp
!
      density_scale1=unit_length/density_scale_cgs
      density_scale=1./density_scale1
      call get_cv1(cv1); cv=1./cv1
      call get_cp1(cp1); cp=1./cp1
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        do k=1,nghost
          if (j==irho .or. j==ilnrho) then
            if (ldensity_nolog) then
              f(:,:,k,j)=f(:,:,n1,j)*exp(-(z(n1)-z(k))*density_scale1)
            else
              f(:,:,k,j)=f(:,:,n1,j) - (z(n1)-z(k))*density_scale
            endif
          else if (j==iss) then
            if (ldensity_nolog) then
              f(:,:,n1-k,j)=f(:,:,n1,j)+(cp-cv)*&
                  (log(f(:,:,n1,j-1))-log(f(:,:,n1-k,j-1)))+&
                  cv*log((z(n1)-z(n1-k))*density_scale+1.)
            else
              f(:,:,n1-k,j)=f(:,:,n1,j)+(cp-cv)*&
                  (f(:,:,n1,j-1)-f(:,:,n1-k,j-1))+&
                  cv*log((z(n1)-z(n1-k))*density_scale+1.)
            endif
          else
            call fatal_error('bc_ism','only for irho, ilnrho, iuz or iss')
          endif
        enddo
!
      case ('top')               ! top boundary
        do k=1,nghost
          if (j==irho .or. j==ilnrho) then
            if (ldensity_nolog) then
              f(:,:,n2+k,j)=f(:,:,n2,j)*exp(-(z(n2+k)-z(n2))*density_scale1)
            else
              f(:,:,n2+k,j)=f(:,:,n2,j) - (z(n2+k)-z(n2))*density_scale1
            endif
          else if (j==iss) then
            if (ldensity_nolog) then
              f(:,:,n2+k,j)=f(:,:,n2,j)+(cp-cv)*&
                  (log(f(:,:,n2,j-1))-log(f(:,:,n2+k,j-1)))+&
                  cv*log((z(n2+k)-z(n2))*density_scale+1.)
            else
              f(:,:,n2+k,j)=f(:,:,n2,j)+(cp-cv)*&
                  (f(:,:,n2,j-1)-f(:,:,n2+k,j-1))+&
                  cv*log((z(n2+k)-z(n2))*density_scale+1.)
            endif
          else
            call fatal_error('bc_ism','only for irho, ilnrho, iuz or iss')
          endif
        enddo
!
      case default
        print*, "bc_ism ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_ism
!***********************************************************************
    subroutine set_consistent_density_boundary(f,dirn,boundtype,tb,rhob,lsuccess)
!
!  This subroutine checks, if the density paramters like type, topbot
!  and boundary value are set consistently with eg. the initial condition.
!
!  26-jun-12/dhruba+joern: coded
!
!  dirn       =     direction                    : 'x','y','z'
!  boundtype  =     type of boundary condition   : 'set','a',...
!  tb         =     top or bottom boundary       : 'top','bot'
!  rhob       =     value at the boundary        : 4.04, 8.35, 10.1
!  lsuccess   =     switch, if it was successful : .true., .false.
!
!  At the moment only the x-direction is implemented
!
      real, dimension (:,:,:,:) :: f
      real, intent(in) :: rhob
      character (len=bclen), intent(in) :: boundtype,tb,dirn
      logical, intent(out) :: lsuccess
!
      character (len=bclen) :: btyp
      logical :: lconsistent=.true.
      real :: boundrho
!
! check for consistency
!
      if (ldensity_nolog) then
        boundrho=rhob
      else
        boundrho=log(rhob)
      endif
      btyp=trim(boundtype)
      select case (dirn)
        case ('x')
          select case (tb)
            case('bot')
              if ((btyp/=bcx12(ilnrho,1)) .or. (rhob/=fbcx(ilnrho,1))) then
                lconsistent=.false.
                bcx12(ilnrho,1)=btyp
                fbcx(ilnrho,1)=boundrho
                call boundconds_x(f,ilnrho,ilnrho)
                if (lroot) print*,'boundcond: density in x at the bottom set to: ', &
                    bcx12(ilnrho,1),', with the value ',fbcx(ilnrho,1)
              endif
            case('top')
              if ((btyp/=bcx12(ilnrho,2)) .or. (rhob/=fbcx(ilnrho,2))) then
                lconsistent=.false.
                bcx12(ilnrho,2)=btyp
                fbcx(ilnrho,2)=boundrho
                call boundconds_x(f,ilnrho,ilnrho)
                if (lroot) print*,'boundcond: density in x at the top set to: ', &
                    bcx12(ilnrho,2),', with the value ',fbcx(ilnrho,2)
              endif
            case default
              call fatal_error('set_consistent_density_boundary','topbot does not match any, aborting')
          endselect
        case ('y')
          call fatal_error('set_consistent_density_boundary','y direction not implemented yet')
        case ('z')
          call fatal_error('set_consistent_density_boundary','z direction not implemented yet')
        case default
          call fatal_error('set_consistent_density_boundary','you have to choose either x,y or z direction')
      endselect
      lsuccess=.true.
!
! density set consistently at the boundary.
!
    endsubroutine set_consistent_density_boundary
!***********************************************************************
    subroutine set_consistent_vel_boundary(f,dirn,boundtype,tb,comp,lsuccess)
!
!  This subroutine checks, if the velocity paramters like type and  topbot
!  are set consistently with eg. the initial condition.
!
!  14-sep-12/joern: coded, adapted from subroutine set_consistent_density_boundary
!
!  dirn       =     direction                    : 'x','y','z'
!  boundtype  =     type of boundary condition   : 'set','a',...
!  tb         =     top or bottom boundary       : 'top','bot'
!  comp       =     component of the velocity    : 'x','y','z'
!  lsuccess   =     switch, if it was successful : .true., .false.
!
!  At the moment only the x-direction is implemented
!
      real, dimension (:,:,:,:) :: f
      character (len=bclen), intent(in) :: boundtype,tb,dirn,comp
      logical, intent(out) :: lsuccess
!
      character (len=bclen) :: btyp
      logical :: lconsistent=.true.
!
! check for consistency
!
      btyp=trim(boundtype)
      select case (dirn)
        case ('x')
          select case (tb)
            case('bot')
              select case (comp)
                case('x')
                  if (btyp/=bcx12(iux,1)) then
                    lconsistent=.false.
                    bcx12(iux,1)=btyp
                    call boundconds_x(f,iux,iux)
                    if (lroot) print*,'boundcond: x velocity in x at the bottom set to: ',bcx12(iux,1)
                  endif
                case('y')
                  if (btyp/=bcx12(iuy,1)) then
                    lconsistent=.false.
                    bcx12(iuy,1)=btyp
                    call boundconds_x(f,iuy,iuy)
                    if (lroot) print*,'boundcond: y velocity in x at the bottom set to: ',bcx12(iuy,1)
                  endif
                case('z')
                  if (btyp/=bcx12(iuz,1)) then
                    lconsistent=.false.
                    bcx12(iuz,1)=btyp
                    call boundconds_x(f,iuz,iuz)
                    if (lroot) print*,'boundcond: z velocity in x at the bottom set to: ',bcx12(iuz,1)
                  endif
                case default
                  call fatal_error('set_consistent_vel_boundary','component does not match any, aborting')
              endselect
            case('top')
              select case (comp)
                case('x')
                  if (btyp/=bcx12(iux,2)) then
                    lconsistent=.false.
                    bcx12(iux,2)=btyp
                    call boundconds_x(f,iux,iux)
                    if (lroot) print*,'boundcond: x velocity in x at the top set to: ',bcx12(iux,2)
                  endif
                case('y')
                  if (btyp/=bcx12(iuy,2)) then
                    lconsistent=.false.
                    bcx12(iuy,2)=btyp
                    call boundconds_x(f,iuy,iuy)
                    if (lroot) print*,'boundcond: y velocity in x at the top set to: ',bcx12(iuy,2)
                  endif
                case('z')
                  if (btyp/=bcx12(iuz,2)) then
                    lconsistent=.false.
                    bcx12(iuz,2)=btyp
                    call boundconds_x(f,iuz,iuz)
                    if (lroot) print*,'boundcond: z velocity in x at the top set to: ',bcx12(iuz,2)
                  endif
                case default
                  call fatal_error('set_consistent_vel_boundary','component does not match any, aborting')
              endselect
            case default
              call fatal_error('set_consistent_vel_boundary','topbot does not match any, aborting')
          endselect
        case ('y')
          call fatal_error('set_consistent_velovity_boundary','y direction not implemented yet')
        case ('z')
          call fatal_error('set_consistent_vel_boundary','z direction not implemented yet')
        case default
          call fatal_error('set_consistent_vel_boundary','you have to choose either x,y or z direction')
      endselect
      lsuccess=.true.
!
! velocity set consistently at the boundary.
!
    endsubroutine set_consistent_vel_boundary
!***********************************************************************
    subroutine set_periodic_boundcond_on_aux(f,ivar)
!
! sets periodic boundary condition on auxiliar variables
!
      real, dimension (:,:,:,:) :: f
      integer :: ivar
!
      call bc_per_x(f,'top',ivar); call bc_per_x(f,'bot',ivar)
      call bc_per_y(f,'top',ivar); call bc_per_y(f,'bot',ivar)
      call bc_per_z(f,'top',ivar); call bc_per_z(f,'bot',ivar)
!
    endsubroutine set_periodic_boundcond_on_aux
!***********************************************************************
    subroutine tayler_expansion(f,topbot,j)
!
      character (len=bclen) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j
!
      select case (topbot)
      case ('top')
        f(l2+1,:,:,j) = + 4.*f(l2,:,:,j)   - 6.*f(l2-1,:,:,j) &
                        + 4.*f(l2-2,:,:,j) -    f(l2-3,:,:,j)
        f(l2+2,:,:,j) = +10.*f(l2,:,:,j)   -20.*f(l2-1,:,:,j) &
                        +15.*f(l2-2,:,:,j) - 4.*f(l2-3,:,:,j)
        f(l2+3,:,:,j) = +20.*f(l2,:,:,j)   -45.*f(l2-1,:,:,j) &
                        +36.*f(l2-2,:,:,j) -10.*f(l2-3,:,:,j)
      case ('bot')
        f(l1-1,:,:,j) = + 4.*f(l1,:,:,j)   - 6.*f(l1+1,:,:,j) &
                        + 4.*f(l1+2,:,:,j) -    f(l1+3,:,:,j)
        f(l1-2,:,:,j) = +10.*f(l1,:,:,j)   -20.*f(l1+1,:,:,j) &
                        +15.*f(l1+2,:,:,j) - 4.*f(l1+3,:,:,j)
        f(l1-3,:,:,j) = +20.*f(l1,:,:,j)   -45.*f(l1+1,:,:,j) &
                        +36.*f(l1+2,:,:,j) -10.*f(l1+3,:,:,j)
      endselect
!
    endsubroutine tayler_expansion
!***********************************************************************
    subroutine copy_BCs(isrc,itarg,num)
!
!  14-apr-15/MR: coded
!
      integer, intent(IN) :: isrc,itarg,num

      integer :: iet, ies

      if (isrc/=0) then

        ies = isrc +num-1
        iet = itarg+num-1

        if (nxgrid>1) then
          bcx (itarg:iet) = bcx (isrc:ies)
          bcx12(itarg:iet,:) = bcx12(isrc:ies,:)
        endif

        if (nygrid>1) then
          bcy (itarg:iet) = bcy (isrc:ies)
          bcy12(itarg:iet,:) = bcy12(isrc:ies,:)
        endif

        if (nzgrid>1) then
          bcz (itarg:iet) = bcz (isrc:ies)
          bcz12(itarg:iet,:) = bcz12(isrc:ies,:)
        endif

      endif
!
    endsubroutine copy_BCs
!***********************************************************************
endmodule Boundcond
