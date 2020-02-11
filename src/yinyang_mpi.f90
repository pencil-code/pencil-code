! $Id$
!
! MODULE_DOC: This module contains Yin-Yang related types and functions
! MODULE_DOC: which are incompatible with FORTRAN 95.
!
!***************************************************************
!
module Yinyang_mpi
!
  use Yinyang
  use Mpicomm
  use Cparam
!
  implicit none

  include 'yinyang_mpi.h'

  private
!
! Variables for z-averages.
!
  integer, parameter :: n_contrib_procs=4   ! supposed to be 3
  integer, dimension(2,n_contrib_procs) :: thrange_gap=0
  integer, dimension(n_contrib_procs) :: yinprocs=-1
  type(ind_coeffs), dimension(n_contrib_procs) :: indweights_zaver_gap
!
  integer, dimension(2) :: thrange_cap=0, capzprocs=-1
  type(ind_coeffs) :: indweights_zaver_cap
!
  integer, dimension(:,:,:), allocatable :: thranges_gap, thranges_cap
  integer :: offset_cap=0, caproot=-1, nycap=0
!
  contains

!*******************************************************************
    subroutine zsum_yy(fnamexy,iname,m,n,a)

      real, dimension(:,:,:), intent(INOUT) :: fnamexy
      integer,                intent(IN)    :: iname,m,n
      real, dimension(:),     intent(IN)    :: a
!
      integer :: i,j,ith

        if (thrange_gap(1,1)>0) then
!
!  collect contributions to the extended phi lines in the gap.
!
          do i=1,3
!
!  Such lines may come from at most 3 different Yin procs.
!
            if (yinprocs(i)==-1) exit

            do j=1,size(indweights_zaver_gap(i)%inds,3)
              ith=indweights_zaver_gap(i)%inds(m,n,j)
              if (ith==0) exit
              fnamexy(iname,:,ith)=fnamexy(iname,:,ith)+a*indweights_zaver_gap(i)%coeffs(m,n,j)
            enddo
          enddo
        endif
!
!  Collect contributions to the phi lines in the cap(s).
!
        if (thrange_cap(1)>0) then
          do j=1,size(indweights_zaver_cap%inds,3)
            ith=indweights_zaver_cap%inds(m,n,j)
            if (ith==0) exit
            fnamexy(iname,:,ith)=fnamexy(iname,:,ith)+a*indweights_zaver_cap%coeffs(m,n,j)
          enddo
        endif

    endsubroutine zsum_yy
!*******************************************************************
    subroutine reduce_zsum(arrm,temp)
!
!  Calculate z-sums (still depending on x and y).
!
!  25-apr-16/MR: outsourced from zaverages_xy
!
      use General, only: find_proc
      use Mpicomm, only: mpireduce_sum, mpisend_real, mpirecv_real, mpiwait, mpibarrier
      use Cdata

      real, dimension(:,:,:), intent(IN) :: arrm
      real, dimension(:,:,:), intent(OUT):: temp

      real, dimension(:,:,:), allocatable :: buffer
!
      integer :: iprocy, iprocz, iproc_yin, iproc_yang, i, offset, len, source, irequest
      integer, dimension(2) :: rng
      integer, dimension(3) :: sz,requests
!
!  Communicate over all processors along z beams.
!  The result is only present on the z-root processors (of Yin grid).
!  On Yang grid, procs (ipx,0,0) and (ipx,nprocy-1,nprocz-1) collect data 
!  for Yin-phi coordinate lines in polar caps.
!
        sz=(/size(arrm,1),size(arrm,2),size(arrm,3)/)
!
!  Summation of local fnamexy arrays (of Yin grid) into fsumxy of z root
!  processors of Yin grid.
!
        if (.not.lyang) call mpireduce_sum(arrm,temp,sz,idir=3)

        if (lyinyang) then
          if (lyang) then
!
!  If processor in Yang grid contributes to any phi coordinate line in gap send
!  data:    
! 
            if (thrange_gap(1,1)>0) then

              offset=0; irequest=0
              do i=1,3              ! These lines can come from at most 3 z root procs Yin grid.
                iproc_yin=yinprocs(i)
                if (iproc_yin==-1) exit
                rng=thrange_gap(:,i)-(thrange_gap(1,i)-1)+offset
                len=rng(2)-rng(1)+1
!
!if (iproc_yin==0 .or. iproc_yin==1) print'(a,3i3,3e12.5)', 'GAP: SEND:
!iproc_world, iproc_yin, len=', iproc_world, iproc_yin, &
!len,
!sum(fnamexy(:,:,rng(1):rng(2))),maxval(fnamexy(:,:,rng(1):rng(2))),minval(fnamexy(:,:,rng(1):rng(2)))
!  Local sum sent to relevant Yin zroot proc (non-blocking).
!
                irequest=irequest+1
!print*, 'iproc_world, iproc_yin,fnamexy=', iproc_world, iproc_yin,
!maxval(fnamexy(:,rng(1):rng(2),:)), minval(fnamexy(:,rng(1):rng(2),:))
                call mpisend_real(arrm(:,:,rng(1):rng(2)),(/sz(1),sz(2),len/), &
                                  iproc_yin,iproc_world,comm=MPI_COMM_WORLD,nonblock=requests(irequest))
                offset=offset+len   !  Update, as contributions for the <=3 procs are stacked in fnamexy.
              enddo
              do i=1,irequest
                call mpiwait(requests(irequest))
              enddo
            endif

            if (lcaproot) then                     ! If proc is cap collector,
              do iprocz=capzprocs(1),capzprocs(2)  ! iterate over all Yang procs
                do iprocy=0,nprocy-1               ! in cap.
                  rng=thranges_cap(:,iprocy,iprocz)
!print*, 'CAP: RECV: iproc, source, rng=', iproc, find_proc(ipx,iprocy,iprocz), rng
                  if (rng(1)==0) cycle
!
!  If processor in Yang grid contributes to any phi coordinate line in cap:
! 
                  source=find_proc(ipx,iprocy,iprocz)
!if (iproc_world==149) print*, 'source, rng:', source, rng
!rng, rng(2)-rng(1)+1
                  if (source==iproc) then          ! no self-communication
                    temp(:,:,rng(1):rng(2))=temp(:,:,rng(1):rng(2))+arrm(:,:,1+offset_cap:rng(2)-rng(1)+1+offset_cap)
                  else                             ! accumulate contribs from all cap procs.            
                    len=rng(2)-rng(1)+1
                    allocate(buffer(sz(1),sz(2),len))
!
!  Receive data from relevant Yang proc (blocking).
!
!print*, 'CAP: RECV: iproc, source, rng=', iproc, source, len
!rng, rng(2)-rng(1)+1
                    call mpirecv_real(buffer,(/sz(1),sz(2),len/),source,source)
                    temp(:,:,rng(1):rng(2)) = temp(:,:,rng(1):rng(2))+buffer
                    deallocate(buffer)
                  endif
                enddo
              enddo
            elseif (thrange_cap(1)>0) then
!
!  If proc is not cap collector, send to it if necessary (blocking).
!
              rng=thrange_cap-(thrange_cap(1)-1)+offset_cap
!print*, 'CAP:SEND,iproc, caproot=', iproc, caproot, rng(2)-rng(1)+1
              call mpisend_real(arrm(:,:,rng(1):rng(2)),(/sz(1),sz(2),rng(2)-rng(1)+1/), & 
                                caproot,iproc)
            endif

          elseif (lfirst_proc_z) then
!
!  Root processors of z beams in Yin grid accumulate contributions from
!  all Yang grid processors in gap.
!
            do iprocz=nprocz/3-1,2*nprocz/3
              do iprocy=0,nprocy-1
!
!  Range of theta of phi coordinate lines to which processor (ipx,iprocy,iprocz)
!  contributes.
! 
                rng=thranges_gap(:,iprocy,iprocz)
                if (rng(1)>0) then                                        ! if range non-empty
                  len=rng(2)-rng(1)+1
                  allocate(buffer(sz(1),sz(2),len))
                  iproc_yang=find_proc(ipx,iprocy,iprocz)+ncpus           ! world rank of source proc
                  call mpirecv_real(buffer,(/sz(1),sz(2),len/),iproc_yang,iproc_yang,comm=MPI_COMM_WORLD)
!if (iproc_world==0 .or. iproc_world==1) print'(a,3i3,3e12.5)', 'GAP: RECV:
!iproc_world, iproc_yang, len=', iproc_world, &
!iproc_yang, len, sum(buffer),maxval(buffer), minval(buffer)
!print*, 'iproc_world, iproc_yang,buffer=', iproc_world, iproc_yang,
!maxval(buffer), minval(buffer)
                  temp(:,:,rng(1):rng(2)) = temp(:,:,rng(1):rng(2))+buffer
                  deallocate(buffer)
                endif

              enddo
            enddo
          endif

          call mpibarrier
        endif
!
    endsubroutine reduce_zsum
!***********************************************************************
    subroutine initialize_zaver_yy(nlines,nycap_)
!
!  Initializations for calculation of z-sums on Yin-Yang grid.
!  Only valid for equidistant grid in y and z.
!  Returns total number of Yin-phi coordinate lines which cross the executing
!  proc.
!  Can be lines in the "gap" of the Yin grid phi<= Pi/4, phi>=7Pi/4
!  or lines which lie completely within Yang grid, then in the "polar caps"
!  theta<=Pi/4, theta>=3Pi/4. Note that executing proc can have a share of both
!  types.
!
!   11-mar-16/MR: coded
!
      use General, only: indgen, yy_transform_strip_other, find_proc, itoa
      use Mpicomm, only: mpireduce_sum_int, mpisend_int, mpirecv_int, mpiwait, mpibarrier, MPI_COMM_GRID, MPI_COMM_YZPLANE
      use Messages, only: warning
      use Cdata

      integer, intent(INOUT):: nlines
      integer, intent(OUT)  :: nycap_

      type(ind_coeffs) :: intcoeffs
      real, dimension(:,:,:), allocatable :: thphprime
      real, dimension(:), allocatable :: yloc, zloc
      integer, dimension(:), allocatable :: requests
      integer :: nok, noks, noks_all, nyl, nzgap, iprocy, iprocz, ifound, request, &
                 n_interproc_gaps, n_interproc_gaps_all, &
                 iproc_yin, iproc_yang, newlines, offset, npoints_tot, irequest, iyl, iyu, igyl, igyu
      integer, dimension(2) :: rng
      logical :: lwith_ghosts
      logical, save :: lcalled=.false.
return !!!
      if (lcalled) then
        nycap_=nycap
        return                               ! routine has already been called
      else
        lcalled=.true.
      endif

      nzgap=floor(2.*xyz0(3)/dz)-1           ! number of points in gap.
      nzgrid_eff=nzgrid+nzgap                ! total number of points on full Yin-phi coordinate line.
                                             ! could be reduced in polar caps.

      if (lyang) then

        npoints_tot=nlines*nprocy*nzgap        ! total number of points on lines continued into gap
        lwith_ghosts = nlines==my
        noks=0; n_interproc_gaps=0

        if (ipz>=nprocz/3-1 .and. ipz<=2*nprocz/3) then
!
!  These Yang procs may see phi coordinate lines which continue from Yin grid into
!  the gap.
!
          allocate(thphprime(2,nlines,nzgap),yloc(nlines),zloc(nzgap))
          zloc=xyz1(3)+indgen(nzgap)*dz        ! valid only for equidistant grids!

          nlines=0
          ifound=0; offset=0
          igyl=1
!
!  Loop over y-beam of processors with iprocz=0 in Yin grid.
!
          do iprocy=0,nprocy-1
!  
            if (lwith_ghosts) then
              if (iprocy==0) then
                iyl=m1; iyu=my
                igyu=igyl+ny+nghost-1
              elseif (iprocy==nprocy-1) then
                iyl=1; iyu=m2
                igyu=igyl+ny+nghost-1
              else
                iyl=1; iyu=my
                igyu=igyl+my-1
              endif
            else
              iyl=1; iyu=ny
              igyu=igyl+ny-1
            endif

            yloc(iyl:iyu)=ygrid(igyl:igyu)
!if (ipz==2) print*, 'iprocy, yloc=', iprocy, yloc(iyl),yloc(iyu)            
            if (lwith_ghosts) then
              igyl=igyu-2
              if (iprocy==0) yloc(1:nghost)=xyz0(2)-indgen(nghost)*dy     ! valid only for equidistant grids!
              if (iprocy==nprocy-1) yloc(iyu:)=xyz1(2)+indgen(nghost)*dy  ! valid only for equidistant grids!
            else
              igyl=igyu+1
            endif

            iproc_yin=find_proc(ipx,iprocy,0)         ! root proc of z-beam (ipx,iprocy)
!
!  yloc x zloc: strip of Yin-phi coordinate lines from Yin-proc iproc_yin.
!
            call yy_transform_strip_other(yloc,zloc,thphprime)
            nok=prep_interp(thphprime,intcoeffs,BILIN,ngap=n_interproc_gaps,th_range=thrange_gap(:,ifound+1))
!
!  nok: number of points of the strip claimed by executing proc; 
!  intcoeffs: interpolation data for these points; thrange_gap: y-range of
!  claimed points,
!  indices refer to local numbering of proc iproc_yin
!            
            if (nok>0) then
print'(a,4(1x,i4))', 'iproc_world,iproc_yin,nok,n_interproc_gap=', iproc_world,iproc_yin,nok,n_interproc_gaps

              ifound=ifound+1
              newlines=thrange_gap(2,ifound)-thrange_gap(1,ifound)+1
              nlines=nlines+newlines         ! accumulate detected lines.
              yinprocs(ifound)=iproc_yin     ! store Yin-proc from which detected lines emanate.
              noks=noks+nok                  ! accumulate claimed points.
!
!  Derive from interpolation data to which phi-coordinate lines a point (m,n)
!  contributes with which weight.
!
              call coeffs_to_weights(intcoeffs,indweights_zaver_gap(ifound))
              where (indweights_zaver_gap(ifound)%inds>0) &
                indweights_zaver_gap(ifound)%inds=indweights_zaver_gap(ifound)%inds-(thrange_gap(1,ifound)-1)+offset
!
! Entries in indweights_zaver_gap%inds now refer to number of line in local
! fnamexy.    
!
!print*, 'GAP: iproc_world, ifound, thrange_gap=', iproc_world, ifound,
!thrange_gap(:,ifound)                     
!if (ipz==nprocz/3.and.ipy==0)
!print*,'iproc_world,ifound,indweights_zaver_gap=', indweights_zaver_gap(ifound)

              rng=thrange_gap(:,ifound)
              offset=offset+newlines         ! offset of line number for proc iproc_yin in local fnamexy.

            else
              rng=(/0,0/)
            endif
!if (rng(1)/=0) print*, 'SEND: iproc_yin, iproc_world, rng=', iproc_yin, iproc_world, rng
!  Tell z-root proc of Yin grid, which of its phi coordinate lines are detected
!  within executing proc (maybe none).
!
            call mpisend_int(rng,2,iproc_yin,iproc_world,MPI_COMM_WORLD)

            if (ifound>3) stop               ! lines from at most 3 Yin procs expected.

          enddo

          offset_cap=nlines                  ! total number of detected gap lines = offset for cap lines.
          deallocate(thphprime,yloc,zloc)
!print*, 'GAP: iproc_world,yinprocs=', iproc_world,yinprocs
        endif

        call mpireduce_sum_int(noks, noks_all,comm=MPI_COMM_YZPLANE) ! sum up number of claimed gap points over Yang grid.
        call mpireduce_sum_int(n_interproc_gaps, n_interproc_gaps_all,comm=MPI_COMM_YZPLANE) ! sum up number of claimed points in interprocessor gaps over Yang grid.
!if (iproc==0) print*, 'GAP: iproc_world,noks_all,n_interproc_gaps_all,npoints_tot=', iproc_world, noks_all, n_interproc_gaps_all, npoints_tot
        if (iproc==0.and.noks_all-n_interproc_gaps_all/4<npoints_tot) &
          call warning('initialize_zaver_yy',  &
                       'only '//trim(itoa(noks_all))//' points in Yin grid gap claimed by Yang procs (goal: ' &
                       //trim(itoa(npoints_tot))//')',ncpus)

        call mpi_barrier(MPI_COMM_GRID)

        nok=0; n_interproc_gaps=0
        if (ipz<=nprocz/3 .or. ipz>=2*nprocz/3-1) then
!
!  These procs may see Yin-phi coordinate lines which lie completely inside Yang
!  grid that is, in the polar caps.
!
          nycap=floor(xyz0(2)/dy)-1           ! number of Yin-phi coordinate lines in polar cap
                                              ! could be reduced for bigger pole distance
                                              ! of closest line.
          if (lwith_ghosts) nycap=nycap+nghost

          allocate(thphprime(2,nycap,nzgrid_eff),yloc(nycap),zloc(nzgrid_eff))

          zloc=indgen(nzgrid_eff)*dz
!if (iproc==0) print*, 'zloc=', zloc(1), zloc(nzgrid_eff)
          if (ipz<=nprocz/3) then
!
! Southern cap
!
            yloc=xyz1(2)+indgen(nycap)*dy
            if (lwith_ghosts) yloc=yloc-dy*nghost

            caproot=find_proc(ipx,0,0)                 ! "cap collector" has lowest rank in layer ipx.
            capzprocs=(/0,nprocz/3/)                   ! range of z ranks in cap.      
!if (iproc==caproot)print*, 'South: yloc=', yloc

          elseif (ipz>=2*nprocz/3-1) then
!
! Northern cap
!
            yloc=xyz0(2)-(nycap+1-indgen(nycap))*dy
            if (lwith_ghosts) yloc=yloc+dy*nghost

            caproot=find_proc(ipx,nprocy-1,nprocz-1)   ! "cap collector" has highest rank in layer ipx.
            capzprocs=(/2*nprocz/3-1,nprocz-1/)        ! range of z ranks in cap.     
!if (iproc==caproot) print*, 'North: yloc=', yloc
          endif
!if (iproc==caproot) print*, 'iproc, yloc=', iproc, yloc(1), yloc(nycap)
!
!  yloc x zloc: strip of all Yin-phi coordinate lines in cap.
!
          call yy_transform_strip_other(yloc,zloc,thphprime)
          nok=prep_interp(thphprime,intcoeffs,BILIN,ngap=n_interproc_gaps,th_range=thrange_cap)
!
!  nok: number of points of the strip claimed by executing proc; 
!  intcoeffs: interpolation data for these points; thrange_cap: y-range of
!  claimed points,
!  indices refer to bunch of all Yin-phi coordinate lines in cap.
!            
          if (nok>0) then

            nlines=nlines+thrange_cap(2)-thrange_cap(1)+1     ! accumulate detected lines.            
!
!  Derive from interpolation data to which phi-coordinate lines a point (m,n)
!  contributes with which weight.
!
            call coeffs_to_weights(intcoeffs,indweights_zaver_cap)
!print*, 'CAP: iproc_world, thrange_cap=', iproc_world, thrange_cap, &
!minval(indweights_zaver_cap%inds), maxval(indweights_zaver_cap%inds)
            where (indweights_zaver_cap%inds>0) &
              indweights_zaver_cap%inds=indweights_zaver_cap%inds-(thrange_cap(1)-1)+offset_cap
!
! Entries in indweights_zaver_cap%inds now refer to number of line in local
! fnamexy.    
!
          endif

          if (iproc==caproot) then      ! if proc is cap collector

            lcaproot=.true.
            allocate(thranges_cap(2,0:nprocy-1,capzprocs(1):capzprocs(2)))
!if (iproc==caproot) print*,'iproc_world, nycap=', iproc_world, nycap
            allocate(requests(nprocy*(capzprocs(2)-capzprocs(1)+1)-1))
!
!  Cap collector receives from each Yang proc in cap, which of the cap lines
!  are detected
!  within proc (maybe none) and stores this in thranges_cap..
!
            irequest=1
            do iprocz=capzprocs(1),capzprocs(2)
              do iprocy=0,nprocy-1
                ipp=find_proc(ipx,iprocy,iprocz)
                if (iproc==ipp) then            ! no self-communication
                  thranges_cap(:,iprocy,iprocz)=thrange_cap
                else
                  call mpirecv_int(thranges_cap(:,iprocy,iprocz),2,ipp,ipp, &
                                   nonblock=requests(irequest))
                  irequest=irequest+1
                endif
              enddo
            enddo

            do irequest=1,size(requests)
              call mpiwait(requests(irequest))
            enddo
          else
!
!  Tell cap collector, which cap lines are detected
!  within executing proc (maybe none).
!
            call mpisend_int(thrange_cap,2,caproot,iproc)
          endif
!print*, 'iproc_world, offset_cap,nlines=', iproc_world, offset_cap,nlines
        endif

        call mpireduce_sum_int(nok, noks_all,comm=MPI_COMM_YZPLANE)    ! sum up number of claimed cap points (both caps) over all Yang procs.
        call mpireduce_sum_int(n_interproc_gaps, n_interproc_gaps_all,comm=MPI_COMM_YZPLANE) ! sum up number of claimed points in interprocessor gaps over Yang grid.
!if (iproc==0) print*, 'CAP: iproc_world,noks_all,n_interproc_gaps_all,npoints_tot=', iproc_world, noks_all, n_interproc_gaps_all, 2*nycap*nzgrid_eff
        if (iproc==0.and.noks_all-n_interproc_gaps_all/4<2*nycap*nzgrid_eff) &
          call warning('initialize_zaver_yy',  &
                       'only '//trim(itoa(noks_all))//' points in polar caps claimed by Yang procs (goal: ' &
                       //trim(itoa(2*nycap*nzgrid_eff))//')',ncpus)  ! ncpus is the global rank of the Yang grid root processor

        nycap_=nycap

      elseif (lfirst_proc_z) then
!
!  Yin z-beam root proc receives from each Yang procs in gap, which of the gap
!  lines are detected
!  within proc (maybe none) and stores this in thranges_gap..
!
        allocate(thranges_gap(2,0:nprocy-1,nprocz/3-1:2*nprocz/3))
        allocate(requests(nprocy*(nprocz/3+2)))
        irequest=1
        do iprocz=nprocz/3-1,2*nprocz/3
          do iprocy=0,nprocy-1
            iproc_yang=find_proc(ipx,iprocy,iprocz)+ncpus
!print*, 'RECV: iproc_yang, iproc_world=', iproc_yang, iproc_world
            call mpirecv_int(thranges_gap(:,iprocy,iprocz),2,iproc_yang,iproc_yang, &
                             MPI_COMM_WORLD,nonblock=requests(irequest))
            irequest=irequest+1
          enddo
        enddo

        do irequest=1,size(requests)
          call mpiwait(requests(irequest))
        enddo

      endif

      call mpibarrier

    endsubroutine initialize_zaver_yy
!***********************************************************************
end module
