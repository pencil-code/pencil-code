! $Id$
!
!***********************************************************************
module Posix
!
! Make C function usleep with argument in microseconds available.
!
  use, intrinsic :: iso_c_binding, only: c_int, c_int32_t
  implicit none

  interface
    ! int usleep(useconds_t useconds)
    function c_usleep(useconds) bind(c, name='usleep')
      import :: c_int, c_int32_t
      integer(kind=c_int32_t), value :: useconds
      integer(kind=c_int)            :: c_usleep
    endfunction c_usleep
  endinterface

endmodule Posix
!***********************************************************************
program remesh
!
!  Program to remesh existing run. The remeshed run get remesh_par 
!  times more grid points in each direction. remesh_par should be
!  a multiple of 2.
!
!  The remeshed result is spread to mulx times as many processors
!  in x direction, muly times as many processors
!  in y direction and to mulz times as many processors in z 
!  direction.
! 
!  The destination directory is specified in remesh.in (this way of 
!  doing it seem a little bit dangerous since specifying the wrong
!  destination will cause the loss of the possibly valuable var.dat 
!  files in that directory)
!
!  remesh.csh now generates proc directory tree (if it doesn't exist)
!  and copies param.nml and proc0/seed.dat to all proc subdirectories
!
!  Timings:
!  from 256^3 (64 procs) to 512^3 (128 procs): 14 min on 2Ghz linux
!
!-----------------------------------------------------------------------
!   18-oct-02/nils: coded
!   13-nov-02/nils: adapted to spread data to more processors
!   18-nov-02/axel: proc dir tree written, seed.dat copies to all
!   22-jun-03/axel: added reading and writing Lx,Ly,Lz (now needed)
!   25-mar-06/axel: added dx_1, dx_tilde, etc, for uniform mesh *only*
!   05-dec-09/simon: added spreading in mulx processors in x direction
!
  use Common
  use Mpicomm
  use Syscalls, only: system_cmd
  use General   !, only: get_from_nml_str,get_from_nml_real,get_from_nml_log,convert_nml
  use Messages, only: memory_usage
  use Posix
!
  implicit none
!
!  define parameters
!
  character (len=130) :: file,destination,file_new,varfile='var.dat'
  character (len=130) :: file2,dimfile,dimfile_loc,dimfile2_loc,gridfile
  character (len=4)   :: overwrite='yes'
  character (len=6)   :: ch,chx,chy,chz
  character (len=1)   :: prec
  integer :: i,j,k,itx=1,ity=1,itz=1
  integer :: kk,jj,ii,cpu
  integer :: dummy,ipxx,ipyy,ipzz
  integer :: mxout_grid,myout_grid,mzout_grid
  integer :: addx, addy, addz
  integer :: iproc_new, cpu_count=1
  integer, dimension (mprocs) :: cpu_global

  real, dimension(:,:,:,:), allocatable :: acoll, a
  real, dimension(mxcoll) :: xcoll
  real, dimension(mycoll) :: ycoll
  real, dimension(mzcoll) :: zcoll

! SC: added axtra dimension in x direction to 'f' array
  real, dimension (mmx_grid,mmy_grid,mmz_grid,mvar) :: f
  real, dimension (:,:,:,:,:), allocatable :: ff
  real, dimension (mmx_grid) :: rx,rdx_1,rdx_tilde
  real, dimension (mmy_grid) :: ry,rdy_1,rdy_tilde
  real, dimension (mmz_grid) :: rz,rdz_1,rdz_tilde
  real, dimension (mmx,mprocs) :: rrx
  real, dimension (mmy,mprocs) :: rry
  real, dimension (mmz,mprocs) :: rrz
  real :: tdummy, t_sp
  real :: facx, facy, facz, x0loc, y0loc, z0loc, w1, w2,dang,dang_old
  integer :: cpu_local,iyy,icpu,out_size,io_len
  integer :: counx, couny, counz, xstart,xstop, ystart,ystop, zstart, zstop
  integer :: nprocxx, nprocyy, nproczz
  logical :: lexist, lshort, lstop=.false.,lok
  integer :: idx, ifxa, ifxe, idy, ifya, ifye, idz, ifza, ifze, iv, icpu0, icpu1, nprocs_rem, ires
  integer :: i1x, i2x, i1y, i2y, i1z, i2z, isx, isy, isz, srcproc, idpx, idpy, idpz, stat
  character(LEN=128) :: clperi, msg
  integer, dimension(3) :: layout_src, layout_dst, layout_rem

  integer(kind=ikind8) :: mmw=mmx*mmy*mmz, mmw_grid=mmx_grid*mmy_grid*mmz_grid, mwcoll=mxcoll*mycoll*mzcoll

  call mpicomm_init_min

  if (lroot) print*,'Per process memory for arrays f:',mmw_grid*mvar,', ff:',mmw*mvar*mprocs, &
                    ', acoll:',mvar*mwcoll,' (Multiply with number of bytes per element!)'

  if (divx<=0) then
    if (lroot) print*, 'Error: divx<=0 not allowed, set it to 1!'
    lstop=.true.
  endif
  if (divy<=0) then
    if (lroot) print*, 'Error: divy<=0 not allowed, set it to 1!'
    lstop=.true.
  endif
  if (divz<=0) then
    if (lroot) print*, 'Error: divz<=0 not allowed, set it to 1!'
    lstop=.true.
  endif

  if (mod(mulx,divx) == 0 .and. divx /= 1 ) then
    if (lroot) print*, 'mulx divisible by divx! Redefine mulx='//trim(itoa(mulx/divx))//', divx=1!'
    lstop=.true.
  endif
  if (mod(muly,divy) == 0 .and. divy /= 1) then
    if (lroot) print*, 'muly divisible by divy! Redefine muly='//trim(itoa(muly/divy))//', divy=1!'
    lstop=.true.
  endif
  if (mod(mulz,divz) == 0 .and. divz /= 1 ) then
    if (lroot) print*, 'mulz divisible by divz! Redefine mulz='//trim(itoa(mulz/divz))//', divz=1!'
    lstop=.true.
  endif

  if ( (all((/mulx,muly,mulz/)==1) .and. all((/divx,divy,divz/)==1)) .and. &
      all((/remesh_parx,remesh_pary,remesh_parz/)==1 )) then
    if (lroot) print*, 'No remeshing needed!'
    lstop=.true.
  endif

  if (lroot) then 
    clperi=get_from_nml_str('LPERI',trim(datadir)//'/param2.nml',lvec=.true.)
    if (clperi=='') then
      print*, 'lperi could not be identified!'
      lstop=.true.
    endif
  endif

  if (lstop) then
    call mpifinalize_min
    stop
  endif

  call mpibcast(clperi)
  call convert_nml(trim(clperi),lperi)
!
  if (lroot) then 
    lyinyang=get_from_nml_log('LYINYANG',lexist,trim(datadir)//'/param.nml')
    if (.not.lexist) then
      print*, 'lyinyang could not be identified!'
      lyinyang=.false.
    endif
  endif
  call mpibcast(lyinyang)
!
  if (lyinyang) then
!
!  Read rel_dang for Yin-Yang grid.
!
    if (lroot) then 
      rel_dang=get_from_nml_real('REL_DANG',lexist,trim(datadir)//'/param2.nml')
      if (.not.lexist) then
        print*, 'rel_dang could not be identified!'
        call mpifinalize_min
        stop
      endif
    endif
    call mpibcast(rel_dang)
    dang_old=rel_dang*min(1./max(1,(nygrid-1)),3./max(1,(nzgrid-1)))*0.5*pi

    if (remesh_pary/=remesh_parz) then
      if (lroot) print*,'remesh_pary/=remesh_parz - this is not recommendable for Yin-Yang grids!'
      dang=dang_old/sqrt(float(remesh_pary*remesh_parz))
    else
      dang=dang_old/remesh_pary
    endif

  endif
!
!  Read input parameters from remesh.in
!
  if (lroot) then 
    open(1,FILE='remesh.in',FORM='formatted')
    read(1,'(a)') destination
    close(1)
!
!  Read input parameters from varfile.in
!
    inquire(FILE='varfile.in',EXIST=lexist)
    if (lexist) then
      open(1,FILE='varfile.in',FORM='formatted')
      read(1,'(a)') varfile
      close(1)
    endif
  endif
  call mpibcast(destination)
  call mpibcast(varfile)
!
  if (lroot) then
    print*,'mx,my,mz,mvar=',mx,my,mz,mvar
!
!  Print out parameters of conversion
!
    print*,'Will increase processor numbers in [xyz] by mul[xyz]=',mulx,muly,mulz
    print*,'Will decrease processor numbers in [xyz] by div[xyz]=',divx,divy,divz
    print*,'Will remesh by factor remesh_par[xyz]=',remesh_parx,remesh_pary,remesh_parz
    print*,'destination='//destination
    print*,'varfile='//varfile
!
  endif
!
!  Confirm remeshing
!
  if (nprocs==1) then
    write(*,*)
    write(*,*) 'Do you really want to overwrite ', trim(destination),'?'
    write(*,'(a$)') ' (Answer yes or no): '
    read(*,*,end=9) overwrite
  else
    overwrite='yes'
  endif
9 if ( overwrite == 'yes') then
!
! Read global dim.dat in order to find precision
!
    if (lroot) then
      call safe_character_assign(dimfile,trim(datadir)//'/dim.dat')
      if (ip<8) print*,'Reading '//trim(dimfile)
      open(1,FILE=dimfile,FORM='formatted')
      read(1,*) dummy
      read(1,*) prec
      close(1)
    endif
    call mpibcast(prec)
!
!  Determine proc layouts
!
   layout_src=(/nprocx,nprocy,nprocz/)
   if (any(mod((/mulx*layout_src(1),muly*layout_src(2),mulz*layout_src(3)/),(/divx,divy,divz/)) /= 0)) then
     if (lroot) print*, 'div[xyz] results in non-integer destination proc number!'
     call mpifinalize_min
     stop
   endif

   layout_dst=(/mulx*layout_src(1)/divx,muly*layout_src(2)/divy,mulz*layout_src(3)/divz/)
   layout_rem=layout_src/(/divx,divy,divz/)    ! proc layout for remesh.x
   nprocs_rem=product(layout_rem)              ! number of processors for remesh.x
   if (lroot) then
     print*, 'Processor layouts of source     :', layout_src
     print*, '                     destination:', layout_dst
     print*, '                     remesh prog:', layout_rem
     flush(6)
   endif
   if (nprocs>1 .and. nprocs/=nprocs_rem) then
     if (lroot) print*, 'Inappropriate number of ranks: nprocs/=nprocs_rem', nprocs,nprocs_rem
     call mpifinalize_min
     stop
   endif
!
!  Determine new size of global array
!
    mxout_grid=remesh_parx*nxgrid+2*nghost
    myout_grid=remesh_pary*nygrid+2*nghost
    mzout_grid=remesh_parz*nzgrid+2*nghost
!
!  Determine new size of local array
!
    nprocxx=nprocx*mulx/divx
    nprocyy=nprocy*muly/divy
    nproczz=nprocz*mulz/divz
!
!  Write size of global array to destination/data/dim.dat
!
    if (lroot) then
      call safe_character_assign(dimfile, trim(destination)//'/'//trim(datadir)//'/dim.dat')
      if (ip<8) print*,'Writing '//trim(dimfile)
      open(1,file=dimfile)
      write(1,'(6i7)') mxout_grid,myout_grid,mzout_grid,mvar,maux,mglobal
      write(1,'(a)') prec
      write(1,'(3i3)') nghost, nghost, nghost
!
! SC: Added iprocz_slowest = 1 in order to solve an issue with the reading 
!  of the remeshed dim.dat file.
!  This should work fine for all cases.
!
      write(1,'(4i3)') nprocxx, nprocyy, nproczz, 1
      close(1)
    endif
!
! Factor for d[xyz].
!
    if (remesh_parx/=1.) then
      if (lperi(1)) then
        facx=1./remesh_parx
      else
        facx=(nxgrid-1.)/(remesh_parx*nxgrid-1.)    ! because dx=Lx/(nxgrid-1) !
      endif
    endif
    if (remesh_pary/=1.) then
      if (lperi(2)) then
        facy=1./remesh_pary
      else
        facy=(nygrid-1.)/(remesh_pary*nygrid-1.)
      endif
    endif
    if (remesh_parz/=1.) then
      if (lperi(3)) then
        facz=1./remesh_parz
      else
        facz=(nzgrid-1.)/(remesh_parz*nzgrid-1.)
      endif
    endif

    if (lroot) then
      print*,'mx,my,mz,mvar=',mx,my,mz,mvar
      print*,'lshear=',lshear
    endif

    if (nprocs==1) then
      icpu0=0; icpu1=nprocs_rem-1
    else
      icpu0=iproc; icpu1=iproc
    endif
!
    call mpibarrier
!       
!  Loop over number of CPU's in original run
!
yinyang_loop: &
    do iyy=0,0   !ncpus,ncpus   !suppressed yinyang-loop for now
    do icpu=icpu0,icpu1
!
! Serialize reading of chunks by the processes.
!
      if (iproc>0) call mpirecv_logical(lok,iproc-1,iproc-1,comm=MPI_COMM_WORLD)

      call find_proc_coords_general(icpu,layout_rem(1),layout_rem(2),layout_rem(3),ipx,ipy,ipz)
      !ipx = modulo(icpu, nprocx)
      !ipy = modulo(icpu/nprocx, nprocy)
      !ipz = icpu/(nprocx*nprocy)
      cpu=icpu+iyy

      allocate(acoll(mxcoll,mycoll,mzcoll,mvar),stat=stat)
      if (stat/=0) print*, 'allocation of acoll fails, iproc=', iproc

      i1z=1; i2z=mz
      do isz=0,divz-1
      i1y=1; i2y=my;
      do isy=0,divy-1
      i1x=1; i2x=mx;
      do isx=0,divx-1

      srcproc=find_proc_general(ipx*divx+isx,ipy*divy+isy,ipz*divz+isz,layout_src(1),layout_src(2),layout_src(3))
!
!  Create name of varfile (this is typically var.dat) 
!
      ip=10
      call chn(srcproc,ch)
      call safe_character_assign(file,trim(datadir)//'/proc'//trim(ch)//'/'//trim(varfile))
      if (lroot.and.ip<8) print*,'Reading '//trim(file)
!
! Check for existence of varfile.
!
      inquire(FILE=trim(file),EXIST=lexist)
      if (lexist) then
!
!  File found for iyy=ncpus, i.e. iproc=ncpus --> a Yin-Yang grid supposed.
!
        if (lyinyang.and.iyy==ncpus.and..not.lyang) then
          lyang=.true.
          if (lroot) print*, 'This run is treated as a Yin-Yang one!'
          if (muly/=mulz) then
            if (lroot) print*, 'muly/=mulz --> nprocz=3*nprocy violated for new layout!'
            call mpifinalize_min
            stop
          endif
        endif
      elseif (iyy==0 .or. lyang) then
        print *, 'Varfile for proc=', cpu, 'not found!'
        stop
      else
!
!  File not found for proc=ncpus: If Yin-Yang run stop, 
!                                 otherwise leave loop.
        if (lyinyang) then
          print *, 'Varfile for proc=', cpu, 'not found!'
          stop
        else
          exit
        endif
      endif

      !print*,'Remeshing processor',cpu+1,'of',ncpus+iyy,'processors'
      if (lroot.and.lyang) print*, '(Yang grid)'
!
!  Read varfile (this is typically var.dat).
!
      open(1,file=file,form='unformatted',action='read')
!
! Possibility to jump here from below
!
!print*, 'i1x,i2x,i1y,i2y,i1z,i2z=', i1x,i2x,i1y,i2y,i1z,i2z
      read(1,iostat=stat,iomsg=msg) acoll(i1x:i2x,i1y:i2y,i1z:i2z,:)
      if (stat/=0) print*,'Error reading file "',trim(file),'"; iproc=',iproc,'; Message: '//trim(msg)
!
!  try first to read with deltay, but if that fails then
!  go to the next possibility without deltay.
!
      lshort=.not.lshear
      if (lshear) then
        read(1,err=991) t_sp,x,y,z,dx,dy,dz,deltay
        if (lroot) print*,'read deltay=',deltay
        goto 993
!
!  this is the end of the program.
!  The following lines are only to be accessed.
!  Here we try to read data without deltay.
!
991     continue
        lshort=.true.
        backspace(1)
      endif
!
!  try again if lshort
!
      if (lshort) then
        read(1,err=994,end=994) t_sp,x,y,z,dx,dy,dz
        goto 993
994     continue
        if (lroot) print*,'lshear,t_sp=',lshear,t_sp
        !backspace(1)
        !read(1) t_sp !,x,y,z !,dx,dy,dz
      endif
!
993   continue
      t=t_sp
      close(1)
!
!  Read grid.dat (needed?)
!
      call safe_character_assign(file,trim(datadir)//'/proc'//trim(ch)//'/grid.dat')
      if (ip<8) print*,'Reading '//trim(file)
      open(1,file=file,form='unformatted',action='read')
      read(1) tdummy,x,y,z,dx,dy,dz
      read(1) dx,dy,dz
      read(1) Lx,Ly,Lz
      read(1) dx_1,dy_1,dz_1
      close(1)

      xcoll(i1x:i2x)=x
      ycoll(i1y:i2y)=y
      zcoll(i1z:i2z)=z

      i1x=i1x+nx; i2x=i2x+nx; 
      enddo
      i1y=i1y+ny; i2y=i2y+ny; 
      enddo
      i1z=i1z+nz; i2z=i2z+nz;
      enddo    ! end collection loop

      if (iproc<nprocs-1) then
        ires=c_usleep(400000)         ! 0.4 sec appropr. for LUMI
        call mpisend_logical(lok,iproc+1,iproc,comm=MPI_COMM_WORLD)
      endif

      !goto 113         ! for testing: neither remeshes nor writes any data then.
!
!  Finding global processor number for new run
!
      cpu_count=1
      do counz=0,mulz-1
        do couny=0,muly-1
          do counx=0,mulx-1
            !iproc_new &
            !    =ipz*nprocy*nprocx*mulz*muly*mulx+counz*muly*mulx*nprocx*nprocy &
            !    +ipy*nprocx*muly*mulx+couny*mulx*nprocx &
            !    +ipx*mulx+counx
!!$            iproc_new=ipz*(nprocy*mulz*muly)+ipy*muly &
!!$                +couny+nprocy*muly*counz
            iproc_new = find_proc_general(mulx*ipx+counx,muly*ipy+couny,mulz*ipz+counz,layout_dst(1),layout_dst(2),layout_dst(3))
            cpu_global(cpu_count)=iproc_new
            cpu_count=cpu_count+1
          enddo
        enddo
      enddo
!print*, 'cpu, cpu_global=', cpu, cpu_global
!
!  Calculate new grid.
!  Should now also work for *coarsening* the grid.
!
      if (remesh_parx/=1.) then
        if (lperi(1).or.ipx>0) then
          x0loc=xcoll(l1)-dx/2.
          dx=dx*facx
          x0loc=x0loc+dx/2.
        else
          x0loc=xcoll(l1)
          dx=dx*facx
        endif
!
!  Correspondingly for the inverse mesh spacing 
!  Do this currently only for uniform meshes.
!  TB checked!
!
        rdx_1=dx_1(1)/facx
        rdx_tilde=dx_tilde(1)/facx

        do i=ll1,ll2
          rx(i)=x0loc+(i-ll1)*dx
        enddo
!
!  Finding the ghost points in x.
!
        do i=1,nghost
          rx(i)=rx(l1)-(nghost+1-i)*dx
          rx(ll2+i)=rx(ll2)+i*dx
        enddo
      else
        rx(1:mx)=xcoll(1:mx); rdx_1(1:mx)=dx_1; rdx_tilde(1:mx)=dx_tilde
      endif
!     
      if (remesh_pary/=1.) then

        if (lperi(2).or.ipy>0) then
          y0loc=ycoll(m1)-dy/2.
          dy=dy*facy
          y0loc=y0loc+dy/2.
        else
          y0loc=ycoll(m1)
          dy=dy*facy
        endif
      
        rdy_1=dy_1(1)/facy
        rdy_tilde=dy_tilde(1)/facy

        do j=mm1,mm2
          ry(j)=y0loc+(j-mm1)*dy
        enddo
!
!  Finding the ghost points in y.
!
        do i=1,nghost
          ry(i)=ry(m1)-(nghost+1-i)*dy
          ry(mm2+i)=ry(mm2)+i*dy
        enddo
      else
        ry(1:my)=ycoll(1:my); rdy_1(1:my)=dy_1; rdy_tilde(1:my)=dy_tilde
      endif
!
      if (remesh_parz/=1.) then

        if (lperi(3).or.ipz>0) then
          z0loc=zcoll(n1)-dz/2.
          dz=dz*facz
          z0loc=z0loc+dz/2.
        else
          z0loc=zcoll(n1)
          dz=dz*facz
        endif
      
        rdz_1=dz_1(1)/facz
        rdz_tilde=dz_tilde(1)/facz

        do k=nn1,nn2
          rz(k)=z0loc+(k-nn1)*dz
        enddo
!
!  Finding the ghost points
!
        do i=1,nghost
          rz(i)=rz(n1)-(nghost+1-i)*dz
          rz(nn2+i)=rz(nn2)+i*dz
        enddo
      else
        rz(1:mz)=zcoll(1:mz); rdz_1(1:mz)=dz_1; rdz_tilde(1:mz)=dz_tilde
      endif
!
!  Interpolating f-array to increased number of mesh points if only one
!  direction is remeshed.
!  Yet to be tested.
!
      if (.false.) then   !m[xyz] -> m[xyz]coll !!!
      if (remesh_parx>1 .and. remesh_pary==1 .and. remesh_parz==1) then
        do i=1,mmx_grid
          do ii=2,mx
            if (rx(i)<=xcoll(ii).and.rx(i)>=xcoll(ii-1)) exit
          enddo
          w2=(rx(i)-xcoll(ii-1))/dx; w1=(xcoll(ii)-rx(i))/dx
          f(i,1:my,1:mz,:)=a(ii-1,1:my,1:mz,:)*w1+a(ii,1:my,1:mz,:)*w2
        enddo
      endif

      if (remesh_parx==1 .and. remesh_pary>1 .and. remesh_parz==1) then
        do i=1,mmy_grid
          do ii=2,my
            if (ry(i)<=ycoll(ii).and.ry(i)>=ycoll(ii-1)) exit
          enddo
          w2=(ry(i)-ycoll(ii-1))/dy; w1=(ycoll(ii)-ry(i))/dy
          f(1:mx,i,1:mz,:)=a(1:mx,ii-1,1:mz,:)*w1+a(1:mx,ii,1:mz,:)*w2
        enddo
      endif

      if (remesh_parx==1 .and. remesh_pary==1 .and. remesh_parz>1) then
        do i=1,mmz_grid
          do ii=2,mz
            if (rz(i)<=zcoll(ii).and.rz(i)>=zcoll(ii-1)) exit
          enddo
          w2=(rz(i)-zcoll(ii-1))/dz; w1=(zcoll(ii)-rz(i))/dz
          f(1:mx,1:my,i,:)=a(1:mx,1:my,ii-1,:)*w1+a(1:mx,1:my,ii,:)*w2
        enddo
      endif
      endif
!print*, 'rz=', cpu, rz
!
!  Scattering f-array to increased number of mesh points. 
!
      itx=remesh_parx
      do i=0,mmx_grid-(2-remesh_parx),remesh_parx
        do ii=0,remesh_parx-1
          ity=remesh_pary
          do j=0,mmy_grid-(2-remesh_pary),remesh_pary
            do jj=0,remesh_pary-1
              itz=remesh_parz
              do k=0,mmz_grid-(2-remesh_parz),remesh_parz
                do kk=0,remesh_parz-1
!
! Tripling or doubling
!
                  if (remesh_parx == 3) then
                    addx=ii 
                    if (itx == 3)    addx=1
                    if (itx == 2) then
                      if (addx == 0) addx=1
                    endif
                    if (itx == mxcoll-2) then
                      if (addx == 2) addx=1
                    endif
                    if (itx == mxcoll-1) addx=0   
                  elseif (remesh_parx == 2) then
                    addx=ii 
                    if (itx == 2)    addx=1
                    if (itx == mxcoll-1) addx=0   
                  else
                    addx=1
                  endif

                  if (remesh_pary == 3) then
                    addy=jj 
                    if (ity == 3)    addy=1
                    if (ity == 2) then
                      if (addy == 0) addy=1
                    endif
                    if (ity == mycoll-2) then
                      if (addy == 2) addy=1
                    endif
                    if (ity == mycoll-1) addy=0   
                    addz=kk 
                    if (itz == 3)    addz=1
                    if (itz == 2) then
                      if (addz == 0) addz=1
                    endif
                    if (itz == mzcoll-2) then
                      if (addz == 2) addz=1
                    endif
                    if (itz == mzcoll-1) addz=0  
                  elseif (remesh_parz == 2) then
                    addy=jj
                    addz=kk
                    if (ity == 2)    addy=1
                    if (ity == mycoll-1) addy=0
                    if (itz == 2)    addz=1
                    if (itz == mzcoll-1) addz=0
                  else
                    addy=1
                    addz=1
                  endif
!
! SC: added 'i+addx' to 'f' array
!
!if (i+addx>mmx_grid .or. j+addy>mmy_grid .or. k+addz>mmz_grid) print*, 'cpu, itx,ity,itz=', cpu, itx,ity,itz
                  f(i+addx,j+addy,k+addz,:)=acoll(itx,ity,itz,:)
                enddo
                itz=itz+1
              enddo
            enddo
            ity=ity+1
          enddo
        enddo
!
! SC: moved the end of the itx loop here
!
        itx=itx+1
      enddo
!   
!if (cpu<4) print*, 'cpu,f=', cpu,f(:,8,8,2)
!
! SC: the spreading into different processors is now done after the remeshing
!     in all dimensions.
!
! Spreading results to different processors.
!
      deallocate(acoll)
      allocate(ff(mmx,mmy,mmz,mvar,mprocs),stat=stat)
      if (stat/=0) print*, 'allocation of ff fails, iproc=', iproc

      if (any((/mulx,muly,mulz/)/=1)) then
        do counx=1,mulx
          do couny=1,muly
            do counz=1,mulz
!                cpu_local=(couny-1)+muly*(counz-1)+1
              cpu_local=(counx-1)+mulx*(couny-1)+mulx*muly*(counz-1)+1
              xstart=1+(counx-1)*nnx
              xstop=counx*nnx+2*nghost
              ystart=1+(couny-1)*nny
              ystop=couny*nny+2*nghost
              zstart=1+(counz-1)*nnz
              zstop=counz*nnz+2*nghost
!
!NILS: This is not OK if mulx>1....should be fixed
! SC: @NILS: is it OK now?
! SC: changed 'ff' array from 'ff(i+addx,...)' to 'ff(:,...)'
! SC: added to 'f' array 'xstart:xstop'
              ff(:,:,:,:,cpu_local)=f(xstart:xstop,ystart:ystop,zstart:zstop,:)
!--             if (itx == 2) then
              rrx(:,cpu_local)=rx(xstart:xstop)
              rry(:,cpu_local)=ry(ystart:ystop)
              rrz(:,cpu_local)=rz(zstart:zstop)
!--             endif
            enddo
          enddo
        enddo
      else
        ff(:mmx_grid,:mmy_grid,:mmz_grid,:mvar,1)=f
        rrx(:mmx_grid,1)=rx; rry(:mmy_grid,1)=ry; rrz(:mmz_grid,1)=rz
      endif
!
! SC: former end of the itx loop
!
!  Smoothing data if any of the remesh_par* is unequal to 1.
!
      if (any((/remesh_parx,remesh_pary,remesh_parz/)/=1.)) then
        do j=1,mprocs
          do i=1,mvar
            call rmwig(ff(:,:,:,:,j),i,1.,.false.)
          enddo
        enddo
      endif
!
!  Write new var.dat 
!
      inquire(IOLENGTH=io_len) t_sp
      out_size=nnx*nny*nnz*io_len
      do i=1,mprocs
        if (lastaroth) then
          call find_proc_coords_general(cpu_global(i),layout_dst(1),layout_dst(2),layout_dst(3),idpx,idpy,idpz)
          call chn(idpx*nnx,chx); call chn(idpy*nny,chy); call chn(idpz*nnz,chz);
          call safe_character_assign(file_new,trim(destination)//'/'//trim(datadir)//'/allprocs/field-')
          call safe_character_assign(file2,'-segment-'//trim(chx)//'-'//trim(chy)//'-'//trim(chz)//'.mesh')
          do j=1,mvar
            !call chn(j,ch)
!print*, trim(file_new)//trim(fields(j))//trim(file2)
            open(91,file=trim(file_new)//trim(fields(j))//trim(file2),form='unformatted',access='direct',recl=out_size)
            write(91,rec=1) ff(nghost+1:nghost+nnx,nghost+1:nghost+nny,nghost+1:nghost+nnz,j,i)
            close(91)
          enddo
        else
          call chn(cpu_global(i),ch)
          call safe_character_assign(file_new,trim(datadir)//'/proc'//trim(ch)//'/var.dat')
          call safe_character_assign(file2,trim(destination)//'/'//trim(file_new))
          if (ip<8) print*,'Writing '//trim(file2)
!print*, 'cpu, cpu_global, rrx=', cpu, cpu_global(i), rrx(:,i)
          open(91,file=file2,form='unformatted')
          write(91) ff(:,:,:,:,i)
          if (lshear) then
            write(91) t_sp,rrx(:,i),rry(:,i),rrz(:,i),dx,dy,dz,deltay
            if (lroot) print*,'wrote deltay=',deltay
          else
            write(91) t_sp,rrx(:,i),rry(:,i),rrz(:,i),dx,dy,dz
          endif
          close(91)
        endif
      enddo
      deallocate(ff)
!
!  Write dim.dat for new local processor(s)
!
      do i=1,mprocs

        iproc_new=cpu_global(i)
        ipxx = modulo(iproc_new,nprocxx)
        ipyy = modulo(iproc_new/nprocxx, nprocyy)
        ipzz = iproc_new/(nprocxx*nprocyy)
        call chn(iproc_new,ch)
        call safe_character_assign(dimfile_loc,trim(datadir)//'/proc'//trim(ch)//'/dim.dat')
        call safe_character_assign(dimfile2_loc,trim(destination)//'/'//trim(dimfile_loc))
        if (ip<8) print*,'Writing ',dimfile2_loc
        open(1,file=dimfile2_loc)
        write(1,'(6i7)') mmx,mmy,mmz,mvar,maux,mglobal
        write(1,'(a)') prec
        write(1,'(3i3)') nghost, nghost, nghost
! SC: Added iprocz_slowest = 1 in order to solve an issue with the reading 
!  of the remeshed dim.dat file.
!  This should work fine for all cases.
        write(1,'(4i3)') ipxx, ipyy, ipzz, 1
        close(1)
      enddo
!
!  Write grid data for new local processor(s)
!
      do i=1,mprocs
        call chn(cpu_global(i),ch)
        call safe_character_assign(gridfile,trim(destination)// &
             '/'//trim(datadir)//'/proc'//trim(ch)//'/grid.dat')
        if (ip<8) print*,'Writing ',gridfile
        open(1,FILE=gridfile,FORM='unformatted')
!if (lroot) print*, '1st rec gridfile:', i, t_sp, size(rrx(:,i)), size(rry(:,i)), size(rrz(:,i)), dx, dy,dz
        write(1) t_sp,rrx(:,i),rry(:,i),rrz(:,i),dx,dy,dz
        write(1) dx,dy,dz
        write(1) Lx,Ly,Lz
        write(1) rdx_1,rdy_1,rdz_1
        write(1) rdx_tilde,rdy_tilde,rdz_tilde
        close(1)
      enddo
!
113   continue
    enddo
    enddo yinyang_loop

    if (lroot) print*, 'Remeshing successfully finished.'

  elseif (overwrite == 'no') then
    print*,'OK, -maybe some other time.'
  else
    print*,'You must answer yes or no!'
  endif

  call mpibarrier
  call memory_usage

  call mpifinalize_min
!
endprogram remesh
!***********************************************************************
    subroutine rmwig_1d(f,df,ivar,awig_input,idir)
!
!  Remove small scale oscillations (`wiggles') in the direction idir from
!  the ivar component of f (normally lnrho).
!  Unlike the version in the actual code, this one does not need
!  updating of the ghost zones, because remeshing works on data
!  where the ghost zones are already updated.
!
!  30-Aug-02/wolf: coded
!  12-nov-02/nils: adapted from equ.f90
!
      use Cdata
      use Common
!
      real, dimension (mmx,mmy,mmz,mvar) :: f,df
      real, dimension (nnx) :: tmp
      real :: awig_input
      integer :: ivar,idir
!
!  do loop over y and z
!  set indices and check whether communication must now be completed
!
      lfirstpoint=.true.        ! true for very first m-n loop
      do m=mmm1,mmm2
        do n=nnn1,nnn2
          call der6(f,ivar,tmp,idir,.true.)
          df(lll1:lll2,m,n,ivar) = 1./64.*tmp
        enddo
      enddo
!
!  Not necessary to do this in a (cache-efficient) loop updating in
!  timestep, since this routine will be applied only infrequently.
!
      f(lll1:lll2,mmm1:mmm2,nnn1:nnn2,ivar) = f(lll1:lll2,mmm1:mmm2,nnn1:nnn2,ivar) + &
           awig_input*df(lll1:lll2,mmm1:mmm2,nnn1:nnn2,ivar)
!
    endsubroutine rmwig_1d
!***********************************************************************
    subroutine rmwig(f,ivar,awig_input,explog)
!
!  30-Aug-02/wolf: coded
!  12-nov-02/nils: adapted from equ.f90
!
      use Cdata
      use Common
!
      real, dimension (mmx,mmy,mmz,mvar) :: f,df
      logical :: explog
      integer :: ivar
      real :: awig_input
!
!  Check whether we want to smooth the actual variable f, or exp(f)
!  The latter can be useful if the variable is lnrho or lncc.
!
      if (explog) then
        f(:,:,:,ivar)=exp(f(:,:,:,ivar))
        if (lroot) print*,'RMWIG: turn f into exp(f), ivar=',ivar
      endif
!
!  x, y, and z directions
!
      call rmwig_1d(f,df,ivar,awig_input,1)
      call rmwig_1d(f,df,ivar,awig_input,2)
      call rmwig_1d(f,df,ivar,awig_input,3)
!
!  Revert back to f if we have been working on exp(f)
!
      if (explog) then
        f(lll1:lll2,mmm1:mmm2,nnn1:nnn2,ivar)=alog(f(lll1:lll2,mmm1:mmm2,nnn1:nnn2,ivar))
        if (lroot) print*,'RMWIG: turn f back into alog(f), ivar=',ivar
      endif
!
    endsubroutine rmwig
!***********************************************************************
    subroutine der6(f,k,df,j,ignoredx)
!
!  Calculate 6th derivative of a scalar, get scalar
!    Used for hyperdiffusion that affects small wave numbers as little as
!  possible (useful for density).
!    The optional flag IGNOREDX is useful for numerical purposes, where
!  you want to affect the Nyquist scale in each direction, independent of
!  the ratios dx:dy:dz.
!
!   8-jul-02/wolf: coded
!  12-nov-02/nils: adapted from deriv.f90
!
      use Cdata
      use Common
!
      real, dimension (mmx,mmy,mmz,mvar) :: f
      real, dimension (nnx) :: df
      real :: fac
      integer :: j,k
      logical :: ignoredx
      logical :: igndx
!
      intent(in)  :: f,k,j,ignoredx
      intent(out) :: df
!
      if (ignoredx) then
        igndx = ignoredx
      else
        igndx = .false.
      endif
!
      if (j==1) then
        if (nxgrid/=1) then
          if (igndx) then; fac=1.; else; fac=1./dx**6; endif
          df=fac*(-20.* f(lll1:lll2,m,n,k) &
                  +15.*(f(lll1+1:lll2+1,m,n,k)+f(lll1-1:lll2-1,m,n,k)) &
                  - 6.*(f(lll1+2:lll2+2,m,n,k)+f(lll1-2:lll2-2,m,n,k)) &
                  +    (f(lll1+3:lll2+3,m,n,k)+f(lll1-3:lll2-3,m,n,k)))
        else
          df=0.
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          if (igndx) then; fac=1.; else; fac=1./dy**6; endif
          df=fac*(-20.* f(lll1:lll2,m  ,n,k) &
                  +15.*(f(lll1:lll2,m+1,n,k)+f(lll1:lll2,m-1,n,k)) &
                  - 6.*(f(lll1:lll2,m+2,n,k)+f(lll1:lll2,m-2,n,k)) &
                  +    (f(lll1:lll2,m+3,n,k)+f(lll1:lll2,m-3,n,k)))
        else
          df=0.
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          if (igndx) then; fac=1.; else; fac=1./dz**6; endif
          df=fac*(-20.* f(lll1:lll2,m,n  ,k) &
                  +15.*(f(lll1:lll2,m,n+1,k)+f(lll1:lll2,m,n-1,k)) &
                  - 6.*(f(lll1:lll2,m,n+2,k)+f(lll1:lll2,m,n-2,k)) &
                  +    (f(lll1:lll2,m,n+3,k)+f(lll1:lll2,m,n-3,k)))
        else
          df=0.
        endif
      endif
!
    endsubroutine der6
!***********************************************************************
    subroutine chn(cpu,ch)
!      
!  12-nov-02/nils: adapted from sub.f90
!
      character (len=6)  :: ch
      integer :: cpu
!
      ch='    '
      if (cpu.lt.0) stop 'chn: lt1'
      if (cpu.lt.10) then
        write(ch(1:1),'(i1)') cpu
      elseif (cpu.lt.100) then
        write(ch(1:2),'(i2)') cpu
      elseif (cpu.lt.1000) then
        write(ch(1:3),'(i3)') cpu
      elseif (cpu.lt.10000) then
        write(ch(1:4),'(i4)') cpu
      elseif (cpu.lt.100000) then
        write(ch(1:5),'(i5)') cpu
      elseif (cpu.lt.1000000) then
        write(ch(1:6),'(i6)') cpu
      else
        print*,'cpu=',cpu
        stop "ch: cpu too large"
      endif
!
    endsubroutine chn
!***********************************************************************
    subroutine safe_character_assign(dest,src)
!
!  Do character string assignement with check against overflow
!
!  08-oct-02/tony: coded
!  25-oct-02/axel: added tag in output to give name of routine
!
      character (len=*), intent(in):: src
      character (len=*), intent(inout):: dest
      integer :: dest_len, src_len
!
      dest_len = len(dest)
      src_len = len(src)
!
      if (dest_len<src_len) then
        print *, "safe_character_assign: ", &
             "RUNTIME ERROR: FORCED STRING TRUNCATION WHEN ASSIGNING '" &
             //src//"' to '"//dest//"'"
        dest=src(1:dest_len)
      else
        dest=src
      endif
!
    endsubroutine safe_character_assign
!***********************************************************************
