! $Id: poisson.f90 12460 2009-12-10 15:19:51Z sven.bingert $
!
!  This module solves the Poisson equation 
!    (d^2/dx^2 + d^2/dy^2 + d^2/dz^2 - h) f = RHS(x,y,z)
!  [which for h/=0 could also be called inhomogenous nonuniform Helmholtz
!  equation] for the function f(x,y,z). The difference between 
!  this and poisson.f90 is that here we solve it in a non-periodic 
!  Cartesian grid, expanding the grid to curb the periodicity of the 
!  Fourier transforms. 
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lpoisson=.true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
module Poisson
!
  use Cdata
  use Cparam
  use Fourier
  use Messages
!
  implicit none
!
  real :: kmax=0.0
  logical :: lrazor_thin=.false., lsemispectral=.false., lklimit_shear=.false.
  logical :: lexpand_grid=.false.
!
  include '../poisson.h'
!
  namelist /poisson_init_pars/ &
      lsemispectral, kmax, lrazor_thin, lklimit_shear, lexpand_grid
!
  namelist /poisson_run_pars/ &
      lsemispectral, kmax, lrazor_thin, lklimit_shear, lexpand_grid
!
!  Variables of the expanded grid are global to this file. 
!
  real :: r2_ext, r2_int
  real, dimension(2*nx)     :: xc,kkx_fft
  real, dimension(2*ny)     :: yc
  real, dimension(2*nygrid) :: kky_fft
  integer :: nnx,nny,iroot,nslice
  real    :: xc0,yc0,dxc1,dyc1
!
  contains
!***********************************************************************
    subroutine initialize_poisson()
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  18-oct-07/anders: adapted
!
      if (lrazor_thin) then
        if (nzgrid/=1) then
          if (lroot) print*, 'inverse_laplacian_fft: razon-thin approximation only works with nzgrid==1'
          call fatal_error('inverse_laplacian_fft','')
        endif
      endif
!
      r2_ext=r_ext**2
      r2_int=r_int**2
!
      if (lklimit_shear) kmax = kx_nyq/sqrt(2.0)
!
    endsubroutine initialize_poisson
!***********************************************************************
    subroutine inverse_laplacian(f,phi)
!
!  Dispatch solving the Poisson equation to inverse_laplacian_fft
!  or inverse_laplacian_semispectral, based on the boundary conditions
!
!  17-jul-2007/wolf: coded wrapper
!
      use General, only: keep_compiler_quiet
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,ny,nz) :: phi
!
      intent(inout) :: phi
!
      if (lcylindrical_coords) then 
        if (lroot) print*,'You are using cylindrical coordinates. '//&
             'Use poisson_cyl.f90 instead'
        call fatal_error("inverse_laplacian","")
      endif
!
      if (lspherical_coords) then 
        if (lroot) then 
          print*,'There is no poisson solver for spherical '
          print*,'coordinates yet. Please feel free to implement it. '
          print*,'Many people will thank you for that.'
        endif
        call fatal_error("inverse_laplacian","")
      endif
!
      if (lsemispectral) then
        call fatal_error('inverse_laplacian', &
            'this file is just for expand grid')
      else
        if (lexpand_grid) then
          call inverse_laplacian_expandgrid(phi)
        else
          call fatal_error('inverse_laplacian', &
              'this file is just for expand grid') 
        endif
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine inverse_laplacian
!***********************************************************************
    subroutine inverse_laplacian_expandgrid(phi)
!
!  Solve the Poisson equation in a 2D global disk by expanding the 
!  grid to avoid the periodicity of the Fourier transform. The geometry 
!  is as follows: 
!                             ...........
!                            .           . 3
!      -----                 ....-----....
!     |-----|3               .  |     |  . 2
!     |-----|2    --->       ...|.....|...
!     |-----|1               .  |     |  . 1
!      ----- 0               ....-----....
!                            .           . 0
!                             ...........
!
!    Original grid        Auxiliary (expanded) grid
!      
!
!  The density field is copied from the original (small) grid
!  to the auxiliary (large) grid. Outside the borders of the 
!  small grid, the density is set to zero. The Poisson equation
!  is then solved in the large grid. The solution (the gravitational
!  potential) is then copied from the large grid back to the small 
!  one. The concept is simple. Most of this subroutine and its 
!  children deal with the bookkeeping, since overlapping parts of 
!  the small/large grid will, in general, be in different processors. 
!
!  This solution curbs the periodicity of the Fourier transform. The 
!  large grid is periodic, but the small grid is far from the boundary 
!  and therefore little affected by it. The numerical and analytical 
!  solution for a galactic exponential disk agree within three percent 
!  (see Fig. 1 of Lyra et al. 2009, A&A, 493, 1125). 
!
!  25-may-2008/wlad: coded
!
      real, dimension (nx,ny,nz) :: phi
      real, dimension (2*nx,2*ny) :: nphi
      logical, dimension(0:ncpus-1,0:ncpus-1) :: lproc_comm_send
      logical, dimension(0:ncpus-1,0:ncpus-1) :: lproc_comm_recv
      logical :: lfirstcall=.true.
!
!  Identify version.
!
      if (lroot .and. ip<10) call svn_id( &
          '$Id: poisson.f90 12460 2009-12-10 15:19:51Z sven.bingert $')
!
!  Break if lshear or 3D
!
      if (lshear) call fatal_error('inverse_laplacian_expandgrid',&
          'not implemented for external shear')
!
!  Break if nprocx>1
!
      if (nprocx>1) &
           call fatal_error("inverse_laplacian_expandgrid",&
           "Not yet implemented for nprocx > 1")
!
      if (lroot.and.lfirstcall) print*,'Entered inverse_laplacian_expandgrid'
!
!  Shortcuts. The 'iroot' variable is of course for the root processor. 
!  The 'zslice=1' integer is because one can only deal with a 3D phi array. 
!
      iroot=0 ; nslice=1
!
!  Now do the magic
!
      call construct_large_grid
      call copy_density_large_grid(phi,nphi, &
                 lproc_comm_send,lproc_comm_recv)
      call calc_potential(nphi)
      call copy_potential_small_grid(nphi,phi, &
                 lproc_comm_send,lproc_comm_recv)
!      
      if (lfirstcall) lfirstcall=.false.
!
    endsubroutine inverse_laplacian_expandgrid
!***********************************************************************
    subroutine construct_large_grid
!
!  Define the expanded Cartesian axes. The x-axis is serial.
!
!  23-apr-10/wlad: coded
!
      real    :: xcn,ycn,dxc,dyc,Lxn,Lyn
      integer :: i,nnxgrid,nnygrid
!
      nnx=2*nx               ; nny=2*ny
      xcn=2*x(l2)+0.5*dx     ; xc0=-xcn
      ycn=xcn                ; yc0=-ycn
      nnxgrid=2*nxgrid       ; nnygrid=2*nygrid
!
      do i=1,nnx
        xc(i)=1.0*(i-1)        /(nnxgrid-1)*(xcn-xc0)+xc0
      enddo
      do m=1,nny
        yc(m)=1.0*(m-1+ipy*nny)/(nnygrid-1)*(ycn-yc0)+yc0
      enddo      
!
      if (ip<=8) print*,'inverse_laplacian_expandgrid: ', & 
           'x-limits of the expanded grid',xc(1),xc(nnx)
      if (ip<=8) print*,'inverse_laplacian_expandgrid: ', &
           'local y-limits of the expanded grid',yc(1),yc(nny)
!
!  Define the wavelengths of the expanded grid
! 
      dxc=xc(2)-xc(1)  ; dyc=dxc
      dxc1=1/dxc       ; dyc1=1/dyc
      Lxn=2*xc(nnx)    ; Lyn=Lxn
!
      kkx_fft = &
          cshift((/(i-(nnxgrid+1)/2,i=0,nnxgrid-1)/),+(nnxgrid+1)/2)*2*pi/Lxn
      kky_fft = &
          cshift((/(i-(nnygrid+1)/2,i=0,nnygrid-1)/),+(nnygrid+1)/2)*2*pi/Lyn

    endsubroutine construct_large_grid
!***********************************************************************
    subroutine copy_density_large_grid(phi,nphi,&
         lproc_comm_send,lproc_comm_recv)
!
!  Copy the small grid to the large one
!
!  23-apr-10/wlad: coded
!
      use Mpicomm
!
      real, dimension (nx,ny,nz) :: phi
      real, dimension (2*nx,2*ny) :: nphi
      real, dimension (nx,ny,0:ncpus-1) :: phi_send,phi_recv
      logical, dimension(0:ncpus-1) :: lproc_comm_loc
      logical, dimension(0:ncpus-1,0:ncpus-1) :: lproc_comm_send
      logical, dimension(0:ncpus-1,0:ncpus-1) :: lproc_comm_recv
      integer :: iy_serial, iproc_send
      integer :: i,j, mm
      real :: rr2
!
      if (lmpicomm) then 
!
        lproc_comm_loc=.false.
!
!  Define, based on (x,y), the processor iproc_send that the 
!  local processor iproc has to send its density field. Store
!  that in phi_send(:,:,iproc_send).
!
        do m=1,ny
          mm=m+m1-1
          iy_serial=nint((y(mm)-yc0)*dyc1)+1
          iproc_send=(iy_serial-1)/nny
          if (iproc_send/=iproc) then
            if (.not.lproc_comm_loc(iproc_send)) then
              lproc_comm_loc(iproc_send)=.true.
            endif
            phi_send(:,m,iproc_send)=phi(:,m,nslice)
          endif
        enddo
!
!  In parallel we need first to get the processors' send/recv matrices. 
!  Do it by now in the form of (ncpus,ncpus) logicals. In the future we 
!  want to avoid ncpus^2 arrays.
!
        call get_communication_matrix(lproc_comm_loc,&
             lproc_comm_send,lproc_comm_recv,iroot)
!
!  The logicals define all send/recv pairs. Each processor loops through 
!  all processors, looking for logicals set to true. When they are 
!  found, send data to or receive data from that processor. Piece by piece, 
!  each processor goes collecting the "phi_recv" density fields it receives.
!  These will be used to build nphi, the density field in the large grid.  
!
        do j=0,ncpus-1
          if (iproc/=j) then
            if (lproc_comm_send(iproc,j)) then
              if (ip<=8) print*,'poisson:proc ',iproc,' sending to proc',j
              call mpisend_real(phi_send(:,:,j),(/nx,ny/),j,111)
            endif
            if (lproc_comm_recv(iproc,j)) then
              if (ip<=8) print*,'poisson:proc ',iproc,' receiving from proc ',j
              call mpirecv_real(phi_recv(:,:,j),(/nx,ny/),j,111)
            endif
          else
!
!  By construction, no send-to-self exist. The parts of 
!  the small and large that overlap in the same 
!  processor are simply copied from phi to phi_recv. 
!
            phi_recv(:,:,iproc)=phi(:,:,nslice)     
          endif
        enddo
! 
      else 
!
!  Non-parallel. Just copy phi to phi_recv. 
!
        phi_recv(:,:,iroot)=phi(:,:,nslice)
!
      endif
!
!  Build nphi, the density field of the large grid, with 
!  the phi_recv received from different processors. Loop 
!  through the cells of the large grid and assign the density 
!  corresponding to the same (x,y) in the smaller grid, else 
!  zero.  
! 
      nphi=0.0
      do i=1,nnx 
        do m=1,nny
!
!  The mass should be zero outside of the outer radius defined on 
!  the small grid. Use squared quantities to avoid using a square root.
!
          rr2=xc(i)**2+yc(m)**2
          if ( (rr2>r2_ext) .or. (rr2<r2_int) ) then
            nphi(i,m)=0.0
          else
!
!  This is inside the domain of the small grid. Copy the 
!  density corresponding to matching (x,y) positions. 
!
            ix=i-nx/2
            iy=m+nny*ipy-nygrid/2
            j=(iy-1)/ny
            do while (iy>ny)
              iy=iy-ny
            enddo
            nphi(i,m)=phi_recv(ix,iy,j)
          endif
        enddo
      enddo
!
    endsubroutine copy_density_large_grid
!***********************************************************************
    subroutine calc_potential(nphi)
!
!  Now that we have the density field in the large grid, 
!  calculate the potential by solving the Poisson equation. 
!
!  23-apr-10/wlad: coded
!
      integer :: ikx, iky
      real, dimension (2*nx,2*ny) :: nphi,nb1
      real :: k_abs
!
!  The right-hand-side of the Poisson equation is purely real.
!
      nb1 = 0.0
!
!  Forward transform (to k-space).
!
      call fourier_transform_xy_xy_other(nphi,nb1)
!
      do iky=1,nny
        do ikx=1,nnx
          if ((kkx_fft(ikx)==0.0) .and. (kky_fft(iky+ipy*nny)==0.0)) then
            nphi(ikx,iky) = 0.0
            nb1(ikx,iky) = 0.0
          else
            if (.not.lrazor_thin) &
              call fatal_error('inverse_laplacian_expandgrid',&
                  '3d case not implemented yet')
!
!  Razor-thin approximation. Here we solve the equation
!    del2Phi=4*pi*G*Sigma(x,y)*delta(z)
!  The solution at scale k=(kx,ky) is
!    Phi(x,y,z)=-(2*pi*G/|k|)*Sigma(x,y)*exp[i*(kx*x+ky*y)-|k|*|z|]
!
            k_abs = sqrt(kkx_fft(ikx)**2+kky_fft(iky+ipy*nny)**2)
            nphi(ikx,iky) = -0.5*nphi(ikx,iky) / k_abs
            nb1(ikx,iky)  = -0.5*nb1(ikx,iky)  / k_abs
!
!  Limit |k| < kmax
!
            if ((kmax>0.0) .and. (k_abs>=kmax)) then
              nphi(ikx,iky) = 0.0
              nb1(ikx,iky) = 0.0
            endif
          endif
        enddo
      enddo
!
!  Inverse transform (to real space).
!
      call fourier_transform_xy_xy_other(nphi,nb1,linv=.true.)
!
!  The same send/recv logicals define the pairs of processors involved 
!  in this communication. This is obvious, since what is happening now 
!  is that every part of the small grid is just getting back the data 

!
    endsubroutine calc_potential
!***********************************************************************
    subroutine copy_potential_small_grid(nphi,phi,&
         lproc_comm_send,lproc_comm_recv)
!
!  The potential is calculated and stored in nphi. This is a large grid 
!  quantity, though. We need to copy the relevant part of this array to 
!  the smaller grid. It is the opposite of the operations done before. 
!
!  23-apr-10/wlad: coded
!
      use Mpicomm
!
      real, dimension (nx,ny,nz) :: phi
      real, dimension (2*nx,2*ny) :: nphi
      real, dimension (nx,ny,0:ncpus-1) :: phi_send,phi_recv
      logical, dimension(0:ncpus-1,0:ncpus-1) :: lproc_comm_send
      logical, dimension(0:ncpus-1,0:ncpus-1) :: lproc_comm_recv
      integer :: xdo,xup,mdo,mup, j
!
      intent (in) :: lproc_comm_send,lproc_comm_recv
!
!  Start by defining the limits of the small grid in the large grid 
!  axes. 
!
      xdo=nx/2+1 ; xup=3*nx/2
      if (ip<=8) print*,'inverse_laplacian_expandgrid: ', &
           'x-limits of large grid to be copied', xc(xdo),xc(xup)
!
      if (lmpicomm) then 
!
!  The same send/recv logicals define the pairs of processors involved 
!  in this communication. This is obvious, since what is happening now 
!  is that every part of the small grid is just getting back the data 
!  it sent. It was just transformed by the large grid from density into 
!  gravitational potential. 
!
        do j=0,ncpus-1
          if (iproc/=j) then
!
!  All processors in the large grid send their chunks of the potential 
!  to the corresponding processor in the small grid. 
!
            if (lproc_comm_recv(iproc,j)) then
              if (ip<=8) print*,'inverse_laplacian_expandgrid: ', &
                   'large to small grid: sending from proc ',iproc,' to ',j
              mdo=(    j*ny+1+nygrid/2)-iproc*nny
              mup=((j+1)*ny  +nygrid/2)-iproc*nny
              if (ip<=8) print*,'inverse_laplacian_expandgrid: ', &
                   'y-limits of large grid to be copied', yc(mdo),yc(mup)
              phi_send(:,:,j)=nphi(xdo:xup,mdo:mup)
              call mpisend_real(phi_send(:,:,j),(/nx,ny/),j,222)
            endif
!
!  The small grid receives all data and builds its potential array phi. 
!
            if (lproc_comm_send(iproc,j)) then
              if (ip<=8) print*,'inverse_laplacian_expandgrid: ', &
                   'large to small grid: proc ',iproc,' receiving from proc ',j
              call mpirecv_real(phi_recv(:,:,j),(/nx,ny/),j,222)
              phi(:,:,nslice)=phi_recv(:,:,j)
            endif
          else
!
!  Same processor small/large grid. Copy the potential from 
!  nphi (in the large grid) to phi (in the small grid) if (x,y) overlap
!
            if ( (yc(1)<=y(m1+1)) .and. (yc(nny)>=y(m2-1)) ) then
              mdo=(    j*ny+1+nygrid/2)-iproc*nny
              mup=((j+1)*ny  +nygrid/2)-iproc*nny
              if (ip<=8) print*,'inverse_laplacian_expandgrid: ', &
                   'y-limits of large grid to be copied', yc(mdo),yc(mup)
              phi(:,:,nslice)=nphi(xdo:xup,mdo:mup)
            endif
          endif
        enddo
      else 
!
!  Non-parallel. Just copy the central square of nphi in 
!  the large grid to phi in the small grid 
!
        mdo=ny/2+1 ; mup=3*ny/2
        phi(:,:,nslice)=nphi(xdo:xup,mdo:mup)
      endif !if mpicomm
!
!  By some reason n=1,nz needs to be filled even though in the 
!  2D case n=nz=nslice=1
!
      do n=1,nz
        phi(:,:,n)=phi(:,:,nslice)
      enddo
!
    endsubroutine copy_potential_small_grid
!***********************************************************************
    subroutine get_communication_matrix(lproc_comm_loc,&
         lproc_comm_send,lproc_comm_recv,iroot)
!
      use Mpicomm
!
!  Define the send/recv communication matrices
!
!  23-apr-10/wlad: coded
!
      logical, dimension(0:ncpus-1,0:ncpus-1) :: lproc_comm_send
      logical, dimension(0:ncpus-1,0:ncpus-1) :: lproc_comm_recv
      logical, dimension(0:ncpus-1) :: lproc_comm_loc,ltmp
      integer :: i,j,iroot
!
      intent (out) :: lproc_comm_send,lproc_comm_recv
!
      lproc_comm_send=.false.
      lproc_comm_recv=.false.
!
      if (.not.lroot) then 
        if (ip<=8) print*,'poisson:log proc ',iproc,' sending to proc',iroot
        call mpisend_logical(lproc_comm_loc,ncpus,iroot,333)
      else
        do i=0,ncpus-1
          if (i==iroot) then 
            ltmp=lproc_comm_loc
          else
            if (ip<=8) print*,'poisson:log proc ',iproc,' receiving from proc',i
            call mpirecv_logical(ltmp,ncpus,i,333)
          endif
          do j=0,ncpus-1
            if (ltmp(j)) then
              lproc_comm_send(i,j)=.true.
              lproc_comm_recv(j,i)=.true.
            endif
          enddo
        enddo
      endif
!
      call mpibcast_logical(lproc_comm_send,(/ncpus,ncpus/))
      call mpibcast_logical(lproc_comm_recv,(/ncpus,ncpus/))
!
    endsubroutine get_communication_matrix
!***********************************************************************
    subroutine read_poisson_init_pars(unit,iostat)
!
!  Read Poisson init parameters.
!
!  17-oct-2007/anders: coded
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=poisson_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=poisson_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_poisson_init_pars
!***********************************************************************
    subroutine write_poisson_init_pars(unit)
!
!  Write Poisson init parameters.
!
!  17-oct-2007/anders: coded
!
      integer, intent(in) :: unit
!
      write(unit,NML=poisson_init_pars)
!
    endsubroutine write_poisson_init_pars
!***********************************************************************
    subroutine read_poisson_run_pars(unit,iostat)
!
!  Read Poisson run parameters.
!
!  17-oct-2007/anders: coded
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=poisson_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=poisson_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_poisson_run_pars
!***********************************************************************
    subroutine write_poisson_run_pars(unit)
!
!  Write Poisson run parameters.
!
!  17-oct-2007/anders: coded
!
      integer, intent(in) :: unit
!
      write(unit,NML=poisson_run_pars)
!
    endsubroutine write_poisson_run_pars
!***********************************************************************
    subroutine inverse_laplacian_semispectral(f,phi)
!
!  Solve the Poisson equation by Fourier transforming on a periodic grid.
!
!  15-may-2006/anders+jeff: dummy
!
      use General, only: keep_compiler_quiet
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,ny,nz) :: phi
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(phi)
      call keep_compiler_quiet(f)
!
    endsubroutine inverse_laplacian_semispectral
!***********************************************************************
endmodule Poisson
