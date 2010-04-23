! $Id: poisson.f90 12460 2009-12-10 15:19:51Z sven.bingert $
!
!  This module solves the Poisson equation
!    (d^2/dx^2 + d^2/dy^2 + d^2/dz^2 - h) f = RHS(x,y,z)
!  [which for h/=0 could also be called inhomogenous nonuniform Helmholtz
!  equation] for the function f(x,y,z).
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
  include 'poisson.h'
!
  namelist /poisson_init_pars/ &
      lsemispectral, kmax, lrazor_thin, lklimit_shear, lexpand_grid
!
  namelist /poisson_run_pars/ &
      lsemispectral, kmax, lrazor_thin, lklimit_shear, lexpand_grid
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
!  Limit the wavenumber to the maximum circular region that is always available
!  in k-space. The radial wavenumber kx changes with time due to shear as
!
!    kx = kx0+qshear*Omega*t*ky
!
!  Considering the available (kx,ky) space, it turns slowly from a square to a
!  parallellogram (the hole for kx,ky<1 is ignored here):
!
!       - - - -                  - - - -
!      |       |               /       /
!      |       |    ->       /       /
!      |       |           /       /
!      |       |         /       /
!       - - - -          - - - -
!     
!  To make the gravity force isotropic at small scales one can limit k to
!  the largest circular region that is present in both the square and the
!  parallellogram. The circle has radius kmax=kNy/sqrt(2). See Gammie (2001).
!
      if (lklimit_shear) kmax = kx_ny/sqrt(2.0)
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
      use Sub, only: keep_compiler_quiet
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
!  grid to avoid the periodicity of the Fourier transform
!
!  Exclude shear and no-razor thin
!
!  25-may-2008/wlad: coded
!
      use Mpicomm
!
      real, dimension (nx,ny,nz) :: phi
      real, dimension (2*nx,2*ny) :: nphi,nb1
!
      real, dimension (nx,ny,0:ncpus-1) :: phi_send,phi_recv
!
      real, dimension(2*nx)     :: xc,kkx_fft
      real, dimension(2*ny)     :: yc
      real, dimension(2*nygrid) :: kky_fft
!      
      logical, dimension(0:ncpus-1) :: lproc_comm_loc,ltmp
      logical, dimension(0:ncpus-1,0:ncpus-1) :: lproc_comm_send
      logical, dimension(0:ncpus-1,0:ncpus-1) :: lproc_comm_recv
!
      real :: rr,k_abs
      integer :: ikx, iky
!
      real    :: x0,xn,y0,yn,dxc,dyc,dxc1,dyc1,Lxn,Lyn
!
      integer :: i,j,m,mm,iproc_send,iy_serial
      integer :: nnx,nny,nnghost
      integer :: nnxgrid,nnygrid
!
      integer :: xdo,xup,iroot,mdo,mup
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
      if (nprocx.gt.1) &
           call fatal_error("inverse_laplacian_expandgrid",&
           "Not yet implemented for nprocx > 1")
!
      if (lroot.and.lfirstcall) print*,'Entered inverse_laplacian_expandgrid'
!
      iroot=0
      nnghost=1
!
!  Define the expanded Cartesian axes
! 
      nnx=2*nx               ; nny=2*ny
      xn=2*x(l2)+0.5*dx      ; x0=-xn
      yn=xn                  ; y0=-yn
      nnxgrid=2*nxgrid       ; nnygrid=2*nygrid
!
      do i=1,nnx
        xc(i)=1.0*(i-1)        /(nnxgrid-1)*(xn-x0)+x0
      enddo
      do m=1,nny
        yc(m)=1.0*(m-1+ipy*nny)/(nnygrid-1)*(yn-y0)+y0
      enddo      
! 
      dxc=xc(2)-xc(1)  ; dyc=dxc
      dxc1=1/dxc       ; dyc1=1/dyc
!
      Lxn=2*xc(nnx)    ; Lyn=Lxn
!
      kkx_fft = &
          cshift((/(i-(nnxgrid+1)/2,i=0,nnxgrid-1)/),+(nnxgrid+1)/2)*2*pi/Lxn
      kky_fft = &
          cshift((/(i-(nnygrid+1)/2,i=0,nnygrid-1)/),+(nnygrid+1)/2)*2*pi/Lyn
!
!  Prepare the sending/receiving between different processors
!
      nnghost=1
      if (lmpicomm) then 
        lproc_comm_loc=.false.
        lproc_comm_send=.false.
        lproc_comm_recv=.false.
        do m=1,ny
          mm=m+m1-1
          iy_serial=nint((y(mm)-y0)*dyc1)+1
          iproc_send=(iy_serial-1)/nny
          if (iproc_send/=iproc) then
            if (.not.lproc_comm_loc(iproc_send)) then
              lproc_comm_loc(iproc_send)=.true.
            endif
            phi_send(:,m,iproc_send)=phi(:,m,nnghost)
          endif
        enddo
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
!  Send to those who need, get from those who send
!
        do j=0,ncpus-1
          if (lproc_comm_send(iproc,j)) then
            if (ip<=8) print*,'poisson:proc ',iproc,' sending to proc',j
            call mpisend_real(phi_send(:,:,j),(/nx,ny/),j,111)
          endif
          if (lproc_comm_recv(iproc,j)) then
            if (ip<=8) print*,'poisson:proc ',iproc,' receiving from proc ',j
            call mpirecv_real(phi_recv(:,:,j),(/nx,ny/),j,111)
          endif
          if (j==iproc) phi_recv(:,:,iproc)=phi(:,:,nnghost)     
        enddo
      else
        phi_recv(:,:,0)=phi(:,:,nnghost)
      endif
!
      nphi=0.0
      do i=1,nnx
        do m=1,nny
!  Zero mass outside of the domain
          rr=sqrt(xc(i)**2+yc(m)**2)
          if ( (rr>r_ext) .or. (rr<r_int) ) then
            nphi(i,m)=0.0
          else
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
!  Copy back to small grid
!
      xdo=nx/2+1 ; xup=3*nx/2
      if (ip<=8) print*,'x-limits of large grid to be copied', &
           xc(xdo),xc(xup)
!
      if (lmpicomm) then 
!
        do j=0,ncpus-1
          if (iproc/=j) then
            !this processor received info, now it has to send back
            if (lproc_comm_recv(iproc,j)) then
              if (ip<=8) print*,'copy grid: sending from proc ',iproc,' to ',j
              mdo=(    j*ny+1+nygrid/2)-iproc*nny
              mup=((j+1)*ny  +nygrid/2)-iproc*nny
              if (ip<=8) print*,'y-limits of large grid to be copied', &
                   yc(mdo),yc(mup)
              phi_send(:,:,j)=nphi(xdo:xup,mdo:mup)
              call mpisend_real(phi_send(:,:,j),(/nx,ny/),j,222)
            endif
            !this one sent, now it receives
            if (lproc_comm_send(iproc,j)) then
              if (ip<=8) print*,'copy grid: proc ',iproc, &
                   ' receiving from proc ',j
              call mpirecv_real(phi_recv(:,:,j),(/nx,ny/),j,222)
              phi(:,:,nnghost)=phi_recv(:,:,j)
            endif
          else
            !only if the big processor encompasses the little
            if ( (yc(1)<=y(m1+1)) .and. (yc(nny)>=y(m2-1)) ) then
              !plus minus 1 just to avoid rounding errors
              mdo=(    j*ny+1+nygrid/2)-iproc*nny
              mup=((j+1)*ny  +nygrid/2)-iproc*nny
              if (ip<=8) print*,'mdo,mup',iproc,j,mdo,mup
              phi(:,:,nnghost)=nphi(xdo:xup,mdo:mup)
            endif
          endif
        enddo
      else
        mdo=ny/2+1 ; mup=3*ny/2
        phi(:,:,nnghost)=nphi(xdo:xup,mdo:mup)
      endif
!
      do n=1,nz
        phi(:,:,n)=phi(:,:,nnghost)
      enddo
!
      if (lfirstcall) lfirstcall=.false.
!
    endsubroutine inverse_laplacian_expandgrid
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
    endsubroutine read_Poisson_run_pars
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
      use Sub, only: keep_compiler_quiet
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
