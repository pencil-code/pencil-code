! $Id$
!
! This module contains (more or less) a replica of poisson.f90, implementing the version of poisson
! equation solving that I have been building (the method outlined in Clement Baruteau's thesis).
!
! The poisson equation on a cylindrical grid is inverted, giving the potential (or the radial and
! azimuthal components of gravitational acceleration) as integrals over cylindrical coordinates. The
! integrals may be replaced by convolution products, by switching to radial coordinate u=ln(r/rmin),
! spacing radial points logarithmically (such that u is evenly spaced), and making the smoothing
! factor proportional to the radius. Given that these conditions are met, the resulting convolution
! product is computed in Fourier space and transformed to physical space at expense O(nlogn) over a
! grid with n points, compared to O(n^2) for the integral.
!
! In order to avoid forcing the periodicity of the Fourier space representation on the final result,
! a grid with twice as many gridpoints in every non-periodic dimension must be used for the Fourier
! analysis. In cylindrical coordinates, this means there must be twice as many points in the radial
! direction (and in the z direction for 3d computations) as there are holding the region of physical
! interest.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lpoisson=.true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
! COMMUNITCATED AUXILIARIES 0
!
!***************************************************************
module Poisson
!
  use Cdata
  use Cparam
  use Fourier
  use Messages
  use General, only: keep_compiler_quiet,linspace,meshgrid
!
  implicit none
!
  include 'poisson.h'
!
  real, pointer :: Gnewton_ptr
  real :: Gnewton
!
  namelist /poisson_init_pars/ Gnewton
!
  namelist /poisson_run_pars/ Gnewton
!
  logical :: luse_fourier_transform = .false.
!
! Variables to be used
!
  real :: innerRadius
  real :: outerRadius
  real :: uMax
  real :: phiMin
  real :: phiMax
  real :: du,dphi
  real, dimension(2*nx,ny) :: u2d,phi2d
  real, dimension(2*nx,ny) :: sigma
  real, dimension(2*nx,ny) :: sr,sr_fft_real,sr_fft_imag
  real, dimension(2*nx,ny) :: kr,kr_fft_real,kr_fft_imag
  real, dimension(2*nx,ny) :: sphi,sphi_fft_real,sphi_fft_imag
  real, dimension(2*nx,ny) :: kphi,kphi_fft_real,kphi_fft_imag
  real, dimension(2*nx,ny) :: gr,gphi
  real, dimension(nx) :: r1d
  real :: B,B2
  integer :: counter=1
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
      use SharedVariables, only: get_shared_variable
!
      call check_setup()
!
!  Dimensionality
!
      call decide_fourier_routine()
!
!  Mask values from elsewhere
!
      innerRadius=x(l1)
      outerRadius=x(l2)
      uMax=log(outerRadius/innerRadius)
      phiMin=y(m1)
      phiMax=y(m2)
!
!  Assign the pointer storing the shared variable for the gravitional constant
!  to another variable to dereference it.
!
      call get_shared_variable('gravitational_const',Gnewton_ptr)
      Gnewton=Gnewton_ptr
!
!  Coordinates u,phi (Fourier grid) and r,x,y (physical grid) are generated
!
      call generate_coordinates()
!
!  Kernals for each integral (gr,gphi,V)
!
      call generate_kernals()
!
    endsubroutine initialize_poisson
!***********************************************************************
    subroutine check_setup
!
!  This routine is only called from initialize_poisson in
!  poisson_logspirals.f90.
!
!  It's sole purpose is to check various initial configuration
!  parameters, to ensure that the poisson equation solver will
!  function properly.
!
!  28-sep-2016/vince: coded
!
      if(lshift_origin(2).or..not.lshift_origin_lower(2)) then
         if (lroot) then
            print*, 'initialize_poisson: logspirals only works when phi=[0,2pi)'
            print*, 'You currently have an origin-shifting configuration for the phi'
            print*, 'coordinate that, because of the way the coordinates are generated' 
            print*, 'by grid.f90, will result in an incompatible range.Please set'
            print*, 'lshift_origin(2) to "F" and lshift_origin_lower(2) to "T" in start.in'
         endif
         call fatal_error('initialize_poisson','')
      endif
      if(xyz0(2)/=0.0) then
         if (lroot) then
            print*, 'initialize_poisson: logspirals only works when phi=[0,2pi)'
            print*, 'You currently have a different range specified for the phi'
            print*, 'coordinate. Please set xyz0(2) to 0.0 in you start.in file.'
            print*, 'Additionally, it should be noted that this method is very'
            print*, 'sensitive to the periodicity of phi, and as such the maximum'
            print*, 'phi value should be specified to machine precision; setting'
            print*, 'xyz1(2) to 6.28318530717958647688 does the trick.'
         endif
         call fatal_error('initialize_poisson','')
      endif
      if (nzgrid/=1) then
        if (lroot) print*, 'initialize_poisson: logspirals only works with nzgrid==1'
        call fatal_error('initialize_poisson','')
      endif
      if (nprocx/=1) then
        if (lroot) print*, 'initialize_poisson: logspirals only works with nprocx==1'
        call fatal_error('initialize_poisson','')
      endif
    endsubroutine check_setup
!***********************************************************************
    subroutine inverse_laplacian(phi)
!
!  Dispatch solving the Poisson equation to inverse_laplacian_fft
!  or inverse_laplacian_semispectral, based on the boundary conditions
!
!  17-jul-2007/wolf: coded wrapper
!
      real, dimension(nx,ny,nz), intent(inout) :: phi
!
      if (lcylindrical_coords) then
         call inverse_laplacian_logradial_fft(phi)
      else if (lspherical_coords) then
         if (lroot) then
          print*,'There is no poisson solver for spherical '
          print*,'coordinates yet. Please feel free to implement it. '
          print*,'Many people will thank you for that.'
          call fatal_error("inverse_laplacian","")
        endif
      else if (lcartesian_coords) then
        if (lroot) then
          print*,'Use poisson.f90 or other poisson solver for cartesian coordinates '
          call fatal_error("inverse_laplacian","")
        endif
      endif
!
    endsubroutine inverse_laplacian
!***********************************************************************
    subroutine inverse_laplacian_semispectral(phi)
!
!  dummy subroutine
!      
!  19-dec-2006/anders: coded
!
      real, dimension (nx,ny,nz) :: phi
!
      call keep_compiler_quiet(phi)
!
    endsubroutine inverse_laplacian_semispectral
!***********************************************************************
    subroutine decide_fourier_routine
!
! Decide, based on the dimensionality and on the geometry
! of the grid, which fourier transform routine is to be
! used. "fourier_transform" and "fourier_tranform_shear"
! are functional only without x-parallelization, and for
! cubic, 2D square, and 1D domains only. Everything else
! should use fft_parallel instead.
!
! 05-dec-2011/wlad: coded
!
      logical :: lxpresent,lypresent,lzpresent
      logical :: l3D,l2Dxy,l2Dyz,l2Dxz
      logical :: l1Dx,l1Dy,l1Dz
!
! Check the dimensionality and store it in logicals
!
      lxpresent=(nxgrid/=1)
      lypresent=(nygrid/=1)
      lzpresent=(nzgrid/=1)
!
      l3D  =     lxpresent.and.     lypresent.and.     lzpresent
      l2Dxy=     lxpresent.and.     lypresent.and..not.lzpresent
      l2Dyz=.not.lxpresent.and.     lypresent.and.     lzpresent
      l2Dxz=     lxpresent.and..not.lypresent.and.     lzpresent
      l1Dx =     lxpresent.and..not.lypresent.and..not.lzpresent
      l1Dy =.not.lxpresent.and.     lypresent.and..not.lzpresent
      l1Dz =.not.lxpresent.and..not.lypresent.and.     lzpresent
!
      if (ldebug.and.lroot) then
        if (l3D)   print*,"This is a 3D simulation"
        if (l2Dxy) print*,"This is a 2D xy simulation"
        if (l2Dyz) print*,"This is a 2D yz simulation"
        if (l2Dxz) print*,"This is a 2D xz simulation"
        if (l1Dx)  print*,"This is a 1D x simulation"
        if (l1Dy)  print*,"This is a 1D y simulation"
        if (l1Dz)  print*,"This is a 1D z simulation"
      endif
!
! The subroutine "fourier_transform" should only be used
! for 1D, square 2D or cubic domains without x-parallelization.
! Everything else uses fft_parallel.
!
      luse_fourier_transform=(nprocx==1.and.&
           (l1dx                                       .or.&
            l1dy                                       .or.&
            l1dz                                       .or.&
            (l2dxy.and.nxgrid==nygrid)                 .or.&
            (l2dxz.and.nxgrid==nzgrid)                 .or.&
            (l2dyz.and.nygrid==nzgrid)                 .or.&
            (l3d.and.nxgrid==nygrid.and.nygrid==nzgrid)))
!
    endsubroutine decide_fourier_routine
!***********************************************************************
    subroutine read_poisson_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=poisson_init_pars, IOSTAT=iostat)
!
    endsubroutine read_poisson_init_pars
!***********************************************************************
    subroutine write_poisson_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=poisson_init_pars)
!
    endsubroutine write_poisson_init_pars
!***********************************************************************
    subroutine read_poisson_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=poisson_run_pars, IOSTAT=iostat)
!
    endsubroutine read_poisson_run_pars
!***********************************************************************
    subroutine write_poisson_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=poisson_run_pars)
!
    endsubroutine write_poisson_run_pars
!***********************************************************************
    subroutine inverse_laplacian_logradial_fft(potential)
!
      !use IO, only: output_globals
!
! Solve the Poisson equation by Fourier transforming on a periodic grid,
! in cylindrical coordinates with a log spaced radial coordinate.
!
! 18-aug-2016/vince: coded
!
      real, dimension(nx,ny,nz), intent(inout) :: potential
!
! The density will be passed in via potential, but needs to be put on the fourier grid
! with 0's in the extra cells
!
      call generate_fourier_density(potential)
!
! Mass fields for each integral (gr,gphi,V)
!
      call generate_massfields()
!
! Mass fields and kernals are used to compute quatities of interest
!
      call compute_acceleration()
!
      potential = 0.0
!
    endsubroutine inverse_laplacian_logradial_fft
!***********************************************************************
    subroutine inverse_laplacian_fft_z(phi)
!
!  15-may-2006/anders+jeff: dummy
!
      real, dimension(nx,ny,nz), intent(in) :: phi
!
      call keep_compiler_quiet(phi)
!
    endsubroutine inverse_laplacian_fft_z
!***********************************************************************
    subroutine generate_coordinates()
!
! Generates all coordinates that will be used, sets the smoothing factor (function of the 
! radial spacing). Coordinates needed for the method are constructed by reading in existing
! coordinates.
!
      real, dimension(ny) :: phi1d
      real, dimension(2*nx) :: u1d
!
      r1d=x(l1:l2)
      u1d(:nx)=log(r1d/xyz0(1))
      du=u1d(2)-u1d(1)
!
! The u coordinate will need to be extended into the Fourier cells by one unit, for the
! purposes of enforcing particular boundary conditions on the Kernal over the Fourier cells
! later. All u values beyond the first cell of the Fourier grid are inconsequential, and
! may be set to 0.0
!
      u1d(nx+1)=u1d(nx)+du
      u1d(nx+2:)=0.0
!
      phi1d=y(m1:m2)
      dphi=dy
!
      B=.01*du
      B2=B**2
!
      call meshgrid(u1d,phi1d,u2d,phi2d)
!
    endsubroutine generate_coordinates
!***********************************************************************
    subroutine generate_fourier_density(potential)
! 
! Put the density (passed in on the physical grid) onto the extended
! Fourier grid, with 0's in the extra Fourier cells.
! 
      real, dimension(nx,ny,nz), intent(in) :: potential
!
      sigma(nx+1:,:)=0.0
      sigma(:nx,:)=potential(:,:,1)
!
    endsubroutine generate_fourier_density
!***********************************************************************
    subroutine generate_kernals()
!
! Calculates the kernals for the three relevant integrals (functions
! only of the geometry of the disk) for the kernal, periodic boundary
! conditions are enforced beyond the edge of the physical disk.
!
      real, dimension(2*nx,ny) :: u_k,krNumerator,kphiNumerator,kernalDenominator
!
      u_k(:nx,:)=u2d(:nx,:)
      u_k(nx+1:,:)=-u2d(nx+1:2:-1,:)
!
      krNumerator=1+B2-exp(-u_k)*cos(phi2d)
      kphiNumerator=sin(phi2d)
      kernalDenominator=(2.*(cosh(u_k)-cos(phi2d))+B2*exp(u_k))**(1.5)
      kr=krNumerator/kernalDenominator
      kphi=kphiNumerator/kernalDenominator
!
      kr_fft_real = kr
      kr_fft_imag = 0.0
      kphi_fft_real = kphi
      kphi_fft_imag = 0.0
      call fft_xy_parallel_2D_x_extended(kr_fft_real,kr_fft_imag,.false.,.true.)
      call fft_xy_parallel_2D_x_extended(kphi_fft_real,kphi_fft_imag,.false.,.true.)
!
    endsubroutine generate_kernals
!***********************************************************************
    subroutine generate_massfields()
!
! Calculates the mass fields (density distributions weighted by
! coordinate locations) for the three relevant integrals. These are
! dependent on the mass distribution, and will therefore change at
! each timestep.
!
      sr=sigma*exp(u2d/2)
      sphi=sigma*exp(3*u2d/2)
!
      sr_fft_real = sr
      sr_fft_imag = 0.0
      sphi_fft_real = sphi
      sphi_fft_imag = 0.0
      call fft_xy_parallel_2D_x_extended(sr_fft_real,sr_fft_imag,.false.,.true.)
      call fft_xy_parallel_2D_x_extended(sphi_fft_real,sphi_fft_imag,.false.,.true.)
!
    endsubroutine generate_massfields
!***********************************************************************
    subroutine compute_acceleration()
!
! Uses kernals and mass fields for the two acceleration integrals to
! calculate gravitational accelerations.
!
      real, dimension(2*nx,ny) :: gr_convolution,gr_factor,gphi_convolution,gphi_factor
      real, dimension(2*nx,ny) :: gr_conv_real,gr_conv_imag
      real, dimension(2*nx,ny) :: gphi_conv_real,gphi_conv_imag
      real, parameter          :: normalization_factor=1./(2.*nxgrid*nygrid)
!
      gr_conv_real = kr_fft_real*sr_fft_real - kr_fft_imag*sr_fft_imag
      gr_conv_imag = kr_fft_real*sr_fft_imag + kr_fft_imag*sr_fft_real
      gphi_conv_real = kphi_fft_real*sphi_fft_real - kphi_fft_imag*sphi_fft_imag
      gphi_conv_imag = kphi_fft_real*sphi_fft_imag + kphi_fft_imag*sphi_fft_real
      call fft_xy_parallel_2D_x_extended(gr_conv_real,gr_conv_imag,.true.,.true.)
      call fft_xy_parallel_2D_x_extended(gphi_conv_real,gphi_conv_imag,.true.,.true.)
      gr_convolution   = normalization_factor*gr_conv_real
      gphi_convolution = normalization_factor*gphi_conv_real
!
      gr_factor=-Gnewton*exp(-u2d/2)*du*dphi
      gphi_factor=-Gnewton*exp(-3*u2d/2)*du*dphi
!
      gr=gr_factor*gr_convolution + Gnewton*sigma*du*dphi/B
      gphi=gphi_factor*gphi_convolution
!
    endsubroutine compute_acceleration
!***********************************************************************
    subroutine fft_xy_parallel_2D_x_extended(a_re,a_im,linv,lneed_im,shift_y)
!
!  Subroutine to do FFT of distributed 2D data in the x- and y-direction,
!  on a grid which is twice the size of the physical grid in the x-diction.
!  For x- and/or y-parallelization the calculation will be done under
!  MPI in parallel on all processors of the corresponding xy-plane.
!  nx is restricted to be an integer multiple of nprocy.
!  ny is restricted to be an integer multiple of nprocx.
!  linv indicates forward (=false, default) or backward (=true) transform.
!  You can set lneed_im=false if the imaginary data can be disregarded.
!  Attention: input data will be overwritten.
!
!  09-feb-2017/vince: adapted from fft_xy_parallel_2D
!
      use Mpicomm, only: transp,transp_other,remap_to_pencil_xy, transp_pencil_xy, unmap_from_pencil_xy
!
      complex, dimension (2*nxgrid) :: ax
      complex, dimension (nygrid) :: ay
      complex, dimension (nzgrid) :: az
      real, dimension (4*(2*nxgrid)+15) :: wsavex
      real, dimension (4*nygrid+15) :: wsavey
      real, dimension (4*nzgrid+15) :: wsavez
!
      real, dimension (2*nx,ny), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im
      real, dimension (2*nxgrid), optional :: shift_y
!
      integer, parameter :: pnx=2*nxgrid, pny=nygrid/nprocxy ! pencil shaped data sizes
      integer, parameter :: tnx=nygrid, tny=2*nxgrid/nprocxy ! pencil shaped transposed data sizes
      real, dimension (:,:), allocatable :: p_re, p_im   ! data in pencil shape
      real, dimension (:,:), allocatable :: t_re, t_im   ! data in transposed pencil shape
      real, dimension (tny) :: deltay_x
      real, dimension (tny) :: dshift_y
!
      integer :: l, m, stat, x_offset
      integer(8) :: nxgrid_fourier,nx_fourier
      logical :: lforward, lcompute_im, lshift
!
      nxgrid_fourier = 2*nxgrid
      nx_fourier = 2*nx
!
      lforward = .true.
      if (present (linv)) lforward = .not. linv
!
      lcompute_im = .true.
      if (present (lneed_im)) lcompute_im = lneed_im
!
      lshift = .false.
      if (present (shift_y)) lshift = .true.
!
!
      ! Check for degenerate cases.
      if (nxgrid == 1) then
        if (lshift) then
          call fft_y_parallel (a_re(1,:), a_im(1,:), .not. lforward, lcompute_im, shift_y(1))
        else
          call fft_y_parallel (a_re(1,:), a_im(1,:), .not. lforward, lcompute_im)
        endif
        return
      endif
      if (nygrid == 1) then
        call fft_x_parallel (a_re(:,1), a_im(:,1), .not. lforward, lcompute_im)
        return
      endif
!
      if (mod (nxgrid, nprocxy) /= 0) &
          call fatal_error ('fft_xy_parallel_2D', 'nxgrid needs to be an integer multiple of nprocx*nprocy', lfirst_proc_xy)
      if (mod (nygrid, nprocxy) /= 0) &
          call fatal_error ('fft_xy_parallel_2D', 'nygrid needs to be an integer multiple of nprocx*nprocy', lfirst_proc_xy)
!
      ! Allocate memory for large arrays.
      allocate (p_re(pnx,pny), p_im(pnx,pny), t_re(tnx,tny), t_im(tnx,tny), stat=stat)
      if (stat > 0) call fatal_error ('fft_xy_parallel_2D', 'Could not allocate memory for p/t', .true.)
!
      if (lshear) then
        x_offset = 1 + (ipx+ipy*nprocx)*tny
        deltay_x = -deltay * (xgrid(x_offset:x_offset+tny-1) - (x0+Lx/2))/Lx
      endif
!
      if (lshift) then
        x_offset = 1 + (ipx+ipy*nprocx)*tny
        dshift_y = shift_y(x_offset:x_offset+tny-1)
      endif
!
      call cffti (nxgrid_fourier, wsavex)
      call cffti (nygrid, wsavey)
!
      if (lforward) then
!
!  Forward FFT:
!
        ! Remap the data we need into pencil shape.
        call remap_to_pencil_xy (a_re, p_re)
        if (lcompute_im) then
          call remap_to_pencil_xy (a_im, p_im)
        else
          p_im = 0.0
        endif
!
        call transp_pencil_xy (p_re, t_re)
        call transp_pencil_xy (p_im, t_im)
!
        ! Transform y-direction.
        do l = 1, tny
          ay = cmplx (t_re(:,l), t_im(:,l))
          call cfftf(nygrid, ay, wsavey)
          if (lshear) ay = ay * exp (cmplx (0, ky_fft * deltay_x(l)))
          if (lshift) ay = ay * exp (cmplx (0,-ky_fft * dshift_y(l)))
          t_re(:,l) = real (ay)
          t_im(:,l) = aimag (ay)
        enddo
!
        call transp_pencil_xy (t_re, p_re)
        call transp_pencil_xy (t_im, p_im)
!
        ! Transform x-direction.
        do m = 1, pny
          ax = cmplx (p_re(:,m), p_im(:,m))
          call cfftf (nxgrid_fourier, ax, wsavex)
          p_re(:,m) = real (ax)
          p_im(:,m) = aimag (ax)
        enddo
!
        ! Unmap the results back to normal shape.
        call unmap_from_pencil_xy (p_re, a_re)
        call unmap_from_pencil_xy (p_im, a_im)
!
        ! Apply normalization factor to fourier coefficients.
        ! Since I am using this internally only for fftconvolve,
        ! which has odd normalization requirements, I will
        ! not normalize here, but rather in fftconvolve

!
      else
!
!  Inverse FFT:
!
        ! Remap the data we need into transposed pencil shape.
        call remap_to_pencil_xy (a_re, p_re)
        call remap_to_pencil_xy (a_im, p_im)
!
        do m = 1, pny
          ! Transform x-direction back.
          ax = cmplx (p_re(:,m), p_im(:,m))
          call cfftb (nxgrid_fourier, ax, wsavex)
          p_re(:,m) = real (ax)
          p_im(:,m) = aimag (ax)
        enddo
!
        call transp_pencil_xy (p_re, t_re)
        call transp_pencil_xy (p_im, t_im)
!
        do l = 1, tny
          ! Transform y-direction back.
          ay = cmplx (t_re(:,l), t_im(:,l))
          if (lshear) ay = ay * exp (cmplx (0, -ky_fft*deltay_x(l)))
          call cfftb (nygrid, ay, wsavey)
          t_re(:,l) = real (ay)
          if (lcompute_im) t_im(:,l) = aimag (ay)
        enddo
!
        call transp_pencil_xy (t_re, p_re)
        if (lcompute_im) call transp_pencil_xy (t_im, p_im)
!
        ! Unmap the results back to normal shape.
        call unmap_from_pencil_xy (p_re, a_re)
        if (lcompute_im) call unmap_from_pencil_xy (p_im, a_im)
!
      endif
!
      deallocate (p_re, p_im, t_re, t_im)
!
    endsubroutine fft_xy_parallel_2D_x_extended
!***********************************************************************
    subroutine get_acceleration(acceleration)
!
      real, dimension(nx,ny,nz,3), intent(out) :: acceleration           !should I (CAN I?) make this allocatable?
      integer :: n
!
      do n=1,nz
        acceleration(:,:,n,1) = gr(1:nx,1:ny)
        acceleration(:,:,n,2) = gphi(1:nx,1:ny)
      enddo
      acceleration(:,:,:,3)=0.
!
    endsubroutine get_acceleration
!***********************************************************************
endmodule Poisson
