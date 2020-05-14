! $Id$
!
! This module contains (more or less) a replica of poisson.f90, implementing the version of poisson
! equation solving that I have been building (the method outlined in Clement Baruteau's thesis).
!
! The derivatives of the poisson equation on a cylindrical grid is inverted, giving the radial and
! azimuthal components of gravitational acceleration as integrals over cylindrical coordinates. The
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
  use Mpicomm, only: mpisend_real,mpirecv_real
!
  implicit none
!
  include 'poisson.h'
!
  real, pointer :: Gnewton_ptr
  real :: Gnewton
  real :: Bsmooth=0.01
!
  namelist /poisson_init_pars/ Gnewton,Bsmooth
!
  namelist /poisson_run_pars/ Gnewton,Bsmooth
!
  logical :: luse_fourier_transform = .false.
!
  real                             :: innerradius,du,dphi,B,B2
  real, dimension(2*nx,ny)         :: u2d,phi2d
  real, dimension(2*nxgrid,nygrid) :: u2d_global,phi2d_global
!
  real, dimension(2*nx,ny) :: sigma
  real, dimension(2*nx,ny) :: sr,sr_fft_real,sr_fft_imag
  real, dimension(2*nx,ny) :: kr,kr_fft_real,kr_fft_imag
  real, dimension(2*nx,ny) :: sphi,sphi_fft_real,sphi_fft_imag
  real, dimension(2*nx,ny) :: kphi,kphi_fft_real,kphi_fft_imag
  real, dimension(nx,ny)   :: gr,gphi
!
  integer :: ipx_fourier,ipy_fourier,iproc_fourier
  integer :: ixlower_fourier,ixupper_fourier
  integer :: iylower_fourier,iyupper_fourier
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
      innerradius=xyz0(1)
!
!  Fourier processor grid coordinates
!
      call physical_to_fourier_procs(ipx,ipy,ipx_fourier,ipy_fourier)
      call physical_to_fourier_proc(iproc,iproc_fourier)
      ixlower_fourier = (2*nx*ipx_fourier + 1)
      ixupper_fourier = 2*nx*(ipx_fourier + 1)
      iylower_fourier = (ny*ipy_fourier + 1)
      iyupper_fourier = ny*(ipy_fourier + 1)
!
!  Assigning the pointer storing the shared variable for the gravitional constant
!  to another variable dereferences it in fortran
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
      !print*, 'leaving initialize_poisson'
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
      if (nprocx /= 1 .and. mod(nprocx,2) /= 0) then
        if (lroot) print*, 'initialize_poisson: logspirals only works with even nprocx'
        call fatal_error('initialize_poisson','')
      endif
    endsubroutine check_setup
!***********************************************************************
    subroutine physical_to_fourier_procs(ipx_physical,ipy_physical,ipx_fourier,ipy_fourier)
!
!  This routine is only called from initialize_poisson in
!  poisson_logspirals.f90.
!
!  It converts from the physical to the fourier proc grid
!  (in a way detailed in Baruteau thesis...add some ascii
!   are represntation)
!
!  22-may-2017/vince: coded
!
      integer, intent(in)  :: ipx_physical,ipy_physical
      integer, intent(out) :: ipx_fourier ,ipy_fourier
!
      if(mod(ipx_physical,2) == 0) then
         ipx_fourier = ipx_physical/2
      else
         ipx_fourier = (nprocx + ipx_physical - 1)/2
      endif
      ipy_fourier = ipy_physical   !same number of processors along x, twice as many points each
!
    endsubroutine physical_to_fourier_procs
!***********************************************************************
    subroutine fourier_to_physical_procs(ipx_fourier,ipy_fourier,ipx_physical,ipy_physical)
!
!  This routine is only called from initialize_poisson in
!  poisson_logspirals.f90.
!
!  It converts from the physical to the fourier proc grid
!  (in a way detailed in Baruteau thesis...add some ascii
!   are represntation)
!
!  22-may-2017/vince: coded
!
      integer, intent(in)  :: ipx_fourier ,ipy_fourier
      integer, intent(out) :: ipx_physical,ipy_physical
!
      if(ipx_fourier*2 < nprocx) then
         ipx_physical = ipx_fourier*2
      else
         ipx_physical = ipx_fourier*2 - nprocx + 1
      endif
      ipy_physical = ipy_fourier   !same number of processors along x, twice as many points each
!
    endsubroutine fourier_to_physical_procs
!***********************************************************************
    subroutine physical_to_fourier_proc(iproc_physical,iproc_fourier)
!
!  This routine is only called from initialize_poisson in
!  poisson_logspirals.f90.
!
!  It converts from the physical to the fourier proc grid
!  (in a way detailed in Baruteau thesis...add some ascii
!   are represntation)
!
!  22-may-2017/vince: coded
!
      integer, intent(in)  :: iproc_physical
      integer, intent(out) :: iproc_fourier
      integer              :: tmpx_physical,tmpy_physical
      integer              :: tmpx_fourier,tmpy_fourier
!
!  iproc = iprocx + nprocx*iprocy
!  iprocx < nprocx (zero indexed)
!
      tmpx_physical = mod(iproc_physical,nprocx)
      tmpy_physical = (iproc_physical - tmpx_physical)/nprocx
!
      call physical_to_fourier_procs(tmpx_physical,tmpy_physical,tmpx_fourier,tmpy_fourier)
!
      iproc_fourier = tmpx_fourier + nprocx*tmpy_fourier
!
    endsubroutine physical_to_fourier_proc
!***********************************************************************
    subroutine fourier_to_physical_proc(iproc_fourier,iproc_physical)
!
!  This routine is only called from initialize_poisson in
!  poisson_logspirals.f90.
!
!  It converts from the physical to the fourier proc grid
!  (in a way detailed in Baruteau thesis...add some ascii
!   are represntation)
!
!  22-may-2017/vince: coded
!
      integer, intent(in)  :: iproc_fourier
      integer, intent(out) :: iproc_physical
      integer              :: tmpx_fourier,tmpy_fourier
      integer              :: tmpx_physical,tmpy_physical
!
!  iproc = iprocx + nprocx*iprocy
!  iprocx < nprocx (zero indexed)
!
      tmpx_fourier = mod(iproc_fourier,nprocx)
      tmpy_fourier = (iproc_fourier - tmpx_fourier)/nprocx
!
      call fourier_to_physical_procs(tmpx_fourier,tmpy_fourier,tmpx_physical,tmpy_physical)
!
      iproc_physical = tmpx_physical + nprocx*tmpy_physical
!
    endsubroutine fourier_to_physical_proc
!***********************************************************************
    subroutine generate_coordinates()
!
!  Generates all coordinates that will be used, sets the smoothing factor (function of the 
!  radial spacing). Coordinates needed for the method are constructed by reading in existing
!  coordinates.
!
!  Because the order in which the processors divide the domain of the fourier grid is different
!  from that in which they divide the domain of the physical grid, care must be taken to make
!  sure that any quantities which will eventually be in the fourier transform contain the portion
!  of the domain that corresponds to their FOURIER grid position, NOT their physical grid position.
!
!  For this reason, the global coordinates are reconstructed, and the appropriate domain
!  decomposition is sliced out of the global coordinates based on the fourier processor grid
!  position (given by ixlower_fourier, ixupper_fourier, etc)
!
      real, dimension(ny)       :: phi1d
      real, dimension(2*nx)     :: u1d
      integer                   :: partner,ipx_loop,ipy_loop
      integer                   :: ixlower,ixupper,iylower,iyupper
      integer, parameter        :: RADIAL_COORDINATE_TAG=1285
      integer, parameter        :: AZIMUTHAL_COORDINATE_TAG=1286
      real, dimension(nxgrid)   :: r1d_global
      real, dimension(nygrid)   :: phi1d_global
      real, dimension(2*nxgrid) :: u1d_global
!
      real, dimension(nx) :: send_buffer_radial,recv_buffer_radial
      real, dimension(ny) :: send_buffer_azimuthal,recv_buffer_azimuthal
!
!  Must exchange local coordinates with all procs. Can be done along x and y separately (to avoid communication
!  of redundant information)
!
      do ipx_loop = 0, nprocx-1
         partner = ipx_loop + nprocx*ipy
         ixlower = (nx*ipx_loop + 1)
         ixupper = nx*(ipx_loop + 1)
         if (iproc == partner) then
            ! data is local
            r1d_global(ixlower:ixupper) = x(l1:l2)
         else
            ! communicate with partner
            send_buffer_radial = x(l1:l2)
            if (iproc > partner) then ! above diagonal: send first, receive then
               call mpisend_real(send_buffer_radial,nx,partner, &
                    RADIAL_COORDINATE_TAG)
               call mpirecv_real(recv_buffer_radial,nx,partner, &
                    RADIAL_COORDINATE_TAG)
               
            else                      ! below diagonal: receive first, send then
               call mpirecv_real(recv_buffer_radial,nx,partner, &
                    RADIAL_COORDINATE_TAG)
               call mpisend_real(send_buffer_radial,nx,partner, &
                    RADIAL_COORDINATE_TAG)
            endif
            r1d_global(ixlower:ixupper) = recv_buffer_radial
         endif
      enddo
!
      do ipy_loop = 0, nprocy-1
        partner = ipx + nprocx*ipy_loop
        iylower = (ny*ipy_loop + 1)
        iyupper = ny*(ipy_loop + 1)
        if (iproc == partner) then
          ! data is local
          phi1d_global(iylower:iyupper) = y(m1:m2)
        else
          ! communicate with partner
          send_buffer_azimuthal = y(m1:m2)
          if (iproc > partner) then ! above diagonal: send first, receive then
            call mpisend_real(send_buffer_azimuthal,ny,partner, &
                 AZIMUTHAL_COORDINATE_TAG)
            call mpirecv_real(recv_buffer_azimuthal,ny,partner, &
                 AZIMUTHAL_COORDINATE_TAG)
          else                      ! below diagonal: receive first, send then
            call mpirecv_real(recv_buffer_azimuthal,ny,partner, &
                 AZIMUTHAL_COORDINATE_TAG)
            call mpisend_real(send_buffer_azimuthal,ny,partner, &
                 AZIMUTHAL_COORDINATE_TAG)
          endif
          phi1d_global(iylower:iyupper) = recv_buffer_azimuthal
        endif
      enddo
!
!  1D global coordinates, and various and sundry coordinate related quantities
!
      u1d_global(:nxgrid) = log(r1d_global/innerradius)
      du                  = u1d_global(2)-u1d_global(1)
      dphi = dy
!
      B  = Bsmooth*du
      B2 = B**2
!
!  For the definition of the kernal given by Baruteau, the u coordinate must be defined one cell beyond the boundary
!  of the physical grid. Its values after that are inconsequatial, and are zeroed
!
      u1d_global(nxgrid+1)  = u1d_global(nxgrid) + du
      u1d_global(nxgrid+2:) = 0.0
!
!  2D global coordinates
!
      call meshgrid(u1d_global,phi1d_global,u2d_global,phi2d_global)
!
!  Slicing processor sized chunks from global coordinates, matching the domain decomposition on the FOURIER grid
!
      u2d   =   u2d_global(ixlower_fourier:ixupper_fourier,iylower_fourier:iyupper_fourier)
      phi2d = phi2d_global(ixlower_fourier:ixupper_fourier,iylower_fourier:iyupper_fourier)
!
    endsubroutine generate_coordinates
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
          print*,'There is no logspiral poisson solver for spherical '
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
! Mass fields for each integral
!
      call generate_massfields()
!
! Mass fields and kernals are used to compute quatities of interest
!
      call compute_acceleration()
!
      potential = 0.0 ! Should I pick a more ridiculous number? Something that will be
                      ! more noticeable if it gets anywhere it shouldn't be?
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
    subroutine inverse_laplacian_z_2nd_neumann(f)
!
!  15-may-2006/anders+jeff: dummy
!
      real, dimension(:,:,:,:), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine inverse_laplacian_z_2nd_neumann
!***********************************************************************
    subroutine generate_fourier_density(potential)
! 
! Put the density (passed in on the physical grid) onto the extended
! Fourier grid, with 0's in the extra Fourier cells.
!
      real, dimension(nx,ny,nz), intent(in) :: potential
!
      real, dimension(nx,ny)              :: send_recv_buffer
      integer, dimension(2)               :: send_recv_size
      integer                             :: target_proc
      integer, parameter                  :: DENSITY_TAG=1828
!
      if(nprocx == 1) then
!
         sigma(:nx,:)   = potential(:,:,1)
         sigma(nx+1:,:) = 0.0
!
      else
!
         send_recv_size(1) = nx
         send_recv_size(2) = ny
!
         if(ipx_fourier < nprocx/2) then
!
!  If this is true, then this processor will handle the portion of the Fourier grid with density consisting of
!  it's own part of the density on the physical grid plus the density from the next radial processor out.
!
            sigma(:nx,:) = potential(:,:,1)
            target_proc  = ipx + 1 + ipy*nprocx
            call mpirecv_real(send_recv_buffer,send_recv_size,target_proc, &
                 DENSITY_TAG)
            sigma(nx+1:,:) = send_recv_buffer
!
         else
!
!  This processor will handle the portion of the Fourier grid with 0's in the density. It must send it's portion
!  of the density to the processor one radial processor interior to it, on the physical grid (see above)
!
            sigma = 0.0
            send_recv_buffer = potential(:,:,1)
            target_proc = ipx - 1 + ipy*nprocx
            call mpisend_real(send_recv_buffer,send_recv_size,target_proc, &
                 DENSITY_TAG)
!
         end if
!
      end if
!
    endsubroutine generate_fourier_density
!***********************************************************************
    subroutine generate_kernals()
!
!  Calculates the kernals for the three relevant integrals. These are functions
!  only of the coordinate geometry of the disk, and will therefore be constant.
!
      real, dimension(2*nxgrid,nygrid) :: uk_global
      real, dimension(2*nx,ny)         :: uk,krNumerator,kphiNumerator,kernalDenominator
!
      uk_global(:nxgrid,:)   =  u2d_global(:nxgrid,:)
      uk_global(nxgrid+1:,:) = -u2d_global(nxgrid+1:2:-1,:)
!
      uk = uk_global(ixlower_fourier:ixupper_fourier,iylower_fourier:iyupper_fourier)
!
      krNumerator       = 1. + B2 - exp(-uk)*cos(phi2d)
      kphiNumerator     = sin(phi2d)
      kernalDenominator = (2.*(cosh(uk) - cos(phi2d)) + B2*exp(uk))**(1.5)
!
      kr   = krNumerator/kernalDenominator
      kphi = kphiNumerator/kernalDenominator
!
      kr_fft_real = kr
      kr_fft_imag = 0.0
      kphi_fft_real = kphi
      kphi_fft_imag = 0.0
      call fft_fouriergrid(kr_fft_real,kr_fft_imag,.false.,.true.)
      call fft_fouriergrid(kphi_fft_real,kphi_fft_imag,.false.,.true.)
!
    endsubroutine generate_kernals
!***********************************************************************
    subroutine generate_massfields()
!
!  Calculates the mass fields (density distributions weighted exponentially by
!  coordinate locations) for the three relevant integrals. These are dependent
!  on the distribution of mass, and will therefore change at each timestep.
!
      sr   = sigma*exp(u2d/2)
      sphi = sigma*exp(3*u2d/2)
!
      sr_fft_real = sr
      sr_fft_imag = 0.0
      sphi_fft_real = sphi
      sphi_fft_imag = 0.0
      call fft_fouriergrid(sr_fft_real,sr_fft_imag,.false.,.true.)
      call fft_fouriergrid(sphi_fft_real,sphi_fft_imag,.false.,.true.)
!
    endsubroutine generate_massfields
!***********************************************************************
    subroutine compute_acceleration()
!
! Uses kernals and mass fields for the two acceleration integrals to
! calculate gravitational accelerations.
!
      real, dimension(2*nx,ny) :: gr_convolution,gphi_convolution
      real, dimension(2*nx,ny) :: gr_conv_real,gr_conv_imag
      real, dimension(2*nx,ny) :: gphi_conv_real,gphi_conv_imag
      real, dimension(2*nx,ny) :: gr_fourier,gphi_fourier,gr_factor,gphi_factor
      real, parameter          :: normalization_factor=1./(2.*nxgrid*nygrid)
!
      real, dimension(nx,ny) :: send_recv_buffer_radial,send_recv_buffer_azimuthal
      integer, dimension(2)  :: send_recv_size
      integer                :: target_proc
      integer, parameter     :: RADIAL_ACC_TAG=1283
      integer, parameter     :: AZIMUTHAL_ACC_TAG=1284
!
      send_recv_size(1) = nx
      send_recv_size(2) = ny
!
      gr_conv_real = kr_fft_real*sr_fft_real - kr_fft_imag*sr_fft_imag
      gr_conv_imag = kr_fft_real*sr_fft_imag + kr_fft_imag*sr_fft_real
      gphi_conv_real = kphi_fft_real*sphi_fft_real - kphi_fft_imag*sphi_fft_imag
      gphi_conv_imag = kphi_fft_real*sphi_fft_imag + kphi_fft_imag*sphi_fft_real
      call fft_fouriergrid(gr_conv_real,gr_conv_imag,.true.,.true.)
      call fft_fouriergrid(gphi_conv_real,gphi_conv_imag,.true.,.true.)
      gr_convolution   = normalization_factor*gr_conv_real
      gphi_convolution = normalization_factor*gphi_conv_real
!
      gr_factor   = -Gnewton*exp(-u2d/2)*du*dphi
      gphi_factor = -Gnewton*exp(-3*u2d/2)*du*dphi
!
      gr_fourier   = gr_factor*gr_convolution + Gnewton*sigma*du*dphi/B
      gphi_fourier = gphi_factor*gphi_convolution
!
      if(nprocx == 1) then
         gr   = gr_fourier(:nx,:)
         gphi = gphi_fourier(:nx,:)
      else
         if(ipx_fourier < nprocx/2) then
!
!  This processor previously received the density from the processor immediately radially exterior
!  to it on the physical grid, and must therefore communicate its results back to that processor.
!
            gr                         = gr_fourier(:nx,:)
            send_recv_buffer_radial    = gr_fourier(nx+1:,:)
            gphi                       = gphi_fourier(:nx,:)
            send_recv_buffer_azimuthal = gphi_fourier(nx+1:,:)
            target_proc = ipx + 1 + ipy*nprocx
            call mpisend_real(send_recv_buffer_radial   ,send_recv_size,target_proc, &
                 RADIAL_ACC_TAG   )
            call mpisend_real(send_recv_buffer_azimuthal,send_recv_size,target_proc, &
                 AZIMUTHAL_ACC_TAG)
         else
!
!  This processor sent its density to another, and must now collect the resulting acceleration on
!  its portion of the physical grid domain.
!
            target_proc = ipx - 1 + ipy*nprocx
            call mpirecv_real(send_recv_buffer_radial   ,send_recv_size,target_proc, &
                 RADIAL_ACC_TAG   )
            call mpirecv_real(send_recv_buffer_azimuthal,send_recv_size,target_proc, &
                 AZIMUTHAL_ACC_TAG)
            gr   = send_recv_buffer_radial
            gphi = send_recv_buffer_azimuthal
         end if
      end if
!
    endsubroutine compute_acceleration
!***********************************************************************
    subroutine fft_fouriergrid(a_re,a_im,linv,lneed_im,shift_y)
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
      complex, dimension (2*nxgrid)     :: ax
      complex, dimension (nygrid)       :: ay
      complex, dimension (nzgrid)       :: az
      real, dimension (4*(2*nxgrid)+15) :: wsavex
      real, dimension (4*nygrid+15)     :: wsavey
      real, dimension (4*nzgrid+15)     :: wsavez
!
      real, dimension (2*nx,ny), intent(inout) :: a_re, a_im
!
      logical, optional, intent(in)        :: linv, lneed_im
      real, dimension (2*nxgrid), optional :: shift_y
!
      integer, parameter                 :: pnx=2*nxgrid, pny=nygrid/nprocxy ! pencil shaped data sizes
      integer, parameter                 :: tnx=nygrid, tny=2*nxgrid/nprocxy ! pencil shaped transposed data sizes
      real, dimension (:,:), allocatable :: p_re, p_im   ! data in pencil shape
      real, dimension (:,:), allocatable :: t_re, t_im   ! data in transposed pencil shape
      real, dimension (tny)              :: deltay_x
      real, dimension (tny)              :: dshift_y
!
      integer :: l, m, stat, x_offset
      integer :: nxgrid_fourier,nx_fourier
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
      if (mod (nxgrid_fourier, nprocxy) /= 0) &
          call fatal_error ('fft_xy_parallel_2D', 'nxgrid_fourier needs to be an integer multiple of nprocx*nprocy', lfirst_proc_xy)
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
        call remap_to_pencil_fouriergrid (a_re, p_re)
        if (lcompute_im) then
          call remap_to_pencil_fouriergrid (a_im, p_im)
        else
          p_im = 0.0
        endif
!
        call transp_pencil_fouriergrid (p_re, t_re)
        call transp_pencil_fouriergrid (p_im, t_im)
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
        call transp_pencil_fouriergrid (t_re, p_re)
        call transp_pencil_fouriergrid (t_im, p_im)
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
        call unmap_from_pencil_fouriergrid (p_re, a_re)
        call unmap_from_pencil_fouriergrid (p_im, a_im)
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
        call remap_to_pencil_fouriergrid (a_re, p_re)
        call remap_to_pencil_fouriergrid (a_im, p_im)
!
        do m = 1, pny
          ! Transform x-direction back.
          ax = cmplx (p_re(:,m), p_im(:,m))
          call cfftb (nxgrid_fourier, ax, wsavex)
          p_re(:,m) = real (ax)
          p_im(:,m) = aimag (ax)
        enddo
!
        call transp_pencil_fouriergrid (p_re, t_re)
        call transp_pencil_fouriergrid (p_im, t_im)
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
        call transp_pencil_fouriergrid (t_re, p_re)
        if (lcompute_im) call transp_pencil_fouriergrid (t_im, p_im)
!
        ! Unmap the results back to normal shape.
        call unmap_from_pencil_fouriergrid (p_re, a_re)
        if (lcompute_im) call unmap_from_pencil_fouriergrid (p_im, a_im)
!
      endif
!
      deallocate (p_re, p_im, t_re, t_im)
!
    endsubroutine fft_fouriergrid
!***********************************************************************
    subroutine remap_to_pencil_fouriergrid(in, out)
!
!  Remaps data distributed on several processors into pencil shape.
!  This routine remaps 2D arrays in x and y only for nprocx>1.
!
!   24-apr-2017/vince: adapted from remap_to_pencil_xy_2D
!
      real, dimension(:,:), intent(in)  :: in
      real, dimension(:,:), intent(out) :: out
!
      integer, parameter    :: inx=2*nx, iny=ny
      integer, parameter    :: onx=2*nxgrid, ony=ny/nprocx
      integer, parameter    :: bnx=2*nx, bny=ny/nprocx ! transfer box sizes
      integer               :: ibox,partner,partner_fourier, alloc_err
      integer, dimension(2) :: send_recv_size
      integer, parameter    :: ltag=104, utag=105
!
      real, dimension(:,:), allocatable   :: send_buf, recv_buf
!
      if (nprocx == 1) then
         out = in
         return
      endif
!
      if (mod (ny, nprocx) /= 0) &
          call stop_fatal ('remap_to_pencil_xy_2D: ny needs to be an integer multiple of nprocx', lfirst_proc_xy)
!
      if ((size (in, 1) /= inx) .or. ((size (in, 2) /= iny))) &
          call stop_fatal ('remap_to_pencil_xy_2D: input array size mismatch /= 2*nx,ny', lfirst_proc_xy)
      if ((size (out, 1) /= onx) .or. ((size (out, 2) /= ony))) &
          call stop_fatal ('remap_to_pencil_xy_2D: output array size mismatch /= 2*nxgrid,ny/nprocx', lfirst_proc_xy)
!
      allocate (send_buf(bnx,bny), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('remap_to_pencil_xy_2D: not enough memory for send_buf!', .true.)
      allocate (recv_buf(bnx,bny), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('remap_to_pencil_xy_2D: not enough memory for recv_buf!', .true.)
      send_recv_size(1) = bnx
      send_recv_size(2) = bny
!
!  These loop traversals need to happen on the FOURIER grid, meaning I will need to translate
!  fourier processor grid positions into physical processor grid positions
!
      do ibox = 0, nprocx-1
        partner_fourier = ipy_fourier*nprocx + ibox
        call fourier_to_physical_proc(partner_fourier,partner)
        if (iproc_fourier == partner_fourier) then
          ! data is local
          out(bnx*ibox+1:bnx*(ibox+1),:) = in(:,bny*ibox+1:bny*(ibox+1))
        else
          ! communicate with partner
          send_buf = in(:,bny*ibox+1:bny*(ibox+1))
          if (iproc_fourier > partner_fourier) then ! above diagonal: send first, receive then
            call mpisend_real(send_buf,send_recv_size,partner,utag)
            call mpirecv_real(recv_buf,send_recv_size,partner,ltag)
          else                      ! below diagonal: receive first, send then
            call mpirecv_real(recv_buf,send_recv_size,partner,utag)
            call mpisend_real(send_buf,send_recv_size,partner,ltag)
          endif
          out(bnx*ibox+1:bnx*(ibox+1),:) = recv_buf
        endif
      enddo
!
      deallocate (send_buf, recv_buf)
!
    endsubroutine remap_to_pencil_fouriergrid
!***********************************************************************
    subroutine unmap_from_pencil_fouriergrid(in, out)
!
!  Unmaps pencil shaped 2D data distributed on several processors back to normal shape.
!  This routine is the inverse of the remap function for nprocx>1.
!
!   4-jul-2010/Bourdin.KIS: coded
!
      real, dimension(:,:), intent(in)  :: in
      real, dimension(:,:), intent(out) :: out
!
      integer, parameter    :: inx=2*nxgrid, iny=ny/nprocx
      integer, parameter    :: onx=2*nx, ony=ny
      integer, parameter    :: bnx=2*nx, bny=ny/nprocx ! transfer box sizes
      integer               :: ibox, partner,partner_fourier, alloc_err
      integer, dimension(2) :: send_recv_size
      integer, parameter    :: ltag=106, utag=107
!
      real, dimension(:,:), allocatable   :: send_buf, recv_buf
!
      if (nprocx == 1) then
        out = in
        return
      endif
!
      if (mod (ny, nprocx) /= 0) &
          call stop_fatal ('unmap_from_pencil_xy_2D: ny needs to be an integer multiple of nprocx', lfirst_proc_xy)
!
      if ((size (in, 1) /= inx) .or. ((size (in, 2) /= iny))) &
          call stop_fatal ('unmap_from_pencil_xy_2D: input array size mismatch /= 2*nxgrid,ny/nprocx', lfirst_proc_xy)
      if ((size (out, 1) /= onx) .or. ((size (out, 2) /= ony))) &
          call stop_fatal ('unmap_from_pencil_xy_2D: output array size mismatch /= 2*nx,ny', lfirst_proc_xy)
!
      allocate (send_buf(bnx,bny), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('unmap_from_pencil_xy_2D: not enough memory for send_buf!', .true.)
      allocate (recv_buf(bnx,bny), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('unmap_from_pencil_xy_2D: not enough memory for recv_buf!', .true.)
      send_recv_size(1) = bnx
      send_recv_size(2) = bny
!
      do ibox = 0, nprocx-1
        partner_fourier = ipy_fourier*nprocx + ibox
        call fourier_to_physical_proc(partner_fourier,partner)
        if (iproc_fourier == partner_fourier) then
          ! data is local
          out(:,bny*ibox+1:bny*(ibox+1)) = in(bnx*ibox+1:bnx*(ibox+1),:)
        else
          ! communicate with partner
          send_buf = in(bnx*ibox+1:bnx*(ibox+1),:)
          if (iproc_fourier > partner_fourier) then ! above diagonal: send first, receive then
            call mpisend_real(send_buf,send_recv_size,partner,utag)
            call mpirecv_real(recv_buf,send_recv_size,partner,ltag)
          else                      ! below diagonal: receive first, send then
            call mpirecv_real(recv_buf,send_recv_size,partner,utag)
            call mpisend_real(send_buf,send_recv_size,partner,ltag)
          endif
          out(:,bny*ibox+1:bny*(ibox+1)) = recv_buf
        endif
      enddo
!
      deallocate (send_buf, recv_buf)
!
    endsubroutine unmap_from_pencil_fouriergrid
!***********************************************************************
    subroutine transp_pencil_fouriergrid(in, out)
!
!  Transpose 2D data distributed on several processors.
!  This routine transposes arrays in x and y only.
!  The data must be mapped in pencil shape, especially for nprocx>1.
!
!   4-jul-2010/Bourdin.KIS: coded, adapted parts of transp_xy
!  21-jun-2013/Bourdin.KIS: reworked, parellized MPI communication
!
      use General, only: count_bits
!
      real, dimension(:,:), intent(in)  :: in
      real, dimension(:,:), intent(out) :: out
!
      integer               :: inx, iny, onx, ony ! sizes of in and out arrays
      integer               :: bnx, bny! destination box sizes and number of elements
      integer, dimension(2) :: send_recv_size
      integer               :: num_it ! number of iterations for box-data transfer
      integer               :: it, ibox, partner,partner_fourier, alloc_err
      integer, parameter    :: ltag=108, utag=109
!
      real, dimension(:,:), allocatable   :: send_buf, recv_buf
!
!
      inx = size (in, 1)
      iny = size (in, 2)
      onx = size (out, 1)
      ony = size (out, 2)
      bnx = onx/nprocxy
      bny = ony
!
      if (mod (onx, nprocxy) /= 0) &
          call stop_fatal ('transp_pencil_xy_2D: onx needs to be an integer multiple of nprocxy', lfirst_proc_xy)
!
      if ((inx /= bny*nprocxy) .or. (iny /= bnx)) &
          call stop_fatal ('transp_pencil_xy_2D: input array has unmatching size', lfirst_proc_xy)
      if ((onx /= bnx*nprocxy) .or. (ony /= bny)) &
          call stop_fatal ('transp_pencil_xy_2D: output array has unmatching size', lfirst_proc_xy)
!
      allocate (send_buf(bnx,bny), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('transp_pencil_xy_2D: not enough memory for send_buf!', .true.)
      allocate (recv_buf(bnx,bny), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('transp_pencil_xy_2D: not enough memory for recv_buf!', .true.)
      send_recv_size(1) = bnx
      send_recv_size(2) = bny
!
      num_it = 2**count_bits (nprocxy-1)
      do it = 0, num_it-1
        ibox = IEOR (iproc_fourier, it)
        partner_fourier = ibox
        call fourier_to_physical_proc(partner_fourier,partner)
        if (ibox >= nprocxy) cycle
        if (iproc_fourier == partner_fourier) then
          ! data is local
          out(bnx*ibox+1:bnx*(ibox+1),:) = transpose (in(bny*ibox+1:bny*(ibox+1),:))
        else
          ! communicate with partner
          send_buf = transpose (in(bny*ibox+1:bny*(ibox+1),:))
          if (iproc_fourier > partner_fourier) then ! above diagonal: send first, receive then
            call mpisend_real(send_buf,send_recv_size,partner,utag)
            call mpirecv_real(recv_buf,send_recv_size,partner,ltag)
          else                      ! below diagonal: receive first, send then
            call mpirecv_real(recv_buf,send_recv_size,partner,utag)
            call mpisend_real(send_buf,send_recv_size,partner,ltag)
          endif
          out(bnx*ibox+1:bnx*(ibox+1),:) = recv_buf
        endif
      enddo
!
      deallocate (send_buf, recv_buf)
!
    endsubroutine transp_pencil_fouriergrid
!***********************************************************************
    subroutine stop_fatal(msg,force)
!
!  Print message and stop. If force, stop without shutting down MPI.
!
!  13-dez-10/Bourdin.KIS: coded
!
      use Mpicomm, only: die_immediately,die_gracefully

      character (len=*) :: msg
      logical, optional :: force
!
      logical :: immediately
!
      immediately = .false.
      if (present (force)) immediately = force
!
      if (lroot .or. immediately) write(*,'(A,A)') 'STOPPED FATAL: ', msg
!
      if (immediately) then
        call die_immediately
      else
        call die_gracefully
      endif
!
    endsubroutine stop_fatal
!***********************************************************************
    subroutine get_acceleration(acceleration)
!
!  Passes the accelerations out of this module for use in selfgravity_logspirals
!
      real, dimension(nx,ny,nz,3), intent(out) :: acceleration
!
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
