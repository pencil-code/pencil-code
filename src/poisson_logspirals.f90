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
! MAUX CONTRIBUTION 3
! COMMUNITCATED AUXILIARIES 3
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
  real, dimension(2*nx,ny) :: u2d,phi2d,sigma,sr,sphi,sv,kr,kphi,kv,gr,gphi
  real, dimension(nx) :: r1d
  real :: B,B2
  integer :: counter=1
!
contains
!***********************************************************************
    subroutine register_poisson()
!
!  Register auxilliary farray variables for poisson_logspirals
!
!  12-dec-16/vince: adapted
!
      use FArrayManager
!
!  Set indices for auxiliary variables
!
      call farray_register_auxiliary('gpotselfx',igpotselfx,communicated=.true.)
      call farray_register_auxiliary('gpotselfy',igpotselfy,communicated=.true.)
      call farray_register_auxiliary('gpotselfz',igpotselfz,communicated=.true.)
!
!  Identify version number (generated automatically by SVN).
!
      if (lroot) call svn_id( &
          "$Id$")
!
    endsubroutine register_selfgravity
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
      if (nzgrid/=1) then
        if (lroot) print*, 'initialize_poisson: logspirals only works with nzgrid==1'
        call fatal_error('initialize_poisson','')
      endif
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
!  Initialize farray contributions
!
      f(:,:,:,igpotselfx:igpotselfz)=0.0
!
!  Boundary condition consistency for gpotself
!
if (ldensity) then
        i = merge(irho, ilnrho, ldensity_nolog)
        if (any(bcx(igpotselfx:igpotselfz)=='p') .and. .not.(bcx(i)=='p')) then
          if (lroot) then
            print*, 'initialize_selfgravity: gpotself has bcx=''p'', but the density is not'
            print*, '                        periodic! (you must set a proper boundary condition'
            print*, '                        for the potential)'
            print*, 'initialize_selfgravity: bcx=', bcx
          endif
          call fatal_error('initialize_selfgravity','')
        endif
        if (any(bcy(igpotselfx:igpotselfz)=='p') .and. .not.(bcy(i)=='p')) then
          if (lroot) then
            print*, 'initialize_selfgravity: gpotself has bcy=''p'', but the density is not'
            print*, '                        periodic! (you must set a proper boundary condition'
            print*, '                        for the potential)'
            print*, 'initialize_selfgravity: bcy=', bcy
          endif
          call fatal_error('initialize_selfgravity','')
        endif
        if (any(bcz(igpotselfx:igpotselfz)=='p') .and. .not.(bcz(i)=='p')) then
          if (lroot) then
            print*, 'initialize_selfgravity: gpotself has bcz=''p'', but the density is not'
            print*, '                        periodic! (you must set a proper boundary condition'
            print*, '                        for the potential)'
            print*, 'initialize_selfgravity: bcz=', bcz
          endif
          call fatal_error('initialize_selfgravity','')
        endif
      endif
!
!     *** potential problems here...ipotself may still be 0, depending on the order in which the
!         register / vs initialize routines are called. Will need to check ***
!
      if (any(bcx(igpotselfx:igpotselfz)=='p') .and. .not.(bcx(ipotself)=='p')) then
        if (lroot) then
          print*, 'initialize_particles_selfgrav: igpotself has bcx=''p'', but the potential is not'
          print*, '                               periodic! (you must set a proper boundary'
          print*, '                               condition for the gradient of the potential)'
          print*, 'initialize_particles_selfgrav: bcx=', bcx
        endif
        call fatal_error('initialize_particles_selfgrav','')
      endif
      if (any(bcy(igpotselfx:igpotselfz)=='p') .and. .not.(bcy(ipotself)=='p')) then
        if (lroot) then
          print*, 'initialize_particles_selfgrav: igpotself has bcy=''p'', but the potential is not'
          print*, '                               periodic! (you must set a proper boundary'
          print*, '                               condition for the gradient of the potential)'
          print*, 'initialize_particles_selfgrav: bcy=', bcy
        endif
        call fatal_error('initialize_particles_selfgrav','')
      endif
      if (any(bcz(igpotselfx:igpotselfz)=='p') .and. .not.(bcz(ipotself)=='p')) then
        if (lroot) then
          print*, 'initialize_particles_selfgrav: igpotself has bcz=''p'', but the potential is not'
          print*, '                               periodic! (you must set a proper boundary'
          print*, '                               condition for the gradient of the potential)'
          print*, 'initialize_particles_selfgrav: bcz=', bcz
        endif
        call fatal_error('initialize_particles_selfgrav','')
      endif
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
!
        if (present(gpotself)) then
          call inverse_laplacian_logradial_fft(phi)
        else
          call fatal_error("inverse_laplacian","poisson_logspirals works with the acceleration only")
        endif
!
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
      call compute_potential(potential)
!      
    endsubroutine inverse_laplacian_logradial_fft
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
      u1d(:nx)=log(r1d/innerRadius)
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
      kv=1/sqrt((1+B2)*exp(u_k)+exp(-u_k)-2*cos(phi2d)) !changed from **-.5 to 1/sqrt
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
      sv = sigma*exp(3*u2d/2)
!
    endsubroutine generate_massfields
!***********************************************************************
    subroutine compute_acceleration()
!
! Uses kernals and mass fields for the two acceleration integrals to
! calculate gravitational accelerations.
!
      real, dimension(2*nx,ny) :: gr_convolution,gr_factor,gr1,gr2,gphi_convolution,gphi_factor
!
      call fftconvolve(kr,sr,gr_convolution)
      call fftconvolve(kphi,sphi,gphi_convolution)
!
      gr_factor=-Gnewton*exp(-u2d/2)*du*dphi
      gphi_factor=-Gnewton*exp(-3*u2d/2)*du*dphi
!
      gr1=gr_factor*gr_convolution
      gr2=Gnewton*sigma*du*dphi/B
      gr=gr1+gr2
!
      gphi=gphi_factor*gphi_convolution
!
      do n=1,nz
        f(l1:l2,m1:m2,n,igpotselfx)=gr(:nx,:)
        f(l1:l2,m1:m2,n,igpotselfy)=gphi(:nx,:)
      enddo
      f(l1:l2,m1:m2,n1:n2,igpotselfz)=0.
!
    endsubroutine compute_acceleration
!***********************************************************************
    subroutine compute_potential(potential)
!
! Uses kernals and mass fields for the potential integral to calculate
! potential on the grid NOTE: Should actually calculate potential from
! the accelerations, since the log radial convolution product produces
! potentials which are known to violate Newton's first law, but produces accelerations
! that are free of this behavior (See Towards Predictive Scenarios of
! Planetary Migration, by Clement Baruteau)
!
      real, dimension(nx,ny,nz), intent(inout) :: potential
      real, dimension(2*nx,ny) :: v_convolution,v_factor
!
      call fftconvolve(kv,sv,v_convolution)
      v_factor=-Gnewton*innerRadius*exp(-u2d/2)*du*dphi
!
      do n=1,nz
        potential(:,:,n)=v_factor(:nx,:)*v_convolution(:nx,:)
      enddo
!
    endsubroutine compute_potential
!***********************************************************************
    subroutine fftconvolve(array1,array2,convolution)
!
!  convolution product using the pendil code fft methods to do the ffts
!
      real, intent(in), dimension(:,:) :: array1
      real, intent(in), dimension(:,:) :: array2
      real, intent(out), dimension(:,:) :: convolution
      complex, dimension(size(convolution,1),size(convolution,2)) :: convolution_fourier
      real, dimension(size(convolution,1),size(convolution,2)) :: array1_fourier_real,array1_fourier_imaginary
      real, dimension(size(convolution,1),size(convolution,2)) :: array2_fourier_real,array2_fourier_imaginary
      real, dimension(size(convolution,1),size(convolution,2)) :: convolution_fourier_real
      real, dimension(size(convolution,1),size(convolution,2)) :: convolution_fourier_imaginary
      integer :: nrow,ncol
      logical :: array1_convolution_diffshape,array2_convolution_diffshape
!
      array1_convolution_diffshape=(size(array1,1)/=size(convolution,1)).or.(size(array1,2)/=size(convolution,2))
      array2_convolution_diffshape=(size(array2,1)/=size(convolution,1)).or.(size(array2,2)/=size(convolution,2))
      if (array1_convolution_diffshape.or.array2_convolution_diffshape) then
         print*,'fftconvolve: input and output arrays must be the same shape as each other.'
         call fatal_error('fftconvolve','')
         stop
      else
!
!  The input arrays are purely real. Here, I copy the arrays into the ones holding the
!  real components for the transform, and fill the arrays holding the imaginary
!  components with 0's
!
         nrow=size(convolution,1)
         ncol=size(convolution,2)
         array1_fourier_real=array1
         array1_fourier_imaginary=0.0
         array2_fourier_real=array2
         array2_fourier_imaginary=0.0
      end if
!
!  Transforming the arrays
!  NOTE: the pencil code normalized Fourier transforms on the way forward, so each
!        array will pick up a factor of 1/(nrow*ncol) here. This will be important later
!
      call fourier_transform_other(array1_fourier_real,array1_fourier_imaginary,.false.)
      call fourier_transform_other(array2_fourier_real,array2_fourier_imaginary,.false.)
!
!  Now, the convolution product must be calculated by multiplying the transformed arrays.
!  They must be multiplied as complex numbers; to save space and function calls, the real
!  and imaginary parts of the product are calculated separately from the real and imaginary
!  parts of the two arrays.
!  NOTE: the convolution product will have a factor of 1/(nrow*ncol)^2, since each array will
!        have picked up a separate normalization factor from the forward transform. This is
!        not desirable: as the convolution product is the result of forward transforming and
!        subserquently inverse transforming, the final product should have a single normalization
!        factor. To this end, a factor of (nrow*ncol) is added to the components of the product
!        in Fourier space, such that the Fourier product has a single normalization factor.
!
      convolution_fourier_real=(nrow*ncol)* &
          (array1_fourier_real*array2_fourier_real-array1_fourier_imaginary*array2_fourier_imaginary)
      convolution_fourier_imaginary=(nrow*ncol)* &
          (array1_fourier_real*array2_fourier_imaginary+array1_fourier_imaginary*array2_fourier_real)
!
!  Finally, the inverse transform moves the convolution product to ordinary space.
!
      call fourier_transform_other(convolution_fourier_real,convolution_fourier_imaginary,.true.)
!
      convolution=convolution_fourier_real
!
    endsubroutine fftconvolve
!***********************************************************************
endmodule Poisson
