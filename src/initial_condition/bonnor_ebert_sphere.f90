! $Id$
!
! Initial condition for a Bonnor-Ebert sphere
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!***************************************************************
!
!
!!!!!!!
!! Eva 27/05/2022: testing the ic generation
!!!!!!!

module InitialCondition
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages

!
  implicit none
!
  include '../initial_condition.h'
!
! Bonnor Ebert sphere initial condition
! xmin and xmax are the min and max values of the dimensionles xi parameter
! rhoc is the central density of the core in code units, 
! cs the sound speed in code units
! nsteps is the desired resolution in the radial direction
!
  double precision :: xmin=1.e-5, xmax=10.d0, rhoc=10.,cs=1.
  integer          :: nsteps = 100
  
  namelist /initial_condition_pars/ &
       xmin,xmaxi,nsteps,rhoc,cs
!
  contains
!***********************************************************************
    subroutine register_initial_condition()
!
!  Register variables associated with this module; likely none.
!
!  26-jan-12/joern: coded
!
      if (lroot) call svn_id( &
           "$Id")
!
    endsubroutine register_initial_condition
!***********************************************************************
    subroutine initialize_initial_condition(f)
!
!  Initialize any module variables which are parameter dependent.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine initial_condition_uu(f)
!
!  Initialize the velocity field.
!
!  Velocity is zero everywhere 
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      integer :: iy, iz
!
      do iy=m1,m2;do iz=n1,n2
         f(:,iy,iz,iuu:iuu+2)=0.d0
      enddo; enddo
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uu
!***********************************************************************
    subroutine initial_condition_lnrho(f)
!
!  Initialize logarithmic density. init_lnrho 
!  will take care of converting it to linear 
!  density if you use ldensity_nolog
!
      use SharedVariables
      use Sub, only: interp1
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f      
      real :: lnrho_r, rgrid
! 
      integer :: ix,iy,iz
    
! This call will return an array of rho as a function of radial distance r (in code units) 
! The xmax and xmin give the extent of the Bonnor Ebert sphere in dimensionless units

      call bonnor_ebert_sphere(xmin,xmax,nsteps,rr,rho)
! eva: not sure about xyz limits and where to get the xyz arrays
      do ix=1,mx
        do iy=1,my
           do ix=1,mx
             rgrid = sqrt(x(ix)**2+y(iy)**2+z(iz)**2)
             if (r_mn.le.xmax*sqrt(4.*pi*rhoc)/cs) then
               ! linear interpolation of rho to rgrid position
               lnrho_r=interp1(rr,rho,nsteps,rgrid)
             else
               lnrho_r=rho(nsteps)
             endif
             f(l,m,n,ilnrho)=lnrho_r
           enddo
        enddo
      enddo

    call keep_compiler_quiet(f)

    endsubroutine initial_condition_lnrho
!********************************************************************
    subroutine initial_condition_aa(f)
!
!  Initialize magnetic vector potential
!
!  eva: not sure what to put here 
!

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f      
!
      call keep_compiler_quiet(f)
!      
    endsubroutine initial_condition_aa
!***********************************************************************
    subroutine initial_condition_ss(f)
!
!  Initialize entropy.
!
!  eva: not sure what to put here: External pressure should match rho*cs**2. at
!  the edge of the BE sphere 
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f      
!      
    endsubroutine initial_condition_ss
!***********************************************************************
    subroutine read_initial_condition_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=initial_condition_pars, IOSTAT=iostat)
!
    endsubroutine read_initial_condition_pars
!***********************************************************************
    subroutine write_initial_condition_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=initial_condition_pars)
!
    endsubroutine write_initial_condition_pars
!***********************************************************************
    subroutine bonnor_ebert_sphere(xmin,xmax,nsteps,rr,rho)

    !============================================================================
    ! This routine integrates the self-similar equations for the density profile
    ! of an isothermal sphere, contained in Bonnor (1956)
    ! The integration will be done using a fourth-order runge-kutta.
    !==========================================================================

    implicit none

!   rho here is actually rho/rho_c = exp(-ys(1))
!   it should be rescaled to whatever physical units we are using

    real, allocatable, intent(out) :: rr(:), rho(:)

    integer             :: i 
    integer, intent(in) :: nsteps
    double precision, intent(in) :: xmax, xmin
    integer, parameter  :: nfunct=2                     !Number of functions to integrate
    double precision    :: y0, dx, xs, mass
    double precision    :: ys(nfunct), dys(nfunct)

    ! Initializations
    if (.not.allocated(rho)) allocate(rho (1:nsteps))
    if (.not.allocated(rr)) allocate(rr (1:nsteps))
    dx = ((xmax-xmin)/real(nsteps))

    xs     = xmin

    ! The initial values of ys for given x0 and alpha:
    call initial_condition(xs,ys)

    ! Actual integration
    do i=1,nsteps
      ! Here we need to turn ksi (xs) into r and potential (ys(1)) into
      ! rho/rho_c: rho/rho_c = exp(-ys(1))
      rr(i)   = sqrt(4.*pi*rhoc)*xs/cs
      rho(i) = rho_c*exp(-ys(1))
      call lane_emden(xs,ys,dys)      ! derivatives     
      call rk4(ys, dys, nfunct, xs, dx, ys)
      xs   = xs + dx
    enddo

    end subroutine bonnor_ebert_sphere 
!***********************************************************************
 subroutine rk4(y, dydx, n, x, h, yout)

 implicit none

 ! Given the values for the variables y(1:n) and their derivatives
 ! dydx(1:n) known at x, use the fourth-order Runge-Kutta method
 ! to advance the solution over an interval h and return the incremented 
 ! variables as yout(1:n), which need not be a distinct array from y.
 ! The user supplies the subroutine derivs(x,y,dydx), which returns the
 ! derivatives dydx at x. 

 integer, parameter :: nmax = 10  ! maximum number of functions 
 integer            :: n
 integer            :: i
 double precision   :: h, x, dydx(n), y(n), yout(n)

 double precision   :: h6, hh, xh, dym(nmax), dyt(nmax), yt(nmax)

 hh = h * 0.5
 h6 = h/6.
 xh = x+hh

 do i=1,n
   yt(i) = y(i) + hh*dydx(i) !First step
 enddo

 call lane_emden(xh,yt,dyt)      !Second step
 do i=1,n
   yt(i) = y(i) + hh*dyt(i)
 enddo

 call lane_emden(xh,yt,dym)      !Third step
 do i=1,n
   yt(i)  = y(i)+h*dym(i)
   dym(i) = dyt(i)+dym(i)
 enddo

 call lane_emden(x+h,yt,dyt)        !Fourth step
 do i=1,n                    !Accumulate increments with proper weights
   yout(i) = y(i) + h6*(dydx(i) + dyt(i) + 2.*dym(i))
 enddo

 return

 end subroutine rk4

!==========================================================================
 subroutine lane_emden(xx,yy,dyydxx)

 implicit none
 ! This routine is needed by the Runge-Kutta to provide the
 ! right-hand side of the equations.  

 integer, parameter :: nf=2
 double precision   :: yy(nf), xx, dyydxx(nf)

 ! y(1) = phi(ksi)
 ! y(2) = dphi(ksi)/d ksi

 dyydxx(1) = yy(2)                     ! derivative of the potential -> will
give the mass parameter
 dyydxx(2) = exp(-yy(1)) - 2*yy(2)/xx   ! potential -> will give the density

 end subroutine lane_emden
!==========================================================================
 subroutine initial_condition(x,y)

 implicit none

 ! Calculates the ksi->0 boundary condition for the Bonnor-Ebert sphere 

 double precision :: x, y(2)

 y(1) = 0.0
 y(2) = 0.0

 end subroutine initial_condition
!==========================================================================

!
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from noinitial_condition.f90 for any    **
!**  InitialCondition routines not implemented in this file        **
!**                                                                **
    include '../initial_condition_dummies.inc'
!********************************************************************
endmodule InitialCondition
