! $Id$
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
  real :: eta0=0.0, k_eta, x0_drop, y0_drop
  real :: Omega_SB=1.0,gamma_parameter=1.0
  real :: v_jet_peak,sigma_jet,r_jet_cent
  real :: Omega_vortex=1.
! Vortex parameters: 
  real :: xv0,yv0, bv,rm,vm,rout  
  integer :: nagrid=10000
!
  character (len=labellen), dimension(ninit) :: init_shallow_density='nothing'
  character (len=labellen), dimension(ninit) :: init_shallow_hydro='nothing'
!
  namelist /initial_condition_pars/  eta0, k_eta, x0_drop, y0_drop, Omega_SB, &
       init_shallow_density,init_shallow_hydro,gamma_parameter,v_jet_peak, &
       sigma_jet,r_jet_cent,Omega_vortex,xv0,yv0,bv,rm,vm,rout,nagrid
!
  real :: planetary_radius=impossible ! derives from gamma_parameter and Omega
!
  contains
!***********************************************************************
    subroutine register_initial_condition()
!
!  Register variables associated with this module; likely none.
!
!  07-may-09/wlad: coded
!
      if (lroot) call svn_id( &
           "$Id$")
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
!  gamma_parameter = Omega/a**2 where a is the planetary radius
!
      planetary_radius = sqrt(Omega_SB/gamma_parameter)
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine initial_condition_lnrho(f)
!
!  Initial condition given by 
!
!     h = eta + Lb
!     rho is g*eta
!
      use General, only: notanumber
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (nx) :: eta,r2
      integer :: j,i

      real, dimension(nx) :: rr,xx,ww
      real :: xs,ys,yy
      real :: b,rm,Ro,Bu
      real :: phi

!
      do j=1,ninit
!
        select case (init_shallow_density(j))
!
        case('linear-zero')
          f(l1:l2,m1:m2,n1:n2,ilnrho)=-impossible 
!
        case('solid-body')   
!
          do n=n1,n2
            do m=m1,m2
!
! Equilibrium -> gh = gh0 + 3/2*Omega**2 * r**2 - Omega**2/(4*a**2) r**4
!
               r2 = x(l1:l2)**2 + y(m)**2
               eta = eta0 + Omega_SB**2*r2 * (1.5 - 0.25*gamma_parameter/Omega_SB * r2)
               f(l1:l2,m,n,ilnrho) = log(eta)
            enddo
          enddo
!
        case('gaussian-blob')   
!
!  with eta=a*exp(-k*(x-x0)**2) 
!
          do n=n1,n2
            do m=m1,m2
              eta = eta0 * exp(-k_eta * ( &  
                   (x(l1:l2)-x0_drop)**2 + (y(m)-y0_drop)**2 & 
                   ))
              f(l1:l2,m,n,ilnrho) = log(eta)
            enddo
          enddo
!
        case('Li-Ingersoll')
          do n=n1,n2
             do m=m1,m2
               xx = x(l1:l2)-xv0
               yy = y(  m  )-yv0
               rr = sqrt(xx**2 + yy**2)
!
               do i=l1,l2
                 call calc_phi(rr(i-l1+1),phi)
                 if (ldensity_linearstart) then
                    f(i,m,n,irho) = eta0 - phi
                 else
                    call fatal_error("initial_condition_lnrho",&
                         "switch ldensity_linearstart=T in density_init_pars")
                 endif
               enddo
!               
             enddo
          enddo
!
       endselect
!
      enddo
!
    endsubroutine initial_condition_lnrho
!***********************************************************************
    subroutine calc_phi(r,intgr)
!
      real :: Lr,dr,vr,tmp,rr
      integer :: ir,alloc_err
      real, intent(in) :: r
      real, allocatable, dimension(:) :: a
      real, intent(out) :: intgr
!
      if (nagrid >= 1) then
        allocate(a(nagrid),stat=alloc_err)
        if (alloc_err > 0) call fatal_error('calc_phi','Could not allocate memory',.true.)
      else
        if (lroot) print*,"nagrid=",nagrid
        call fatal_error("calc_phi",&
              "give nagrid a positive integer greater than 1")
      endif
!
      Lr = rout - r
      dr = Lr/nagrid
      do ir=1,nagrid
        rr = r + (ir-1)*dr
        vr = vm * (rr/rm) * exp(1/bv * (1-(rr/rm)**bv))
        tmp = vr**2/rr + 2*Omega_SB*vr 
        a(ir) = tmp*dr
      enddo
!      
      intgr = 0.5*(a(1)+a(nagrid-1)) + sum(a(2:nagrid-1))
!      
      deallocate(a)
!
    endsubroutine calc_phi
!***********************************************************************    
    subroutine initial_condition_uu(f)
!
!  Initial condition given by 
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(nx) :: rr,xx,vr
      integer :: j
      real :: yy
!
      do j=1,ninit

        select case (init_shallow_hydro(j))
!
        case('solid-body')
          do n=n1,n2
            do m=m1,m2
              f(l1:l2,m,n,iux) = -Omega_SB * y(  m  )
              f(l1:l2,m,n,iuy) =  Omega_SB * x(l1:l2)
            enddo
          enddo
!
! Add Gaussian jet at the polar region (Gaussian in radius)
! v_phi = v_jet_peak * exp( - (rr(i)/ sigma_jet)**2),
! where v_jet_peak is the peak azimuthal velocity, rr is the distance
! from the jet center to the grid point, and sigma_jet is the jet width.
!
        case('zonal-jet')
          do n=n1,n2
            do m=m1,m2
              rr = sqrt(x(l1:l2)**2 + y(m)**2)
              f(l1:l2,m,n,iux) = -v_jet_peak * &
                                exp(- (((rr - r_jet_cent) / sigma_jet)**2)) * &
                                y(m) / rr
                f(l1:l2,m,n,iuy) = v_jet_peak * &
                                exp(- (((rr - r_jet_cent) / sigma_jet)**2)) * &
                                x(l1:l2) / rr
            enddo
          enddo
!
       case('vortex','point-vortex')
          do n=n1,n2
             do m=m1,m2
               xx = x(l1:l2)-xv0
               yy = y(  m  )-yv0
               rr = sqrt(xx**2 + yy**2)
               f(l1:l2,m,n,iux) = -Omega_vortex/(2*pi) * yy / rr**2
               f(l1:l2,m,n,iuy) =  Omega_vortex/(2*pi) * xx / rr**2
             enddo
          enddo
!
        case('Li-Ingersoll')
          do n=n1,n2
             do m=m1,m2
               xx = x(l1:l2)-xv0
               yy = y(  m  )-yv0
               rr = sqrt(xx**2 + yy**2)
!
               vr = vm * (rr/rm) * exp(1/bv * (1-(rr/rm)**bv))
!
               f(l1:l2,m,n,iux) = - vr * yy/rr
               f(l1:l2,m,n,iuy) =   vr * xx/rr
               
             enddo
          enddo
!
        endselect
      enddo
!
    endsubroutine initial_condition_uu
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
