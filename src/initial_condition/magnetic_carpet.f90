!  Initial condition (density, magnetic field, velocity) 
!  for a particular configuration of magnetic carpets.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!***************************************************************
module InitialCondition
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet  
  use Mpicomm
  use Messages
  use Boundcond ! for the core boundary communication
!
  implicit none
!
  include '../initial_condition.h'
!
! ampl = amplitude of the magnetic field
!
  real :: ampl = 1.0
  character (len=labellen) :: config='parasitic_polarities'
!
  namelist /initial_condition_pars/ &
    ampl, config
!
  contains
!***********************************************************************
  subroutine register_initial_condition()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  30-april-15/simon: coded
!
  endsubroutine register_initial_condition
!***********************************************************************
  subroutine initial_condition_uu(f)
!
!  Initialize the velocity field.
!
!  07-may-09/wlad: coded
!
    real, dimension (mx,my,mz,mfarray), intent(inout) :: f
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
!  07-may-09/wlad: coded
!
    real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
    call keep_compiler_quiet(f)
!
  endsubroutine initial_condition_lnrho
!***********************************************************************
  subroutine initial_condition_aa(f)
!
!  Initialize the magnetic vector potential.
!
!  30-april-15/simon: coded
!
!  Magnetic carpet consisting of three spots.
!
!  Created 2015-04-30 by Simon Candelaresi (Iomsn)
!
    use Mpicomm, only: stop_it
    use Poisson
    use Sub
!    
    real, dimension (mx,my,mz,mfarray) :: f
!        
!   The next variables are used for the uncurling.
    integer :: l, j, ju, k
    real, dimension (nx,ny,nz,3) :: jj, tmpJ  ! This is phi for poisson.f90
!
!   clear the magnetic field to zero
    f(:,:,:,iax:iaz) = 0.
!
    if (config == 'homogeneous') then
        do m = 1, my, 1
            do l = 1, mx, 1
                f(l,m,:,iax) = -ampl*y(m)/2
                f(l,m,:,iay) = ampl*x(l)/2
            enddo
        enddo
    endif
    
    if (config == 'parasitic_polarities') then
        do n = 1, mz, 1
            do m = 1, my, 1
                do l = 1, mx, 1
                    f(l,m,n,iax) = ampl * (&
                        2. * ((x(l) - 10.) ** 2 + y(m) ** 2 + (z(n) + 0.85) ** 2)&
                        ** (-3.0 / 2.) * y(m) + 2. * ((x(l) + 8.) ** 2 + y(m) **&
                        2 + (z(n) + 0.85) ** 2) ** (-3.0 / 2.) * y(m) + 2. * ((x(l)&
                        + 10.) ** 2 + y(m) ** 2 + (z(n) + 0.85) ** 2) ** (-3.0&
                        / 2.) * y(m) + 2. * ((x(l) + 6.) ** 2 + y(m) ** 2 + (z(n)&
                        + 0.85) ** 2) ** (-3.0 / 2.) * y(m) + 2. * (x(l) ** 2 + y(m)&
                        ** 2 + (z(n) + 0.85) ** 2) ** (-3.0 / 2.) * y(m) + 2. *&
                        ((x(l) + 2.) ** 2 + y(m) ** 2 + (z(n) + 0.85) ** 2) ** (-3.0&
                        / 2.) * y(m) + 2. * ((x(l) - 2.) ** 2 + y(m) ** 2 + (z(n)&
                        + 0.85) ** 2) ** (-3.0 / 2.) * y(m) + 2. * ((x(l) - 8.)&
                        ** 2 + (y(m) - 8.) ** 2 + (z(n) + 0.85) ** 2) ** (-3.0 /&
                        2.) * (y(m) - 8.) + 2. * ((x(l) - 6.) ** 2 + (y(m) - 8.)&
                        ** 2 + (z(n) + 0.85) ** 2) ** (-3.0 / 2.) * (y(m) - 8.) +&
                        2. * ((x(l) - 10.) ** 2 + (y(m) - 8.) ** 2 + (z(n) + 0.85)&
                        ** 2) ** (-3.0 / 2.) * (y(m) - 8.) + 2. * ((x(l) - 8.) **&
                        2 + (y(m) + 8.) ** 2 + (z(n) + 0.85) ** 2) ** (-3.0 / 2.)&
                        * (y(m) + 8.) + 2. * ((x(l) - 6.) ** 2 + (y(m) + 8.) ** 2&
                        + (z(n) + 0.85) ** 2) ** (-3.0 / 2.) * (y(m) + 8.) + 2. *&
                        ((x(l) - 10.) ** 2 + (y(m) + 8.) ** 2 + (z(n) + 0.85) **&
                        2) ** (-3.0 / 2.) * (y(m) + 8.) + 2. * ((x(l) + 8.) ** 2 +&
                        (y(m) - 8.) ** 2 + (z(n) + 0.85) ** 2) ** (-3.0 / 2.) * (y(m)&
                        - 8.) + 2. * ((x(l) + 10.) ** 2 + (y(m) - 8.) ** 2 +&
                        (z(n) + 0.85) ** 2) ** (-3.0 / 2.) * (y(m) - 8.) + 2. * ((x(l)&
                        + 6.) ** 2 + (y(m) - 8.) ** 2 + (z(n) + 0.85) ** 2) **&
                        (-3.0 / 2.) * (y(m) - 8.) + 2. * ((x(l) + 8.) ** 2 + (y(m)&
                        + 8.) ** 2 + (z(n) + 0.85) ** 2) ** (-3.0 / 2.) * (y(m)&
                        + 8.) + 2. * ((x(l) + 10.) ** 2 + (y(m) + 8.) ** 2 + (z(n)&
                        + 0.85) ** 2) ** (-3.0 / 2.) * (y(m) + 8.) + 2. * ((x(l)&
                        + 6.) ** 2 + (y(m) + 8.) ** 2 + (z(n) + 0.85) ** 2) ** (-3.0&
                        / 2.) * (y(m) + 8.) + 2. * (x(l) ** 2 + (y(m) + 8.) **&
                        2 + (z(n) + 0.85) ** 2) ** (-3.0 / 2.) * (y(m) + 8.) + 2.&
                        * ((x(l) + 2.) ** 2 + (y(m) + 8.) ** 2 + (z(n) + 0.85) **&
                        2) ** (-3.0 / 2.) * (y(m) + 8.) + 2. * ((x(l) - 2.) ** 2&
                        + (y(m) + 8.) ** 2 + (z(n) + 0.85) ** 2) ** (-3.0 / 2.) *&
                        (y(m) + 8.) + 2. * (x(l) ** 2 + (y(m) - 8.) ** 2 + (z(n) +&
                        0.85) ** 2) ** (-3.0 / 2.) * (y(m) - 8.) + 2. * ((x(l) +&
                        2.) ** 2 + (y(m) - 8.) ** 2 + (z(n) + 0.85) ** 2) ** (-3.0&
                        / 2.) * (y(m) - 8.) + 2. * ((x(l) - 2.) ** 2 + (y(m) - 8.)&
                        ** 2 + (z(n) + 0.85) ** 2) ** (-3.0 / 2.) * (y(m) - 8.)&
                        + 2. * ((x(l) - 8.) ** 2 + y(m) ** 2 + (z(n) + 0.85) ** 2)&
                        ** (-3.0 / 2.) * y(m) + 2. * ((x(l) - 6.) ** 2 + y(m) **&
                        2 + (z(n) + 0.85) ** 2) ** (-3.0 / 2.) * y(m))
                    f(l,m,n,iay) = ampl * (&
                        -2. * ((x(l) - 10.) ** 2 + y(m) ** 2 + (z(n) + 0.85) **&
                        2) ** (-3.0 / 2.) * (x(l) - 10.) - 2. * ((x(l) + 8.) ** 2&
                        + y(m) ** 2 + (z(n) + 0.85) ** 2) ** (-3.0 / 2.) * (x(l) +&
                        8.) - 2. * ((x(l) + 10.) ** 2 + y(m) ** 2 + (z(n) + 0.85)&
                        ** 2) ** (-3.0 / 2.) * (x(l) + 10.) - 2. * ((x(l) + 6.) **&
                        2 + y(m) ** 2 + (z(n) + 0.85) ** 2) ** (-3.0 / 2.) * (x(l)&
                        + 6.) - 2. * (x(l) ** 2 + y(m) ** 2 + (z(n) + 0.85) **&
                        2) ** (-3.0 / 2.) * x(l) - 2. * ((x(l) + 2.) ** 2 + y(m) **&
                        2 + (z(n) + 0.85) ** 2) ** (-3.0 / 2.) * (x(l) + 2.) - 2.&
                        * ((x(l) - 2.) ** 2 + y(m) ** 2 + (z(n) + 0.85) ** 2) **&
                        (-3.0 / 2.) * (x(l) - 2.) - 2. * ((x(l) - 8.) ** 2 + (y(m)&
                        - 8.) ** 2 + (z(n) + 0.85) ** 2) ** (-3.0 / 2.) * (x(l)&
                        - 8.) - 2. * ((x(l) - 6.) ** 2 + (y(m) - 8.) ** 2 + (z(n)&
                        + 0.85) ** 2) ** (-3.0 / 2.) * (x(l) - 6.) - 2. * ((x(l) -&
                        10.) ** 2 + (y(m) - 8.) ** 2 + (z(n) + 0.85) ** 2) ** (-3.0&
                        / 2.) * (x(l) - 10.) - 2. * ((x(l) - 8.) ** 2 + (y(m) +&
                        8.) ** 2 + (z(n) + 0.85) ** 2) ** (-3.0 / 2.) * (x(l) - 8.)&
                        - 2. * ((x(l) - 6.) ** 2 + (y(m) + 8.) ** 2 + (z(n) + 0.85)&
                        ** 2) ** (-3.0 / 2.) * (x(l) - 6.) - 2. * ((x(l) - 10.)&
                        ** 2 + (y(m) + 8.) ** 2 + (z(n) + 0.85) ** 2) ** (-3.0&
                        / 2.) * (x(l) - 10.) - 2. * ((x(l) + 8.) ** 2 + (y(m) - 8.)&
                        ** 2 + (z(n) + 0.85) ** 2) ** (-3.0 / 2.) * (x(l) + 8.)&
                        - 2. * ((x(l) + 10.) ** 2 + (y(m) - 8.) ** 2 + (z(n) + 0.85)&
                        ** 2) ** (-3.0 / 2.) * (x(l) + 10.) - 2. * ((x(l) + 6.)&
                        ** 2 + (y(m) - 8.) ** 2 + (z(n) + 0.85) ** 2) ** (-3.0 /&
                        2.) * (x(l) + 6.) - 2. * ((x(l) + 8.) ** 2 + (y(m) + 8.) **&
                        2 + (z(n) + 0.85) ** 2) ** (-3.0 / 2.) * (x(l) + 8.) - 2.&
                        * ((x(l) + 10.) ** 2 + (y(m) + 8.) ** 2 + (z(n) + 0.85)&
                        ** 2) ** (-3.0 / 2.) * (x(l) + 10.) - 2. * ((x(l) + 6.) **&
                        2 + (y(m) + 8.) ** 2 + (z(n) + 0.85) ** 2) ** (-3.0 / 2.)&
                        * (x(l) + 6.) - 2. * (x(l) ** 2 + (y(m) + 8.) ** 2 + (z(n)&
                        + 0.85) ** 2) ** (-3.0 / 2.) * x(l) - 2. * ((x(l) + 2.)&
                        ** 2 + (y(m) + 8.) ** 2 + (z(n) + 0.85) ** 2) ** (-3.0 / 2.)&
                        * (x(l) + 2.) - 2. * ((x(l) - 2.) ** 2 + (y(m) + 8.) **&
                        2 + (z(n) + 0.85) ** 2) ** (-3.0 / 2.) * (x(l) - 2.) - 2.&
                        * (x(l) ** 2 + (y(m) - 8.) ** 2 + (z(n) + 0.85) ** 2) **&
                        (-3.0 / 2.) * x(l) - 2. * ((x(l) + 2.) ** 2 + (y(m) - 8.)&
                        ** 2 + (z(n) + 0.85) ** 2) ** (-3.0 / 2.) * (x(l) + 2.) -&
                        2. * ((x(l) - 2.) ** 2 + (y(m) - 8.) ** 2 + (z(n) + 0.85)&
                        ** 2) ** (-3.0 / 2.) * (x(l) - 2.) - 2. * ((x(l) - 8.) **&
                        2 + y(m) ** 2 + (z(n) + 0.85) ** 2) ** (-3.0 / 2.) * (x(l)&
                        - 8.) - 2. * ((x(l) - 6.) ** 2 + y(m) ** 2 + (z(n) + 0.85)&
                        ** 2) ** (-3.0 / 2.) * (x(l) - 6.) + x(l))
                    f(l,m,n,iaz) = 0
                enddo          
            enddo          
        enddo
    elseif (config == 'dominant_polarities') then
        do n = 1, mz, 1
            do m = 1, my, 1
                do l = 1, mx, 1
                    f(l,m,n,iax) = ampl * (&
                        -1. * ((x(l) - 10.) ** 2 + ((y(m) + 8) ** 2) + (z(n) + 0.5)&
                        ** 2) ** (-3.0 / 2.) * (y(m) + 8) - 1. * ((x(l) - 8.)&
                        ** 2 + ((y(m) - 8) ** 2) + (z(n) + 0.5) ** 2) ** (-3.0 / 2.)&
                        * (y(m) - 8) - 1. * ((x(l) - 6.) ** 2 + ((y(m) - 8) **&
                        2) + (z(n) + 0.5) ** 2) ** (-3.0 / 2.) * (y(m) - 8) - 1. *&
                        ((x(l) - 10.) ** 2 + ((y(m) - 8) ** 2) + (z(n) + 0.5) **&
                        2) ** (-3.0 / 2.) * (y(m) - 8) - 1. * (x(l) ** 2 + (y(m) **&
                        2) + (z(n) + 0.5) ** 2) ** (-3.0 / 2.) * (y(m)) - 1. * ((x(l)&
                        + 2.) ** 2 + (y(m) ** 2) + (z(n) + 0.5) ** 2) ** (-3.0&
                        / 2.) * (y(m)) - 1. * ((x(l) - 2.) ** 2 + (y(m) ** 2) +&
                        (z(n) + 0.5) ** 2) ** (-3.0 / 2.) * (y(m)) - 1. * ((x(l)&
                        + 8.) ** 2 + (y(m) ** 2) + (z(n) + 0.5) ** 2) ** (-3.0 / 2.)&
                        * (y(m)) - 1. * ((x(l) + 10.) ** 2 + (y(m) ** 2) + (z(n)&
                        + 0.5) ** 2) ** (-3.0 / 2.) * (y(m)) - 1. * ((x(l) + 6.)&
                        ** 2 + (y(m) ** 2) + (z(n) + 0.5) ** 2) ** (-3.0 / 2.) *&
                        (y(m)) - 1. * ((x(l) + 8.) ** 2 + ((y(m) - 8) ** 2) + (z(n)&
                        + 0.5) ** 2) ** (-3.0 / 2.) * (y(m) - 8) - 1. * ((x(l) +&
                        10.) ** 2 + ((y(m) - 8) ** 2) + (z(n) + 0.5) ** 2) ** (-3.0&
                        / 2.) * (y(m) - 8) - 1. * ((x(l) + 6.) ** 2 + ((y(m) -&
                        8) ** 2) + (z(n) + 0.5) ** 2) ** (-3.0 / 2.) * (y(m) - 8)&
                        - 1. * (x(l) ** 2 + ((y(m) + 8) ** 2) + (z(n) + 0.5) ** 2)&
                        ** (-3.0 / 2.) * (y(m) + 8) - 1. * ((x(l) + 2.) ** 2 + ((y(m)&
                        + 8) ** 2) + (z(n) + 0.5) ** 2) ** (-3.0 / 2.) * (y(m)&
                        + 8) - 1. * ((x(l) - 2.) ** 2 + ((y(m) + 8) ** 2) + (z(n)&
                        + 0.5) ** 2) ** (-3.0 / 2.) * (y(m) + 8) - 1. * (x(l) **&
                        2 + ((y(m) - 8) ** 2) + (z(n) + 0.5) ** 2) ** (-3.0 / 2.)&
                        * (y(m) - 8) - 1. * ((x(l) + 2.) ** 2 + ((y(m) - 8) ** 2)&
                        + (z(n) + 0.5) ** 2) ** (-3.0 / 2.) * (y(m) - 8) - 1. * ((x(l)&
                        - 2.) ** 2 + ((y(m) - 8) ** 2) + (z(n) + 0.5) ** 2)&
                        ** (-3.0 / 2.) * (y(m) - 8) - 1. * ((x(l) - 8.) ** 2 + (y(m)&
                        ** 2) + (z(n) + 0.5) ** 2) ** (-3.0 / 2.) * (y(m)) - 1.&
                        * ((x(l) - 6.) ** 2 + (y(m) ** 2) + (z(n) + 0.5) ** 2) **&
                        (-3.0 / 2.) * (y(m)) - 1. * ((x(l) - 10.) ** 2 + (y(m) **&
                        2) + (z(n) + 0.5) ** 2) ** (-3.0 / 2.) * (y(m)) - 1. * ((x(l)&
                        + 8.) ** 2 + ((y(m) + 8) ** 2) + (z(n) + 0.5) ** 2) **&
                        (-3.0 / 2.) * (y(m) + 8) - 1. * ((x(l) + 10.) ** 2 + ((y(m)&
                        + 8) ** 2) + (z(n) + 0.5) ** 2) ** (-3.0 / 2.) * (y(m)&
                        + 8) - 1. * ((x(l) + 6.) ** 2 + ((y(m) + 8) ** 2) + (z(n)&
                        + 0.5) ** 2) ** (-3.0 / 2.) * (y(m) + 8) - 1. * ((x(l) -&
                        8.) ** 2 + ((y(m) + 8) ** 2) + (z(n) + 0.5) ** 2) ** (-3.0&
                        / 2.) * (y(m) + 8) - 1. * ((x(l) - 6.) ** 2 + ((y(m) + 8)&
                        ** 2) + (z(n) + 0.5) ** 2) ** (-3.0 / 2.) * (y(m) + 8))
                    f(l,m,n,iay) = ampl * (&
                        1. * (x(l) ** 2 + ((y(m) - 8) ** 2) + (z(n) + 0.5) ** 2)&
                        ** (-3.0 / 2.) * x(l) + 1. * ((x(l) + 2.) ** 2 + ((y(m) -&
                        8) ** 2) + (z(n) + 0.5) ** 2) ** (-3.0 / 2.) * (x(l) + 2.)&
                        + 1. * ((x(l) - 2.) ** 2 + ((y(m) - 8) ** 2) + (z(n) + 0.5)&
                        ** 2) ** (-3.0 / 2.) * (x(l) - 2.) + 1. * ((x(l) + 8.)&
                        ** 2 + ((y(m) + 8) ** 2) + (z(n) + 0.5) ** 2) ** (-3.0 /&
                        2.) * (x(l) + 8.) + 1. * ((x(l) + 10.) ** 2 + ((y(m) + 8)&
                        ** 2) + (z(n) + 0.5) ** 2) ** (-3.0 / 2.) * (x(l) + 10.) +&
                        1. * ((x(l) + 6.) ** 2 + ((y(m) + 8) ** 2) + (z(n) + 0.5)&
                        ** 2) ** (-3.0 / 2.) * (x(l) + 6.) + 1. * ((x(l) - 8.) **&
                        2 + ((y(m) + 8) ** 2) + (z(n) + 0.5) ** 2) ** (-3.0 / 2.)&
                        * (x(l) - 8.) + 1. * ((x(l) - 6.) ** 2 + ((y(m) + 8) ** 2)&
                        + (z(n) + 0.5) ** 2) ** (-3.0 / 2.) * (x(l) - 6.) + 1. *&
                        ((x(l) - 10.) ** 2 + ((y(m) + 8) ** 2) + (z(n) + 0.5) **&
                        2) ** (-3.0 / 2.) * (x(l) - 10.) + 1. * ((x(l) - 8.) ** 2&
                        + ((y(m) - 8) ** 2) + (z(n) + 0.5) ** 2) ** (-3.0 / 2.) *&
                        (x(l) - 8.) + 1. * ((x(l) - 6.) ** 2 + ((y(m) - 8) ** 2) +&
                        (z(n) + 0.5) ** 2) ** (-3.0 / 2.) * (x(l) - 6.) + 1. * ((x(l)&
                        - 10.) ** 2 + ((y(m) - 8) ** 2) + (z(n) + 0.5) ** 2)&
                        ** (-3.0 / 2.) * (x(l) - 10.) + 1. * (x(l) ** 2 + (y(m) **&
                        2) + (z(n) + 0.5) ** 2) ** (-3.0 / 2.) * x(l) + 1. * ((x(l)&
                        + 2.) ** 2 + (y(m) ** 2) + (z(n) + 0.5) ** 2) ** (-3.0&
                        / 2.) * (x(l) + 2.) + 1. * ((x(l) - 2.) ** 2 + (y(m) ** 2)&
                        + (z(n) + 0.5) ** 2) ** (-3.0 / 2.) * (x(l) - 2.) + 1. *&
                        ((x(l) + 8.) ** 2 + (y(m) ** 2) + (z(n) + 0.5) ** 2) ** (-3.0&
                        / 2.) * (x(l) + 8.) + 1. * ((x(l) + 10.) ** 2 + (y(m)&
                        ** 2) + (z(n) + 0.5) ** 2) ** (-3.0 / 2.) * (x(l) + 10.) +&
                        1. * ((x(l) + 6.) ** 2 + (y(m) ** 2) + (z(n) + 0.5) ** 2)&
                        ** (-3.0 / 2.) * (x(l) + 6.) + 1. * ((x(l) - 8.) ** 2 + (y(m)&
                        ** 2) + (z(n) + 0.5) ** 2) ** (-3.0 / 2.) * (x(l) - 8.)&
                        + 1. * ((x(l) - 6.) ** 2 + (y(m) ** 2) + (z(n) + 0.5) **&
                        2) ** (-3.0 / 2.) * (x(l) - 6.) + 1. * ((x(l) - 10.) **&
                        2 + (y(m) ** 2) + (z(n) + 0.5) ** 2) ** (-3.0 / 2.) * (x(l)&
                        - 10.) + 1. * ((x(l) + 8.) ** 2 + ((y(m) - 8) ** 2) + (z(n)&
                        + 0.5) ** 2) ** (-3.0 / 2.) * (x(l) + 8.) + 1. * ((x(l)&
                        + 10.) ** 2 + ((y(m) - 8) ** 2) + (z(n) + 0.5) ** 2) **&
                        (-3.0 / 2.) * (x(l) + 10.) + 1. * ((x(l) + 6.) ** 2 + ((y(m)&
                        - 8) ** 2) + (z(n) + 0.5) ** 2) ** (-3.0 / 2.) * (x(l)&
                        + 6.) + 1. * (x(l) ** 2 + ((y(m) + 8) ** 2) + (z(n) + 0.5)&
                        ** 2) ** (-3.0 / 2.) * x(l) + 1. * ((x(l) + 2.) ** 2 + ((y(m)&
                        + 8) ** 2) + (z(n) + 0.5) ** 2) ** (-3.0 / 2.) * (x(l)&
                        + 2.) + 1. * ((x(l) - 2.) ** 2 + ((y(m) + 8) ** 2) + (z(n)&
                        + 0.5) ** 2) ** (-3.0 / 2.) * (x(l) - 2.) + x(l))
                    f(l,m,n,iaz) = 0
                enddo          
            enddo          
        enddo
    endif
! !     Transform the magnetic field into a vector potential
! !
! !     communicate the core boundaries for taking the curl
!       call MPI_BARRIER(MPI_comm_world, ierr)
!       call boundconds_x(f)
!       call boundconds_y(f)
!       call boundconds_z(f)
!       call initiate_isendrcv_bdry(f)
!       call finalize_isendrcv_bdry(f)
!       call MPI_BARRIER(MPI_comm_world, ierr)
! 
! !     Compute curl(B) = J for the Poisson solver
!       do m=m1,m2
!         do n=n1,n2
!           call curl(f,iaa,jj(:,m-nghost,n-nghost,:))
!         enddo
!       enddo
!       tmpJ = -jj
! !
! !     Use the Poisson solver to solve \nabla^2 A = -J for A
!       do j=1,3
!         call inverse_laplacian(tmpJ(:,:,:,j))
!       enddo
! !
! !     Overwrite the f-array with the correct vector potential A
!       do j=1,3
!         ju=iaa-1+j
!         f(l1:l2,m1:m2,n1:n2,ju) = tmpJ(:,:,:,j)
!       enddo
! !
! !     Add a background field to the braid
!       do l=1,mx
!         do m=1,my
!           f(l,m,:,iax) = f(l,m,:,iax) - y(m)*B_bkg/2.
!           f(l,m,:,iay) = f(l,m,:,iay) + x(l)*B_bkg/2.
!         enddo
!       enddo
!
! !     communicate the core boundaries for taking the curl
!       call MPI_BARRIER(MPI_comm_world, ierr)
!       call boundconds_x(f)
!       call initiate_isendrcv_bdry(f)
!       call finalize_isendrcv_bdry(f)
!       call boundconds_y(f)
!       call boundconds_z(f)
!       call MPI_BARRIER(MPI_comm_world, ierr)
! !
!
  endsubroutine initial_condition_aa
!***********************************************************************
    subroutine read_initial_condition_pars(iostat)
!
      use File_io, only: get_unit
!
      integer, intent(out) :: iostat
      include "../parallel_unit.h"
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
!***********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from noinitial_condition.f90 for any    **
!**  InitialCondition routines not implemented in this file        **
!**                                                                **
  include '../initial_condition_dummies.inc'
!********************************************************************
endmodule InitialCondition
