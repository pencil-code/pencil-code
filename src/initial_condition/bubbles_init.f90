! $Id$
!
! Initial condition for bubbles in a stratified medium.
! Set bubble size, bubble density, surrounding medium density,
! temperatures, stratification and magnetic field.
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
!! NB: This is a very early stage of the initial condition. DO NOT USE IT!!!
!!!!!!!

module InitialCondition
!
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet

!
  implicit none
!
  include '../initial_condition.h'
!
  ! r_b = radius of the bubble
  ! [xyz]_b = position of the bubble
  ! rho_b = density of the bubble
  ! T_b = temperature of the bubble
  ! rho_m_0 = density of the medium at z = z_0
  ! z_0 > 0
  ! T_0 = temperature of the medium at z = z_0
  ! gamma = adiabatic index
  ! passive_scalar = determines if a passive scalar at the bubble should be set
  ! n_smooth = exponent of the smoothing function for the vector potential
  ! k_aa = wave vector of the initial magnetic vector potential in terms of the bubble radius
  ! ampl = amplitude of the initial magnetic vector potential
  
  real :: r_b = 1.0
  real :: x_b = 0, y_b = 0, z_b = 0
  real :: rho_b = 0.1, T_b = 5.0
  real :: rho_m_0 = 1.0, T_0 = 1.0
  real :: z_0
  integer :: passive_scalar = 0, n_smooth = 2
  real :: k_aa = 1., ampl = 1., asym_factor = 1., sigma_b = 1.
  
  namelist /initial_condition_pars/ &
      r_b, x_b, y_b, z_b, rho_b, T_b, rho_m_0, T_0, z_0, passive_scalar, k_aa, ampl, asym_factor, sigma_b, n_smooth
!
  contains
!***********************************************************************
    subroutine register_initial_condition()
!
!  Register variables associated with this module; likely none.
!
!  20-may-11/simon: coded
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
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine initial_condition_uu(f)
!
!  Initialize the velocity field.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
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
!  07-sep-11/simon: coded
!
      use SharedVariables
      use EquationOfState
      
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      integer :: l, m, n
      real :: log_rho_b, log_T_b, log_rho_m_0, log_T_0
      real :: r, z_tilde, rho_ex, T_ex  ! working variables for smoothing
!
    ! initialize the density of the medium and the bubble
    log_rho_b = log(rho_b)
    log_T_b = log(T_b)
    log_rho_m_0 = log(rho_m_0)
    log_T_0 = log(T_0)
    
    write(*,*) 'init_bubble: k_aa = ', k_aa
    do n = n1, n2, 1
      do m = m1, m2, 1
        do l = l1, l2, 1
          ! check if this point lies in the bubble
          if (((x(l) - x_b)**2 + (y(m) - y_b)**2 + (z(n) - z_b)**2) .le. r_b**2) then
            ! reduce the shocks at the start of the simulation by reducing the
            ! density and temperature gradient
!             r = sqrt((x(l) - x_b)**2 + (y(m) - y_b)**2 + (z(n) - z_b)**2)
!             z_tilde = z_b - (z(n)-z_b)/r*r_b
            
!             rho_ex = exp(log_rho_m_0 + log(1-(gamma-1)/gamma*z_tilde/z_0) / (gamma-1))
!             f(l,m,n,ilnrho) = log(rho_b * (1 - (r/r_b)**2) + rho_ex*(r/r_b)**2)
!             
!             T_ex = exp(log_T_0 + log(1-(gamma-1)/gamma*z_tilde/z_0))
!             f(l,m,n,ilnTT) = log(T_b * (1 - (r/r_b)**2) + T_ex*(r/r_b)**2)

            f(l,m,n,ilnrho) = log_rho_b

!             if (passive_scalar == 1) then
!                 f(l,m,n,ilncc) = 1.
!             endif
            
          else
            f(l,m,n,ilnrho) = log_rho_m_0 + log(1-(gamma-1)/gamma*z(n)/z_0) / (gamma-1)
          endif
        enddo
      enddo
    enddo    

    endsubroutine initial_condition_lnrho
!********************************************************************
    subroutine initial_condition_aa(f)
!
!  Initialize entropy.
!
!  07-sep-11/simon: coded
!
      use SharedVariables
      use EquationOfState

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f      
      integer :: l, m, n
      real :: log_T_b, log_T_0
!
    write(*,*) 'ini_bubble: aa'
    ! initialize the temperature of the medium and the bubble
    do n = n1, n2, 1
      do m = m1, m2, 1
        do l = l1, l2, 1
          ! check if this point lies in the bubble
          if (((x(l) - x_b)**2 + (y(m) - y_b)**2 + (z(n) - z_b)**2) .le. r_b**2) then
            f(l,m,n,iax) = ampl * (cos((y(m)-y_b)*k_aa/r_b) + sin((z(n)-z_b)*k_aa/r_b))
            f(l,m,n,iay) = ampl * (cos((z(n)-z_b)*k_aa/r_b) + sin((x(l)-x_b)*k_aa/r_b))
            f(l,m,n,iaz) = ampl * (cos((x(l)-x_b)*k_aa/r_b) + sin((y(m)-y_b)*k_aa/r_b))
            f(l,m,n,iax:iaz) = f(l,m,n,iax:iaz) * &
                (1-((x(l)-x_b)**n_smooth+(y(m)-y_b)**n_smooth+(z(n)-z_b)**n_smooth)/r_b**n_smooth)
          endif
        enddo
      enddo
    enddo
!      
    endsubroutine initial_condition_aa
!***********************************************************************
    subroutine initial_condition_ss(f)
!
!  Initialize entropy.
!
!  07-sep-11/simon: coded
!
      use SharedVariables
      use EquationOfState

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f      
      integer :: l, m, n
      real :: log_T_b, log_T_0
!
    ! initialize the temperature of the medium and the bubble
    log_T_b = log(T_b)
    log_T_0 = log(T_0)
    do n = n1, n2, 1
      do m = m1, m2, 1
        do l = l1, l2, 1
          ! check if this point lies in the bubble
          if (((x(l) - x_b)**2 + (y(m) - y_b)**2 + (z(n) - z_b)**2) .le. r_b**2) then
            f(l,m,n,ilnTT) = log_T_b
          else
            f(l,m,n,ilnTT) = log_T_0 + log(1-(gamma-1)/gamma*z(n)/z_0)
          endif
        enddo
      enddo
    enddo
!      
    endsubroutine initial_condition_ss
!***********************************************************************
    subroutine initial_condition_lncc(f)
!
!  Initialize entropy.
!
!  07-sep-11/simon: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f      
      integer :: l, m, n
!
    write(*,*) 'init cc'
    ! initialize the temperature of the medium and the bubble
    if (passive_scalar == 1) then
        do n = n1, n2, 1
            do m = m1, m2, 1
                do l = l1, l2, 1
                ! check if this point lies in the bubble
                if (((x(l) - x_b)**2 + (y(m) - y_b)**2 + (z(n) - z_b)**2) .le. r_b**2) then
                    f(l,m,n,ilncc) = 1.
                endif
                enddo
            enddo
        enddo
    endif
!      
    endsubroutine initial_condition_lncc
!***********************************************************************
    subroutine read_initial_condition_pars(unit,iostat)
!
!  07-may-09/wlad: coded
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=initial_condition_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=initial_condition_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_initial_condition_pars
!***********************************************************************
    subroutine write_initial_condition_pars(unit)
!     
      integer, intent(in) :: unit
!
      write(unit,NML=initial_condition_pars)
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
