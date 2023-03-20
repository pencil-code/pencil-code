! $Id$
!
!  This module provide a way for users to specify custom initial
!  conditions.
!
!  The module provides a set of standard hooks into the Pencil Code
!  and currently allows the following customizations:
!
!   Description                               | Relevant function call
!  ------------------------------------------------------------------------
!   Initial condition registration            | register_initial_condition
!     (pre parameter read)                    |
!   Initial condition initialization          | initialize_initial_condition
!     (post parameter read)                   |
!                                             |
!   Initial condition for momentum            | initial_condition_uu
!   Initial condition for density             | initial_condition_lnrho
!   Initial condition for entropy             | initial_condition_ss
!   Initial condition for magnetic potential  | initial_condition_aa
!                                             |
!   Initial condition for all in one call     | initial_condition_all
!     (called last)                           |
!
!   And a similar subroutine for each module with an "init_XXX" call.
!   The subroutines are organized IN THE SAME ORDER THAT THEY ARE CALLED.
!   First uu, then lnrho, then ss, then aa, and so on.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!***************************************************************
!
! HOW TO USE THIS FILE
! --------------------
!
! Change the line above to
!    linitial_condition = .true.
! to enable use of custom initial conditions.
!
! The rest of this file may be used as a template for your own initial
! conditions. Simply fill out the prototypes for the features you want
! to use.
!
! Save the file with a meaningful name, e.g. mhs_equilibrium.f90, and
! place it in the $PENCIL_HOME/src/initial_condition directory. This
! path has been created to allow users to optionally check their
! contributions in to the Pencil Code SVN repository. This may be
! useful if you are working on/using an initial condition with
! somebody else or may require some assistance from one from the main
! Pencil Code team. HOWEVER, less general initial conditions should
! not go here (see below).
!
! You can also place initial condition files directly in the run
! directory. Simply create the folder 'initial_condition' at the same
! level as the *.in files and place an initial condition file there.
! With pc_setupsrc this file is linked automatically into the local
! src directory. This is the preferred method for initial conditions
! that are not very general.
!
! To use your additional initial condition code, edit the
! Makefile.local in the src directory under the run directory in which
! you wish to use your initial condition. Add a line that says e.g.
!
!    INITIAL_CONDITION =   initial_condition/mhs_equilibrium
!
! Here mhs_equilibrium is replaced by the filename of your new file,
! not including the .f90 extension.
!
! This module is based on Tony's special module.
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
  real :: l0_d, l_d1, psipn0, psipn1, psi0, psi1, lnrho0_c
  real :: m0, r0_d=40.0, rs=1.0, apara, h_d, dsteep, ngamma, &
          rho0_c=1.0e-5, cs0_c=1.0, Tc
!
  namelist /initial_condition_pars/ m0, r0_d, rs, apara, & 
  h_d, rho0_c, dsteep, ngamma, cs0_c, Tc, l0_d
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
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine initial_condition_uu(f)
!
!  Initialize the velocity field.
!
!  07-may-09/wlad: coded

!  /mayank

      use Sub,            only: get_radial_distance
      use SharedVariables, only: get_shared_variable
      
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx) :: rr_sph,rr_cyl,g_r
      real, dimension (mx) :: vphi_d, l_d
      real :: m00
      real, pointer :: g0

      call get_shared_variable('g0',g0)
      m00=(g0/G_Newton)
      l0_d=sqrt(g0*r0_d**3.0)/(r0_d-rs)

      do m=m1, m2;do n=n1, n2
          call get_radial_distance(rr_sph,rr_cyl) 

! Specific angular momentum

      l_d=l0_d*((rr_cyl(l1:l2))/r0_d)**apara

! azimuthal velocity. The value of loop variable 'n' has been used as the z coordinate.
      
      vphi_d=(1.0-tanh(abs(z(n)/h_d))**dsteep)*tanh(abs(rr_sph/rs)**dsteep)*l_d/(rr_cyl)
      f(:,m,n,iuz) =  vphi_d     
      enddo;enddo

!
    endsubroutine initial_condition_uu
!***********************************************************************
    subroutine initial_condition_lnrho(f)
!
!  Initialize logarithmic density. init_lnrho will take care of
!  converting it to linear density if you use ldensity_nolog.
!
!  07-may-09/wlad: coded

!/mayank

      use EquationOfState, only: get_cp1, cs0, cs20, cs2bot, cs2top, rho0, lnrho0, &
                             gamma, gamma1, gamma_m1
      use FArrayManager,   only: farray_use_global
      use Sub,             only: get_radial_distance
      use SharedVariables, only: get_shared_variable

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx) :: rr_sph,rr_cyl
      real, dimension (nx) :: lnrho_d, lnrho_c, psipn, psi, l_d
      real, pointer :: g0
      real :: cp1, m00

! Some constants

      lnrho0_c=log(rho0_c)
      l0_d=sqrt(G_Newton*m0*r0_d**3.0)/(r0_d-rs)
      l_d1=l0_d*(2.0*rs/r0_d)**apara
      psipn0=-(G_Newton*m0)/(r0_d-rs)
      psipn1=-G_Newton*m0/rs
      psi0=psipn0+1.0/(2.0-2.0*apara)*(l0_d/(r0_d))**2.0
      psi1=psipn1+1.0/(2.0-2.0*apara)*(l_d1/(2.0*rs))**2.0
!    
      call get_cp1(cp1) 
      call get_shared_variable('g0',g0)
      m00=(g0/G_Newton)
      do n=n1,n2
        do m=m1,m2
          call get_radial_distance(rr_sph,rr_cyl)
!
! Specific angular momentum
!        
          l_d=l0_d*((rr_cyl(l1:l2))/r0_d)**apara
!       
!Pseudo-Newtonian potential
!        
        psipn=-g0/sqrt((rr_sph(l1:l2)-rs)**2+(1e-3*rs)**2)

! Gravitational Potential for the system
        
        psi=psipn+1.0/(2.0-2.0*apara)*(l_d/(sqrt(rr_cyl(l1:l2)**2+(1e-3*rs)**2)))**2.0

!  Disk Density
        
        lnrho_d=ngamma*log(abs(1.0-gamma*(psi-psi0)/(cs0*(ngamma+1.0))))+lnrho0
          print*, m, n, minval(-gamma*(psi-psi0)/(cs0*(ngamma+1.))),maxval(-gamma*(psi-psi0)/(cs0*(ngamma+1.)))

! Corona Density
        
        lnrho_c=lnrho0_c-(psipn-psipn1)*gamma/cs0_c**2
        f(l1:l2,m,n,ilnrho) = lnrho_d+lnrho_c
        
        enddo
      enddo

!
    endsubroutine initial_condition_lnrho
!***********************************************************************
    subroutine initial_condition_ss(f)
!
!  Initialize entropy.
!
!  07-may-09/wlad: coded

!/mayank

      use FArrayManager,   only: farray_use_global
      use EquationOfState, only: get_cp1, cs0, cs20, cs2bot, cs2top, rho0, lnrho0, &
                             gamma, gamma1, gamma_m1
      
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx) :: lnT_d, lnT_c
      real, dimension (mx) :: lnrho_d
      real :: cp1

      call get_cp1(cp1)
      do n=n1,n2
        do m=m1,m2
        lnrho_d= f(:,m,n,ilnrho)

! Disk Temperature

      lnT_d=(lnrho_d-lnrho0)/ngamma+log(cs0/(gamma_m1/cp1))
      f(:,m,n,ilnTT) = lnT_d
      
      end do
     end do
! Corona Temperature is constant. Defined in initial_condition_pars as 'Tc'.

!
    endsubroutine initial_condition_ss
!***********************************************************************
    subroutine initial_condition_aa(f)
!
!  Initialize the magnetic vector potential.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_aa
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
!********************************************************************!
endmodule InitialCondition
