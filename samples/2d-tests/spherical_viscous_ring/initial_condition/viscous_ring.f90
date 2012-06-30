!  $Id: mhs_equilibrium.f90 14134 2010-06-16 18:21:01Z wladimir.lyra $
!
!  Initial condition for spherical viscous ring, according 
!  to the test of Frederic Masset, 
!
!    2D 1/2 thick off-centered Keplerian viscous ring spread. 
!    http://www.maths.qmul.ac.uk/~masset/hd/tests.html
!
!  The original files can be found at at:  
!
!   http://www.maths.qmul.ac.uk/~masset/hd/tests/vkoffring2d.ini
!   http://www.maths.qmul.ac.uk/~masset/hd/tests/vkoffring2d.pot
!
!  07-may-09/wlad: adapted from noinitial_condition.f90
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
  use Messages
!
  implicit none
!
  include '../initial_condition.h'
!
  real :: nu, cs20
  real :: time0=0.018
  real :: sigmaz=0.3
!
  namelist /initial_condition_pars/ nu, cs20
!
  contains
!***********************************************************************
    subroutine register_initial_condition()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  07-oct-09/wlad: coded
!
!  Identify CVS/SVN version information.
!
      if (lroot) call svn_id( &
           "$Id: mhs_equilibrium.f90 14134 2010-06-16 18:21:01Z wladimir.lyra $")
!
    endsubroutine register_initial_condition
!***********************************************************************
    subroutine initial_condition_uu(f)
!
!  Initialize the velocity field.
!
!  07-may-09/wlad: coded
!
      use Sub,            only: get_radial_distance
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (nx) :: rr_sph,rr_cyl,urad,uphi
      real, dimension (nx) :: omega,pressure_correction
      real, dimension (nx,3) :: uu
!
      do n=n1,n2
        do m=m1,m2
!
          call get_radial_distance(rr_sph,rr_cyl)
!
          pressure_correction = 2*cs20*(rr_cyl-1)/(rr_cyl*time0)
!
          omega=sqrt(1./rr_cyl**3 - pressure_correction)
!              
          uphi = rr_cyl*omega
          urad = -1.5*nu/rr_cyl+6.0*nu*(rr_cyl-1.0)/time0
!
          if (lspherical_coords) then 
            uu(:,1) = urad*sinth(m)
            uu(:,2) = urad/rr_sph*costh(m)
            uu(:,3) = uphi
          elseif (lcylindrical_coords) then 
            uu(:,1) = urad
            uu(:,2) = uphi
            uu(:,3) = 0.
          endif
!    
          f(l1:l2,m,n,iux:iuz) = uu
!
        enddo
      enddo
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
      use Sub, only: get_radial_distance
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx) :: rr_sph,rr_cyl,z_mn
      real, dimension (mx) :: tmp,logrhor,logrhoz
!
!  Pencilize the density allocation.
!
      do n=1,mz
        do m=1,my
!
          call get_radial_distance(rr_sph,rr_cyl)
          if (lspherical_coords) then 
            z_mn = rr_sph*costh(m)
          elseif (lcylindrical_coords) then
            z_mn = z(n)
          endif
!
          tmp = 1./(2*pi*sqrt(pi*time0)*rr_cyl**0.75)
          logrhor = -1./time0 *(rr_cyl-1)**2
          logrhoz = -1./sigmaz*(z_mn  -1)**2
!
          f(:,m,n,ilnrho) = log(tmp) + logrhor + logrhoz
!
        enddo
      enddo
!
    endsubroutine initial_condition_lnrho
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
