module DensityMethods
!
!  11-mar-15/MR:  Created to avoid circular dependencies with EquationOfState.
!
  use Cparam

  include 'density_methods.h'

  real, pointer :: rho0, lnrho0

  contains
!***********************************************************************
    subroutine initialize_density_methods

      use SharedVariables, only: get_shared_variable
       
      call get_shared_variable('rho0',rho0, caller='initialize_density_methods')
      call get_shared_variable('lnrho0',lnrho0, caller='initialize_density_methods')

    endsubroutine initialize_density_methods
!***********************************************************************
    function getrho_s(f,lf)

      real                :: getrho_s
      real,    intent(in) :: f
      integer, intent(in) :: lf        ! here dummy parameter only

      getrho_s = rho0

    endfunction getrho_s
!***********************************************************************
    subroutine getrho_1d(f,rho)

      real, dimension(mx), intent(in) :: f
      real, dimension(nx), intent(out):: rho

      rho=rho0

    endsubroutine getrho_1d
!***********************************************************************
    subroutine getlnrho_1d_x(f,lnrho)

      real, dimension(mx), intent(in) :: f
      real, dimension(nx), intent(out):: lnrho

      lnrho=lnrho0

    endsubroutine getlnrho_1d_x
!***********************************************************************
    subroutine getlnrho_1d_y(f,lnrho,topbot)

      real, dimension(my), intent(in) :: f
      real, dimension(ny), intent(out):: lnrho
      integer,             intent(in) :: topbot

      lnrho=lnrho0

    endsubroutine getlnrho_1d_y
!***********************************************************************
    subroutine getlnrho_2dxy(f,lnrho)

      real, dimension(mx,my), intent(in) :: f
      real, dimension(mx,my), intent(out):: lnrho

      lnrho=lnrho0

    endsubroutine getlnrho_2dxy
!***********************************************************************
    subroutine getlnrho_2dyz(f,lnrho,ix)

      real, dimension(my,mz), intent(in) :: f
      real, dimension(my,mz), intent(out):: lnrho
      integer,                intent(in) :: ix

      lnrho=lnrho0

    endsubroutine getlnrho_2dyz
!***********************************************************************
    subroutine getrho_2dxy(f,rho)

      real, dimension(mx,my), intent(in) :: f
      real, dimension(mx,my), intent(out):: rho

      rho=rho0

    endsubroutine getrho_2dxy
!***********************************************************************
    subroutine getrho_2dyz(f,rho,ix)

      real, dimension(my,mz), intent(in) :: f
      real, dimension(my,mz), intent(out):: rho
      integer,                intent(in) :: ix

      rho=rho0

    endsubroutine getrho_2dyz 
!***********************************************************************
    subroutine getdlnrho_z(f,in,dlnrho)

      integer,                       intent(in) :: in
      real, dimension(mx,my,-in:in), intent(in) :: f
      real, dimension(mx,my),        intent(out):: dlnrho

      dlnrho = 0.
!
    endsubroutine getdlnrho_z
!***********************************************************************
    subroutine getdlnrho_x(f,il,dx2,dlnrho,ix)

      integer,                       intent(in) :: il,ix
      real, dimension(-il:il,my,mz), intent(in) :: f
      real,                          intent(in) :: dx2
      real, dimension(my,mz),        intent(out):: dlnrho

      dlnrho = 0.

    endsubroutine getdlnrho_x
!***********************************************************************
endmodule DensityMethods
