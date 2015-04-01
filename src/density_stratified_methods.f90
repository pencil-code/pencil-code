module DensityMethods
!
!  11-mar-15/MR:  Created to avoid circular dependencies with EquationOfState.
!                 Yet to be adapted to density_stratified.
  use Cparam
  use Cdata

  implicit none

  include 'density_methods.h'  
!
  public :: putrho, putlnrho

  interface putrho
    module procedure putrho_s
    module procedure putrho_v
  endinterface
!
  real, dimension(mz) :: rho0z

  contains
!***********************************************************************
    subroutine initialize_density_methods
!
      rho0z = 0.0
      if (lstratz) call get_stratz(z, rho0z)
!
    endsubroutine initialize_density_methods
!***********************************************************************
    function getrho_s(f,iz)

      real                :: getrho_s
      real,    intent(in) :: f
      integer, intent(in) :: iz

      getrho_s = rho0z(iz)*(1.+f)

    endfunction getrho_s
!***********************************************************************
    subroutine getrho_1d(f,rho)

      real, dimension(mx), intent(in) :: f
      real, dimension(nx), intent(out):: rho

      rho=f(l1:l2)

    endsubroutine getrho_1d
!***********************************************************************
    subroutine getlnrho_1d_x(f,lnrho)

      real, dimension(mx), intent(in) :: f
      real, dimension(nx), intent(out):: lnrho

      lnrho=log(f(l1:l2))

    endsubroutine getlnrho_1d_x
!***********************************************************************
    subroutine getlnrho_1d_y(f,ix,lnrho)

      real, dimension(my), intent(in) :: f
      real, dimension(ny), intent(out):: lnrho
      integer,             intent(in) :: ix

      lnrho=log(f(m1:m2))

    endsubroutine getlnrho_1d_y
!***********************************************************************
    subroutine getlnrho_2d(f,lnrho)

      real, dimension(:,:), intent(in) :: f
      real, dimension(:,:), intent(out):: lnrho

      lnrho=log(f)

    endsubroutine getlnrho_2d
!***********************************************************************
    subroutine getlnrho_2dyz(f,ix,lnrho)

      real, dimension(my,mz), intent(in) :: f
      real, dimension(my,mz), intent(out):: lnrho
      integer,                intent(in) :: ix

      lnrho=log(f)

    endsubroutine getlnrho_2dyz
!***********************************************************************
    subroutine getrho_2d(f,rho)

      real, dimension(:,:), intent(in) :: f
      real, dimension(:,:), intent(out):: rho

      rho=f

    endsubroutine getrho_2d
!***********************************************************************
    subroutine getrho_2dyz(f,ix,rho)

      real, dimension(my,mz), intent(in) :: f
      real, dimension(my,mz), intent(out):: rho
      integer,                intent(in) :: ix

      rho=f

    endsubroutine getrho_2dyz 
!***********************************************************************
    subroutine putrho_v(f,rho)

      real, dimension(mx), intent(out):: f
      real, dimension(nx), intent(in) :: rho
!
      f(l1:l2)=rho

    endsubroutine putrho_v
!***********************************************************************
    subroutine putrho_s(f,rho)

      real, dimension(mx,my), intent(out):: f
      real,                   intent(in) :: rho
!
      integer :: m

      f(l1:l2,:)=rho

    endsubroutine putrho_s
!***********************************************************************
    subroutine putlnrho(f,lnrho)

      real, dimension(mx), intent(out):: f
      real, dimension(nx), intent(in) :: lnrho
!
      f(l1:l2)=exp(lnrho)
      
    endsubroutine putlnrho
!***********************************************************************
    subroutine getdlnrho_z(f,in,dlnrho)

      integer,                       intent(in) :: in
      real, dimension(mx,my,-in:in), intent(in) :: f
      real, dimension(mx,my),        intent(out):: dlnrho

      dlnrho = (f(:,:,in)-f(:,:,-in))/f(:,:,0)
!
    endsubroutine getdlnrho_z
!***********************************************************************
    subroutine getdlnrho_y(f,im,dlnrho)

      integer,                       intent(in) :: im
      real, dimension(mx,-im:im,mz), intent(in) :: f
      real, dimension(mx,mz),        intent(out):: dlnrho
!
      call fatal_error('getdlnrho_y', 'not implemented')
!
    endsubroutine getdlnrho_y
!***********************************************************************
    subroutine getdlnrho_x(f,il,ix,rho,dlnrho)

      integer,                       intent(in) :: il,ix
      real, dimension(-il:il,my,mz), intent(in) :: f
      real, dimension(my,mz),        intent(in) :: rho
      real, dimension(my,mz),        intent(out):: dlnrho

      dlnrho = (f(il,:,:)-f(-il,:,:))/rho

    endsubroutine getdlnrho_x
!***********************************************************************
endmodule DensityMethods
