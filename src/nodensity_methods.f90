module DensityMethods
!
!  11-mar-15/MR:  Created to avoid circular dependencies with EquationOfState.
!
  use Cparam
  use General, only: keep_compiler_quiet
  use Messages, only: fatal_error

  include 'density_methods.h'
!
  interface getrho
    module procedure getrho_1d
    module procedure getrho_2dyz
    module procedure getrho_2d
  endinterface
!
  interface getlnrho
    module procedure getlnrho_1d_x
    module procedure getlnrho_1d_y
    module procedure getlnrho_2dyz
    module procedure getlnrho_2d
  endinterface

  interface putlnrho
    module procedure putlnrho_s
    module procedure putlnrho_v
  endinterface

 real, pointer :: rho0, lnrho0

  contains
!***********************************************************************
    subroutine initialize_density_methods

      use SharedVariables, only: get_shared_variable
      use Messages, only: warning
      use Cdata, only: lstratz

      if (lstratz) &
        call warning('initialize_density_methods', &
                     'density methods not yet implemented for density_stratified')
       
      call get_shared_variable('rho0',rho0,caller='initialize_density_methods')
      call get_shared_variable('lnrho0',lnrho0)

    endsubroutine initialize_density_methods
!***********************************************************************
    function getrho_s(f,irf)

      real                :: getrho_s
      real,    intent(in) :: f
      integer, intent(in) :: irf        ! here dummy parameter only

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
    subroutine getlnrho_1d_y(f,ix,lnrho)

      real, dimension(my), intent(in) :: f
      real, dimension(ny), intent(out):: lnrho
      integer,             intent(in) :: ix

      lnrho=lnrho0

    endsubroutine getlnrho_1d_y
!***********************************************************************
    subroutine getlnrho_2d(f,lnrho)

      real, dimension(:,:), intent(in) :: f
      real, dimension(:,:), intent(out):: lnrho

      lnrho=lnrho0

    endsubroutine getlnrho_2d
!***********************************************************************
    subroutine getlnrho_2dyz(f,ix,lnrho)

      real, dimension(my,mz), intent(in) :: f
      real, dimension(my,mz), intent(out):: lnrho
      integer,                intent(in) :: ix

      lnrho=lnrho0

    endsubroutine getlnrho_2dyz
!***********************************************************************
    subroutine getrho_2d(f,rho)

      real, dimension(:,:), intent(in) :: f
      real, dimension(:,:), intent(out):: rho

      rho=rho0

    endsubroutine getrho_2d
!***********************************************************************
    subroutine getrho_2dyz(f,ix,rho)

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
    subroutine getdlnrho_x(f,il,ix,rho,dlnrho)

      integer,                       intent(in) :: il,ix
      real, dimension(-il:il,my,mz), intent(in) :: f
      real, dimension(my,mz),        intent(in) :: rho
      real, dimension(my,mz),        intent(out):: dlnrho

      dlnrho = 0.

    endsubroutine getdlnrho_x
!***********************************************************************
    subroutine getdlnrho_y(f,im,dlnrho)

      integer,                       intent(in) :: im
      real, dimension(mx,-im:im,mz), intent(in) :: f
      real, dimension(mx,mz),        intent(out):: dlnrho

      dlnrho = 0.

    endsubroutine getdlnrho_y
!***********************************************************************
   subroutine putrho(f,rho)

      real, dimension(mx), intent(out):: f
      real, dimension(nx), intent(in) :: rho

      call fatal_error('putrho', 'not implemented in nodensity.')
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(rho)
!
    endsubroutine putrho
!***********************************************************************
    subroutine putlnrho_v(f,lnrho)

      real, dimension(mx), intent(out):: f
      real, dimension(nx), intent(in) :: lnrho

      call fatal_error('putlnrho', 'not implemented in nodensity.')
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lnrho)
!
    endsubroutine putlnrho_v
!***********************************************************************
    subroutine putlnrho_s(f,lnrho)

      real, dimension(mx,my), intent(out):: f
      real,                   intent(in ):: lnrho

      call fatal_error('putlnrho', 'not implemented in nodensity.')
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lnrho)
!
    endsubroutine putlnrho_s
!***********************************************************************
endmodule DensityMethods
