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
  interface getrho1
    module procedure getrho1_1d
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
    subroutine getrho1_1d(f,rho1)
!
!  Fetches inverse of density.
!
!   4-oct-17/MR: derived from getrho_1d.
!

      real, dimension(mx), intent(in) :: f
      real, dimension(nx), intent(out):: rho1

      rho1=1./rho0
      call keep_compiler_quiet(f)

    endsubroutine getrho1_1d
!***********************************************************************
    function getrho_s(f,irf)

      real                :: getrho_s
      real,    intent(in) :: f
      integer, intent(in) :: irf        ! here dummy parameter only

      getrho_s = rho0
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(irf)

    endfunction getrho_s
!***********************************************************************
    subroutine getrho_1d(f,rho)

      real, dimension(mx), intent(in) :: f
      real, dimension(nx), intent(out):: rho

      rho=rho0
      call keep_compiler_quiet(f)

    endsubroutine getrho_1d
!***********************************************************************
    subroutine getlnrho_1d_x(f,lnrho)

      real, dimension(mx), intent(in) :: f
      real, dimension(nx), intent(out):: lnrho

      lnrho=lnrho0
      call keep_compiler_quiet(f)

    endsubroutine getlnrho_1d_x
!***********************************************************************
    subroutine getlnrho_1d_y(f,ix,lnrho)

      real, dimension(my), intent(in) :: f
      real, dimension(ny), intent(out):: lnrho
      integer,             intent(in) :: ix

      lnrho=lnrho0
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(ix)

    endsubroutine getlnrho_1d_y
!***********************************************************************
    subroutine getlnrho_2d(f,lnrho)

      real, dimension(:,:), intent(in) :: f
      real, dimension(:,:), intent(out):: lnrho

      lnrho=lnrho0
      call keep_compiler_quiet(f)

    endsubroutine getlnrho_2d
!***********************************************************************
    subroutine getlnrho_2dyz(f,ix,lnrho)

      real, dimension(my,mz), intent(in) :: f
      real, dimension(my,mz), intent(out):: lnrho
      integer,                intent(in) :: ix

      lnrho=lnrho0
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(ix)

    endsubroutine getlnrho_2dyz
!***********************************************************************
    subroutine getrho_2d(f,rho)

      real, dimension(:,:), intent(in) :: f
      real, dimension(:,:), intent(out):: rho

      rho=rho0
      call keep_compiler_quiet(f)

    endsubroutine getrho_2d
!***********************************************************************
    subroutine getrho_2dyz(f,ix,rho)

      real, dimension(my,mz), intent(in) :: f
      real, dimension(my,mz), intent(out):: rho
      integer,                intent(in) :: ix

      rho=rho0
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(ix)

    endsubroutine getrho_2dyz 
!***********************************************************************
    subroutine getdlnrho_x(f,rl,il,rho,dlnrho)

      integer,                   intent(in) :: rl,il
      real, dimension(mx,my,mz), intent(in) :: f
      real, dimension(my,mz),    intent(in) :: rho
      real, dimension(my,mz),    intent(out):: dlnrho

      dlnrho = 0.
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(rho)
      call keep_compiler_quiet(rl,il)

    endsubroutine getdlnrho_x
!***********************************************************************
    subroutine getdlnrho_y(f,rm,im,dlnrho)

      integer,                   intent(in) :: rm,im
      real, dimension(mx,my,mz), intent(in) :: f
      real, dimension(mx,mz),    intent(out):: dlnrho

      dlnrho = 0.
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(rm,im)

    endsubroutine getdlnrho_y
!***********************************************************************
    subroutine getdlnrho_z(f,rn,in,dlnrho)

      integer,                   intent(in) :: rn,in
      real, dimension(mx,my,mz), intent(in) :: f
      real, dimension(mx,my),    intent(out):: dlnrho

      dlnrho = 0.
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(rn,in)
!
    endsubroutine getdlnrho_z
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
    subroutine getderlnrho_z(f,iz,derlnrho)
!
!  Evaluates derlnrho as d_z ln(rho) for all x,y at z-position iz.
!
!  30-sep-16/MR: coded
!
      use Deriv, only: der

      integer,                             intent(in) :: iz
      real, dimension(:,:,:,:),            intent(in) :: f
      real, dimension(size(f,1),size(f,2)),intent(out):: derlnrho

      call fatal_error('getderlnrho_z', 'not implemented in nodensity.')
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(iz)
      call keep_compiler_quiet(derlnrho)
!
    endsubroutine getderlnrho_z
!***********************************************************************
endmodule DensityMethods
