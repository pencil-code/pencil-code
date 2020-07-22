module DensityMethods
!
!  Dummy module for density_stratified, with which all density queries
!  should be done in EquationOfState.
!
  use Cparam
  use General, only: keep_compiler_quiet
  use Messages, only: fatal_error
!
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
!
  interface putlnrho
    module procedure putlnrho_s
    module procedure putlnrho_v
  endinterface

 !real, pointer :: rho0, lnrho0

  contains
!***********************************************************************
    subroutine initialize_density_methods
!
!  03-apr-15/ccyang: dummy.
!
      use SharedVariables, only: get_shared_variable
      use Messages, only: warning
      use Cdata, only: lstratz
!
      if (lstratz) &
        call warning('initialize_density_methods', &
                     'density methods not yet implemented for density_stratified')
!
      !call get_shared_variable('rho0',rho0,caller='initialize_density_methods')
      !call get_shared_variable('lnrho0',lnrho0)
!
    endsubroutine initialize_density_methods
!***********************************************************************
    real function getrho_s(f, irf)
!
!  03-apr-15/ccyang: dummy.
!
      real, intent(in) :: f
      integer, intent(in) :: irf        ! here dummy parameter only
!
      call fatal_error('getrho_s', 'not implemented.')
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(irf)
!
    endfunction getrho_s
!***********************************************************************
    subroutine getrho_1d(f, rho)
!
!  03-apr-15/ccyang: dummy.
!
      real, dimension(mx), intent(in) :: f
      real, dimension(nx), intent(out):: rho
!
      call fatal_error('getrho_1d', 'not implemented.')
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(rho)
!
    endsubroutine getrho_1d
!***********************************************************************
    subroutine getrho1_1d(f,rho1)
!
!  Fetches inverse of density.
!
!   4-oct-17/MR: derived from getrho_1d.
!
      real, dimension(mx), intent(in) :: f
      real, dimension(nx), intent(out):: rho1
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(rho1)
!
    endsubroutine getrho1_1d
!***********************************************************************           
    subroutine getlnrho_1d_x(f, lnrho)
!
!  03-apr-15/ccyang: dummy.
!
      real, dimension(mx), intent(in) :: f
      real, dimension(nx), intent(out):: lnrho
!
      call fatal_error('getlnrho_1d_x', 'not implemented.')
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lnrho)
!
    endsubroutine getlnrho_1d_x
!***********************************************************************
    subroutine getlnrho_1d_y(f, ix, lnrho)
!
!  03-apr-15/ccyang: dummy.
!
      real, dimension(my), intent(in) :: f
      real, dimension(ny), intent(out):: lnrho
      integer, intent(in) :: ix
!
      call fatal_error('getlnrho_1d_y', 'not implemented.')
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(ix)
      call keep_compiler_quiet(lnrho)
!
    endsubroutine getlnrho_1d_y
!***********************************************************************
    subroutine getlnrho_2d(f, lnrho)
!
!  03-apr-15/ccyang: dummy.
!
      real, dimension(:,:), intent(in) :: f
      real, dimension(:,:), intent(out):: lnrho
!
      call fatal_error('getlnrho_2d', 'not implemented.')
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lnrho)
!
    endsubroutine getlnrho_2d
!***********************************************************************
    subroutine getlnrho_2dyz(f, ix, lnrho)
!
!  03-apr-15/ccyang: dummy.
!
      real, dimension(my,mz), intent(in) :: f
      real, dimension(my,mz), intent(out):: lnrho
      integer, intent(in) :: ix
!
      call fatal_error('getlnrho_2dyz', 'not implemented.')
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(ix)
      call keep_compiler_quiet(lnrho)
!
    endsubroutine getlnrho_2dyz
!***********************************************************************
    subroutine getrho_2d(f, rho)
!
!  03-apr-15/ccyang: dummy.
!
      real, dimension(:,:), intent(in) :: f
      real, dimension(:,:), intent(out):: rho
!
      call fatal_error('getrho_2d', 'not implemented.')
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(rho)
!
    endsubroutine getrho_2d
!***********************************************************************
    subroutine getrho_2dyz(f, ix, rho)
!
!  03-apr-15/ccyang: dummy.
!
      real, dimension(my,mz), intent(in) :: f
      real, dimension(my,mz), intent(out):: rho
      integer, intent(in) :: ix
!
      call fatal_error('getrho_2dyz', 'not implemented.')
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(ix)
      call keep_compiler_quiet(rho)
!
    endsubroutine getrho_2dyz 
!***********************************************************************
    subroutine getdlnrho_x(f, rl, il, rho, dlnrho)
!
!  03-apr-15/ccyang: dummy.
!
      integer, intent(in) :: rl, il
      real, dimension(mx,my,mz), intent(in) :: f
      real, dimension(my,mz), intent(in) :: rho
      real, dimension(my,mz), intent(out):: dlnrho
!
      call fatal_error('getdlnrho_x', 'not implemented.')
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(rl,il)
      call keep_compiler_quiet(rho)
      call keep_compiler_quiet(dlnrho)
!
    endsubroutine getdlnrho_x
!***********************************************************************
    subroutine getdlnrho_y(f, rm, im, dlnrho)
!
!  03-apr-15/ccyang: dummy.

      integer, intent(in) :: rm, im
      real, dimension(mx,my,mz), intent(in) :: f
      real, dimension(mx,mz), intent(out):: dlnrho
!
      call fatal_error('getdlnrho_y', 'not implemented.')
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(rm,im)
      call keep_compiler_quiet(dlnrho)
!
    endsubroutine getdlnrho_y
!***********************************************************************
    subroutine getdlnrho_z(f, rn, in, dlnrho)
!
!  03-apr-15/ccyang: dummy.
!
      integer, intent(in) :: rn, in
      real, dimension(mx,my,mz), intent(in) :: f
      real, dimension(mx,my), intent(out):: dlnrho
!
      call fatal_error('getdlnrho_z', 'not implemented.')
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(rn,in)
      call keep_compiler_quiet(dlnrho)
!
    endsubroutine getdlnrho_z
!***********************************************************************
    subroutine getderlnrho_z(f,iz,derlnrho)
!
!  16-jul-20/wlyra: dummy.
!
      integer,                             intent(in) :: iz
      real, dimension(:,:,:,:),            intent(in) :: f
      real, dimension(size(f,1),size(f,2)),intent(out):: derlnrho
      !integer, intent(in) :: rn, in
      !real, dimension(mx,my,mz), intent(in) :: f
      !real, dimension(mx,my), intent(out):: dlnrho
!
      call fatal_error('getderlnrho_z', 'not implemented.')
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(iz)
      call keep_compiler_quiet(derlnrho)
!
    endsubroutine getderlnrho_z
!***********************************************************************
   subroutine putrho(f,rho)

      real, dimension(mx), intent(out):: f
      real, dimension(nx), intent(in) :: rho

      call fatal_error('putrho', 'not implemented.')
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(rho)
!
    endsubroutine putrho
!***********************************************************************
    subroutine putlnrho(f,lnrho)

      real, dimension(mx), intent(out):: f
      real, dimension(nx), intent(in) :: lnrho

      call fatal_error('putlnrho', 'not implemented.')
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lnrho)
!
    endsubroutine putlnrho
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
