module DensityMethods
!
!  Dummy module for density_stratified, with which all density queries
!  should be done in EquationOfState.
!
  use Cparam
  use General, only: keep_compiler_quiet
  use Messages
!
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
!
  contains
!***********************************************************************
    subroutine initialize_density_methods
!
!  03-apr-15/ccyang: dummy.
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
    subroutine getdlnrho_z(f, in, dlnrho)
!
!  03-apr-15/ccyang: dummy.
!
      integer, intent(in) :: in
      real, dimension(mx,my,-in:in), intent(in) :: f
      real, dimension(mx,my), intent(out):: dlnrho
!
      call fatal_error('getdlnrho_z', 'not implemented.')
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(in)
      call keep_compiler_quiet(dlnrho)
!
    endsubroutine getdlnrho_z
!***********************************************************************
    subroutine getdlnrho_x(f, il, ix, rho, dlnrho)
!
!  03-apr-15/ccyang: dummy.
!
      integer, intent(in) :: il, ix
      real, dimension(-il:il,my,mz), intent(in) :: f
      real, dimension(my,mz), intent(in) :: rho
      real, dimension(my,mz), intent(out):: dlnrho
!
      call fatal_error('getdlnrho_x', 'not implemented.')
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(il)
      call keep_compiler_quiet(ix)
      call keep_compiler_quiet(rho)
      call keep_compiler_quiet(dlnrho)
!
    endsubroutine getdlnrho_x
!***********************************************************************
    subroutine getdlnrho_y(f, im, dlnrho)
!
!  03-apr-15/ccyang: dummy.

      integer, intent(in) :: im
      real, dimension(mx,-im:im,mz), intent(in) :: f
      real, dimension(mx,mz), intent(out):: dlnrho
!
      call fatal_error('getdlnrho_y', 'not implemented.')
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(im)
      call keep_compiler_quiet(dlnrho)
!
    endsubroutine getdlnrho_y
!***********************************************************************
endmodule DensityMethods
