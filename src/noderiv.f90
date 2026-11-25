! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM integer, parameter :: nghost = 0
! nghost=0 creates trouble, as then mx=nx, creating conflicts in case constructs etc.
!
!***************************************************************
module Deriv
!
  use Cparam
  use Quiet
!
  implicit none
!
  private
!
  include 'deriv.h'
!
  contains
!
!***********************************************************************
    subroutine initialize_deriv
!
!  Dummy routine
!
    endsubroutine initialize_deriv
!***********************************************************************
    subroutine calc_coeffs_1( grid, coeffs )
!
!  dummy
! 
!      use Messages, only: fatal_error
      use Cdata,only: lroot
!
      !real, dimension(-2:3), intent(in ) :: grid
      !real, dimension(-3:3), intent(out) :: coeffs
      real, dimension(1:0), intent(in ) :: grid
      real, dimension(-0:0), intent(out) :: coeffs
!
      if (lroot) print*,'calc_coeffs_1 is not evaluated'
!--   call fatal_error("calc_coeffs_1","not coded for deriv_2nd")
      call keep_compiler_quiet(grid)
      call keep_compiler_quiet(coeffs)
! 
    endsubroutine calc_coeffs_1
!***********************************************************************
    subroutine der_main(f, k, df, j, ignoredx)
!
!  Dummy routine
!
      real, contiguous,dimension(:,:,:,:), intent(in) :: f
      real, dimension(nx), intent(in) :: df
      integer, intent(in) :: j, k
      logical, intent(in), optional :: ignoredx

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(k,j)
      call keep_compiler_quiet(ignoredx)
!
    endsubroutine der_main
!***********************************************************************
    subroutine der_x(f,df)
!
!  Dummy routine
!
      real, dimension (mx), intent(in)  :: f
      real, dimension (nx), intent(in) :: df

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
!
    endsubroutine der_x
!***********************************************************************
    subroutine der2_x(f,df2)
!
!  Dummy routine
!
      real, dimension (mx), intent(in)  :: f
      real, dimension (nx), intent(in) :: df2

      call keep_compiler_quiet(f,df2)
!
    endsubroutine der2_x
!***********************************************************************
    subroutine der_z(f,df)
!
!  Dummy routine
!
      real, dimension (mz), intent(in)  :: f
      real, dimension (nz), intent(in) :: df

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
!
    endsubroutine der_z
!***********************************************************************
    subroutine der2_z(f,df2)
!
!  Dummy routine
!
      real, dimension (mz), intent(in)  :: f
      real, dimension (nz), intent(in) :: df2

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df2)
!
    endsubroutine der2_z
!***********************************************************************
    subroutine der_other(f,df,j)
!
!  Dummy routine
!
      real, dimension (mx,my,mz) :: f
      real, dimension (:) :: df
      integer :: j

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(j)
!
    endsubroutine der_other
!***********************************************************************
    subroutine der_pencil(j,pencil,df)
!
!  Dummy routine
!
      real, dimension (:) :: pencil,df
      integer :: j
      intent(in)  :: df, j, pencil

      call keep_compiler_quiet(j)
      call keep_compiler_quiet(pencil)
      call keep_compiler_quiet(df)
!
    endsubroutine der_pencil
!***********************************************************************
  subroutine distr_der(arr,idir,der,order)
!
!  Dummy routine
!
    real, dimension(:,:) :: arr, der
    integer :: idir
    integer, optional :: order

     call keep_compiler_quiet(arr)
     call keep_compiler_quiet(idir)
     call keep_compiler_quiet(der)
     call keep_compiler_quiet(order)
!
  endsubroutine distr_der
!***********************************************************************
    subroutine der2_main(f,k,df2,j,lwo_line_elem)
!
!  Dummy routine
!
      real, contiguous,dimension(:,:,:,:) :: f
      real, dimension (nx) :: df2
      integer :: j,k
      logical, optional :: lwo_line_elem
!
      intent(in) :: f,df2,k,j,lwo_line_elem

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(k)
      call keep_compiler_quiet(df2)
      call keep_compiler_quiet(j)
      call keep_compiler_quiet(lwo_line_elem)
!
    endsubroutine der2_main
!***********************************************************************
    subroutine der2_other(f,df2,j)
!
!  Dummy routine
!
      real, dimension (mx,my,mz) :: f
      real, dimension (nx) :: df2
      integer :: j
!
      intent(in) :: f,df2,j
      
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df2)
      call keep_compiler_quiet(j)
!
    endsubroutine der2_other
!***********************************************************************
    subroutine der2_pencil(j,pencil,df2)
!
!  Dummy routine
!
      real, dimension (:) :: pencil,df2
      integer :: j
      intent(in)  :: j, pencil, df2

      call keep_compiler_quiet(j)
      call keep_compiler_quiet(pencil)
      call keep_compiler_quiet(df2)
!
    endsubroutine der2_pencil
!***********************************************************************
    subroutine der3(f,k,df,j,ignoredx)
!
!  Dummy routine
!
      real, contiguous,dimension(:,:,:,:) :: f
      real, dimension (nx) :: df
      integer :: j,k
      logical, optional :: ignoredx
      intent(in)  :: f,k,df,j,ignoredx

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(k)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(j)
      call keep_compiler_quiet(ignoredx)
!
    endsubroutine der3
!***********************************************************************
    subroutine der4(f,k,df,j,ignoredx,upwind)
!
!  Dummy routine
!
      real, contiguous,dimension(:,:,:,:) :: f
      real, dimension (nx) :: df
      integer :: j,k
      logical, optional :: ignoredx,upwind
      intent(in)  :: f,k,j,ignoredx,upwind

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(k)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(j)
      call keep_compiler_quiet(ignoredx)
      call keep_compiler_quiet(upwind)
!
    endsubroutine der4
!***********************************************************************
    subroutine der5(f,k,df,j,ignoredx)
!
!  Dummy routine
!
      real, contiguous,dimension(:,:,:,:) :: f
      real, dimension (nx) :: df
      integer :: j,k
      logical, optional :: ignoredx
      intent(in)  :: f,k,df,j,ignoredx

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(k)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(j)
      call keep_compiler_quiet(ignoredx)
!
    endsubroutine der5
!***********************************************************************
    subroutine der6_main(f,k,df,j,ignoredx,upwind)
!
!  Dummy routine
!
      real, contiguous,dimension(:,:,:,:) :: f
      real, dimension (nx) :: df
      integer :: j,k
      logical, optional :: ignoredx,upwind
      intent(in)  :: f,k,df,j,ignoredx,upwind

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(k)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(j)
      call keep_compiler_quiet(ignoredx)
      call keep_compiler_quiet(upwind)
!
    endsubroutine der6_main
!***********************************************************************
    subroutine der2_minmod(f,j,delfk,delfkp1,delfkm1,k)
!
!  Dummy routine
!
      intent(in) :: f,k,j
      intent(in) :: delfk,delfkp1,delfkm1
      real, contiguous,dimension(:,:,:,:) :: f
      real, dimension (nx) :: delfk,delfkp1,delfkm1
      integer :: j,k
      
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(j)
      call keep_compiler_quiet(delfk)
      call keep_compiler_quiet(delfkp1)
      call keep_compiler_quiet(delfkm1)
      call keep_compiler_quiet(k)
!
    endsubroutine der2_minmod
!***********************************************************************
!   real function minmod(a,b,c)
!
!  Dummy routine
!
!     real :: a,b,c
!
!   endfunction minmod
!***********************************************************************
    subroutine der6_other(f,df,j,ignoredx,upwind)
!
!  Dummy routine
!
      real, dimension (mx,my,mz) :: f
      real, dimension (nx) :: df
      integer :: j
      logical, optional :: ignoredx,upwind
      intent(in)  :: f,df,j,ignoredx,upwind

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(j)
      call keep_compiler_quiet(ignoredx)
      call keep_compiler_quiet(upwind)
!
    endsubroutine der6_other
!***********************************************************************
    subroutine der6_pencil(j,pencil,df6,ignoredx,upwind)
!
!  Calculate 6th derivative of any x, y or z pencil.
!
!  20-jul-20/wlyra: adapted from der2_pencil
!
      real, dimension (:) :: pencil,df6
      integer :: j
      logical, optional :: ignoredx,upwind

      call keep_compiler_quiet(j)
      call keep_compiler_quiet(pencil)
      call keep_compiler_quiet(df6)
      call keep_compiler_quiet(ignoredx)
      call keep_compiler_quiet(upwind)
!
    endsubroutine der6_pencil
!***********************************************************************
    real function der5_single(f,j,dc1)
!
!  Dummy routine
!
      real, dimension(:),  intent(in) :: f, dc1
      integer           ,  intent(in) :: j
!
      der5_single=0.

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(j)
      call keep_compiler_quiet(dc1)
!
    endfunction der5_single
!***********************************************************************
    subroutine derij_main(f,k,df,i,j,lwo_line_elem)
!
!  Dummy routine
!
      real, contiguous,dimension(:,:,:,:) :: f
      real, dimension (nx) :: df
      integer :: i,j,k
      logical, optional :: lwo_line_elem
      intent(in) :: f,k,df,i,j

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(k)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(i)
      call keep_compiler_quiet(j)
      call keep_compiler_quiet(lwo_line_elem)
!
    endsubroutine derij_main
!***********************************************************************
    subroutine derij_other(f,df,i,j)
!
!  Dummy routine
!
      real, dimension (mx,my,mz) :: f
      real, dimension (nx) :: df
      integer :: i,j
      intent(in)  :: f,df,i,j

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(i)
      call keep_compiler_quiet(j)
!
    endsubroutine derij_other
!***********************************************************************
    subroutine der5i1j(f,k,df,i,j)
!
!  Dummy routine
!
      real, contiguous,dimension(:,:,:,:) :: f
      real, dimension (nx) :: df
      integer :: i,j,k
      intent(in) :: f,k,df,i,j

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(k)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(i)
      call keep_compiler_quiet(j)
!
    endsubroutine der5i1j
!***********************************************************************
    subroutine der4i2j(f,k,df,i,j)
!
!  Dummy routine
!
      real, contiguous,dimension(:,:,:,:) :: f
      real, dimension (nx) :: df
      integer :: i,j,k
      intent(in) :: f,k,df,i,j

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(k)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(i)
      call keep_compiler_quiet(j)
!
    endsubroutine der4i2j
!***********************************************************************
    subroutine der2i2j2k(f,k,df)
!      
!  Dummy routine
!
      real, contiguous,dimension(:,:,:,:),intent(in) :: f
      integer,intent(in) :: k
      real, dimension(nx), intent(in) :: df

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(k)
      call keep_compiler_quiet(df)
!
    endsubroutine der2i2j2k
!***********************************************************************
    subroutine der3i3j(f,k,df,i,j)
!
!  Dummy routine
!
      real, contiguous,dimension(:,:,:,:) :: f
      real, dimension (nx), intent(in) :: df
      integer, intent(in) :: k,i,j

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(k)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(i)
      call keep_compiler_quiet(j)
!
    endsubroutine der3i3j
!***********************************************************************          
    subroutine der3i2j1k(f,ik,df,i,j,k)
!
!  Dummy routine
!
      real, contiguous,dimension(:,:,:,:) :: f
      real, dimension (nx), intent(in) :: df
      integer, intent(in) :: ik,i,j,k
      
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(ik)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(i)
      call keep_compiler_quiet(j)
      call keep_compiler_quiet(k)
!
    endsubroutine der3i2j1k
!***********************************************************************
    subroutine der4i1j1k(f,ik,df,i,j,k)
!
!  Dummy routine
!
      real, contiguous,dimension(:,:,:,:) :: f
      real, dimension (nx), intent(in) :: df
      integer, intent(in) :: ik,i,j,k

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(ik)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(i)
      call keep_compiler_quiet(j)
      call keep_compiler_quiet(k)
!
    endsubroutine der4i1j1k
!***********************************************************************
    subroutine der_upwind1st(f,uu,k,df,j)
!
!  Dummy routine
!
      real, contiguous,dimension(:,:,:,:) :: f
      real, dimension (nx,3) :: uu
      real, dimension (nx) :: df
      integer :: j,k
      intent(in)  :: f,uu,k,j
      intent(in) :: df

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(uu)
      call keep_compiler_quiet(k)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(j)
!
    endsubroutine der_upwind1st
!***********************************************************************
    subroutine der_onesided_4_slice_main(f,sgn,k,df,pos,j)
!
!  Dummy routine
!
      real, contiguous,dimension(:,:,:,:) :: f
      real, dimension (:,:) :: df
      integer :: pos,k,sgn,j
      intent(in)  :: f,k,pos,sgn,j
      intent(in) :: df

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(sgn)
      call keep_compiler_quiet(k)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(pos)
      call keep_compiler_quiet(j)
!
    endsubroutine der_onesided_4_slice_main
!***********************************************************************
    subroutine der_onesided_4_slice_main_pt(f,sgn,k,df,lll,mmm,nnn,j)
!
!  Dummy routine
!
      real, contiguous,dimension(:,:,:,:) :: f
      real  :: df
      integer :: lll,mmm,nnn,k,sgn,j
      intent(in)  :: f,k,df,lll,mmm,nnn,sgn,j

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(sgn)
      call keep_compiler_quiet(k)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(lll)
      call keep_compiler_quiet(mmm)
      call keep_compiler_quiet(nnn)
      call keep_compiler_quiet(j)
!
    endsubroutine der_onesided_4_slice_main_pt
!***********************************************************************
    subroutine der_onesided_4_slice_other(f,sgn,df,pos,j)
!
!  Dummy routine
!
      real, dimension (mx,my,mz) :: f
      real, dimension (:,:) :: df
      integer :: pos,sgn,j
      intent(in)  :: f,df,pos,sgn,j

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(sgn)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(pos)
      call keep_compiler_quiet(j)
!
    endsubroutine der_onesided_4_slice_other
!***********************************************************************
    subroutine der_onesided_4_slice_other_pt(f,sgn,df,lll,mmm,nnn,j)
!
!  Dummy routine
!
      real, dimension (mx,my,mz) :: f
      real :: df
      integer :: lll,mmm,nnn,sgn,j
      intent(in)  :: f,df,lll,mmm,nnn,sgn,j

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(sgn)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(lll)
      call keep_compiler_quiet(mmm)
      call keep_compiler_quiet(nnn)
      call keep_compiler_quiet(j)
!
    endsubroutine der_onesided_4_slice_other_pt
!***********************************************************************
    subroutine finalize_deriv
!
!  Dummy
!
    endsubroutine finalize_deriv
!***********************************************************************
    subroutine deri_3d_inds(f,df,inds,j,lignored,lnometric)
!
!  dummy routine for compatibility
!
      real, dimension (mx,my,mz)          :: f
      real, dimension (nx)                :: df
      integer                             :: j
      logical,                   optional :: lignored, lnometric
      integer, dimension(nx)              :: inds
!
      intent(in)  :: f,df,j,inds,lignored,lnometric

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(inds)
      call keep_compiler_quiet(j)
      call keep_compiler_quiet(lignored)
      call keep_compiler_quiet(lnometric)
!
    endsubroutine deri_3d_inds
!************************************************************************
    logical function heatflux_deriv_x(f, inh, fac, topbot)
!
!   dummy routine
!
      real, contiguous,dimension(:,:,:,:), intent(IN):: f
      real, dimension(my,mz)           , intent(IN):: inh
      real                             , intent(IN):: fac
      integer                          , intent(IN):: topbot
!
      heatflux_deriv_x = .true.

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(inh)
      call keep_compiler_quiet(fac)
      call keep_compiler_quiet(topbot)
!
    endfunction heatflux_deriv_x
!***********************************************************************
    subroutine bval_from_neumann_scl(f,topbot,j,idir,val)
!
!  Dummy routine
!
      real, contiguous, dimension(:,:,:,:) :: f
      integer, intent(IN) :: topbot
      integer :: j,idir
      real :: val

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
      call keep_compiler_quiet(j)
      call keep_compiler_quiet(idir)
      call keep_compiler_quiet(val)
!
    endsubroutine bval_from_neumann_scl
!***********************************************************************
    subroutine bval_from_3rd_scl(f,topbot,j,idir,val)
!
!  Dummy routine
!
      real, contiguous, dimension(:,:,:,:) :: f
      integer, intent(IN) :: topbot
      integer :: j,idir
      real :: val
      
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
      call keep_compiler_quiet(j)
      call keep_compiler_quiet(idir)
      call keep_compiler_quiet(val)
!
    endsubroutine bval_from_3rd_scl
!***********************************************************************
    subroutine bval_from_4th_scl(f,topbot,j,idir,val)
!
!  Dummy routine
!
      real, contiguous, dimension(:,:,:,:) :: f
      integer, intent(IN) :: topbot
      integer :: j,idir
      real :: val

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
      call keep_compiler_quiet(j)
      call keep_compiler_quiet(idir)
      call keep_compiler_quiet(val)
!
    endsubroutine bval_from_4th_scl
!***********************************************************************
    subroutine set_ghosts_for_onesided_ders(f,topbot,j,idir,l2nd)
!
!  Dummy routine
!
      real, contiguous, dimension(:,:,:,:) :: f
      integer, intent(IN) :: topbot
      integer :: j,idir
      logical, optional :: l2nd

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
      call keep_compiler_quiet(j)
      call keep_compiler_quiet(idir)
      call keep_compiler_quiet(l2nd)
!
    endsubroutine set_ghosts_for_onesided_ders
!***********************************************************************
    subroutine bval_from_neumann_arr(f,topbot,j,idir,val)
!
!  Dummy routine
!
      real, contiguous, dimension(:,:,:,:) :: f
      integer, intent(IN) :: topbot
      integer :: j,idir
      real, dimension(:,:) :: val

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
      call keep_compiler_quiet(j)
      call keep_compiler_quiet(idir)
      call keep_compiler_quiet(val)
!
    endsubroutine bval_from_neumann_arr
!***********************************************************************
    subroutine bval_from_3rd_arr(f,topbot,j,idir,val)
!
!  Dummy routine
!
      real, contiguous, dimension(:,:,:,:) :: f
      integer, intent(IN) :: topbot
      integer :: j,idir
      real, dimension(:,:) :: val

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
      call keep_compiler_quiet(j)
      call keep_compiler_quiet(idir)
      call keep_compiler_quiet(val)
!
    endsubroutine bval_from_3rd_arr
!***********************************************************************
    subroutine bval_from_4th_arr(f,topbot,j,idir,val)
!
!  Dummy routine
!
      real, contiguous, dimension(:,:,:,:) :: f
      integer, intent(IN) :: topbot
      integer :: j,idir
      real, dimension(:,:) :: val

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
      call keep_compiler_quiet(j)
      call keep_compiler_quiet(idir)
      call keep_compiler_quiet(val)
!
    endsubroutine bval_from_4th_arr
!***********************************************************************
 endmodule Deriv
