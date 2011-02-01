! $Id$
!
!  This module calculates the chemistry contribution to df in the case
!  where chemistry is solved separately using the LSODE solver
!
module LsodeForChemistry
!
  use Cdata
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  private
!
  public :: pde_chemistry, lsode_fc
!
  contains
!***********************************************************************
    include 'pencil_init.inc' ! defines subroutine initialize_pencils()
!***********************************************************************
    subroutine pde_chemistry(f,df,p)
!
!  29-nov-10/julien: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine pde_chemistry
!***********************************************************************
   subroutine lsode_fc (t1,t2,f,df)
! 
!  29-nov-10/Julien: coded, interface for the LSODE routine   
!   
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df  
      real :: t1, t2
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(t1,t2)
!
   end subroutine lsode_fc
!***********************************************************************
endmodule LsodeForChemistry
