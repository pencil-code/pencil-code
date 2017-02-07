!***********************************************************************
!  FROM NOMPICOMM.F90
!***********************************************************************
!  ROUTINE
!    initiate_isendrcv_bdry_ogrid
!    finalize_isendrcv_bdry_ogrid
!    isendrcv_bdry_x_ogrid

module Solid_Cells_Mpicomm

  use Cparam, only: ALWAYS_FALSE

  implicit none

  contains 
!***********************************************************************
    subroutine initiate_isendrcv_bdry(f)
!
!  For one processor, use periodic boundary conditions.
!  Dummy
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid) :: f
!
      if (ALWAYS_FALSE) print*, f
!
    endsubroutine initiate_isendrcv_bdry
!***********************************************************************
    subroutine finalize_isendrcv_bdry(f)
!
!  Apply boundary conditions.
!  Dummy
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid) :: f
!
      if (ALWAYS_FALSE) print*, f
!
    endsubroutine finalize_isendrcv_bdry
!***********************************************************************
    subroutine isendrcv_bdry_x(f)
!
!  Dummy
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid) :: f
!
      if (ALWAYS_FALSE) print *, f
!
    endsubroutine isendrcv_bdry_x
!***********************************************************************
end module Solid_Cells_Mpicomm
