module Solid_Cells_Mpicomm

  use Cparam

  implicit none

  include 'solid_cells_mpi.h'

  contains 
!***********************************************************************
    subroutine initiate_isendrcv_bdry(f)
!
!  For one processor, use periodic boundary conditions.
!  Dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
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
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (ALWAYS_FALSE) print*, f
!
    endsubroutine finalize_isendrcv_bdry
!***********************************************************************
    subroutine initiate_isendrcv_bdry_ogrid(f)
!
!  For one processor, use periodic boundary conditions.
!  Dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (ALWAYS_FALSE) print*, f
!
    endsubroutine initiate_isendrcv_bdry_ogrid
!***********************************************************************
    subroutine finalize_isendrcv_bdry_ogrid(f)
!
!  Apply boundary conditions.
!  Dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (ALWAYS_FALSE) print*, f
!
    endsubroutine finalize_isendrcv_bdry_ogrid
!***********************************************************************
    subroutine isendrcv_bdry_x(f)
!
!  Dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (ALWAYS_FALSE) print *, f
!
    endsubroutine isendrcv_bdry_x
!***********************************************************************
end module Solid_Cells_Mpicomm
