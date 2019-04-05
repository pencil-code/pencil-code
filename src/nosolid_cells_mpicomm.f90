module Solid_Cells_Mpicomm

  use Cparam

  implicit none

  include 'solid_cells_mpi.h'

  contains 
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
    subroutine finalize_isend_init_interpol(ireq1D,ireq2D,nreq1D,nreq2D)
!
!  Dummy
!
      integer :: nreq1D, nreq2D
      integer, dimension(nreq1D) :: ireq1D
      integer, dimension(nreq2D) :: ireq2D
!
      if (ALWAYS_FALSE) print *, nreq1D,nreq2D,ireq1D,ireq2D
!
    endsubroutine finalize_isend_init_interpol
!***********************************************************************
    subroutine initialize_mpicomm_ogrid(lf)
      logical, intent(in) :: lf
!
!  Dummy
!
      if(ALWAYS_FALSE) print*, lf
!
    endsubroutine initialize_mpicomm_ogrid
!***********************************************************************
    subroutine cyclic_parallel_y(a,b,c,alpha,beta,r,x,n)
!
!  Dummy
!
      integer, intent(in) :: n
      real, dimension(n) :: a,b,c,r
      real, dimension(n) :: x
      real :: alpha,beta

      if(ALWAYS_FALSE) print*, a,b,c,alpha,beta,r,x,n

    endsubroutine cyclic_parallel_y
!***********************************************************************
    subroutine initiate_isendrcv_bdry_filter(f_og,Hsize)
!
      real, dimension (:,:,:,:) ::  f_og
      integer :: Hsize
!
      if(ALWAYS_FALSE) print*, f_og,Hsize

    endsubroutine initiate_isendrcv_bdry_filter
!***********************************************************************
    subroutine finalize_isendrcv_bdry_filter(f_Hlox,f_Hupx,f_Hloy,f_Hupy,Hsize)
!
      integer, intent(in) :: Hsize
      real, dimension (:,:,:,:) ::  f_Hlox,f_Hupx,f_Hloy,f_Hupy

      if(ALWAYS_FALSE) print*, f_Hlox,f_Hupx,f_Hloy,f_Hupy,Hsize
    endsubroutine finalize_isendrcv_bdry_filter
!***********************************************************************
    subroutine tridag_parallel_x(a,b,c,r,u,n)
!
!
      integer, intent(in) :: n
      real, dimension(n) :: a,b,c,r,u
!
      if(ALWAYS_FALSE) print*, n,a,b,c,r,u
    endsubroutine tridag_parallel_x
!***********************************************************************
end module Solid_Cells_Mpicomm
