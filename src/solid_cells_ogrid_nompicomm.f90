!***********************************************************************
!  FROM NOMPICOMM.F90
!***********************************************************************
!  ROUTINE

module Solid_Cells_Mpicomm

  use Cparam
  use General, only: keep_compiler_quiet

  implicit none

  include 'solid_cells_mpi.h'

  integer, parameter :: nx_ogrid=nxgrid_ogrid/nprocx,ny_ogrid=nygrid_ogrid/nprocy,nz_ogrid=nzgrid_ogrid/nprocz
  integer, parameter :: mx_ogrid=nx_ogrid+2*nghost,l1_ogrid=1+nghost,l2_ogrid=mx_ogrid-nghost
  integer, parameter :: my_ogrid=ny_ogrid+2*nghost,m1_ogrid=1+nghost,m2_ogrid=my_ogrid-nghost
  integer, parameter :: mz_ogrid=nz_ogrid+2*nghost,n1_ogrid=1+nghost,n2_ogrid=mz_ogrid-nghost

  contains 
!***********************************************************************
    subroutine initiate_isendrcv_bdry_ogrid(f)
!
!  For one processor, use periodic boundary conditions.
!  Dummy
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray) :: f
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
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray) :: f
!
      if (ALWAYS_FALSE) print*, f
!
    endsubroutine finalize_isendrcv_bdry_ogrid
!***********************************************************************
    subroutine isendrcv_bdry_x_ogrid(f)
!
!  Dummy
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray) :: f
!
      if (ALWAYS_FALSE) print *, f
!
    endsubroutine isendrcv_bdry_x_ogrid
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
    subroutine initialize_mpicomm_ogrid(a)
!
!  Dumy
!
      logical :: a

      call keep_compiler_quiet(a)

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
    subroutine tridag_parallel_x(a,b,c,r,u,n)
!
!  Dummy
!
      integer, intent(in) :: n
      real, dimension(n), intent(in) :: a,b,c,r
      real, dimension(n), intent(out) :: u

      call keep_compiler_quiet(a)
      call keep_compiler_quiet(b)
      call keep_compiler_quiet(c)
      call keep_compiler_quiet(r)
      call keep_compiler_quiet(u)
      call keep_compiler_quiet(n)

    endsubroutine tridag_parallel_x
!***********************************************************************
    subroutine tridag_parallel_y
    endsubroutine tridag_parallel_y
!***********************************************************************
    subroutine tridag_parallel_z
    endsubroutine tridag_parallel_z
!***********************************************************************
    subroutine finalize_isendrcv_bdry_filter(f_Hlox,f_Hupx,f_Hloy,f_Hupy,Hsize)
!
!  Dummy
!
      integer, intent(in) :: Hsize
      real, dimension (Hsize,my_ogrid,nz_ogrid,mvar) ::  f_Hlox,f_Hupx
      real, dimension (mx_ogrid,Hsize,nz_ogrid,mvar) ::  f_Hloy,f_Hupy
      intent(inout) :: f_Hlox,f_Hupx,f_Hloy,f_Hupy
      
      call keep_compiler_quiet(Hsize)
      call keep_compiler_quiet(f_Hlox)
      call keep_compiler_quiet(f_Hupx)
      call keep_compiler_quiet(f_Hloy)
      call keep_compiler_quiet(f_Hupy)
!
    endsubroutine finalize_isendrcv_bdry_filter
!***********************************************************************
    subroutine initiate_isendrcv_bdry_filter(f_og,Hsize)
!
!  Dummy
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray) :: f_og
      integer :: Hsize
      intent(in) :: f_og,Hsize

      call keep_compiler_quiet(f_og)
      call keep_compiler_quiet(Hsize)
    endsubroutine initiate_isendrcv_bdry_filter
!***********************************************************************
!***********************************************************************
!***********************************************************************


end module Solid_Cells_Mpicomm
