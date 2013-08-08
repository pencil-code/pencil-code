! This is a pseudo code that describes the structure of a possible
! implementation of a GPU pencil code.
!
! f, df, are all "GPU" memory now (on-board, not in chip).  They need
! to be copied to "CPU" memory for I/O, initial conditions, etc
!
! p is still the temporary array but it is 3D with dimension
! 
!   (pnx + 6, pny + 6, 7)
!
! where pnx and pny are limited by the size of shared memory.  On
! Tesla it is given by
!
!   (pnx + 6) * (pny + 6) * 7 * number_of_variable * 
!     sizeof(float) < 48K
!=====================================================================
subroutine time_step(f,df) ! Because p is in the shared memory, it is
                           ! not known outside the kernel.

  ! pnx * pny is the number of threads that we will run in parallel.
  ! For good performance, it should be larger than 256.
  dimBlock = dim3(pnx, pny, 1)
  ! nx, ny, and nz are the actual dimension of the computational grid.
  ! pnx and pny are usually smaller than nx and ny.  Therefore, we
  ! need to launch many blocks to cover the whole domain.
  dimGrid  = dim3((nx-1) / pnx + 1, (nx-1) / pny + 1, 1 )

  ...
  ! Kernel launching
  call pde<<<dimGrid, dimBlcok>>>(f,df)
  ...

endsubroutine time_step

!=====================================================================
subroutine pde(f,df)

  intent(inout)  :: f
  intent(out)    :: df
  real,device,dimension(pnx+6,pny+6,7),shared :: p
  real,device :: l%ugu,...

  ! Calculate the global array indices according to the thread and
  ! block indices.
  i = blockIdx%y * blockDim%y + threadIdx%y
  j = blockIdx%z * blockDim%z + threadIdx%z

  ! Preload global memory into shared memory (cache); extra codes are
  ! needed to load the ghost cells
  p%uu = f(i,j,-2:3,iux:iuz)
  ...

  ! Start rolling the cache down along the z-direction
  n_loop: do n = 1, ...

    ! Load one layer of data; extra codes are needed to load the ghost
    ! cells
    p%uu=f(i,j,n+3,iux:iuz)
    ...

    ! Compute derivatives using the shared memory; local variables are
    ! in register or local memory
    l%ugu=...
    ...

    ! Update output values
    df(i,j,n,iux:iuz)=df(i,j,n,iux:iuz)-p%ugu
    ...

  enddo n_loop

endsubroutine pde
