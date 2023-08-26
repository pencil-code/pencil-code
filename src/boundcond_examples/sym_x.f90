!
!  Example of offloaded boundcond --- Symmetric boundcond
!  p_farray_in is a pointer to fields stored on GPU through Astaroth
!
!***********************************************************************
    subroutine bc_sym_x(f,ivar1_opt,ivar2_opt)
!
!  Envelope for being called from C code.
!
!$    use omp_lib
      use GPU
      use Cparam
      real, dimension (mx,my,mz,mfarray) :: f
      integer, optional :: ivar1_opt, ivar2_opt
      integer :: i, j
      integer :: num_devices, nteams, nthreads
      logical :: initial_device
      integer :: ivar1, ivar2
      ivar1=1; ivar2=min(mcom,size(f,4))
  !$omp target has_device_addr(p_farray_in) map(ivar1,ivar2) 
  !$omp distribute teams parallel do collapse(2)
      do j=ivar1,ivar2
        do i=1,3
            p_farray_in(l2+i,:,:,j) = p_farray_in(l2-i,:,:,j)
            p_farray_in(l1-i,:,:,j)  = p_farray_in(l1+i,:,:,j)
        enddo
      enddo
  !$omp end distribute teams parallel do
  !$omp end target

    endsubroutine bc_sym_x 
!***********************************************************************