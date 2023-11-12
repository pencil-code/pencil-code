subroutine copyin_farray_alloc
!$omp parallel copyin( &
!$omp f  &
!$omp ,df  &
!$omp )
!$omp end parallel
end subroutine copyin_farray_alloc
