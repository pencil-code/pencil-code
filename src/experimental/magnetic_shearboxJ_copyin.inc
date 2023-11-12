subroutine copyin_magnetic
!$omp parallel copyin( &
!$omp uxbb  &
!$omp ,nbvec  &
!$omp )
!$omp end parallel
end subroutine copyin_magnetic
