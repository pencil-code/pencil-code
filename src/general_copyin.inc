subroutine copyin_general
!$omp parallel copyin( &
!$omp rstate  &
!$omp ,rstate2  &
!$omp )
!$omp end parallel
end subroutine copyin_general
