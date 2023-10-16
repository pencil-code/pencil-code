subroutine copyin_magnetic
!$omp parallel copyin( &
!$omp nbvec  &
!$omp ,eta_total  &
!$omp ,uxbb  &
!$omp )
!$omp end parallel
end subroutine copyin_magnetic
