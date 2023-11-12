subroutine copyin_dustvelocity
!$omp parallel copyin( &
!$omp scolld  &
!$omp ,md  &
!$omp ,mi  &
!$omp )
!$omp end parallel
end subroutine copyin_dustvelocity
