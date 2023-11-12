subroutine copyin_special
!$omp parallel copyin( &
!$omp fcoriolis  &
!$omp ,fgravity  &
!$omp )
!$omp end parallel
end subroutine copyin_special
