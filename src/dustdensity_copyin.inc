subroutine copyin_dustdensity
!$omp parallel copyin( &
!$omp dkern  &
!$omp ,lcalcdkern  &
!$omp )
!$omp end parallel
end subroutine copyin_dustdensity
