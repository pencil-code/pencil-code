subroutine copyin_weno_transport
!$omp parallel copyin( &
!$omp f  &
!$omp ,df  &
!$omp )
!$omp end parallel
end subroutine copyin_weno_transport
