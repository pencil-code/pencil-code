subroutine copyin_special
!$omp parallel copyin( &
!$omp bavg  &
!$omp ,uavg  &
!$omp ,rhoavg  &
!$omp )
!$omp end parallel
end subroutine copyin_special
