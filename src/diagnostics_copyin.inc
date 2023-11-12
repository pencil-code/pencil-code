subroutine copyin_diagnostics
!$omp parallel copyin( &
!$omp phiavg_profile  &
!$omp ,phiavg_norm  &
!$omp )
!$omp end parallel
end subroutine copyin_diagnostics
