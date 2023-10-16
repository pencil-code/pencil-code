subroutine copyin_chemistry
!$omp parallel copyin( &
!$omp rhs_y_full  &
!$omp ,net_react_m  &
!$omp ,net_react_p  &
!$omp )
!$omp end parallel
end subroutine copyin_chemistry
