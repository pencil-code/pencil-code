subroutine copyin_messages
!$omp parallel copyin( &
!$omp errors  &
!$omp ,fatal_errors  &
!$omp ,scaller  &
!$omp ,message_stored  &
!$omp )
!$omp end parallel
end subroutine copyin_messages
