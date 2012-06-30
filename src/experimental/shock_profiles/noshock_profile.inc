!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
!
!***********************************************************************
  subroutine shock_calc_body(f)
    use Cdata
    use Mpicomm, only: stop_it
    real, dimension (mx,my,mz,mvar+maux) :: f
     real :: dA,dB,dC
     real :: face_factor1
     real :: face_factor2
call fatal_error('shock_calc_body',& 
'There is no shock profile compiled in. Set SHOCK_PROFILE in Makefile.local.')
    if (ALWAYS_FALSE) print *,f(1,1,1,1)
    endsubroutine shock_calc_body
!***********************************************************************
  subroutine shock_calc_internalboundary(f)
    use Cdata
    use Mpicomm, only: stop_it
    real, dimension (mx,my,mz,mvar+maux) :: f
call fatal_error('shock_calc_body',& 
'There is no shock profile compiled in. Set SHOCK_PROFILE in Makefile.local.')
    if (ALWAYS_FALSE) print *,f(1,1,1,1)
  endsubroutine shock_calc_internalboundary
!***********************************************************************
  subroutine shock_calc_externalboundary(f)
    use Cdata
    use Mpicomm, only: stop_it
    real, dimension (mx,my,mz,mvar+maux) :: f
call fatal_error('shock_calc_body',& 
'There is no shock profile compiled in. Set SHOCK_PROFILE in Makefile.local.')
    if (ALWAYS_FALSE) print *,f(1,1,1,1)
  endsubroutine shock_calc_externalboundary
!***********************************************************************
