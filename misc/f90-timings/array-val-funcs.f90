!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   array-val-funcs.f90   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!  Author: wd (Wolfgang.Dobler@ucalgary.ca)
!!!  Date:   18-Sep-2006
!!!
!!!  Description:
!!!   Compare array-valued functions with subroutines. Only 1d is
!!!   relevant for the Pencil Code.
!!!
!!!  Results on Frenesi (large values are bad), normalized by Subroutine
!!!  time:
!!!   Compiler:               g95      gfortran  openf90
!!!   Version:                0.92     4.4.1     4.2.1
!!!                           Jun 24
!!!                           2009
!!!   ---------------------------------------------------
!!!   Subroutine:             1.0      1.0       1.0
!!!   Array-valued fn:        1.9      2.6       2.1
!!!   Fair array-valued fn:   1.2      1.8       1.1
!!!   Hybrid w/ intrinsic     0.79     0.67      0.66
!!!     vector operations:
!!!
!!!  Same results in absolute times (comparable since tests were done on
!!!  the same machine):
!!!   Compiler:               g95     gfortran  openf90
!!!   Version:                0.92    4.4.1     4.2.1
!!!   --------------------------------------------------
!!!   Subroutine:             20.      6.0       7.5
!!!   Array-valued fn:        37.     16.       16.
!!!   Fair array-valued fn:   26.     11.        8.4
!!!   Hybrid w/ intrinsic     15.     4.0        5.0
!!!     vector operations:



! ---------------------------------------------------------------------- !

program Array_val_funcs

  use Timings
  use Array_valued

  implicit none

  real :: t0,t1,t2,t3,t4

  call init(.false.)

  t0 = mpiwtime()
  call sub1()
  t1 = mpiwtime()
  call sub2()
  t2 = mpiwtime()
  call sub3()
  t3 = mpiwtime()
  call sub4()
  t4 = mpiwtime()

  call report_timings(&
      (/ 'Subroutine                         ', &
         'Array-valued fn                    ', &
         'Fair array-valued fn               ', &
         'Hybrid w/ intrinsic vector operator' /) , &
      (/ t1-t0, t2-t1, t3-t2, t4-t3 /) &
  )

contains

!***********************************************************************
  subroutine sub1()
!
    integer :: i
!
    do i=1,niter
      !
      ! Mutliply first by scal1, then by scal2=1/scal1 to avoid over/underflow
      !
      call multsv_mn(scal_field1, vect_field1, vect_field2)
      call multsv_mn(scal_field2, vect_field2, vect_field1)
    enddo
!
    call dummy_use(vect_field2)
!
  endsubroutine sub1
!***********************************************************************
  subroutine sub2()
!
    integer :: i
!
    do i=1,niter
      !
      ! Mutliply first by scal1, then by scal2=1/scal1 to avoid over/underflow
      !
      vect_field2 = scal_field1 * vect_field1
      vect_field1 = scal_field2 * vect_field2
    enddo
!
    call dummy_use(vect_field2)
!
  endsubroutine sub2
!***********************************************************************
  subroutine sub3()
!
    integer :: i
!
    do i=1,niter
      !
      ! Mutliply first by scal1, then by scal2=1/scal1 to avoid over/underflow
      !
      vect_field1 = scal_field2 * scal_field1 * vect_field1
    enddo
    vect_field2 = vect_field1/scal_field2

    call dummy_use(vect_field2)
!
  endsubroutine sub3
!***********************************************************************
  subroutine sub4()
!
    integer :: i
!
    do i=1,niter
      !
      ! Mutliply first by scal1, then by scal2=1/scal1 to avoid over/underflow
      !
      call multsv_mn(scal_field1*scal_field2, vect_field1, vect_field1)
    enddo
    call multsv_mn(1./scal_field2, vect_field1, vect_field2)

    call dummy_use(vect_field2)
!
  endsubroutine sub4


endprogram Array_val_funcs

!!! End of file array-val-funcs.f90
