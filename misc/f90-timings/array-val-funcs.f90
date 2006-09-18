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
!!!  Results (large values are bad):
!!!   Compiler
!!!       Variant           : absolute  relative
!!!   ----------------------+---------------------
!!!   - g95 -O4
!!!       Subroutine        : 2.6       1.0
!!!       Array-valued-func
!!!         + overloading * : 7.5       2.9

program Array_val_funcs

  use Timings

  implicit none

  interface operator (*)
    module procedure pencil_multiply1
  endinterface

  interface operator (/)
    module procedure pencil_divide1
  endinterface


  real :: t0,t1,t2,t3,t4

  call init()

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
      (/ 'Subroutine', 'Array-valued fn',  'Fair array-valued fn', &
                  'Hybrid w/ intrinsic vector operator'/) , &
      (/ t1-t0       , t2-t1, t3-t2, t4-t3         /) &
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
!***********************************************************************
  subroutine multsv_mn(a,b,c)
!
!  vector multiplied with scalar, gives vector
!   22-nov-01/nils erland: coded
!   10-oct-03/axel: a is now the scalar (now consistent with old routines)
!
!    use Cdata
!
    intent(in)  :: a,b
    intent(out) :: c
!
    real, dimension (nx,3) :: b,c
    real, dimension (nx) :: a
    integer :: i
!
    do i=1,3
      c(:,i)=a*b(:,i)
    enddo
!
  endsubroutine multsv_mn
!***********************************************************************
  function pencil_multiply1(s,v)
!
!  The `*' operator may be extended through this function to allow
!  elementwise multiplication of a `pencil-scalar' with a `pencil-vector'
!
!   6-Sep-05/tobi: coded
!
!    use Cdata

    real, dimension(nx), intent(in) :: s
    real, dimension(nx,3), intent(in) :: v
    real, dimension(nx,3) :: pencil_multiply1

    integer :: i

    do i=1,3
      pencil_multiply1(:,i) = s(:) * v(:,i)
    enddo

  endfunction pencil_multiply1
!***********************************************************************
  function pencil_divide1(v,s)
!
!  The `*' operator may be extended through this function to allow
!  elementwise multiplication of a `pencil-scalar' with a `pencil-vector'
!
!   18-Sep-05/tony: coded
!
!    use Cdata

    real, dimension(nx), intent(in) :: s
    real, dimension(nx,3), intent(in) :: v
    real, dimension(nx,3) :: pencil_divide1

    integer :: i

    do i=1,3
      pencil_divide1(:,i) = v(:,i) / s(:)
    enddo

  endfunction pencil_divide1
!***********************************************************************


endprogram Array_val_funcs

!!! End of file array-val-funcs.f90
