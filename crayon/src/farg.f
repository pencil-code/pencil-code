! farg.f
! Needed with some combinations of g77-copiled mpich and F90 compileri
! to fix the error message
!   undefined reference to `f__xargc'
!
! See
!   http://www.pgroup.com/support/link.htm
! or
!   http://www.clumeq.mcgill.ca/faq.html
!
! Adapted this to work with one underscore appended (ifort 8.1 with
! -us option [on by default])
!
! Usage:
!   ifort -c farg.f
!   ifort *.o farg.o -lmpich -l... -o run.x
!
      integer function mpir_iargc_()
      mpir_iargc_ = iargc()
      return
      end
!
      subroutine mpir_getarg_( i, s )
      integer       i
      character*(*) s
      call getarg(i,s)
      return
      end
