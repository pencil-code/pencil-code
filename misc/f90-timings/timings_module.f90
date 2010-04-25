!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   timings_module.f90   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!  Author: wd (Wolfgang.Dobler@ucalgary.ca)
!!!  Date:   18-Sep-2006
!!!
!!!  Description:
!!!   Utility module for code timings.

module Timings

  implicit none

  include "param.inc"

  interface dummy_use
    module procedure dummy_use_scal
    module procedure dummy_use_vect
  endinterface

  real, dimension(nx)   :: scal_field1,scal_field2,scal_field3
  real, dimension(nx,3) :: vect_field1,vect_field2,vect_field3
  !
  logical               :: false ! always false, but optimizer must not know

  contains

!***********************************************************************
    subroutine init(dummy_flag)
!
! Initialize 1d variables and similar stuff
!
      integer :: i
      real    :: x
      logical :: dummy_flag
!
      if (dummy_flag) then
        write(0,*) &
            'Unexpected: dummy_flag should be false to avoid excessive output'
      endif
      false = dummy_flag
!
! Use chaotic logistic map to initialize scal1 with `random' values
! in [0,1]. Could just as well use random number, I guess, but this is
! definitively portable.
!
      x = sqrt(0.7)
      do i=1,nx
        x = 4*x*(1-x)
        scal_field1(i) = x
      enddo
!
! Set scal2 to inverse of scal1, scal3 to reverse sequence
!
      scal_field2 = 1./ scal_field1
      scal_field3 = scal_field1(nx:1:-1)
!
! Splice scalar values into vectors
!
      vect_field1(:,1) = scal_field1
      vect_field1(:,2) = scal_field2
      vect_field1(:,3) = scal_field3
!
      vect_field2 = 1./ vect_field1
      vect_field3 = cshift(vect_field1,nx/2,DIM=1)
!
    endsubroutine init
!***********************************************************************
    subroutine dummy_use_scal(scal_field)
!
! Call this to make the compiler think we are using the argument
! (although we are not..), in order to avoid optimizing the whole loop
! away.
!
! Scalar version.
!
      real, dimension(nx) :: scal_field
!
      if (false) print*, scal_field(1)
!
    endsubroutine dummy_use_scal
!***********************************************************************
    subroutine dummy_use_vect(vect_field)
!
! Call this to make the compiler think we are using the argument
! (although we are not..), in order to avoid optimizing the whole loop
! away.
!
! Vector version
!
      real, dimension(nx,3) :: vect_field
!
      if (false) print*, vect_field(1,2)
!
    endsubroutine dummy_use_vect
!***********************************************************************
    subroutine report_timings(labels, measured_times)
!
      character(LEN=*), dimension(:)  :: labels
      real, dimension(size(labels,1)) :: measured_times

      real              :: abstime, reltime
      integer           :: ntimes,i
      character(LEN=40) :: fmt1,fmt2
      logical           :: imprecise ! flag for unreliable timings 
!
      ntimes    = size(labels,1)
      imprecise = .false.
!
! Print out results
!
      fmt1 = '(A30,": ",1p,2G10.2,A)'
      fmt2 = '(A30,": ",1p,2G10.2)'
      write(*,fmt1) 'variant        ', 'absolute', 'relative', ' [large is bad]'
      write(*,*) '  ---------------------------+-------------------------------'
     do i=1,ntimes
        abstime = measured_times(i)
        reltime = abstime / measured_times(1)
        write(*,fmt2) labels(i), abstime, reltime
        imprecise = assess_precision(abstime)
      enddo
!
! Warn about possible bad precision
!
      if (imprecise) then
        write(*,*) 'Warning: some values are imprecise;' &
                  // ' consider increasing NITER!'
      endif
!
    endsubroutine report_timings
!***********************************************************************
    logical function assess_precision(abstime)
!
! Set flag IMPRECISE to T if the absolute time indicates impreceise
! results. Assumes that ABSTIME is in seconds and uses some handwaving gut
! numbers (you get the picture...).
!
      real,    intent(in)  :: abstime

      real,    save :: clocktick=0. 
      logical, save :: first_call=.true.
!
! Get and store clocktick value
!
      if (first_call) then
        clocktick = mpiwtick()
        first_call = .false.
      endif
!
      if (      (abstime < 2) &                ! < 2 seconds is unreliable
           .or. (abstime/clocktick < 10)) then ! < 10 ticks means > 10% error
        assess_precision = .true.
      else
        assess_precision = .false.
      endif
!
    endfunction assess_precision
!***********************************************************************
    function mpiwtime()
!
!  Mimic the MPI_WTIME() timer function. On many machines, the
!  implementation through system_clock() will overflow after about 50
!  minutes, so MPI_WTIME() is better.
!
!  Usage:
!    real :: time0,time1,time2
!
!    call init()
!    time0 = mpiwtime()
!    call work1()
!    time1 = mpiwtime()
!    call work2()
!    time2 = mpiwtime()
!
!    print*, 'Work1: , time1-time0
!    print*, 'Work2: , time2-time1
!
!  5-oct-2002/wolf: coded
!
      double precision :: mpiwtime
      integer :: count_rate,time
!
      call system_clock(COUNT_RATE=count_rate)
      call system_clock(COUNT=time)

      if (count_rate /= 0) then
        mpiwtime = (time*1.)/count_rate
      else                      ! occurs with ifc 6.0 after long (> 2h) runs
        mpiwtime = 0
      endif
!
    endfunction mpiwtime
!***********************************************************************
    function mpiwtick()
!
!  Mimic the MPI_WTICK() function for measuring timer resolution.
!
!   5-oct-2002/wolf: coded
!
      double precision :: mpiwtick
      integer :: count_rate
!
      call system_clock(COUNT_RATE=count_rate)
      if (count_rate /= 0) then
        mpiwtick = 1./count_rate
      else                      ! occurs with ifc 6.0 after long (> 2h) runs
        mpiwtick = 0
      endif
!
    endfunction mpiwtick
!***********************************************************************


endmodule Timings

!!! End of file timings_module.f90
