module General

!!! Module with utility subroutines which are used by procedures in Sub
!!! and in Mpicomm.

  implicit none

  contains

!***********************************************************************
    subroutine chn(n,ch)
      character*4 ch
      integer :: n
!
!  make a character out of a number
!  take care of numbers that have less than 4 digits
!  30-sep-97/axel: coded
!
      ch='    '
      if (n.lt.0) stop 'chn: lt1'
      if (n.lt.10) then
        write(ch(1:1),'(i1)') n
      elseif (n.lt.100) then
        write(ch(1:2),'(i2)') n
      elseif (n.lt.1000) then
        write(ch(1:3),'(i3)') n
      elseif (n.lt.10000) then
        write(ch(1:4),'(i4)') n
      else
        print*,'n=',n
        stop "chn: n too large"
      endif
!
    endsubroutine chn
!***********************************************************************
    subroutine chk_time(label,time1,time2)
!
      integer :: time1,time2,count_rate
      character*(*) label
!
!  prints time in seconds
!
      call system_clock(count=time2,count_rate=count_rate)
      print*,label,(time2-time1)/real(count_rate)
      time1=time2
!
    endsubroutine chk_time
!***********************************************************************
!
end module General
