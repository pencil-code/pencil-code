! $Id$
! 
! CPARAM logical, parameter :: ltavg = .false.
!
module Timeavg
!
  use Cdata
  use General, only: keep_compiler_quiet
!
  implicit none
!
  include 'timeavg.h'
!
  integer :: idx_tavg=0         ! just scalar, since unused and no mtavg known
!
  contains
!***********************************************************************
    subroutine initialize_timeavg(a)
!
      real, dimension(mx,my,mz,mfarray) :: a
!
      intent (in) :: a
!
      call keep_compiler_quiet(a)
!
    endsubroutine initialize_timeavg
!***********************************************************************
    subroutine update_timeavgs(a,dt,init)
!
      real, dimension(mx,my,mz,mfarray) :: a
      real :: dt
      logical, optional :: init
!
      intent (in) :: a
!
      call keep_compiler_quiet(a)
      call keep_compiler_quiet(dt)
      if (present(init)) call keep_compiler_quiet(init)
!
    endsubroutine update_timeavgs
!***********************************************************************
    subroutine wsnap_timeavgs(chsnap,enum,flist)
!
      character (len=*) :: chsnap,flist
      logical :: enum
      optional :: flist
!
      call keep_compiler_quiet(chsnap)
      call keep_compiler_quiet(enum)
      if (present(flist)) call keep_compiler_quiet(flist)
!
    endsubroutine wsnap_timeavgs
!***********************************************************************
    subroutine pushpars2c(p_par)

      use Syscalls, only: copy_addr

      integer, parameter :: n_pars=0
      integer(KIND=ikind8), dimension(n_pars) :: p_par

      call keep_compiler_quiet(p_par)

    endsubroutine pushpars2c
!***********************************************************************
endmodule Timeavg
