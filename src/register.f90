! $Id: register.f90,v 1.30 2002-06-01 09:36:38 brandenb Exp $

!!!  A module for setting up the f-array and related variables (`register' the
!!!  entropy, magnetic, etc modules). Didn't know where else to put this:
!!!  Entropy uses Sub and Init must be used by both, Start and Run.


module Register

  implicit none 

  contains

!***********************************************************************
    subroutine initialize()
!
!  Call all initialisation hooks, i.e. initialise MPI and register
!  physics modules.
!
!  6-nov-01/wolf: coded
!
      use Cdata
      use Mpicomm
      use Sub
      use Gravity
      use Hydro
      use Forcing
      use Entropy
      use Magnetic
!
      call mpicomm_init
!
      nvar = 0                  ! to start with
      call register_hydro
      call register_forcing
      call register_ent
      call register_aa
      call register_grav
!
      if (nvar /= mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('Initialize: nvar /= mvar. Fix mvar in cparam.f90')
      endif
!
    endsubroutine initialize
!***********************************************************************
    subroutine rprint_list(lreset)
!
!  read in output times from control file
!
!   3-may-01/axel: coded
!
      use Cdata
      use Hydro
      use Entropy
      use Magnetic
!
      integer :: iname
      logical :: lreset
!
!  read in the list of variables to be printed
!
      open(1,file='print.in')
      do iname=1,mname
        read(1,*,end=99) cname(iname)
      enddo
      close(1)
99    nname=iname-1
      if (lroot) print*,'nname=',nname
!
!  check which variables are set
!
      call rprint_hydro(lreset)
      call rprint_entropy(lreset)
      call rprint_magnetic(lreset)
!
    endsubroutine rprint_list
!***********************************************************************

endmodule Register

!!! End of file register.f90
