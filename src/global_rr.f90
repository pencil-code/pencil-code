module Global

!
!  A module container for additional variables which are globally needed 
!  --- here this is spherical radius (for calculating gravity and heating
!  /cooling functions) 
!  Put into a module, so one can easily switch.
!  NB: These variables use half as much memory as a new variable, so
!  keep their number at minimum.
!
  use Cparam

  implicit none

  real, dimension (mx,my,mz) :: rr

  contains

!***********************************************************************
    subroutine write_global()
!
!  write global variables
!
!  10-jan-02/wolf: coded
!
      use Cdata
      use Mpicomm
      use IO
!
      call output(trim(directory)//'/rr.dat',rr,1)
!
    endsubroutine write_global
!***********************************************************************
    subroutine read_global()
!
!  read global variables
!
!  10-jan-02/wolf: coded
!
      use Cdata
      use Mpicomm
      use IO
!
      call input(trim(directory)//'/rr.dat',rr,1,0)
!
    endsubroutine read_global
!***********************************************************************

endmodule Global
