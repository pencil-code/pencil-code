module Global

!
!  A module container for additional variables which are globally needed 
!  --- here this is the (negative) gravity potential.
!  Put into a module, so one can easily switch.
!  NB: These variables use half as much memory as a new variable, so
!  keep their number at minimum.
!
  use Cparam

  implicit none

  real, dimension (mx,my,mz) :: m_pot

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
      call output(trim(directory)//'/global.dat',m_pot,1)
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
      call input(trim(directory)//'/global.dat',m_pot,1,0)
!
    endsubroutine read_global
!***********************************************************************

endmodule Global
