module Global

!  Additional global variables.
!  Put into a module so one can easily switch.
!  NB: Each of these variables takes half the space as one physical
!    variable f(:,:,:,i), so keep their number at minimum.

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
      use Sub
!
!
      call output(trim(directory)//'/global.dat',m_pot,1)
!
    endsubroutine write_global
!***********************************************************************
    subroutine read_global()
!
!  write global variables
!
!  10-jan-02/wolf: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
!
      call input(trim(directory)//'/global.dat',m_pot,1,0)
!
    endsubroutine read_global
!***********************************************************************

endmodule Global
