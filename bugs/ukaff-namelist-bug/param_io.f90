! $Id$ 

module Param_IO

  use Hydro
  use Magnetic
 
  implicit none 
 
  contains

!***********************************************************************
    subroutine read_runpars()
!
      implicit none
!
!  open namelist file
!

!
!  Ukaff bug workaround: uncomment following line
!
!  namelist /magnetic_run_pars/ eta

      open(1,file='run.in',form='formatted')

      write(0,*) '1.'
      read(1,NML=hydro_run_pars )
      write(0,*) '2.'
      read(1,NML=magnetic_run_pars)
      write(0,*) '3.'
      write(0,*) 'Finished successfully'

      close(1)
!
    endsubroutine read_runpars
!***********************************************************************

endmodule Param_IO


