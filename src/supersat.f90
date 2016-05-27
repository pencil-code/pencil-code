! $Id$
!
! This module is used to solve the passive scalar 
! equation of supersaturation for swarm model. 
!
! **************************************************
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lsupersat = .true.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED cc; cc1
! PENCILS PROVIDED gcc; ugcc

!***************************************************************
module Superstat
!
  use Cdata
  use Cparam
  use Messages
!
  implicit none
!

include 'supersat.h'
!
!  Init parameters.
!
  logical :: nosupersat=.false., reinitialize_cc=.false.
  logical :: reinitialize_lncc=.false.
  character (len=labellen) :: initlncc='impossible'
 !
  namelist /supersat_init_pars/ &
           initlncc
!
!  Run parameters.
!
  real :: supersat_diff=0.0
  real :: supersat_sink=0.0, Rsupersat_sink=0.5
  real, dimension(3) :: gradC0=(/0.0,0.0,0.0/)
  logical :: lsupersat_sink=.false.

  namelist /supersat_run_pars/ &
      lsupersat_sink, Rsupersat_sink, supersat_sink &
      supersat_diff, gradC0
!**********************************
    subroutine register_supersat()
      use FArrayManager
!
      !lpscalar_nolog = .true.
      lsupersat = .true.
!
      call farray_register_pde('cc', icc)
      ilncc = 0                 ! needed for idl
!
!  Identify version number.
!
      if (lroot) call svn_id( &
          "$Id$")
!
    endsubroutine register_supersat
!***********************************************************************
    subroutine initialize_supersat(f)
!
!  Perform any necessary post-parameter read initialization
!  Since the passive scalar is often used for diagnostic purposes
!  one may want to reinitialize it to its initial distribution.
!
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!
!  set to zero and then call the same initial condition
!  that was used in start.csh
!
      if (reinitialize_cc) then
        f(:,:,:,icc)=0.
        call init_lncc(f)
      endif
!
    endsubroutine initialize_supersat
!
!**********************************************************************
    subroutine calc_pencils_supersat(f,p)
!
!  Calculate pscalar Pencils.
!  Most basic pencils should come first, as others may depend on them.
!
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
      
! cc
      if (lpencil(i_cc)) p%cc=f(l1:l2,m,n,icc)
! cc1
      if (lpencil(i_cc1)) p%cc1=1./p%cc
! gcc
      if (lpencil(i_gcc)) then
          call grad(f,icc+i-1,p%gcc)
      endif
! ugcc
      if (lpencil(i_ugcc)) then
          call u_dot_grad(f,icc,p%gcc,p%uu,p%ugcc,UPWIND=lupw_cc)
      endif
!
    endsubroutine calc_pencils_supersat
!***********************************************************************
    subroutine dlncc_dt(f,df,p)
      use Diagnostics
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: diff_op,diff_op2,bump,gcgu
      real :: cc_xyaver
      real :: lam_gradC_fact=1., om_gradC_fact=1., gradC_fact=1.
      integer :: j, k
      integer, parameter :: nxy=nxgrid*nygrid
!
      intent(in)  :: f
      intent(out) :: df
!
      character(len=2) :: id
!  Passive scalar equation.
!
        df(l1:l2,m,n)=df(l1:l2,m,n)-p%ugcc
!
!  Passive scalar sink.
!
        if (supersat_sink) then
          if (Rsupersat_sink==0) then
            bump=supersat_sink
          else
            bump=supersat_sink*exp(-0.5*(x(l1:l2)**2+y(m)**2+z(n)**2)/Rsupersat_sink**2)
          endif
          df(l1:l2,m,n)=df(l1:l2,m,n)-spread(bump,2)*p%cc
        endif
    endsubroutine dlncc_dt
!
!***********************************************************************
    subroutine read_supersat_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=supersat_init_pars, IOSTAT=iostat)
!
    endsubroutine read_supersat_init_pars
!***********************************************************************
    subroutine write_supersat_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=supersat_init_pars)
!
    endsubroutine write_supersat_init_pars
!***********************************************************************
    subroutine read_supersat_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=supersat_run_pars, IOSTAT=iostat)
!
    endsubroutine read_supersat_run_pars
!***********************************************************************
    subroutine write_supersat_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=supersat_run_pars)
!
    endsubroutine write_supersat_run_pars
!
endmodule Supersat
