! $Id: testperturb.f90,v 1.1 2008-03-21 07:55:13 brandenb Exp $

!  test perturbation method

module TestPerturb

  use Sub
  use Cdata
  use Messages

  implicit none

  real :: dummy

  include 'testperturb.h'

  namelist /testperturb_init_pars/ &
      dummy

  namelist /testperturb_run_pars/ &
      ltestperturb

  integer :: idiag_alp11=0

  contains

!***********************************************************************
    subroutine register_testperturb()
!
!  Initialise variables
!
!  2-july-02/nils: coded
!
      use Mpicomm
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_testperturb called twice')
      first = .false.
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: testperturb.f90,v 1.1 2008-03-21 07:55:13 brandenb Exp $")
!
    endsubroutine register_testperturb
!***********************************************************************
    subroutine initialize_testperturb()
!
!  21-nov-02/tony: coded
!  08-jul-04/anders: Sshear calculated whenever qshear /= 0

!  calculate shear flow velocity; if qshear is given then Sshear=-qshear*Omega
!  is calculated. Otherwise Sshear keeps its value from the input list.
!
      if (qshear/=0.0) Sshear=-qshear*Omega
      if (lroot .and. ip<=12) &
          print*,'initialize_testperturb: Sshear,qshear=',Sshear,qshear
!
    endsubroutine initialize_testperturb
!***********************************************************************
    subroutine read_testperturb_init_pars(unit,iostat)
!
!  read initial testperturb parameters
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=testperturb_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=testperturb_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_testperturb_init_pars
!***********************************************************************
    subroutine write_testperturb_init_pars(unit)
!
!  write initial testperturb parameters
!
      integer, intent(in) :: unit
!
      write(unit,NML=testperturb_init_pars)
!
    endsubroutine write_testperturb_init_pars
!***********************************************************************
    subroutine read_testperturb_run_pars(unit,iostat)
!
!  read run testperturb parameters
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=testperturb_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=testperturb_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_testperturb_run_pars
!***********************************************************************
    subroutine write_testperturb_run_pars(unit)
!
!  write run testperturb parameters
!
      integer, intent(in) :: unit
!
      write(unit,NML=testperturb_run_pars)
!
    endsubroutine write_testperturb_run_pars
!***********************************************************************
    subroutine testperturbing(finit,f,df)
!
!  Calculate the response to perturbations after one timestep.
!  The name "testperturbing" is chosen, because a subroutine's name
!  cannot be the same as that of a module (testperturb).
!  We have followed here the convection from the shear module.
!
!  19-mar-08/axel: coded
!
      use Timestep, only: rk_2n
      use Hydro, only: calc_pencils_hydro

      real, dimension (mx,my,mz,mfarray) :: finit,f
      real, dimension (mx,my,mz,mvar) :: df
!
      real, dimension (nx,3) :: uu,bb,uxb
      logical :: headtt_save
      type (pencil_case) :: p
!
!  print identifier
!
      if (headtt.or.ldebug) print*, 'testperturbing'
!
      finit=f
!
!  Time advance
!
          call rk_2n(f,df,p)
!
!  calculate emf for alpha effect (for imposed field)
!  Note that uxbm means <EMF.B0>/B0^2, so it gives already alpha=EMF/B0.
!
      lfirstpoint=.true.
      headtt_save=headtt
      do m=m1,m2
      do n=n1,n2
!       if (idiag_uxbm/=0 .or. idiag_uxbmx/=0 .or. idiag_uxbmy/=0 &
!           .or. idiag_uxbmz/=0) then
!         uxb_dotB0=B_ext(1)*p%uxb(:,1)+B_ext(2)*p%uxb(:,2)+B_ext(3)*p%uxb(:,3)
!         uxb_dotB0=uxb_dotB0*B_ext21
!         if (idiag_uxbm/=0) call sum_mn_name(uxb_dotB0,idiag_uxbm)
!         if (idiag_uxbmx/=0) call sum_mn_name(uxbb(:,1),idiag_uxbmx)
!         if (idiag_uxbmy/=0) call sum_mn_name(uxbb(:,2),idiag_uxbmy)
!         if (idiag_uxbmz/=0) call sum_mn_name(uxbb(:,3),idiag_uxbmz)
          call calc_pencils_hydro(f,p)
          call curl(f,iaa,bb)
          call cross_mn(p%uu,bb,uxb)
          call sum_mn_name(uxb(:,1),idiag_alp11)
        headtt=.false.
        lfirstpoint=.false.
      enddo
      enddo
!
!  reset f to what it was before
!
      f=finit
!
!  Since this routine was called before any regular call to the
!  timestepping routine, we must reset headtt to what it was before.
!
      headtt=headtt_save
!
    endsubroutine testperturbing
!***********************************************************************
    subroutine advance_testperturb(f,df,dt_shear)
!
!  Advance shear distance, deltay, using dt. Using t instead introduces
!  significant errors when nt = t/dt exceeds ~100,000 steps.
!  This formulation works also when Sshear is changed during the run.
!
!  18-aug-02/axel: incorporated from nompicomm.f90
!
      use Cdata
      use Fourier, only: fourier_shift_y
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real :: dt_shear
!
      real, dimension (nx,ny,nz) :: tmp
      real, dimension (nx) :: uy0
      integer :: l, ivar
!
!  Print identifier.
!
      if (headtt.or.ldebug) print*, 'advance_shear: deltay=',deltay
!
    end subroutine advance_testperturb
!***********************************************************************
    subroutine rprint_testperturb(lreset,lwrite)
!
!  reads and registers print parameters relevant to testperturb
!
!  20-mar-08/axel: adapted from shear.f90
!
      use Cdata
      use Sub
!
      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_alp11=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'alp11',idiag_alp11)
      enddo
!
!  write column where which testperturb variable is stored
!
      if (lwr) then
        write(3,*) 'i_alp11=',idiag_alp11
      endif
!
    endsubroutine rprint_testperturb
!***********************************************************************
  endmodule TestPerturb
