! $Id$
!
! This module is used to solve the equation of supersaturation
! for either the Smoluchowski approach or the swarm model.
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
! PENCILS PROVIDED ssat
! PENCILS PROVIDED gssat(3); ugssat
! PENCILS PROVIDED del2ssat
! PENCILS PROVIDED tausupersat
!
!***************************************************************
module Supersat
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
  real :: ssat_const=0., amplssat=0., widthssat=0.
  logical :: nosupersat=.false., reinitialize_ssat=.false.
  character (len=labellen) :: initssat='nothing'
!
  namelist /supersat_init_pars/ &
           initssat, ssat_const, amplssat, widthssat
!
!  Run parameters.
!
  real :: supersat_diff=0.0
  real :: supersat_sink=0.0
  real :: vapor_mixing_ratio_qvs=0.0
  real, dimension(3) :: gradssat0=(/0.0,0.0,0.0/)
  logical :: lsupersat_sink=.false., Rsupersat_sink=.false.
  logical :: lupw_ssat=.false., lcondensation_rate=.false.

  namelist /supersat_run_pars/ &
      lupw_ssat, lsupersat_sink, Rsupersat_sink, supersat_sink, &
      supersat_diff, gradssat0, lcondensation_rate, vapor_mixing_ratio_qvs
!
! Declare index of new variables in f array
!
!  integer :: issat=0
! XY: type of "issat" is defined in "cdata.f90"
!
!  Diagnostics variables
!
  integer :: idiag_ssatrms=0, idiag_ssatmax=0, idiag_ssatmin=0
  integer :: idiag_uxssatm=0, idiag_uyssatm=0, idiag_uzssatm=0
!
  contains
!***********************************************************************
    subroutine register_supersat()
!
!  Initialise the ssat variable and increase nvar accordingly
!
!   3-jun-16/xiangyu: adapted from pscalar_nolog
!
      use FArrayManager
!
      call farray_register_pde('ssat', issat)
!NILS: issat is given a value in farray_register_pde, then it can not be
!NILS: set to zero right below. I have therefore commented out this line.
!      issat = 0                 ! needed for idl
!
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
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (lroot) print*, 'Supersaturation routine'
!
!  set to zero and then call the same initial condition
!  that was used in start.csh
!
      if (reinitialize_ssat) then
        f(:,:,:,issat)=0.
        call init_ssat(f)
      endif
!
    endsubroutine initialize_supersat
!
!***********************************************************************
    subroutine init_ssat(f)
!
!  initialise passive scalar field; called from start.f90
!
      use Sub
      use Initcond
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: mudry1, muvap1
      integer :: l
!
      select case (initssat)
        case ('nothing')
        case ('zero'); f(:,:,:,issat)=0.0
        case ('constant'); f(:,:,:,issat)=ssat_const
        case ('tanhx')
          do m=m1,m2;do n=n1,n2
             f(:,m,n,issat)=ssat_const+amplssat*tanh(x/widthssat)
          enddo;enddo
        case ('tanhz')
          do l=l1,l2; do m=m1,m2
             f(l,m,:,issat)=ssat_const+amplssat*tanh(z/widthssat)
          enddo;enddo
!
!  Catch unknown values.
!
        case default
          call fatal_error('initssat','initssat value not recognised')
       endselect
!
!  modify the density to have lateral pressure equilibrium
!
!print*,'AXEL6', mudry1, muvap1
       mudry1=1./29.
       muvap1=1./18.
       f(:,:,:,ilnrho)=f(:,:,:,ilnrho)-alog((1-f(:,:,:,issat))*mudry1+f(:,:,:,issat)*muvap1)
!XX
!
    endsubroutine init_ssat
!***********************************************************************
    subroutine pencil_criteria_supersat()
!
!  All pencils that the Supersat module depends on are specified here.
!
      integer :: i
!
!  ssat itself always
!
      lpenc_requested(i_ssat)=.true.
      if (ldustdensity) then
!       lpenc_requested(i_nd)=.true.
!       lpenc_requested(i_ad)=.true.
      endif
!
!  background gradient 
!
      do i=1,3
        if (gradssat0(i)/=0.) lpenc_requested(i_uu)=.true.
      enddo
!
      if (lsupersat_sink) then 
        lpenc_requested(i_ssat)=.true.
        lpenc_requested(i_tausupersat)=.true.
      endif
      if (supersat_diff/=0.) lpenc_requested(i_del2ssat)=.true.
 
      lpenc_diagnos(i_ssat)=.true.
    endsubroutine pencil_criteria_supersat
!***********************************************************************
    subroutine pencil_interdep_supersat(lpencil_in)
!  Interdependency among pencils provided by the Pscalar module
!  is specified here.
      logical, dimension(npencils) :: lpencil_in
!
      lpencil_in(i_ssat)=.true.
      lpencil_in(i_ugssat)=.true.
      if (lpencil_in(i_ugssat)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_gssat)=.true.
      endif
      lpencil_in(i_tausupersat)=.true.
    endsubroutine pencil_interdep_supersat
!**********************************************************************
    subroutine calc_pencils_supersat(f,p)
!
!  Calculate supersat Pencils.
!  Most basic pencils should come first, as others may depend on them.
!
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      integer :: j
!
      intent(in) :: f
      intent(inout) :: p
!
! ssat
      if (lpencil(i_ssat)) p%ssat=f(l1:l2,m,n,issat)
!
!  Compute gssat. Add imposed spatially constant gradient of ssat.
!  (This only makes sense for periodic boundary conditions.)
!
      if (lpencil(i_gssat)) then
        call grad(f,issat,p%gssat)
        do j=1,3
          if (gradssat0(j)/=0.) p%gssat(:,j)=p%gssat(:,j)+gradssat0(j)
        enddo
      endif
! ugssat
      if (lpencil(i_ugssat)) then
        call u_dot_grad(f,issat,p%gssat,p%uu,p%ugssat,UPWIND=lupw_ssat)
      endif
! del2ssat
      if (lpencil(i_del2ssat)) then
        call del2(f,issat,p%del2ssat)
      endif
!
    endsubroutine calc_pencils_supersat
!***********************************************************************
    subroutine dssat_dt(f,df,p)
!
!  Active scalar evolution for supersation
!  Calculate dssat/dt=-uu.gssat + supersat_diff*[del2ssat + glnrho.gssat].
!
!  27-may-16/xiangyu: adapted from pscalar_nolog
!   4-sep-16/axel: added more diagnostics
!
      use Diagnostics
!      use Dustdensity
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: diff_op,diff_op2,bump,gcgu
      real, dimension (nx) :: radius_sum, condensation_rate_Cd, qv, superaturation
      real :: ssat_xyaver
      real :: tau=10., A1=5e-4
      real :: lam_gradC_fact=1., om_gradC_fact=1., gradC_fact=1.
      integer, parameter :: nxy=nxgrid*nygrid
      integer :: k
!
      intent(in)  :: f
      intent(out) :: df
!
      character(len=2) :: id
!
!  Identify module and boundary conditions.
!
      if (nosupersat) then
        if (headtt.or.ldebug) print*,'not SOLVED: dssat_dt'
      else
        if (headtt.or.ldebug) print*,'SOLVE dssat_dt'
      endif
      if (headtt) then
        write(id,'(i0)')
        call identify_bcs('ssat'//trim(id),issat)
      endif
!
!  Passive scalar equation.
!
      df(l1:l2,m,n,issat)=df(l1:l2,m,n,issat)-p%ugssat
      if (supersat_diff/=0.) &
          df(l1:l2,m,n,issat)=df(l1:l2,m,n,issat)+supersat_diff*p%del2ssat
!
        ! 1-June-16/XY coded: to be completed, here is only for 1D case 
        if (lsupersat_sink) then
          if (Rsupersat_sink) then
            bump=A1*p%uu(:,1)
            !print*,'tau=',p%tausupersat
          else
            bump=A1*p%uu(:,1)-f(l1:l2,m,n,issat)/tau
            !bump=A1*p%uu(:,1)-f(l1:l2,m,n,issat)/p%tausupersat
            !print*,'tau=',f(l1:l2,m,n,itausupersat)
            !print*,'tau=',p%tausupersat
          endif
        !df(l1:l2,m,n,issat)=df(l1:l2,m,n,issat)-p%ugssat+bump
        df(l1:l2,m,n,issat)=df(l1:l2,m,n,issat)+bump
      endif
!
      if (ldustdensity.and.lcondensation_rate) then
!
!  compute sum of particle radii
!
        qv=f(l1:l2,m,n,issat)
        superaturation=qv/vapor_mixing_ratio_qvs-1.
        radius_sum=0.
!       do k=1,ndustspec
!         radius_sum=radius_sum+p%ad(:,k)*p%nd(:,k)
!         print*,k,radius_sum(1:5)
!       enddo
        condensation_rate_Cd=4.*pi*superaturation*radius_sum
        df(l1:l2,m,n,issat)=df(l1:l2,m,n,issat)-condensation_rate_Cd
        if (ltemperature) df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT)+condensation_rate_Cd/(p%cp*p%TT)
      endif
!
!  Diagnostics
!
      if (ldiagnos) then
        if (idiag_ssatrms/=0) call sum_mn_name(p%ssat**2,idiag_ssatrms,lsqrt=.true.)
        if (idiag_ssatmax/=0) call max_mn_name(p%ssat,idiag_ssatmax)
        if (idiag_ssatmin/=0) call max_mn_name(-p%ssat,idiag_ssatmin,lneg=.true.)
        if (idiag_uxssatm/=0) call sum_mn_name(p%uu(:,1)*p%ssat,idiag_uxssatm)
        if (idiag_uyssatm/=0) call sum_mn_name(p%uu(:,2)*p%ssat,idiag_uyssatm)
        if (idiag_uzssatm/=0) call sum_mn_name(p%uu(:,3)*p%ssat,idiag_uzssatm)
      endif
!
    endsubroutine dssat_dt
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
!***********************************************************************
    subroutine rprint_supersat(lreset,lwrite)
      use Diagnostics
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname
      logical :: lwr
!
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!      
      if (lreset) then
        idiag_ssatrms=0; idiag_ssatmax=0; idiag_ssatmin=0
        idiag_uxssatm=0; idiag_uyssatm=0; idiag_uzssatm=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'ssatrms',idiag_ssatrms)
        call parse_name(iname,cname(iname),cform(iname),'ssatmax',idiag_ssatmax)
        call parse_name(iname,cname(iname),cform(iname),'ssatmin',idiag_ssatmin)
        call parse_name(iname,cname(iname),cform(iname),'uxssatm',idiag_uxssatm)
        call parse_name(iname,cname(iname),cform(iname),'uyssatm',idiag_uyssatm)
        call parse_name(iname,cname(iname),cform(iname),'uzssatm',idiag_uzssatm)
      enddo
!
      if (lwr) then 
        write(3,*) 'issat = ', issat
      endif
!
    endsubroutine rprint_supersat 
!***********************************************************************
endmodule Supersat
