! $Id: newton_cooling.f90,v 1.2 2015/12/11 02:56:37 wlyra Exp $

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 2
!
!***************************************************************

!-------------------------------------------------------------------
!
! HOW TO USE THIS FILE
! --------------------
!
! The rest of this file may be used as a template for your own
! special module.  Lines which are double commented are intended
! as examples of code.  Simply fill out the prototypes for the 
! features you want to use.
!
! Save the file with a meaningful name, eg. geo_kws.f90 and place
! it in the $PENCIL_HOME/src/special directory.  This path has
! been created to allow users ot optionally check their contributions
! in to the Pencil-Code CVS repository.  This may be useful if you
! are working on/using the additional physics with somebodyelse or 
! may require some assistance from one of the main Pencil-Code team.
!
! To use your additional physics code edit the Makefile.local in
! the src directory under the run directory in which you wish to
! use your additional physics.  Add a line with all the module 
! selections to say something like:
!
!    SPECIAL=special/nstar
!
! Where nstar it replaced by the filename of your new module
! upto and not including the .f90
!
!--------------------------------------------------------------------
!
!  This method solves for flux-limited diffusion, as per
!  
!    Levermore & Pomraning 1981, ApJ, 248, 321
!    Kley 1989, A&A, 208, 98
!  
!  11-dec-14/wlad: coded
!
module Special
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include '../special.h'
!
  real :: taucool_floor = impossible
  logical :: laddheatingrate=.true.
!
  namelist /special_init_pars/ taucool_floor
!
  namelist /special_run_pars/ laddheatingrate,taucool_floor
!
  type InternalPencils
     real, dimension(nx)   :: kappa,tau,tau1,taucool,taucool1
  endtype InternalPencils
!
  real, dimension(nx,nz,0:nprocy-1) :: tau_column
  real, dimension(nx,ny,nz) :: dtau
!
  type (InternalPencils) :: q
!
  integer :: ikappar,itau
!
  integer :: idiag_kappam=0,idiag_kappamax=0,idiag_kappamin=0
  integer :: idiag_taum=0,idiag_taumax=0,idiag_taumin=0
  integer :: idiag_taucoolm=0,idiag_taucoolmax=0,idiag_taucoolmin=0
  integer :: idiag_dts=0
!
  contains
!***********************************************************************
    subroutine register_special()
!
!  Configure pre-initialised (i.e. before parameter read) variables 
!  which should be know to be able to evaluate
!
!  14-jul-09/wlad: coded
!
      use Cdata
      use FArrayManager, only: farray_register_auxiliary
!
      if (lroot) call svn_id( &
           "$Id: newton_cooling.f90,v 1.2 2015/12/11 02:56:37 wlyra Exp $")
!
      call farray_register_auxiliary('kappar',ikappar)
      call farray_register_auxiliary('tau',itau)
!
      if (lroot) then
        open(4,file=trim(datadir)//'/index_special.pro',status='replace')
        write(4,*) 'tau    ',itau
        write(4,*) 'kappar ',ikappar
        close(4)
      endif
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  14-jul-09/wlad: coded
!
      use Cdata 
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (lroot) print*,'calling special before boundary'
      call special_before_boundary(f)
      if (lroot) print*,'done calling special before boundary'
!
      call keep_compiler_quiet(f)
!
    endsubroutine  init_special
!***********************************************************************
    subroutine special_before_boundary(f)
!
!  This subroutine calculates the full potential due to the turbulence.
!
!  03-oct-12/wlad: coded
!
      use EquationOfState, only: cs20,rho0,gamma_m1,get_cp1,get_cv1
      use Sub, only: grad,dot
      use General, only: notanumber
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(nx) :: rho,TT
      real :: TT0,rho01,lnTT0,cp1,cv1,kappa_cgs
      integer :: i,j      
!
      call get_cp1(cp1)
      call get_cv1(cv1)

      lnTT0=log(cs20*cp1/gamma_m1)
      TT0=exp(lnTT0)
      rho01=1./rho0
!
      do n=n1,n2; do m=m1,m2
!
        if (ldensity_nolog) then
          rho=f(l1:l2,m,n,irho)
        else
          rho=exp(f(l1:l2,m,n,ilnrho))
        endif
        TT = TT0 * (rho*rho01)**gamma_m1 * exp(f(l1:l2,m,n,iss)*cv1) 
!
! lnrho = log(f(i,m,n,irho))
! lnTT = lnTT0 + cv1*f(i,m,n,iss) + gamma_m1*(lnrho-lnrho0)
! TT_csg = exp(lnTT) * temperature_unit_cgs
!
        do i=1,nx 
           ! spits out opacity in cgs units
           call calc_opacity(TT(i)*unit_temperature,&
                rho(i)*unit_density,&
                kappa_cgs)
!
! kappa_code = kappa_cgs/ unit_kappa
! unit_kappa = 1/(unit_density*unit_length)
!
           f(i+nghost,m,n,ikappar) = kappa_cgs * (unit_density*unit_length)
!
        enddo
!
! Also precompute dtau for optical depth calculation. 
!
        dtau(:,m-nghost,n-nghost) = f(l1:l2,m,n,ikappar)*rho* x(l1:l2)/dy_1(m)
!
      enddo;enddo
!
!  Calculate now opacity
!
      call integrate_optical_depth(f)
!
!  Boundaries for limiter (gradient will be needed)
!
    endsubroutine special_before_boundary
!***********************************************************************
    subroutine integrate_optical_depth(f)
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(my) :: tau_below_local,tau_above_local
      real, dimension(my) :: tau_below,tau_above
      real :: tau_above_sum,tau_below_sum
      integer :: ipy_loc,i
!
      if (lmpicomm)  call calc_column_tau(f)
!
!  For every point, get stuff above and below -- assume meriodional to start with
!
      do n=n1,n2; do i=l1,l2
          !if tht < 90, gas parcel is above midplane. Sum stuff that is above is (lower thetas)
          !if tht > 90, gas parcel is below midplane. Sum stuff that is below it (higher thetas)
          !In fact, for all, take both all that is above and all that is below, take optical depth
          !  as the smaller of the two numbers. 
!
        do m=m1,m2
          tau_above_local(m) = sum(dtau(i-nghost,m1-nghost:m-nghost,n-nghost))
        enddo
        tau_above_sum=0.
        if (lmpicomm) then
          if (ipy/=0) then
            do ipy_loc=0,ipy-1
              tau_above_sum=tau_above_sum+tau_column(i-nghost,n-nghost,ipy_loc)
            enddo
          endif
        endif
        do m=m1,m2
          tau_above(m)=tau_above_sum+tau_above_local(m)
        enddo
!
        do m=m2,m1,-1
          tau_below_local(m) = sum(dtau(i-nghost,m-nghost:m2-nghost,n-nghost))
        enddo
        tau_below_sum=0.
        if (lmpicomm) then
          if (ipy/=nprocy-1) then
            do ipy_loc=ipy+1,nprocy-1
              tau_below_sum=tau_below_sum+tau_column(i-nghost,n-nghost,ipy_loc)
            enddo
          endif
        endif
        do m=m1,m2
          tau_below(m)=tau_below_sum+tau_below_local(m)
        enddo
!
        do m=m1,m2
           !if (f(i,m,n,ikappar) /= impossible) then !T >> T9
              f(i,m,n,itau)=min(tau_below(m),tau_above(m))
           !else
           !   !force it to cool! 
           !   f(i,m,n,itau)=1.
           !endif
        enddo
!
      enddo;enddo
!
    endsubroutine integrate_optical_depth
!***********************************************************************
    subroutine calc_column_tau(f)
!
!  Integrate the density vertically (meridionally): we need column 
!  density for the X-ray and cosmic ray attenuation factor.
!
      use Mpicomm,only: mpibcast_real,mpisend_real,mpirecv_real
!
      real, dimension(mx,my,mz,mfarray) :: f
      integer :: i,n,j
      real, dimension(nx,nz) :: tau_column_local,tmp
      real, dimension(nx,nz,0:nprocy-1) :: tau_column_yroot,tmp3
      integer, parameter :: y1tag=117
      integer, parameter :: y2tag=217
      integer :: partner,py,broadcaster,collector
!
      do n=1,nz; do i=1,nx
        tau_column_local(i,n)=sum(dtau(i,:,n))
      enddo; enddo
!      
!  Collect all y data
! 
      collector = ipx + ipz*nprocxy !root in y
      if (iproc == collector) then
        do py = 0,nprocy-1 
          partner = ipx + py*nprocx + ipz*nprocxy
          if (iproc == partner) then
            !data is local
            tau_column_yroot(:,:,py)=tau_column_local
          else
            call mpirecv_real(tmp,(/nx,nz/),partner,y1tag)
            tau_column_yroot(:,:,py)=tmp
          endif
        enddo
      else
        call mpisend_real(tau_column_local,(/nx,nz/),collector,y1tag)
      endif
!
!  Broadcast the data back to all y-processors
!
      broadcaster=ipx + ipz*nprocxy !root in y
      if (iproc == broadcaster) then
        do py=0,nprocy-1
          partner = ipx + py*nprocx + ipz*nprocxy
          if (iproc == partner) then 
            !data is local
            tau_column=tau_column_yroot
          else
            ! send to partner
             call mpisend_real(tau_column_yroot,(/nx,nz,nprocy/),partner,y2tag)
          endif
        enddo
      else
         call mpirecv_real(tau_column,(/nx,nz,nprocy/),broadcaster,y2tag)
      endif
!
    endsubroutine calc_column_tau
!***********************************************************************
    subroutine calc_cooling_time(f,p,taucool,taucool1)
!
      use EquationOfState, only: gamma1,gamma_m1,get_cp1,get_cv1,cs20,rho0
!
      real, dimension(mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      real, dimension(nx) :: tmp,tau_eff,Rd,OOK1
      real, dimension(nx), intent(out) :: taucool,taucool1
      real :: cp1,cp,cv1,cv
!
      integer :: i
!
      tmp=p%cp**1.5 * gamma1 * sqrt(gamma_m1) / (3.*sigmaSB) * p%rho * p%TT**(-2.5)
      tau_eff = 0.375*q%tau + .25*sqrt(3.) + .25*q%tau1
      Rd=tmp*tau_eff
      OOK1 = (x(l1:l2)*sinth(m))**(1.5)  !inverse Keplerian frequency
      taucool = Rd*OOK1
!
      if (taucool_floor/=impossible) & 
           where (taucool < taucool_floor) taucool=taucool_floor
!
      taucool1 = 1./taucool
!
    endsubroutine calc_cooling_time
!***********************************************************************
    subroutine pencil_criteria_special()
!
      lpenc_requested(i_TT)=.true.
      lpenc_requested(i_TT1)=.true.
      lpenc_requested(i_rho)=.true.
      lpenc_requested(i_cp)=.true.
      lpenc_requested(i_cp1)=.true.
      lpenc_requested(i_cv)=.true.
      lpenc_requested(i_ss)=.true.
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!   14-jul-09/wlad: coded
!
      real, dimension(mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      q%kappa=f(l1:l2,m,n,ikappar)
!
      q%tau=f(l1:l2,m,n,itau)
      q%tau1=1./q%tau
!
      call calc_cooling_time(f,p,q%taucool,q%taucool1)
!
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine calc_opacity(TT,rho,kk)
!
!  Piece-wise opacities from Bell et al. 1997.
!
!  01-aug-11/wlad: coded
!
      use General, only: notanumber
      use Mpicomm, only: stop_it
!
      real, intent(in) :: TT,rho
      real, intent(out) :: kk
      real :: k,a,b,logkk,logk
!
      real, save :: T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11
      logical, save :: lfirstcall=.true.
!
      if (lfirstcall) then
        T1=132. ; T2=170. ; T3=375. ; T4=390.
        T5=580. ; T6=680. ; T7=960. ; T8=1570. ; T9=3730.
        T10=1e4; T11=1e5
        lfirstcall=.false.
      endif
!
      if (TT < 0.0)        call inevitably_fatal_error("calc_opacity", "Negative temperature",.true.)
      if (rho < 0.0)       call inevitably_fatal_error("calc_opacity", "Negative density",.true.)
      if (notanumber(TT))  call inevitably_fatal_error("calc_opacity", "notanumber in temperature",.true.)
      if (notanumber(rho)) call inevitably_fatal_error("calc_opacity", "notanumber in density",.true.)
!
      if (TT <= T1) then
        k=2d-4 ; a=0 ; b= 2.1  ; kk = k * TT**b
      else if ((TT > T1) .and. (TT <= T2)) then
        k=3.   ; a=0 ; b=-0.01 ; kk = k * TT**b
      else if ((TT > T2) .and. (TT <= T3)) then
        k=0.01 ; a=0 ; b= 1.1  ; kk = k * TT**b
      else if ((TT > T3) .and. (TT <= T4)) then
        k=5d4  ; a=0 ; b=-1.5  ; kk = k * TT**b
      else if ((TT > T4) .and. (TT <= T5)) then
        k=0.1  ; a=0 ;  b= 0.7 ; kk = k * TT**b
      else if ((TT > T5) .and. (TT <= T6)) then
        k=2d15 ; a=0 ; b=-5.2  ; kk = k * TT**b
      else if ((TT > T6) .and. (TT <= T7)) then
        k=0.02 ; a=0 ; b= 0.8  ; kk = k * TT**b
      else if ((TT > T7) .and. (TT <= T8)) then
        logk=81.3010 ; a=1. ; b=-24.
        logkk = logk + a*alog10(rho) + b*alog10(TT)
        kk=10**(logkk)
      else if ((TT > T8) .and. (TT <= T9)) then
        k=1d-8 ; a=2./3 ; b=3.
        kk = k * rho**a * TT**b
      else if ((TT > T9) .and. (TT <= T10)) then
         logk=-36. ; a=1./3 ; b=10.
         logkk = logk + a*alog10(rho) + b*alog10(TT)
         kk=10**(logkk)
      else if ((TT > T10) .and. (TT <= T11)) then
         k=1.5d20 ; a=1. ; b=-2.5
         kk = k * rho**a * TT**b
      else if (TT > T11) then 
         kk=0.348
      endif
! 
    endsubroutine calc_opacity
!***********************************************************************
    subroutine read_special_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=special_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=special_init_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=special_init_pars)

    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=special_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=special_run_pars,ERR=99)
      endif
!
99    return
endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
      integer, intent(in) :: unit
!
      write(unit,NML=special_run_pars)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
      use Diagnostics
!
!  reads and registers print parameters relevant to special
!
!   14-jul-09/wlad: coded
!
      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Write information to index.pro
!
      if (lreset) then
        idiag_kappam=0;idiag_kappamax=0;idiag_kappamin=0 
        idiag_taum=0;idiag_taumax=0;idiag_taumin=0 
        idiag_taucoolm=0;idiag_taucoolmax=0;idiag_taucoolmin=0 
        idiag_dts=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'dts',idiag_dts)
!
        call parse_name(iname,cname(iname),cform(iname),'kappam',idiag_kappam)
        call parse_name(iname,cname(iname),cform(iname),'kappamax',idiag_kappamax)
        call parse_name(iname,cname(iname),cform(iname),'kappamin',idiag_kappamin)
!
        call parse_name(iname,cname(iname),cform(iname),'taum',idiag_taum)
        call parse_name(iname,cname(iname),cform(iname),'taumax',idiag_taumax)
        call parse_name(iname,cname(iname),cform(iname),'taumin',idiag_taumin)
!
        call parse_name(iname,cname(iname),cform(iname),'taucoolm',idiag_taucoolm)
        call parse_name(iname,cname(iname),cform(iname),'taucoolmax',idiag_taucoolmax)
        call parse_name(iname,cname(iname),cform(iname),'taucoolmin',idiag_taucoolmin)

      enddo
!
      if (lwr) then
        write(3,*) 'i_dts=',idiag_dts
!
        write(3,*) 'i_kappam=',idiag_kappam
        write(3,*) 'i_kappamax=',idiag_kappamax
        write(3,*) 'i_kappamin=',idiag_kappamin
!
        write(3,*) 'i_taum=',idiag_taum
        write(3,*) 'i_taumax=',idiag_taumax
        write(3,*) 'i_taumin=',idiag_taumin
!
        write(3,*) 'i_taucoolm=',idiag_taucoolm
        write(3,*) 'i_taucoolmax=',idiag_taucoolmax
        write(3,*) 'i_taucoolmin=',idiag_taucoolmin
      endif
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine special_calc_energy(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the 
!   momentum equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
      use Cdata
      use Diagnostics
      use EquationOfState, only: cs20,gamma_m1,gamma1
      use General, only: notanumber
!      
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
      real, dimension(nx) :: rr_cyl,TT_init,heating_rate
      integer :: i
!
!  Modified momentum equation
!
      if (laddheatingrate) then
         rr_cyl=x(l1:l2)*sinth(m)
         TT_init = cs20/(p%cp*gamma_m1)  * r_ref / rr_cyl
         heating_rate = p%cv * (p%TT-TT_init)*q%taucool1
!
         df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) - p%TT1 * heating_rate
         if (lfirst.and.ldt) dt1_max=max(dt1_max,q%taucool1/cdts)
       endif
!
      if (ldiagnos) then 
        if (idiag_kappam/=0)    call sum_mn_name(q%kappa,idiag_kappam)
        if (idiag_kappamax/=0)  call max_mn_name(q%kappa,idiag_kappamax)
        if (idiag_kappamin/=0)  call max_mn_name(-q%kappa,idiag_kappamin,lneg=.true.)
!
        if (idiag_taum/=0)   call sum_mn_name(q%tau,idiag_taum)
        if (idiag_taumax/=0) call max_mn_name(q%tau,idiag_taumax)
        if (idiag_taumin/=0) call max_mn_name(-q%tau,idiag_taumin,lneg=.true.)
!
        if (idiag_taucoolm/=0)   call sum_mn_name(q%taucool,idiag_taucoolm)
        if (idiag_taucoolmax/=0) call max_mn_name(q%taucool,idiag_taucoolmax)
        if (idiag_taucoolmin/=0) call max_mn_name(-q%taucool,idiag_taucoolmin,lneg=.true.)
!
        if (idiag_dts/=0) call max_mn_name(q%taucool1/cdts,idiag_dts,l_dt=.true.)
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_energy
!***********************************************************************
!
!***********************************************************************
!********************************************************************
!
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************
endmodule Special

