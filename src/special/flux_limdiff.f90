! $Id: baroclinic_run.f90 19193 2012-06-30 12:55:46Z wdobler $

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
  type InternalPencils
     real, dimension(nx,3) :: gkappa,glambda,glnkappa,glnlambda,gksi
     real, dimension(nx)   :: divflux,kappa,kappa1,lambda,lambda1,gTTgksi
  endtype InternalPencils
!
  type (InternalPencils) :: q
!
  integer :: ikappar,ilambda
!
  integer :: idiag_kappam=0,idiag_kappamax=0,idiag_kappamin=0
  integer :: idiag_lambdam=0,idiag_lambdamax=0,idiag_lambdamin=0
  integer :: idiag_divfluxm=0,idiag_divflux2m=0
  integer :: idiag_divfluxmax=0,idiag_divfluxmin=0
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
           "$Id: baroclinic_run.f90 19193 2012-06-30 12:55:46Z wdobler $")
!
      call farray_register_auxiliary('kappar',ikappar)
      call farray_register_auxiliary('lambda',ilambda)
!
      if (lroot) then
        open(4,file=trim(datadir)//'/index_special.pro',status='replace')
        write(4,*) 'ikappar=  ',ikappar
        write(4,*) 'ilambda=  ',ilambda
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
      call special_before_boundary(f)
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
      !use Boundcond, only: update_ghosts
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(nx,3) :: grho,gss,glnrho,glnTT
      real, dimension(nx) :: rho,TT,modglnTT,tmp,RR,rho1
      real :: TT0,rho01,lnTT0,cp1,cv1,kappa_cgs
      integer :: i,j      
!
!  Rossland mean opacity and flux limiter.
!
      call get_cp1(cp1)
      call get_cv1(cv1)

      lnTT0=log(cs20*cp1/gamma_m1)
      TT0=exp(lnTT0)
      rho01=1./rho0
!
      do n=n1,n2; do m=m1,m2
!
        rho=f(l1:l2,m,n,irho)
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
          f(i+l1-1,m,n,ikappar) = kappa_cgs * (unit_density*unit_length)
        enddo
      !enddo;enddo
!
!  Calculate now flux limiter.
!
      !do n=n1,n2; do m=m1,m2
!
        call grad(f,irho,grho)
        call grad(f,iss,gss)
        rho1=1./rho
        do j=1,3
          glnrho(:,j)=grho(:,j)*rho1
        enddo
        glnTT = gamma_m1*glnrho + cv1*gss
!
        call dot(glnTT,glnTT,tmp)
        modglnTT=sqrt(tmp)
!
        RR = 4*modglnTT*rho1/f(l1:l2,m,n,ikappar)
!
        f(l1:l2,m,n,ilambda) = (2 + RR)/(RR**2 + 3*RR + 6 )
!
      enddo; enddo
!
!  Boundaries for limiter (gradient will be needed)
!
      call update_ghosts_local(f)
!
    endsubroutine special_before_boundary
!***********************************************************************
    subroutine update_ghosts_local(f)
!
      use Mpicomm, only: initiate_isendrcv_bdry, finalize_isendrcv_bdry
!
      real, dimension(mx,my,mz,mfarray) :: f
      integer :: i
!
!  Do 's' boundary
!
      do i=1,nghost
        !call bc_sym_x(f,+1,'bot',ivar)
        f(l1-i,:,:,ikappar:ilambda)=f(l1+i,:,:,ikappar:ilambda)
        !call bc_sym_x(f,+1,'top',ivar)
        f(l2+i,:,:,ikappar:ilambda)=f(l2-i,:,:,ikappar:ilambda)
      enddo
!
      call initiate_isendrcv_bdry(f,ikappar)
      call finalize_isendrcv_bdry(f,ilambda)
!
      do i=1,nghost
        f(:,m1-i,:,ikappar:ilambda)=f(:,m1+i,:,ikappar:ilambda)
        f(:,m2+i,:,ikappar:ilambda)=f(:,m2-i,:,ikappar:ilambda)
      enddo
!
      do i=1,nghost
        f(:,:,n1-i,ikappar:ilambda)=f(:,:,n1+i,ikappar:ilambda)
        f(:,:,n2+i,ikappar:ilambda)=f(:,:,n2-i,ikappar:ilambda)
      enddo
!
    endsubroutine update_ghosts_local
!***********************************************************************
    subroutine pencil_criteria_special
!
      real, dimension(mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      lpenc_requested(i_gTT)=.true.
      lpenc_requested(i_TT)=.true.
      lpenc_requested(i_TT1)=.true.
      lpenc_requested(i_rho1)=.true.
      lpenc_requested(i_del2TT)=.true.
      lpenc_requested(i_glnTT)=.true.
      lpenc_requested(i_glnrho)=.true.
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!   14-jul-09/wlad: coded
!
      use Sub, only: grad,dot
!
      real, dimension(mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      integer :: j,i
!
      q%kappa=f(l1:l2,m,n,ikappar)
      q%kappa1=1./q%kappa
      call grad(f,ikappar,q%gkappa)
      do j=1,3 
        q%glnkappa(:,j) = q%kappa1*q%gkappa(:,j)
      enddo
!
      q%lambda=f(l1:l2,m,n,ilambda)
      q%lambda1=1./q%lambda
      call grad(f,ilambda,q%glambda)
      do j=1,3 
        q%glnlambda(:,j) = q%lambda1*q%glambda(:,j)
      enddo
!
      q%gksi = q%glnlambda + 3*p%glnTT - q%glnkappa - p%glnrho
!
      call dot(p%gTT,q%gksi,q%gTTgksi)
!      
      q%divflux = -16*sigmaSB*q%lambda*p%TT**3*p%rho1*q%kappa1*(p%del2TT+q%gTTgksi)
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
      real, intent(in) :: TT,rho
      real, intent(out) :: kk
      real :: k,a,b,logkk,logk
!
      real, save :: t1,t2,t3,t4,t5,t6,t7,t8,t9
      logical, save :: lfirstcall=.true.
!
      if (lfirstcall) then
        T1=132. ; T2=170. ; T3=375. ; T4=390.
        T5=580. ; T6=680. ; T7=960. ; T8=1570. ; T9=3730.
        lfirstcall=.false.
      endif
!
      kk=0.
!
      if (TT < 0.0) then
        call fatal_error("calc_opacity", "Negative temperature")
      endif
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
        k=1d33 ! why this?
      else if ((TT > T8) .and. (TT <= T9)) then
        k=1d-8 ; a=2./3 ; b=3.
        kk = k * rho**a * TT**b
      else
        if (lroot) then
          print*,'calc_opacity: density ',rho,' g/cm3'
          print*,'calc_opacity: temperature ',TT,' K. Higher '//&
               'than maximum allowed, ',T9
        endif
        call fatal_error("","")
      endif
!
    endsubroutine calc_opacity
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
!
        idiag_lambdam=0;idiag_lambdamax=0;idiag_lambdamin=0 
!
        idiag_divfluxm=0;idiag_divflux2m=0
        idiag_divfluxmin=0;idiag_divfluxmax=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'kappam',idiag_kappam)
        call parse_name(iname,cname(iname),cform(iname),'kappamax',idiag_kappamax)
        call parse_name(iname,cname(iname),cform(iname),'kappamin',idiag_kappamin)
!
        call parse_name(iname,cname(iname),cform(iname),'lambdam',idiag_lambdam)
        call parse_name(iname,cname(iname),cform(iname),'lambdamax',idiag_lambdamax)
        call parse_name(iname,cname(iname),cform(iname),'lambdamin',idiag_lambdamin)
!
        call parse_name(iname,cname(iname),cform(iname),'divfluxm',idiag_divfluxm)
        call parse_name(iname,cname(iname),cform(iname),'divflux2m',idiag_divflux2m)
        call parse_name(iname,cname(iname),cform(iname),'divfluxmax',idiag_divfluxmax)
        call parse_name(iname,cname(iname),cform(iname),'divfluxmin',idiag_divfluxmin)
      enddo
!
      if (lwr) then
        write(3,*) 'i_kappam=',idiag_kappam
        write(3,*) 'i_kappamax=',idiag_kappamax
        write(3,*) 'i_kappamin=',idiag_kappamin
!
        write(3,*) 'i_lambdam=',idiag_lambdam
        write(3,*) 'i_lambdamax=',idiag_lambdamax
        write(3,*) 'i_lambdamin=',idiag_lambdamin
!
        write(3,*) 'i_divfluxm=',idiag_divfluxm
        write(3,*) 'i_divflux2m=',idiag_divflux2m
        write(3,*) 'i_divfluxmax=',idiag_divfluxmax
        write(3,*) 'i_divfluxmin=',idiag_divfluxmin
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
!      
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
!  Modified momentum equation
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) - q%divflux*p%rho1*p%TT1
!
      if (ldiagnos) then 
        if (idiag_kappam/=0)    call sum_mn_name(q%kappa,idiag_kappam)
        if (idiag_kappamax/=0)  call max_mn_name(q%kappa,idiag_kappamax)
        if (idiag_kappamin/=0)  call max_mn_name(-q%kappa,idiag_kappamin,lneg=.true.)
!
        if (idiag_lambdam/=0)   call sum_mn_name(q%lambda,idiag_lambdam)
        if (idiag_lambdamax/=0) call max_mn_name(q%lambda,idiag_lambdamax)
        if (idiag_lambdamin/=0) call max_mn_name(-q%lambda,idiag_lambdamin,lneg=.true.)
!
        if (idiag_divfluxm/=0)  call sum_mn_name(q%divflux,idiag_divfluxm)
        if (idiag_divflux2m/=0) call sum_mn_name(q%divflux**2,idiag_divflux2m)
        if (idiag_divfluxmax/=0)  call max_mn_name(q%divflux,idiag_divfluxmax)
        if (idiag_divfluxmin/=0)  call max_mn_name(-q%divflux,idiag_divfluxmin,lneg=.true.)
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

