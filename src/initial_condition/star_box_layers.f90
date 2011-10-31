! $Id$
!
!  31-oct-2011/dintrans: coded (moved here my two subroutines from the 
!  entropy module after discussions during the 2011 PC meeting in Toulouse).
!  This module initializes the star-in-a-box setup with a central heating.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!***************************************************************
!
module InitialCondition
!
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
  use EquationOfState, only: gamma_inv, mpoly0, mpoly1, mpoly2, gamma_m1, cs20
!
  implicit none
!
  include 'initial_condition.h'
!
  real, pointer :: hcond0
  real :: wheat,luminosity,r_bcz,widthss,alpha_MLT

!
!!  integer :: dummy
!
!!  namelist /initial_condition_pars/ dummy
!
  contains
!***********************************************************************
    subroutine initial_condition_ss(f)
!
!  Initialize entropy for two superposed polytropes with a central heating
!
!  20-dec-06/dintrans: coded
!  28-nov-07/dintrans: merged with strat_heat_grav
!  05-jan-10/dintrans: now the gravity field is computed in gravity_r
!
      use EquationOfState, only: rho0, lnrho0, get_soundspeed,eoscalc, ilnrho_TT
      use Sub, only: step, interp1, erfunc
      use SharedVariables, only: get_shared_variable
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      integer, parameter   :: nr=100
      integer              :: i,l,iter,ierr
      real, dimension (nr) :: r,lnrho,temp,hcond,g,flux
      real                 :: u,r_mn,lnrho_r,temp_r,cs2,ss,lumi,Rgas
      real                 :: rhotop, rbot,rt_old,rt_new,rhobot,rb_old, &
        rb_new,crit,r_max,hcond1,hcond2
      real, dimension(:), pointer :: star_params
!
!  Get the needed parameters from the entropy module.
!
      call get_shared_variable('hcond0',hcond0,ierr)
      if (ierr/=0) call stop_it(" initial_condition_ss: "//&
           "there was a problem when getting hcond0")
      call get_shared_variable('star_params',star_params,ierr)
      if (ierr/=0) call stop_it(" initial_condition_ss: "//&
           "there was a problem when getting star_params")
      wheat=star_params(1)
      luminosity=star_params(2)
      r_bcz=star_params(3)
      widthss=star_params(4)
      alpha_MLT=star_params(5)
      print*, "star_params=", star_params
!
!  Define the radial grid r=[0,r_max] and needed radial profiles.
!
      if (nzgrid == 1) then
        r_max=sqrt(xyz1(1)**2+xyz1(2)**2)
      else
        r_max=sqrt(xyz1(1)**2+xyz1(2)**2+xyz1(3)**2)
      endif
      do i=1,nr
        r(i)=r_max*float(i-1)/(nr-1)
        u=r(i)/sqrt(2.)/wheat
        if (i==1) then
          flux(1)=0.
        else
          if (nzgrid==1) then
            lumi=luminosity*(1.-exp(-u**2))
            flux(i)=lumi/(2.*pi*r(i))
          else
            lumi=luminosity*(erfunc(u)-2.*u/sqrt(pi)*exp(-u**2))
            flux(i)=lumi/(4.*pi*r(i)**2)
          endif
        endif
      enddo
      Rgas=1.-gamma_inv
      g=-flux*Rgas*(mpoly0+1.)/hcond0
!
!  The bottom density value we want at r=r_bcz, actually given by rho0.
!
      rbot=rho0
      rt_old=0.1*rbot
      rt_new=0.12*rbot
!
!  Need to iterate for rhobot=1.
!  Produce first estimate.
!
      rhotop=rt_old
      call strat_heat(nr,r,flux,g,lnrho,temp,rhotop,rhobot)
      rb_old=rhobot
!
!  Next estimate.
!
      rhotop=rt_new
      call strat_heat(nr,r,flux,g,lnrho,temp,rhotop,rhobot)
      rb_new=rhobot
!
      do 10 iter=1,10
!
!  New estimate.
!
        rhotop=rt_old+(rt_new-rt_old)/(rb_new-rb_old)*(rbot-rb_old)
!
        crit=abs(rhotop/rt_new-1.)
        if (crit<=1e-4) goto 20
        call strat_heat(nr,r,flux,g,lnrho,temp,rhotop,rhobot)
!
!  Update new estimates.
!
        rt_old=rt_new ; rb_old=rb_new
        rt_new=rhotop ; rb_new=rhobot
 10   continue
 20   print*,'- iteration completed: rhotop,crit=',rhotop,crit
!
!  Redefine rho0 and lnrho0 (important for eoscalc!).
!
      rho0=rhotop
      lnrho0=log(rhotop)
      print*,'final rho0, lnrho0=',rho0, lnrho0
!      T0=cs20/gamma_m1
!      print*,'final rho0, lnrho0, T0=',rho0, lnrho0, T0
!
!  Compute rho and ss arrays from interpolations.
!
      do imn=1,ny*nz
        n=nn(imn)
        m=mm(imn)
        do l=l1,l2
          if (nzgrid == 1) then
            r_mn=sqrt(x(l)**2+y(m)**2)
          else
            r_mn=sqrt(x(l)**2+y(m)**2+z(n)**2)
          endif
          lnrho_r=interp1(r,lnrho,nr,r_mn)
          temp_r=interp1(r,temp,nr,r_mn)
          f(l,m,n,ilnrho)=lnrho_r
          call eoscalc(ilnrho_TT,lnrho_r,temp_r,ss=ss)
          f(l,m,n,iss)=ss
        enddo
      enddo
!
      if (lroot) then
        hcond1=(mpoly1+1.)/(mpoly0+1.)
        hcond2=(mpoly2+1.)/(mpoly0+1.)
        hcond=1.+(hcond1-1.)*step(r,r_bcz,-widthss) &
                +(hcond2-1.)*step(r,r_ext,widthss)
        hcond=hcond0*hcond
        print*,'--> writing initial setup to data/proc0/setup.dat'
        open(unit=11,file=trim(directory)//'/setup.dat')
        write(11,'(a1,a5,5a14)') '#','r','rho','ss','cs2','grav','hcond'
        do i=nr,1,-1
          call get_soundspeed(temp(i),cs2)
          call eoscalc(ilnrho_TT,lnrho(i),temp(i),ss=ss)
          write(11,'(f6.3,4e14.5,1pe14.5)') r(i),exp(lnrho(i)),ss, &
          cs2,g(i),hcond(i)
        enddo
        close(11)
      endif
!
    endsubroutine initial_condition_ss
!***********************************************************************
    subroutine strat_heat(nr,r,flux,g,lnrho,temp,rhotop,rhobot)
!
!  Compute the radial stratification for two superposed polytropic
!  layers and a central heating.
!
!  17-jan-07/dintrans: coded
!
    use EquationOfState, only: lnrho0
    use Sub, only: step, erfunc, interp1
!
    integer :: i,nr
    integer, dimension(1) :: i_ext
    real, dimension (nr) :: r,flux,g,lnrho,temp
    real :: dtemp,dlnrho,dr,rhotop,rhobot, &
            polyad,delad,fr_frac,fc_frac,fc,del,Rgas
!
!  Needed for computing a MLT stratification.
!
    polyad=1./gamma_m1
    delad=1.-gamma_inv
    fr_frac=delad*(mpoly0+1.)
    fc_frac=1.-fr_frac
    Rgas=1.-gamma_inv
!
!  Start from surface values for rho and temp.
!
    temp(nr)=cs20/gamma_m1 ; lnrho(nr)=alog(rhotop)
    dr=r(2)
    do i=nr-1,1,-1
      if (r(i+1) > r_ext) then
!
!  Isothermal exterior for r > r_ext.
!
        del=0.
      elseif (r(i+1) > r_bcz) then
!
!  Convection zone for r_bcz < r < r_ext: MLT solution if alpha_MLT/=0.
!
        fc=fc_frac*flux(i+1)
        del=delad+alpha_MLT*(fc/ &
                (exp(lnrho(i+1))*(gamma_m1*temp(i+1))**1.5))**.6666667
      else
!
!  Radiative zone for r < r_bcz.
!
        del=1./(mpoly1+1.)
      endif
      dtemp=-g(i+1)/Rgas*del
      dlnrho=-g(i+1)/(Rgas*temp(i+1))*(1.-del)
      temp(i)=temp(i+1)+dtemp*dr
      lnrho(i)=lnrho(i+1)+dlnrho*dr
    enddo
!
!  Find the value of rhobot at the bottom of convection zone.
!
!    lnrhobot=interp1(r,lnrho,nr,r_bcz)
!    rhobot=exp(lnrhobot)
!  new constraint for the iterative computation of stratification:
!  --> we impose Mtot=1 in r=(0,r_ext)
!
    i_ext=minloc(abs(r-r_ext))
    if (nzgrid == 1) then
      rhobot=sum(exp(lnrho(2:i_ext(1)-1))*r(2:i_ext(1)-1))+ &
        0.5*exp(lnrho(i_ext(1)))*r(i_ext(1))
      rhobot=rhobot*2.*pi*dr
    else
      rhobot=sum(exp(lnrho(2:i_ext(1)-1))*r(2:i_ext(1)-1)**2)+ &
        0.5*exp(lnrho(i_ext(1)))*r(i_ext(1))**2
      rhobot=rhobot*4.*pi*dr
    endif
    print*, 'total mass=', rhobot
!
    endsubroutine strat_heat
!***********************************************************************
!
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from noinitial_condition.f90 for any    **
!**  InitialCondition routines not implemented in this file        **
!**                                                                **
    include 'initial_condition_dummies.inc'
!********************************************************************
endmodule InitialCondition
