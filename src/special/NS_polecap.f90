! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of special_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 3
! COMMUNICATED AUXILIARIES 3
!
!! PENCILS PROVIDED stress_ij(6), hijij, gijij
!***************************************************************
module Special
!
  use Cparam
  use Cdata
  use Initcond
  use General, only: keep_compiler_quiet
  use Messages, only: svn_id, fatal_error
!
  implicit none
!
  include '../special.h'
!
! Declare index of new variables in f array (if any).
!
  integer :: icond=0, icond_par=0, icond_perp=0, icond_hall=0
  real, dimension(nz) :: q_heating=0., glnrho_prof=0.
  real, dimension(mz) :: rho_prof=1., rho1_prof
  real, dimension(mx) :: B_r=0., B_theta=0.
  integer :: l_polecap=0
  real :: hcond0_kramers=0.
  logical :: lpotekhin_cond
!
! run parameters
!
  real :: r_polecap=0., d_heating=0., T_topobs=0., B_dipole=0., &
          rho_bot=0., rho_pow=3.          

  namelist /special_run_pars/ r_polecap, d_heating, T_topobs, B_dipole, &
                              rho_bot, rho_pow, hcond0_kramers, lpotekhin_cond
!
  integer, parameter :: nkramers=1
  real :: unit_cond=1.

  contains
!***********************************************************************
    subroutine initialize_special(f)
!
!  Called after reading parameters, but before the time loop.
!
!  17-dec-18/tony: coded
!
      use General, only: find_index
      use Deriv, only: der_z
      use Sub, only: step

      real, dimension (mx,my,mz,mfarray) :: f
!
      real :: T_bottom

      T_bottom  = fbcz(iTT,1)
      if (r_polecap>=x(l1)) l_polecap = find_index(x(l1:l2), r_polecap)
      q_heating = step(z(n1:n2),d_heating,0.1*d_heating) &
                  *sigmaSB * (T_topobs**4 - T_bottom**4)/d_heating

      B_r = B_dipole*cos(x)
!print*, rho_bot, xyz0(3), Lxyz(3), rho_pow
      rho_prof = rho_bot-.001*((z-xyz0(3))/Lxyz(3))**rho_pow
!print*, 'rhoprof=', rho_prof(1:20), rho_prof(mz-20:mz)
      rho1_prof = 1./rho_prof

      call der_z(alog(rho_prof),glnrho_prof)

      unit_cond = unit_density*unit_length**5/unit_time/unit_temperature

      call keep_compiler_quiet(f)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine register_special()
!
!  Set up indices for variables in special modules.
!
!  6-oct-03/tony: coded
!
      use FArrayManager
!
      if (lroot) call svn_id( &
           "$Id$")
!
!  Set indices for auxiliary variables.
!
      call farray_register_auxiliary('cond',icond, vector=3,communicated=.true.)
      icond_par=icond
      icond_perp=icond+1
      icond_hall=icond+2
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_special
!***********************************************************************
    subroutine pencil_criteria_special()
!
!  All pencils that this special module depends on are specified here.
!
      lpenc_requested(i_TT)=.true.
      lpenc_requested(i_TT1)=.true.
      lpenc_requested(i_glnTT)=.true.
      lpenc_requested(i_r_mn1)=.true.

    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine special_before_boundary(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      integer :: ll, nn
      real, dimension(nx) :: Krad, g2
      double precision :: SIGMA,CKAPPA,QJ,SIGMAT,CKAPPAT,QJT,SIGMAH,CKAPPAH,QJH
      real, parameter :: Zion=26,CMI=56,Zimp=0.1
!
      if (lpotekhin_cond) then

        do nn=n1,n2
!if (nn>=(n1+n2+1)/2.and.nn<=(n1+n2+11)/2) print*, 'rho=', rho_prof(nn)*unit_density
         do ll=l1,mx
!if (nn==n1) print*, 'B=', B_r(ll)*unit_magnetic
          if (ltemperature_nolog) then 
            call CONDCONV(f(ll,m1,nn,ilnTT)*unit_temperature,rho_prof(nn)*unit_density, &
                          B_r(ll)*unit_magnetic,dble(Zion),dble(CMI),dble(Zimp), &
                          SIGMA,CKAPPA,QJ,SIGMAT,CKAPPAT,QJT,SIGMAH,CKAPPAH,QJH)
          else
            call CONDCONV(exp(f(ll,m1,nn,ilnTT))*unit_temperature,rho_prof(nn)*unit_density, &
                          B_r(ll)*unit_magnetic,dble(Zion),dble(CMI),dble(Zimp), &
                          SIGMA,CKAPPA,QJ,SIGMAT,CKAPPAT,QJT,SIGMAH,CKAPPAH,QJH)
          endif
          f(ll,m1,nn,icond_par )=CKAPPA
          f(ll,m1,nn,icond_perp)=CKAPPAT
          f(ll,m1,nn,icond_hall)=CKAPPAH
if (ipx==0.and.nn==100) write(101+ipz,*) ll, x(ll), &
f(ll,m1,nn,ilnTT)*unit_temperature,rho_prof(nn)*unit_density,CKAPPA, CKAPPAT, CKAPPAH
        enddo; enddo

        f(ll,m1,nn,icond_par:icond_hall) = f(ll,m1,nn,icond_par:icond_hall)/unit_cond
      else
        f(ll,m1,nn,icond_par:icond_hall) = 0.
      endif
!
    endsubroutine special_before_boundary
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!  Calculate Special pencils.
!  Most basic pencils should come first, as others may depend on them.

      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine special_calc_energy(f,df,p)

      use Deriv, only: der

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: df,p

      real, dimension(nx) :: thdiff, ckappa, ckappat, heatcap_Fe
      real, dimension(nx) :: gckappa_z, gckappat_x
     
      if (l_polecap>0) then
        heatcap_Fe = 4.4e12*(1.+0.024*rho_prof(n)**(-2./3.)*p%TT)*rho_prof(n)   ! cv*rho
        if (ltemperature_nolog) then
          df(l1:l2,m,n,iTT) = df(l1:l2,m,n,iTT) + q_heating(n-nghost)/heatcap_Fe
        else
          df(l1:l2,m,n,iTT) = df(l1:l2,m,n,iTT) + q_heating(n-nghost)*p%TT1/heatcap_Fe
        endif
      endif

      thdiff=0.
      if (hcond0_kramers>0) call kramers_cond(p,thdiff)

      if (lpotekhin_cond) then

        ckappa=f(l1:l2,m,n,icond_par)
        ckappat=f(l1:l2,m,n,icond_perp)

        call der(f,icond_par,gckappa_z,3)
        call der(f,icond_perp,gckappat_x,1)

        thdiff = thdiff + (ckappa*(p%d2zlnTT + p%glnTT(:,3)**2) + ckappat*(p%d2xlnTT + p%glnTT(:,1)**2) &
                        + (ckappat*p%r_mn1 + gckappat_x)*p%glnTT(:,1) + gckappa_z*p%glnTT(:,3))/rho_prof(n)
        
      endif
      df(l1:l2,m,n,iTT) = df(l1:l2,m,n,iTT) + thdiff/heatcap_Fe
      
      call keep_compiler_quiet(f)

    endsubroutine special_calc_energy
!***********************************************************************
    subroutine kramers_cond(p,thdiff)

      use Sub, only: dot

      type (pencil_case) :: p
      real, dimension(nx) :: Krho1,thdiff

      real, dimension(nx) :: g2

      Krho1 = hcond0_kramers*rho1_prof(n)**(2.*nkramers+1.)*p%TT**(6.5*nkramers)   ! = K/rho
      !Krho1 = hcond0_kramers*exp(-p%lnrho*(2.*nkramers+1.)+p%lnTT*(6.5*nkramers))   ! = K/rho 
      !if(chimax_kramers>0.) &
      !  Krho1 = max(min(Krho1,chimax_kramers/p%cp1),chimin_kramers/p%cp1)
      call dot(-2.*nkramers*glnrho_prof(n-nghost)+(6.5*nkramers+1)*p%glnTT,p%glnTT,g2)
      thdiff = Krho1*(p%del2lnTT+g2)

    endsubroutine kramers_cond
!***********************************************************************
    subroutine special_boundconds(f,bc)

      real, dimension (mx,my,mz,mfarray) :: f
      type (boundary_condition) :: bc

      real, dimension(nx) :: TT,flux,Kheat
      integer :: i

      if (bc%bcname=='bbr'.and.bc%ivar==ilnTT) then
        if (bc%location==-3) then     ! bottom

          if (ltemperature_nolog) then
            TT=f(l1:l2,m1,n1,ilnTT)
            flux=sigmaSB*TT**4
          else
            TT=exp(f(l1:l2,m1,n1,ilnTT))
            flux=sigmaSB*TT**3
          endif

          Kheat = hcond0_kramers*rho1_prof(1)**(2.*nkramers+1.)*TT**(6.5*nkramers)
          if (lpotekhin_cond) Kheat=Kheat+f(l1:l2,m1,n1,icond_par)

          do i=1,nghost
            f(l1:l2,m1,n1-i,ilnTT) = f(l1:l2,m1,n1+i,ilnTT) - flux/Kheat*2*i*dz
          enddo
          bc%done=.true.

        elseif (bc%location==3) then  ! top

          if (ltemperature_nolog) then
            TT=f(l1:l2,m1,n2,ilnTT)
            flux=sigmaSB*TT**4
          else
            TT=exp(f(l1:l2,m1,n2,ilnTT))
            flux=sigmaSB*TT**3
          endif

          Kheat = hcond0_kramers*rho1_prof(nz)**(2.*nkramers+1.)*TT**(6.5*nkramers)
          if (lpotekhin_cond) Kheat=Kheat+f(l1:l2,m1,n2,icond_par)
               
          do i=1,nghost
            f(l1:l2,m1,n2+i,ilnTT) = f(l1:l2,m1,n2-i,ilnTT) + flux/Kheat*2*i*dz
          enddo
          bc%done=.true.

        endif
      endif

    endsubroutine special_boundconds
!***********************************************************************
    subroutine read_special_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=special_run_pars, IOSTAT=iostat)
!
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine read_special_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      iostat=0
! 
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=special_run_pars)
!
    endsubroutine write_special_run_pars
!***********************************************************************
!***********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************
endmodule Special
