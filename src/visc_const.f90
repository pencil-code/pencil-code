! $Id: visc_const.f90,v 1.11 2003-06-16 04:41:11 brandenb Exp $

!  This modules implements viscous heating and diffusion terms
!  here for cases 1) nu constant, 2) mu = rho.nu 3) constant and 

module Viscosity

  use Cparam
  use Cdata
  use Density

  implicit none

  logical :: lvisc_shock=.false.
  character (len=labellen) :: ivisc='nu-const'

  ! input parameters
  integer :: dummy1
  namelist /viscosity_init_pars/ dummy1

  ! run parameters
  namelist /viscosity_run_pars/ nu, ivisc
 
  contains

!***********************************************************************
    subroutine register_viscosity()
!
!  19-nov-02/tony: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_viscosity called twice')
      first = .false.
!
      lviscosity = .true.
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_viscosity: constant viscosity'
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: visc_const.f90,v 1.11 2003-06-16 04:41:11 brandenb Exp $")


! Following test unnecessary as no extra variable is evolved
!
!      if (nvar > mvar) then
!        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
!        call stop_it('Register_viscosity: nvar > mvar')
!      endif
!
    endsubroutine register_viscosity
!***********************************************************************
    subroutine initialize_viscosity()
!
!  20-nov-02/tony: coded

      use Cdata

      if (nu /= 0. .and. (ivisc=='nu-const')) then
         lneed_sij=.true.
         lneed_glnrho=.true.
      endif

    endsubroutine initialize_viscosity

!***********************************************************************
    subroutine calc_viscous_heat(f,df,glnrho,divu,rho1,cs2,TT1)
!
!  calculate viscous heating term for right hand side of entropy equation
!
!  20-nov-02/tony: coded
!
      use Cdata
      use Mpicomm
      use Sub

      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx)   :: rho1,TT1,cs2
      real, dimension (nx)   :: sij2, divu
      real, dimension (nx,3) :: glnrho
!
!  traceless strainmatrix squared
!
      call multm2_mn(sij,sij2)
!
      select case(ivisc)
       case ('simplified', '0')
         if (headtt) print*,'no heating: ivisc=',ivisc
       case('rho_nu-const', '1')
         if (headtt) print*,'viscous heating: ivisc=',ivisc
         df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) + TT1*2.*nu*sij2*rho1
       case('nu-const', '2')
         if (headtt) print*,'viscous heating: ivisc=',ivisc
         df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) + TT1*2.*nu*sij2
       case default
         if (lroot) print*,'ivisc=',trim(ivisc),' -- this could never happen'
         call stop_it("")
      endselect
      if(ip==0) print*,f,cs2,divu,glnrho  !(keep compiler quiet)
    endsubroutine calc_viscous_heat

!***********************************************************************
    subroutine calc_viscous_force(f,df,glnrho,divu,rho1)
!
!  calculate viscous heating term for right hand side of entropy equation
!
!  20-nov-02/tony: coded
!
      use Cdata
      use Mpicomm
      use Sub

      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: glnrho,del2u,del6u,graddivu,fvisc,sglnrho
      real, dimension (nx) :: murho1,rho1,divu
      integer :: i

      intent (in) :: f, glnrho, rho1
      intent (out) :: df
!
!  viscosity operator
!  rho1 is pre-calculated in equ
!
      if (nu /= 0.) then
        select case (ivisc)

        case ('simplified', '0')
          !
          !  viscous force: nu*del2v
          !  -- not physically correct (no momentum conservation), but
          !  numerically easy and in most cases qualitatively OK
          !
          if (headtt) print*,'viscous force: nu*del2v'
          call del2v(f,iuu,del2u)
          fvisc=nu*del2u
          maxdiffus=amax1(maxdiffus,nu)

        case('rho_nu-const', '1')
          !
          !  viscous force: mu/rho*(del2u+graddivu/3)
          !  -- the correct expression for rho*nu=const (=rho0*nu)
          !
          if (headtt) print*,'viscous force: mu/rho*(del2u+graddivu/3)'
          if (.not.ldensity) &
               print*, "ldensity better be .true. for ivisc='rho_nu-const'"
          murho1=(nu*rho0)*rho1  !(=mu/rho)
          call del2v_etc(f,iuu,del2u,GRADDIV=graddivu)
          do i=1,3
            fvisc(:,i)=murho1*(del2u(:,i)+1./3.*graddivu(:,i))
          enddo
          maxdiffus=amax1(maxdiffus,murho1)

        case('nu-const')
          !
          !  viscous force: nu*(del2u+graddivu/3+2S.glnrho)
          !  -- the correct expression for nu=const
          !
          if (headtt) print*,'viscous force: nu*(del2u+graddivu/3+2S.glnrho)'
          call del2v_etc(f,iuu,del2u,GRADDIV=graddivu)
          if(ldensity) then
            call multmv_mn(sij,glnrho,sglnrho)
            fvisc=2*nu*sglnrho+nu*(del2u+1./3.*graddivu)
            maxdiffus=amax1(maxdiffus,nu)
          else
            if(lfirstpoint) &
                 print*,"ldensity better be .true. for ivisc='nu-const'"
          endif

       case ('hyper6')
          !
          !  viscous force: nu*del6v
          !
          if (headtt) print*,'viscous force: nu*del6v'
          call del6v(f,iuu,del6u)
          fvisc=nu*del6u
          maxdiffus=amax1(maxdiffus,nu)

        case default
        !
        !  Catch unknown values
        !
        if (lroot) print*, 'No such such value for ivisc: ', trim(ivisc)
        call stop_it('calc_viscous_forcing')

        endselect

        df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)+fvisc
      else ! (nu=0)
        if (headtt) print*,'no viscous force: (nu=0)'
      endif

      if(ip==0) print*,divu  !(keep compiler quiet)
    end subroutine calc_viscous_force
endmodule Viscosity
