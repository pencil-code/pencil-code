! $Id: visc_const.f90,v 1.26 2004-01-30 14:26:50 dobler Exp $

!  This modules implements viscous heating and diffusion terms
!  here for cases 1) nu constant, 2) mu = rho.nu 3) constant and 

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Viscosity

  use Cparam
  use Cdata
  use Density

  implicit none

  character (len=labellen) :: ivisc='nu-const'
  integer :: i_dtnu=0
  ! dummy logical
  logical :: lvisc_first=.false.

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
      lvisc_shock=.false.
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_viscosity: constant viscosity'
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: visc_const.f90,v 1.26 2004-01-30 14:26:50 dobler Exp $")


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
!*******************************************************************
    subroutine rprint_viscosity(lreset,lwrite)
!
!  Writes ishock to index.pro file
!
!  24-nov-03/tony: adapted from rprint_ionization
!
      use Cdata
      use Sub
! 
      logical :: lreset
      logical, optional :: lwrite
      integer :: iname
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        i_dtnu=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if(lroot.and.ip<14) print*,'rprint_viscosity: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'dtnu',i_dtnu)
      enddo
!
!  write column where which ionization variable is stored
!
      if (present(lwrite)) then
        if (lwrite) then
          write(3,*) 'i_dtnu=',i_dtnu
          write(3,*) 'ihyper=',ihyper
          write(3,*) 'ishock=',ishock
          write(3,*) 'itest=',0
        endif
      endif
!   
      if(ip==0) print*,lreset  !(to keep compiler quiet)
    endsubroutine rprint_viscosity
!!***********************************************************************
    subroutine calc_viscosity(f)
      real, dimension (mx,my,mz,mvar+maux) :: f
      if(ip==0) print*,f  !(to keep compiler quiet)
    endsubroutine calc_viscosity
!!***********************************************************************
    subroutine calc_viscous_heat(f,df,glnrho,divu,rho1,cs2,TT1,shock)
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
      real, dimension (nx)   :: sij2, divu,shock
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
         df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + TT1*2.*nu*sij2*rho1
       case('nu-const', '2')
         if (headtt) print*,'viscous heating: ivisc=',ivisc
         df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + TT1*2.*nu*sij2
       case default
         if (lroot) print*,'ivisc=',trim(ivisc),' -- this could never happen'
         call stop_it("")
      endselect
      if(ip==0) print*,f,cs2,divu,glnrho,shock  !(keep compiler quiet)
    endsubroutine calc_viscous_heat

!***********************************************************************
    subroutine calc_viscous_force(f,df,glnrho,divu,rho1,shock,gshock)
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
      real, dimension (nx,3) :: glnrho,del2u,del6u,graddivu,fvisc,sglnrho,gshock
      real, dimension (nx) :: murho1,rho1,divu,shock
      integer :: i

      intent (in) :: f, glnrho, rho1
      intent (out) :: df,shock,gshock
!
!  viscosity operator
!  rho1 is pre-calculated in equ
!
      shock=0.
      gshock=0.

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
          maxdiffus=max_for_dt(maxdiffus,nu)

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
          maxdiffus=max_for_dt(maxdiffus,murho1)

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
            maxdiffus=max_for_dt(maxdiffus,nu)
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
          maxdiffus=max_for_dt(maxdiffus,nu)
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
!
!  set viscous time step
!
      if (ldiagnos.and.i_dtnu/=0) then
        call max_mn_name(spread(nu,1,nx)/dxmin**2/cdtvDim,i_dtnu,l_dt=.true.)
      endif
!
      if(ip==0) print*,divu  !(keep compiler quiet)
    end subroutine calc_viscous_force
!***********************************************************************

endmodule Viscosity
