! $Id: visc_var.f90,v 1.2 2003-04-28 10:47:29 tarek Exp $

!  This modules implements viscous heating and diffusion terms
!  here for cases 1) nu constant, 2) mu = rho.nu 3) constant and 

module Viscosity

  use Cparam
  use Cdata
  use Density

  implicit none

!  real :: nu=0.
  character (len=labellen) :: ivisc='nu-const'
  real :: nu_var, q_DJO=2., t0_DJO=0., nuf_DJO,ti_DJO=1.,tf_DJO
  real :: pp
  integer :: dummy, i_nu_var

  ! input parameters
  namelist /viscosity_init_pars/ dummy

  ! run parameters
  namelist /viscosity_run_pars/ ivisc,q_DJO,t0_DJO,nu,nuf_DJO,ti_DJO,tf_DJO

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
           "$Id: visc_var.f90,v 1.2 2003-04-28 10:47:29 tarek Exp $")


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
!    subroutine rprint_viscosity(lreset)
!!
!!  reads and registers print parameters relevant for compressible part
!!
!!   3-may-02/axel: coded
!!  27-may-02/axel: added possibility to reset list
!!  19-mar-03/tarek: adapted  fro, rprint_density
!      use Sub
!!
!      integer :: iname
!      logical :: lreset
!!
!!  reset everything in case of reset
!!  (this needs to be consistent with what is defined above!)
!!
!      if (lreset) then
!        i_nu_var=nu
!      endif
!!
!!  iname runs through all possible names that may be listed in print.in
!!
!      if(lroot.and.ip<14) print*,'run through parse list'
!      do iname=1,nname
!        call parse_name(iname,cname(iname),cform(iname),'nu_var',i_nu_var)
!      enddo
!!
!!  write column where which viscosity variables stored
!!
!      write(3,*) 'q_DJO=',q_DJO
!      write(3,*) 't0_DJO=',t0_DJO
!!
!    endsubroutine rprint_viscosity

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

      real, dimension (mx,my,mz,mvar) :: f,df
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
         call stop_it("") !"
      endselect
      if(ip==0) print*,f,cs2,divu,glnrho  !(keep compiler quiet)
    endsubroutine calc_viscous_heat

!***********************************************************************
    subroutine calc_viscous_force(f,df,glnrho,divu,rho1)
!
      use Cdata
      use Mpicomm
      use Sub

      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3) :: glnrho,del2u, graddivu,fvisc,sglnrho
      real, dimension (nx) :: murho1,rho1,divu
      integer :: i

      intent (in) :: f, glnrho, rho1
      intent (out) :: df
!
!  viscosity operator
!  rho1 is pre-calculated in equ
!
!      if (nu /= 0.) then
       
        select case (ivisc)
        
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

        case('nu-DJO')
         !  Ditlivsen et at. show that
         !  E=k^q*psi((t*k^(3+q)/2),nu t^{-(1-q)/(3+q)}), t_code=t-t0
         !  Viscosity is choosen here such that second argument remains constant
         !  with time. q and t0  must be guest apriori. q=3.67 might be a good 
         !  choice.
         
         !  if nu is given in namelist this will correspond to nu at t=1. 
         !  alternativly nu_min=nu(t_max) and t_max can be set in namelist 
         ! 
         !  viscous force: nu*(del2u+graddivu/3+2S.glnrho)
         !  -- the correct expression for nu=const
         !
        

           pp = - (1.-q_DJO)/(3.+q_DJO)  
           if (t.lt.ti_DJO) then 
             nu_var = nu
           else 
             nu_var= nu*((ti_DJO + t0_DJO)/(t + t0_DJO))**pp
           endif 
          tf_DJO = (nu/nuf_DJO)**(1./pp)*(ti_DJO + t0_DJO) - t0_DJO
          if (lroot.and.lfirstpoint.and.lout) call write_viscosity
          if (headtt) print*,'Using DJO variable viscosity. with '
          if (headtt.and.(ip.lt.6)) then   
                print*,'VISC_VAR: nu, nuf,ti,tf,t0,q'
                print*,nu,nuf_DJO,ti_DJO,tf_DJO,t0_DJO,q_DJO
          endif 
          if (headtt) print*,'viscous force: nu_var*(del2u+graddivu/3+2S.glnrho)'
          call del2v_etc(f,iuu,del2u,GRADDIV=graddivu)
          if(ldensity) then
            call multmv_mn(sij,glnrho,sglnrho)
            fvisc=2.*nu_var*sglnrho+nu_var*(del2u+1./3.*graddivu)
            maxdiffus=amax1(maxdiffus,nu_var)
          else
            if(lfirstpoint) &
                 print*,"ldensity better be .true. for ivisc=",ivisc
          endif

        case('nu-DJO-nobulk')
           
          pp = - (1.-q_DJO)/(3.+q_DJO)  
          if (t.lt.ti_DJO) then 
             nu_var = nu
           else 
             nu_var= nu*((ti_DJO + t0_DJO)/(t + t0_DJO))**pp
           endif 
          tf_DJO = (nu/nuf_DJO)**(1./pp)*(ti_DJO + t0_DJO) - t0_DJO
          if (lroot.and.lfirstpoint.and.lout) call write_viscosity
          if (headtt.and.(ip.lt.6)) then   
                print*,'VISC_VAR: nu, nuf,ti,tf,t0,q'
                print*,nu,nuf_DJO,ti_DJO,tf_DJO,t0_DJO,q_DJO
          endif 
          if (headtt) print*,'VARIABLE viscous force: nu_var*del2v'
          call del2v(f,iuu,del2u)
          fvisc=nu_var *del2u

        case default
        !
        !  Catch unknown values
        !
        if (lroot) print*, 'Program compiled with VISCOSITY = visc_var '
        if (lroot) print*, 'No such such value for ivisc: ', trim(ivisc)
        call stop_it('calc_viscous_forcing')
        
      endselect

        df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)+fvisc
!      else ! (nu=0)
!        if (headtt) print*,'no viscous force: (nu=0)'
!      endif

      if(ip==0) print*,divu  !(keep compiler quiet)
    end subroutine calc_viscous_force


    subroutine write_viscosity
     open(1,file=trim(datadir)//'/visc_var.dat',position='append')
     write(1,*) it,t,nu_var
     close(1)
    end subroutine
endmodule Viscosity
