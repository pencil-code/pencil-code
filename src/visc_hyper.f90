! $Id: visc_hyper.f90,v 1.1 2003-12-12 08:04:38 nilshau Exp $

!  This modules implements viscous heating and diffusion terms
!  here for third order hyper viscosity 

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 3
!
!***************************************************************

module Viscosity

  use Cparam
  use Cdata
  use Density

  implicit none

  character (len=labellen) :: ivisc='hyper3'
  real :: maxeffectivenu
  logical :: lvisc_first=.false.

  ! input parameters
  integer :: dummy
  namelist /viscosity_init_pars/ dummy

  ! run parameters
  namelist /viscosity_run_pars/ nu, lvisc_first,ivisc
 
  contains

!***********************************************************************
    subroutine register_viscosity()
!
!  19-nov-02/tony: coded
!  24-nov-03/nils: adapted from visc_shock
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
      lvisc_hyper = .true.
!
      ihyper = mvar + naux + 1
      naux = naux + 3 
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_viscosity: hyper viscosity nvar = ', nvar
        print*, 'ihyper = ', ihyper
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: visc_hyper.f90,v 1.1 2003-12-12 08:04:38 nilshau Exp $")
!
! Check we aren't registering too many auxiliary variables
!
      if (naux > maux) then
        if (lroot) write(0,*) 'naux = ', naux, ', maux = ', maux
        call stop_it('register_viscosity: naux > maux')
      endif
!
!  Writing files for use with IDL
!
      if (naux < maux) aux_var(aux_count)=',hyper $'
      if (naux == maux) aux_var(aux_count)=',hyper'
      aux_count=aux_count+3
      if (lroot) write(15,*) 'hyper = fltarr(mx,my,mz)*one'
!
    endsubroutine register_viscosity
!***********************************************************************
    subroutine initialize_viscosity()
!
!  20-nov-02/tony: coded
!  24-nov-03/nils: adapted from visc_shock
!
       use Cdata
!
       lneed_sij=.false.
       lneed_glnrho=.false.
!
        if (headtt.and.lroot) print*,'viscosity: nu=',nu

    endsubroutine initialize_viscosity
!*******************************************************************
    subroutine rprint_viscosity(lreset,lwrite)
!
!  Writes ihyper to index.pro file
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
!      if (lreset) then
!        i_TTm=0
!      endif
!
!  iname runs through all possible names that may be listed in print.in
!
!      if(lroot.and.ip<14) print*,'rprint_ionization: run through parse list'
!      do iname=1,nname
!        call parse_name(iname,cname(iname),cform(iname),'yHm',i_yHm)
!      enddo
!
!  write column where which ionization variable is stored
!
      if (present(lwrite)) then
        if (lwrite) then
          write(3,*) 'ihyper=',ihyper
          write(3,*) 'ishock=',ishock
        endif
      endif
!   
      if(ip==0) print*,lreset  !(to keep compiler quiet)
    endsubroutine rprint_viscosity
!***********************************************************************
    subroutine calc_viscosity(f)
!
!  calculate viscous heating term for right hand side of entropy equation
!
!  24-nov-03/nils: coded
!
      use IO
      use Mpicomm
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,3) :: tmp
!
      if(headt) print*,'calc_viscosity: dxmin=',dxmin
!
!  calculate grad4 grad div u
!
      if (ivisc .eq. 'hyper3') then
        call divgrad_2nd(f,tmp,iux)
        f(:,:,:,ihyper:ihyper+2)=tmp
        call del2v_2nd(f,tmp,ihyper)
        f(:,:,:,ihyper:ihyper+2)=tmp
        call del2v_2nd(f,tmp,ihyper)
        f(:,:,:,ihyper:ihyper+2)=tmp
      elseif (ivisc .eq. 'hyper2') then
        call divgrad_2nd(f,tmp,iux)
        f(:,:,:,ihyper:ihyper+2)=tmp
        call del2v_2nd(f,tmp,ihyper)
        f(:,:,:,ihyper:ihyper+2)=tmp
      else
        call stop_it('visc_hyper:no such ivisc')          
      endif
!
!  max effective nu is the max of shock viscosity and the ordinary viscosity
!
      !maxeffectivenu = maxval(f(:,:,:,ihyper))*nu
!
    endsubroutine calc_viscosity
!***********************************************************************
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
      real, dimension (nx) :: rho1,TT1,cs2
      real, dimension (nx) :: sij2, divu,shock
      real, dimension (nx,3) :: glnrho

!
!  traceless strain matrix squared
!
      !call multm2_mn(sij,sij2)
!      if (headtt) print*,'viscous heating: ',ivisc

      !df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + TT1 * &
       !    (2.*nu*sij2  & 
        !   + nu_shock * shock * divu**2)

      !maxheating=amax1(maxheating,df(l1:l2,m,n,iss))
!
      !if(ip==0) print*,glnrho,rho1,cs2 !(to keep compiler quiet)
    endsubroutine calc_viscous_heat

!***********************************************************************
    subroutine calc_viscous_force(f,df,glnrho,divu,rho1,shock,gshock)
!
!  calculate viscous force
!
!  24-nov-03/nils: coded
!
      use Cdata
      use Mpicomm
      use Sub

      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: hyper
      real, dimension (nx,3) :: glnrho,del6u,fvisc,gshock,del4u
      real, dimension (nx) :: rho1,divu,shock
!
      intent (out) :: df
!
!
      if (nu /= 0.) then
        if (ivisc .eq. 'hyper3') then
          !  viscous force:nu*(del6u+del4*graddivu/3)
          !  (Assuming rho*nu=const)
          hyper=f(l1:l2,m,n,ihyper:ihyper+2)/3.
          if (headtt) print*,'viscous force: nu*(del6u+del4*graddivu/3)'
          call del6v(f,iuu,del6u)
          fvisc=nu*(del6u+hyper)
          maxdiffus=amax1(maxdiffus,nu)
          df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)+fvisc
        elseif (ivisc .eq. 'hyper2') then
          !  viscous force:nu*(del4u+del2*graddivu/3)
          !  (Assuming rho*nu=const)
          hyper=f(l1:l2,m,n,ihyper:ihyper+2)/3.
          if (headtt) print*,'viscous force: nu*(del4u+del2*graddivu/3)'
          call del4v(f,iuu,del4u)
          fvisc=nu*(del4u+hyper)
          maxdiffus=amax1(maxdiffus,nu)
          df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)+fvisc
        else
          call stop_it('visc_hyper:no such ivisc')  
        endif
      else ! (nu=0)
        if (headtt.and.lroot) print*,'no viscous force: (nu=0)'
      endif
!
      if(ip==0) print*,rho1,divu,shock,gshock !(to keep compiler quiet)
    end subroutine calc_viscous_force
!***********************************************************************
    subroutine derij_2nd(var,tmp,j)
!
!  24-nov-03/nils: coded
!
      use Cdata
!
      real, dimension (mx,my,mz) :: var
      real, dimension (mx,my,mz) :: tmp
      integer :: j
!
      intent (in) :: var,j
      intent (out) :: tmp
!
      tmp=0.

      if (j==1 .and. nxgrid/=1) then
          tmp(     1,:,:) = (- 3.*var(1,:,:) &
                            + 4.*var(2,:,:) &
                            - 1.*var(3,:,:))/(2.*dx) 
          tmp(2:mx-1,:,:) = (- 1.*var(1:mx-2,:,:) &
                            + 1.*var(3:mx  ,:,:))/(2.*dx) 
          tmp(    mx,:,:) = (+ 1.*var(mx-2,:,:) &
                            - 4.*var(mx-1,:,:) &
                            + 3.*var(mx  ,:,:))/(2.*dx) 
      endif
!
      if (j==2 .and. nygrid/=1) then
          tmp(:,     1,:) = (- 3.*var(:,1,:) &
                            + 4.*var(:,2,:) &
                            - 1.*var(:,3,:))/(2.*dy) 
          tmp(:,2:my-1,:) = (- 1.*var(:,1:my-2,:) &
                            + 1.*var(:,3:my  ,:))/(2.*dy) 
          tmp(:,    my,:) = (+ 1.*var(:,my-2,:) &
                            - 4.*var(:,my-1,:) &
                            + 3.*var(:,my  ,:))/(2.*dy) 
      endif
!
      if (j==3 .and. nzgrid/=1) then
          tmp(:,:,     1) = (- 3.*var(:,:,1) &
                            + 4.*var(:,:,2) &
                            - 1.*var(:,:,3))/(2.*dz) 
          tmp(:,:,2:mz-1) = (- 1.*var(:,:,1:mz-2) &
                            + 1.*var(:,:,3:mz  ))/(2.*dz) 
          tmp(:,:,    mz) = (+ 1.*var(:,:,mz-2) &
                            - 4.*var(:,:,mz-1) &
                            + 3.*var(:,:,mz  ))/(2.*dz) 
      endif
!
    end subroutine derij_2nd
!***********************************************************************
    subroutine del2v_2nd(f,del2f,k)
!
!  24-nov-03/nils: adapted from del2v
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,3) :: del2f
      real, dimension (mx,my,mz) :: tmp
      integer :: i,k,k1
!
      intent (in) :: f, k
      intent (out) :: del2f
!
      del2f=0.
!
!  do the del2 diffusion operator
!
      k1=k-1
      do i=1,3
        call del2_2nd(f,tmp,k1+i)
        del2f(:,:,:,i)=tmp
      enddo
!
    end subroutine del2v_2nd
!***********************************************************************
    subroutine del2_2nd(f,del2f,k)
!
!  calculate del2 of a scalar, get scalar
!  24-nov-03/nils: adapted from del2
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: del2f,d2fd
      integer :: i,k,k1
!
      intent (in) :: f, k
      intent (out) :: del2f
!
      k1=k-1
      call der2_2nd(f,d2fd,k,1)
      del2f=d2fd
      call der2_2nd(f,d2fd,k,2)
      del2f=del2f+d2fd
      call der2_2nd(f,d2fd,k,3)
      del2f=del2f+d2fd
!
    endsubroutine del2_2nd
!***********************************************************************
    subroutine der2_2nd(f,der2f,i,j)
!
!  24-nov-03/nils: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: der2f
      integer :: i,j
!
      intent (in) :: f,i,j
      intent (out) :: der2f
!
      der2f=0.
!
      if (j==1 .and. nxgrid/=1) then
        der2f(1     ,:,:) =    (+ 2.*f(1,:,:,i) &
                               - 5.*f(2,:,:,i) &
                               + 4.*f(3,:,:,i) &
                               - 1.*f(4,:,:,i) ) &
                               / (dx**2) 
        der2f(2:mx-1,:,:) =    (+ 1.*f(1:mx-2,:,:,i) &
                               - 2.*f(2:mx-1,:,:,i) &
                               + 1.*f(3:mx  ,:,:,i) ) &
                               / (dx**2) 
        der2f(mx    ,:,:) =    (+ 2.*f(mx  ,:,:,i) &
                               - 5.*f(mx-1,:,:,i) &
                               + 4.*f(mx-2,:,:,i) &
                               - 1.*f(mx-3,:,:,i) ) &
                               / (dx**2) 
      endif
!
     if (j==2 .and. nygrid/=1) then
        der2f(:,1     ,:) =    (+ 2.*f(:,1,:,i) &
                               - 5.*f(:,2,:,i) &
                               + 4.*f(:,3,:,i) &
                               - 1.*f(:,4,:,i) ) &
                               / (dy**2) 
        der2f(:,2:my-1,:) =    (+ 1.*f(:,1:my-2,:,i) &
                               - 2.*f(:,2:my-1,:,i) &
                               + 1.*f(:,3:my  ,:,i) ) &
                               / (dy**2) 
        der2f(:,my    ,:) =    (+ 2.*f(:,my  ,:,i) &
                               - 5.*f(:,my-1,:,i) &
                               + 4.*f(:,my-2,:,i) &
                               - 1.*f(:,my-3,:,i) ) &
                               / (dy**2) 
      endif
!
     if (j==3 .and. nzgrid/=1) then
        der2f(:,:,1     ) =    (+ 2.*f(:,:,1,i) &
                               - 5.*f(:,:,2,i) &
                               + 4.*f(:,:,3,i) &
                               - 1.*f(:,:,4,i) ) &
                               / (dz**2) 
        der2f(:,:,2:mz-1) =    (+ 1.*f(:,:,1:mz-2,i) &
                               - 2.*f(:,:,2:mz-1,i) &
                               + 1.*f(:,:,3:mz  ,i) ) &
                               / (dz**2) 
        der2f(:,:,mz    ) =    (+ 2.*f(:,:,mz  ,i) &
                               - 5.*f(:,:,mz-1,i) &
                               + 4.*f(:,:,mz-2,i) &
                               - 1.*f(:,:,mz-3,i) ) &
                               / (dz**2) 
      endif
!
    endsubroutine der2_2nd
!***********************************************************************
    subroutine divgrad_2nd(f,divgradf,k)
!
!  24-nov-03/nils: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,3) :: divgradf
      real, dimension (mx,my,mz) :: tmp,tmp2
      integer :: k,k1
!
      intent (in) :: f, k
      intent (out) :: divgradf
!
      divgradf=0.
!
      k1=k-1
!
      call der2_2nd(f,tmp,k1+1,1)
      divgradf(:,:,:,1)=tmp
      call derij_2nd(f(:,:,:,k1+2),tmp,1)
      call derij_2nd(tmp,tmp2,2)
      divgradf(:,:,:,1)=divgradf(:,:,:,1)+tmp2
      call derij_2nd(f(:,:,:,k1+3),tmp,1)
      call derij_2nd(tmp,tmp2,3)
      divgradf(:,:,:,1)=divgradf(:,:,:,1)+tmp2
!
      call der2_2nd(f,tmp,k1+2,2)
      divgradf(:,:,:,2)=tmp
      call derij_2nd(f(:,:,:,k1+1),tmp,1)
      call derij_2nd(tmp,tmp2,2)
      divgradf(:,:,:,2)=divgradf(:,:,:,2)+tmp2
      call derij_2nd(f(:,:,:,k1+3),tmp,2)
      call derij_2nd(tmp,tmp2,3)
      divgradf(:,:,:,2)=divgradf(:,:,:,2)+tmp2
!
      call der2_2nd(f,tmp,k1+3,3)
      divgradf(:,:,:,3)=tmp
      call derij_2nd(f(:,:,:,k1+1),tmp,1)
      call derij_2nd(tmp,tmp2,3)
      divgradf(:,:,:,3)=divgradf(:,:,:,3)+tmp2
      call derij_2nd(f(:,:,:,k1+2),tmp,2)
      call derij_2nd(tmp,tmp2,3)
      divgradf(:,:,:,3)=divgradf(:,:,:,3)+tmp2
!
    endsubroutine divgrad_2nd
!***********************************************************************

endmodule Viscosity
