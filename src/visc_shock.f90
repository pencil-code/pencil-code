! $Id: visc_shock.f90,v 1.4 2002-11-26 10:15:57 mee Exp $

!  This modules implements viscous heating and diffusion terms
!  here for shock viscosity nu_total = nu + nu_shock * dx * smooth(max5(-(div u)))) 

module Viscosity

  use Cparam
  use Cdata
  use Density

  implicit none

!  real :: nu=0.
  real :: nu_shock = 0.
  real, dimension(mx,my,mz) :: shock_characteristic
  character (len=labellen) :: ivisc=''
  integer :: icalculated = -1

  ! input parameters
  integer :: dummy
  namelist /viscosity_init_pars/ dummy

  ! run parameters
  namelist /viscosity_run_pars/ nu, nu_shock
 
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
           "$Id: visc_shock.f90,v 1.4 2002-11-26 10:15:57 mee Exp $")


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

 !     if (nu /= 0. .and. (ivisc=='nu-const')) then
         lneed_sij=.true.
         lneed_glnrho=.true.
 !     endif

    endsubroutine initialize_viscosity
!***********************************************************************
    subroutine calc_viscosity(f)
!
!  calculate viscous heating term for right hand side of entropy equation
!
!  23-nov-02/tony: coded
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: tmp

      if (nu_shock /= 0.) then
         call shock_divu(f,shock_characteristic)

         shock_characteristic=amax1(0., -shock_characteristic)
         
         call shock_max5(shock_characteristic,tmp)
         call shock_smooth(tmp,shock_characteristic)

!ajwm Shouldn't be dx in all directions! min(dx,dy,dz) ?         
         shock_characteristic = shock_characteristic * dx**2 
      else
         shock_characteristic = 0.
      end if

      icalculated = it
!ajwm debug onyl line:-
!print *,'nu_effective: max min avg = ',maxval(nu_effective),minval(nu_effective),sum(nu_effective*1.D0)/(mx*my*mz*1.D0)
    endsubroutine calc_viscosity
!***********************************************************************
!ajwm Utility routines - poss need moving elsewhere
    subroutine shock_max5(f,maxf)
!
!  return array maxed with by 2 points either way
!  skipping 1 data point all round
!
!  23-nov-02/tony: coded - from sub.f90 - nearmax
!
      use Cdata
!
      real, dimension (mx,my,mz) :: f
      real, dimension (mx,my,mz) :: maxf
      
!
!  check for degeneracy
!ajwm not correct condition
      if (mx.gt.1) then
         maxf(1     ,:,:) = amax1(f(1     ,:,:),  &
                                  f(2     ,:,:),  &
                                  f(3     ,:,:))
         maxf(2     ,:,:) = amax1(f(1     ,:,:), &
                                  f(2     ,:,:), &
                                  f(3     ,:,:), &
                                  f(4     ,:,:))
         maxf(3:mx-2,:,:) = amax1(f(1:mx-4,:,:), &
                                  f(2:mx-3,:,:), &
                                  f(3:mx-2,:,:), &
                                  f(4:mx-1,:,:), &
                                  f(5:mx  ,:,:))
         maxf(mx-1,:,:)   = amax1(f(mx-3  ,:,:), &
                                  f(mx-2  ,:,:), &
                                  f(mx-1  ,:,:), &
                                  f(mx    ,:,:))
         maxf(mx  ,:,:)   = amax1(f(mx-2  ,:,:), &
                                  f(mx-1  ,:,:), &
                                  f(mx    ,:,:))
      else
        maxf=f
      endif
!
!  check for degeneracy
!ajwm not correct condition
      if (my.gt.1) then
         f(:,1     ,:) = amax1(maxf(:,1     ,:),  &
                               maxf(:,2     ,:),  &
                               maxf(:,3     ,:))
         f(:,2     ,:) = amax1(maxf(:,1     ,:), &
                               maxf(:,2     ,:), &
                               maxf(:,3     ,:), &
                               maxf(:,4     ,:))
         f(:,3:my-2,:) = amax1(maxf(:,1:my-4,:), &
                               maxf(:,2:my-3,:), &
                               maxf(:,3:my-2,:), &
                               maxf(:,4:my-1,:), &
                               maxf(:,5:my  ,:))
         f(:,my-1,:)   = amax1(maxf(:,my-3  ,:), &
                               maxf(:,my-2  ,:), &
                               maxf(:,my-1  ,:), &
                               maxf(:,my    ,:))
         f(:,my  ,:)   = amax1(maxf(:,my-2  ,:), &
                               maxf(:,my-1  ,:), &
                               maxf(:,my    ,:))
      else
        f=maxf
      endif
!
!  check for degeneracy
!ajwm not correct condition
      if (mz.gt.1) then
      maxf(:,:,1     ) = amax1(f(:,:,1     ),  &
                               f(:,:,2     ),  &
                               f(:,:,3     ))
      maxf(:,:,2     ) = amax1(f(:,:,1     ), &
                               f(:,:,2     ), &
                               f(:,:,3     ), &
                               f(:,:,4     ))
      maxf(:,:,3:mz-2) = amax1(f(:,:,1:mz-4), &
                               f(:,:,2:mz-3), &
                               f(:,:,3:mz-2), &
                               f(:,:,4:mz-1), &
                               f(:,:,5:mz  ))
      maxf(:,:,mz-1  ) = amax1(f(:,:,mz-3  ), &
                               f(:,:,mz-2  ), &
                               f(:,:,mz-1  ), &
                               f(:,:,mz    ))
      maxf(:,:,mz    ) = amax1(f(:,:,mz-2  ), &
                               f(:,:,mz-1  ), &
                               f(:,:,mz    ))
      else
        maxf=f
      endif
!
    endsubroutine shock_max5

    subroutine shock_smooth(f,smoothf)
!
!  return array smoothed with by 2 points either way
!  skipping 3 data point all round 
!  i.e. result valid ()
!
!  23-nov-02/tony: coded
!
      use Cdata
!
      real, dimension (mx,my,mz) :: f
      real, dimension (mx,my,mz) :: smoothf
      
!
!  check for degeneracy
!ajwm not correct condition
      if (mx.gt.1) then
         smoothf(1     ,:,:) = 0.25 * (3.*f(1     ,:,:) +  &
                                          f(2     ,:,:))

         smoothf(2:mx-1,:,:) = 0.25 * (   f(1:mx-2,:,:) + &
                                       2.*f(2:mx-1,:,:) + &
                                          f(3:mx  ,:,:))
                                
         smoothf(mx    ,:,:) = 0.25 * (   f(mx-1  ,:,:) + &
                                       3.*f(mx    ,:,:))
      else
        smoothf=f
      endif
!
!  check for degeneracy
!ajwm not correct condition
      if (my.gt.1) then
         f(:,1     ,:) = 0.25 * (3.*smoothf(:,1     ,:) +  &
                                    smoothf(:,2     ,:))

         f(:,2:my-1,:) = 0.25 * (   smoothf(:,1:my-2,:) + &
                                 2.*smoothf(:,2:my-1,:) + &
                                    smoothf(:,3:my  ,:))
                                
         f(:,my    ,:) = 0.25 * (   smoothf(:,my-1  ,:) + &
                                 3.*smoothf(:,my    ,:))
      else
        f=smoothf
      endif
!
!  check for degeneracy
!ajwm not correct condition
      if (mz.gt.1) then
         smoothf(:,:,1     ) = 0.25 * (3.*f(:,:,1     ) +  &
                                          f(:,:,2     ))

         smoothf(:,:,2:mz-1) = 0.25 * (   f(:,:,1:mz-2) + &
                                       2.*f(:,:,2:mz-1) + &
                                          f(:,:,3:mz  ))
                                
         smoothf(:,:,mz    ) = 0.25 * (   f(:,:,mz-1  ) + &
                                       3.*f(:,:,mz    ))
      else
        smoothf=f
      endif

    endsubroutine shock_smooth
!***********************************************************************
    subroutine shock_divu(f,df)
!
!  calculate divergence of a vector U, get scalar
!  accurate to 2nd order, explicit,  centred an left and right biased 
!
!  23-nov-02/tony: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: df
      real :: fac
!
!ajwm If using mx,my,mz do we need degenerate n(xyz)grid=1 cases??
!ajwm Much slower using partial array?
!      fac=1./(2.*dx)
  df(1     ,:,:) =   (  4.*f(2,:,:,iux) &
                      - 3.*f(1,:,:,iux) &
                      -    f(3,:,:,iux) &
                      )/(2.*dx) 
  df(2:mx-1,:,:) =   (f(3:mx,:,:,iux)-f(1:mx-2,:,:,iux))/(2.*dx) 
  df(mx    ,:,:) =   (  3.*f(mx  ,:,:,iux) &
                      - 4.*f(mx-1,:,:,iux) &
                      +    f(mx-2,:,:,iux) &
                      )/(2.*dx)
 

  df(:,1     ,:) = df(:,1     ,:) &
                   + (  4.*f(:,2,:,iuy) &
                      - 3.*f(:,1,:,iuy) &
                      -    f(:,3,:,iuy))/(2.*dy) 
  df(:,2:my-1,:) = df(:,2:my-1,:) &
                   + (f(:,3:my,:,iuy)-f(:,1:my-2,:,iuy))/(2.*dy) 
  df(:,my    ,:) = df(:,my    ,:) & 
                   + (  3.*f(:,my  ,:,iuy) &
                      - 4.*f(:,my-1,:,iuy) &
                      +    f(:,my-2,:,iuy))/(2.*dy)

  df(:,:,1     ) = df(:,:,1     ) &
                   + (  4.*f(:,:,2,iuz) &
                      - 3.*f(:,:,1,iuz) &
                      -    f(:,:,3,iuz))/(2.*dz) 
  df(:,:,2:mz-1) = df(:,:,2:mz-1) & 
                   + (f(:,:,3:mz,iuz)-f(:,:,1:mz-2,iuz))/(2.*dz) 
  df(:,:,mz    ) = df(:,:,mz    ) &   
                   + (  3.*f(:,:,mz  ,iuz) &
                      - 4.*f(:,:,mz-1,iuz) &
                      +    f(:,:,mz-2,iuz))/(2.*dz)

    endsubroutine shock_divu

!***********************************************************************
    subroutine calc_viscous_heat(f,df,rho1,TT1)
!
!  calculate viscous heating term for right hand side of entropy equation
!
!  20-nov-02/tony: coded
!
      use Cdata
      use Mpicomm
      use Sub

      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx) :: rho1,TT1
      real, dimension (nx) :: sij2
      real, dimension (nx,3) :: gshock_characteristic, sgshock_characteristic

      if ( icalculated<it ) call calc_viscosity(f)
 !          call grad(spread(shock_characteristic,4,1),1,gshock_characteristic)
 !           call multmv_mn(sij,spread((nu_shock*shock_characteristic(l1:l2,m,n) + nu),2,3) * glnrho,sglnrho)
 !           call multmv_mn(sij,gshock_characteristic,sgshock_characteristic)

!            fvisc=2*sglnrho+nu*(del2u+1./3.*graddivu)
!            fvisc=fvisc + 2*nu_shock*sgshock_characteristic
!
!  traceless strainmatrix squared
!
      call multm2_mn(sij,sij2)
!
!      select case(ivisc)
!       case ('simplified', '0')
!         if (headtt) print*,'no heating: ivisc=',ivisc
!       case('rho_nu-const', '1')
!         if (headtt) print*,'viscous heating: ivisc=',ivisc
!         df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) + TT1*2.*nu*sij2*rho1
!       case('nu-const', '2')
         if (headtt) print*,'viscous heating: ivisc=',ivisc
         df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) + TT1*2.*nu*sij2 + TT1*nu_shock*(shock_characteristic(l1:l2,m,n)**2)
!       case default
!         if (lroot) print*,'ivisc=',trim(ivisc),' -- this could never happen'
!         call stop_it("")
!      endselect
    endsubroutine calc_viscous_heat

!***********************************************************************
    subroutine calc_viscous_force(f,df,glnrho,rho1)
!
!  calculate viscous heating term for right hand side of entropy equation
!
!  20-nov-02/tony: coded
!
      use Cdata
      use Mpicomm
      use Sub

      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3) :: glnrho, del2u, graddivu, fvisc, sglnrho 
      real, dimension (nx,3) :: gshock_characteristic, sgshock_characteristic
      real, dimension (nx) :: murho1,rho1
      integer :: i

      intent (in) :: f, glnrho, rho1
      intent (out) :: df

      if ( icalculated<it ) call calc_viscosity(f)
!
!  viscosity operator
!  rho1 is pre-calculated in equ
!
      if (nu_shock /= 0.) then

         !
         !  viscous force: nu*(del2u+graddivu/3+2S.glnrho)
         !  -- the correct expression for nu=const
         !
         if (headtt) print*,'viscous force: (nu_total) *(del2u+graddivu/3+2S.glnrho) + 2S.gnu_total)'
         call del2v_etc(f,iuu,del2u,GRADDIV=graddivu)
         if(ldensity) then
            call grad(spread(shock_characteristic,4,1),1,gshock_characteristic)
            call multmv_mn(sij,spread((nu_shock*shock_characteristic(l1:l2,m,n) + nu),2,3) * glnrho,sglnrho)
            call multmv_mn(sij,gshock_characteristic,sgshock_characteristic)

            fvisc=2*sglnrho+nu*(del2u+1./3.*graddivu)
            fvisc=fvisc + 2*nu_shock*sgshock_characteristic
!ajwm what's this???
            maxdiffus=amax1(maxdiffus,nu)
         else
            if(lfirstpoint) &
                 print*,"ldensity better be .true. for shock viscosity"
         endif
         
         df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)+fvisc
      else ! (nu_shock=0)
         if (nu /= 0.) then
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
            
            df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)+fvisc
         else ! (nu=0)
            if (headtt) print*,'no viscous force: (nu=0)'
         endif
      endif

    end subroutine calc_viscous_force
endmodule Viscosity
