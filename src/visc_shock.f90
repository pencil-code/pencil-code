! $Id: visc_shock.f90,v 1.53 2003-12-10 14:47:20 nilshau Exp $

!  This modules implements viscous heating and diffusion terms
!  here for shock viscosity nu_total = nu + nu_shock*dx*smooth(max5(-(div u)))) 

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 2
!
!***************************************************************

module Viscosity

  use Cparam
  use Cdata
  use Density

  implicit none

  real :: nu_shock = 0.
  character (len=labellen) :: ivisc=''
  real :: nutotal_max
  integer :: itest
  logical :: lvisc_first=.false.,lshock_max5=.true.

  ! input parameters
  integer :: dummy
  namelist /viscosity_init_pars/ dummy

  ! run parameters
  namelist /viscosity_run_pars/ nu,nu_shock,lvisc_first, &
       lshock_max5
 
  ! other variables (needs to be consistent with reset list below)
  integer :: i_dtnu=0,i_shockmax=0

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
      lvisc_shock = .true.
!
      ishock = mvar + naux + 1
      itest = mvar + naux + 2
      !naux = naux + 1
      naux = naux + 2
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_viscosity: shock viscosity nvar = ', nvar
        print*, 'ishock = ', ishock
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: visc_shock.f90,v 1.53 2003-12-10 14:47:20 nilshau Exp $")
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
      if (naux < maux) aux_var(aux_count)=',shock $'
      if (naux == maux) aux_var(aux_count)=',shock'
      aux_count=aux_count+1
      if (lroot) write(15,*) 'shock = fltarr(mx,my,mz)*one'
!
    endsubroutine register_viscosity
!***********************************************************************
    subroutine initialize_viscosity()
!
!  20-nov-02/tony: coded
!
       use CData
 !     if (nu /= 0. .and. (ivisc=='nu-const')) then
         lneed_sij=.true.
         lneed_glnrho=.true.
 !     endif

        if (headtt.and.lroot) print*,'viscosity: nu=',nu,', nu_shock=',nu_shock

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
        i_shockmax=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if(lroot.and.ip<14) print*,'rprint_viscosity: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'dtnu',i_dtnu)
        call parse_name(iname,cname(iname),cform(iname),'shockmax',i_shockmax)
      enddo
!
!  write column where which ionization variable is stored
!
      if (present(lwrite)) then
        if (lwrite) then
          write(3,*) 'i_dtnu=',i_dtnu
          write(3,*) 'i_shockmax=',i_dtnu
          write(3,*) 'ihyper=',ihyper
          write(3,*) 'ishock=',ishock
          write(3,*) 'itest=',itest
        endif
      endif
!   
      if(ip==0) print*,lreset  !(to keep compiler quiet)
    endsubroutine rprint_viscosity
!!***********************************************************************
    subroutine calc_viscosity(f)
!
!  calculate viscous heating term for right hand side of entropy equation
!
!  23-nov-02/tony: coded
!
      use IO
!     use Sub, only: grad
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: tmp
!      real, dimension (mx,my,mz,3) :: savef
!      real, dimension (nx,3) :: tmpuu
!      real :: radius 
!
      if(headt) print*,'calc_viscosity: dxmin=',dxmin
!
!  calculate shock viscosity only when nu_shock /=0
!
      if (nu_shock /= 0.) then
!
!  calculate (-divu)_+
!
!savef=f(:,:,:,iux:iuz)
!        radius=dxmin*10.
!        f(:,:,:,ishock)= spread(spread(exp((-x/radius)**2),2,my),3,mz) &
!                        *spread(spread(exp((-y/radius)**2),1,mx),3,mz) &
!                        *spread(spread(exp((-z/radius)**2),1,mx),2,my)
        
!         do n=n1,n2
!         do m=m1,m2
!           call grad(f,ishock,tmpuu)
!           f(l1:l2,m,n,iux:iuz)=tmpuu    
!         enddo
!         enddo

         call shock_divu(f,f(:,:,:,ishock))
         f(:,:,:,ishock)=amax1(0., -f(:,:,:,ishock))
!         f(:,:,:,ishock)=0.
!         f(mx/2,my/2,mz/2,ishock)=10.
!
!  take the max over 5 neighboring points and smooth.
!  Note that this means that we'd need 4 ghost zones, so
!  we use one-sided formulae on processor boundaries.
!  Alternatively, to get the same result with and without MPI
!  you may want to try lshock_max5=.false. (default is .true.)  
!
         if (lshock_max5) then
           call shock_max5(f(:,:,:,ishock),tmp)
         else
           call shock_max3(f(:,:,:,ishock),tmp)
         endif

! Save test data and scale to match the maximum expected result of smoothing 
         f(:,:,:,itest)=tmp  

         !f(:,:,:,itest)=f(:,:,:,ishock)
         call shock_smooth(tmp,f(:,:,:,ishock))
         !f(:,:,:,ishock)=tmp
!
!  scale with dxmin**2
!
         f(:,:,:,ishock) = f(:,:,:,ishock) * dxmin**2 
      else
         f(:,:,:,ishock) = 0.
      end if
!f(:,:,:,iux:iuz)=savef
!
!ajwm debug only line:-
! if (ip=0) call output(trim(directory_snap)//'/shockvisc.dat',f(:,:,:,ishock),1)
    endsubroutine calc_viscosity
!***********************************************************************
!ajwm Utility routines - poss need moving elsewhere
    subroutine shock_max5(f,maxf)
!
!  return array maxed with by 2 points either way
!  skipping 1 data point all round
!
!  23-nov-02/tony: coded - from sub.f90 - nearmax
!  26-may-03/axel: maxf and f where interchanged in y-chunk
!
      use Cdata
!
      real, dimension (mx,my,mz) :: f
      real, dimension (mx,my,mz) :: maxf, tmp
!
!  x-direction, f -> maxf
!  check for degeneracy
!
      if (nxgrid/=1) then
         if (mx.ge.5) then
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
         elseif (mx.eq.4) then
            maxf(1,:,:)=amax1(f(1,:,:),f(2,:,:),f(3,:,:))
            maxf(2,:,:)=amax1(f(1,:,:),f(2,:,:),f(3,:,:),f(4,:,:))
            maxf(3,:,:)=maxf(2,:,:)
            maxf(4,:,:)=amax1(f(2,:,:),f(3,:,:),f(4,:,:))
         elseif (mx.eq.3) then
            maxf(1,:,:)=amax1(f(1,:,:),f(2,:,:),f(3,:,:))
            maxf(2,:,:)=maxf(1,:,:)
            maxf(3,:,:)=maxf(1,:,:)
         elseif (mx.eq.2) then
            maxf(1,:,:)=amax1(f(1,:,:),f(2,:,:))
            maxf(2,:,:)=maxf(1,:,:)
         else
            maxf=f
         endif
      else
         maxf=f
      endif
!
!  y-direction, maxf -> f (swap back again)
!  check for degeneracy
!
      if (nygrid/=1) then
         if (my.ge.5) then
            tmp(:,1     ,:) = amax1(maxf(:,1     ,:),  &
                                  maxf(:,2     ,:),  &
                                  maxf(:,3     ,:))
            tmp(:,2     ,:) = amax1(maxf(:,1     ,:), &
                                  maxf(:,2     ,:), &
                                  maxf(:,3     ,:), &
                                  maxf(:,4     ,:))
            tmp(:,3:my-2,:) = amax1(maxf(:,1:my-4,:), &
                                  maxf(:,2:my-3,:), &
                                  maxf(:,3:my-2,:), &
                                  maxf(:,4:my-1,:), &
                                  maxf(:,5:my  ,:))
            tmp(:,my-1,:)   = amax1(maxf(:,my-3  ,:), &
                                  maxf(:,my-2  ,:), &
                                  maxf(:,my-1  ,:), &
                                  maxf(:,my    ,:))
            tmp(:,my  ,:)   = amax1(maxf(:,my-2  ,:), &
                                  maxf(:,my-1  ,:), &
                                  maxf(:,my    ,:))
         elseif (my.eq.4) then
            tmp(:,1,:)=amax1(maxf(:,1,:),maxf(:,2,:),maxf(:,3,:))
            tmp(:,2,:)=amax1(maxf(:,1,:),maxf(:,2,:),maxf(:,3,:),maxf(:,4,:))
            tmp(:,3,:)=tmp(:,2,:)
            tmp(:,4,:)=amax1(maxf(:,2,:),maxf(:,3,:),maxf(:,4,:))
         elseif (my.eq.3) then
            tmp(:,1,:)=amax1(maxf(:,1,:),maxf(:,2,:),maxf(:,3,:))
            tmp(:,2,:)=f(:,1,:)
            tmp(:,3,:)=f(:,1,:)
         elseif (my.eq.2) then
            tmp(:,1,:)=amax1(maxf(:,1,:),maxf(:,2,:))
            tmp(:,2,:)=tmp(:,1,:)
         else
            tmp=maxf
         endif
      else
         tmp=maxf
      endif
!
!  z-direction, f -> maxf
!  check for degeneracy
!
      if (nzgrid/=1) then
         if (mz.ge.5) then
            maxf(:,:,1     ) = amax1(tmp(:,:,1     ),  &
                                     tmp(:,:,2     ),  &
                                     tmp(:,:,3     ))
            maxf(:,:,2     ) = amax1(tmp(:,:,1     ), &
                                     tmp(:,:,2     ), &
                                     tmp(:,:,3     ), &
                                     tmp(:,:,4     ))
            maxf(:,:,3:mz-2) = amax1(tmp(:,:,1:mz-4), &
                                     tmp(:,:,2:mz-3), &
                                     tmp(:,:,3:mz-2), &
                                     tmp(:,:,4:mz-1), &
                                     tmp(:,:,5:mz  ))
            maxf(:,:,mz-1  ) = amax1(tmp(:,:,mz-3  ), &
                                     tmp(:,:,mz-2  ), &
                                     tmp(:,:,mz-1  ), &
                                     tmp(:,:,mz    ))
            maxf(:,:,mz    ) = amax1(tmp(:,:,mz-2  ), &
                                     tmp(:,:,mz-1  ), &
                                     tmp(:,:,mz    ))
         elseif (mz.eq.4) then
            maxf(:,:,1)=amax1(tmp(:,:,1),tmp(:,:,2),tmp(:,:,3))
            maxf(:,:,2)=amax1(tmp(:,:,1),tmp(:,:,2),tmp(:,:,3),tmp(:,:,4))
            maxf(:,:,3)=maxf(:,:,2)
            maxf(:,:,4)=amax1(tmp(:,:,2),tmp(:,:,3),tmp(:,:,4))
         elseif (mz.eq.3) then
            maxf(:,:,1)=amax1(tmp(:,:,1),tmp(:,:,2),tmp(:,:,3))
            maxf(:,:,2)=maxf(:,:,1)
            maxf(:,:,3)=maxf(:,:,1)
         elseif (mz.eq.2) then
            maxf(:,:,1)=amax1(tmp(:,:,1),tmp(:,:,2))
            maxf(:,:,2)=maxf(:,:,1)
         else
            maxf=tmp
         endif
      else
         maxf=tmp
      endif
!
    endsubroutine shock_max5

!***********************************************************************
!ajwm Utility routines - poss need moving elsewhere
    subroutine shock_max3(f,maxf)
!
!  return array maxed with by 2 points either way
!  skipping 1 data point all round
!
!  29-nov-03/axel: adapted from shock_max5
!
      use Cdata
!
      real, dimension (mx,my,mz) :: f
      real, dimension (mx,my,mz) :: maxf, tmp
!
!  x-direction, f -> maxf
!  check for degeneracy
!
      if (nxgrid/=1) then
         if (mx.ge.3) then
            maxf(1     ,:,:) = amax1(f(1     ,:,:), &
                                     f(2     ,:,:))
            maxf(2:mx-1,:,:) = amax1(f(1:mx-2,:,:), &
                                     f(2:mx-1,:,:), &
                                     f(3:mx  ,:,:))
            maxf(  mx  ,:,:) = amax1(f(  mx-1,:,:), &
                                     f(  mx  ,:,:))
         else
            maxf=f
         endif
      else
         maxf=f
      endif
!
!  y-direction, maxf -> f (swap back again)
!  check for degeneracy
!
      if (nygrid/=1) then
         if (my.ge.3) then
            tmp(:,1     ,:) = amax1(maxf(:,1     ,:),  &
                                    maxf(:,2     ,:))
            tmp(:,2:my-1,:) = amax1(maxf(:,1:my-2,:), &
                                    maxf(:,2:my-1,:), &
                                    maxf(:,3:my  ,:))
            tmp(:,  my  ,:) = amax1(maxf(:,  my-1,:), &
                                    maxf(:,  my  ,:))
         else
            tmp=maxf
         endif
      else
         tmp=maxf
      endif
!
!  z-direction, f -> maxf
!  check for degeneracy
!
      if (nzgrid/=1) then
         if (mz.ge.3) then
            maxf(:,:,1     ) = amax1(tmp(:,:,1     ),  &
                                     tmp(:,:,2     ))
            maxf(:,:,2:mz-1) = amax1(tmp(:,:,1:mz-2), &
                                     tmp(:,:,2:mz-1), &
                                     tmp(:,:,3:mz  ))
            maxf(:,:,mz    ) = amax1(tmp(:,:,  mz-1), &
                                     tmp(:,:,  mz  ))
         endif
      else
         maxf=tmp
      endif
!
    endsubroutine shock_max3

!***********************************************************************
    subroutine shock_smooth(f,smoothf)
!
!  return array smoothed with by 2 points either way
!  skipping 3 data point all round 
!  i.e. result valid ()
!
!  23-nov-02/tony: coded
!
      use Cdata, only: mx,my,mz
!
      real, dimension (mx,my,mz) :: f
      real, dimension (mx,my,mz) :: smoothf, tmp
!
!  check for degeneracy
!
      if (nxgrid/=1) then
         if (mx.ge.3) then
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
      else
         smoothf=f
      endif
!
!  check for degeneracy
!
      if (nygrid/=1) then
         if (my.ge.3) then
            tmp(:,1     ,:) = 0.25 * (3.*smoothf(:,1     ,:) +  &
                                         smoothf(:,2     ,:))

            tmp(:,2:my-1,:) = 0.25 * (   smoothf(:,1:my-2,:) + &
                                      2.*smoothf(:,2:my-1,:) + &
                                         smoothf(:,3:my  ,:))
                                
            tmp(:,my    ,:) = 0.25 * (   smoothf(:,my-1  ,:) + &
                                      3.*smoothf(:,my    ,:))
         else
            tmp=smoothf
         endif
      else
         tmp=smoothf
      endif
!
!  check for degeneracy
!
      if (nzgrid/=1) then
         if (mz.ge.3) then
            smoothf(:,:,1     ) = 0.25 * (3.*tmp(:,:,1     ) +  &
                                             tmp(:,:,2     ))

            smoothf(:,:,2:mz-1) = 0.25 * (   tmp(:,:,1:mz-2) + &
                                          2.*tmp(:,:,2:mz-1) + &
                                             tmp(:,:,3:mz  ))
                                
            smoothf(:,:,mz    ) = 0.25 * (   tmp(:,:,mz-1  ) + &
                                          3.*tmp(:,:,mz    ))
         else
            smoothf=tmp
         endif
      else
         smoothf=tmp
      end if
!
!smoothf=0.
!      do n=2,mz-1
!      do m=2,my-1
!      do l=2,mx-1
!        smoothf(l,m,n) =  8. * (  f(l,m,n) ) +                                                    & 
!                          4. * (  f(l,m-1,n  )+f(l+1,m,n  )+f(l-1,m,n  )+f(l,m+1,n  )             &
!                                + f(l,m,n+1)+f(l,m,n-1)) +                                        & 
!                          2. * (  f(l,m-1,n+1)+f(l+1,m,n+1)+f(l-1,m,n+1)+f(l,m+1,n+1)             &
!                                + f(l+1,m-1,n)+f(l+1,m+1,n)+f(l-1,m-1,n)+f(l-1,m+1,n)             &
!                                + f(l,m-1,n-1)+f(l+1,m,n-1)+f(l-1,m,n-1)+f(l,m+1,n-1) ) +         &
!                               (  f(l+1,m-1,n+1)+f(l+1,m+1,n+1)+f(l-1,m-1,n+1)+f(l-1,m+1,n+1)     &
!                                + f(l+1,m-1,n-1)+f(l+1,m+1,n-1)+f(l-1,m-1,n-1)+f(l-1,m+1,n-1) )  
!      enddo
!      enddo
!      enddo
!      smoothf = smoothf / 64.

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
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: df
!
!ajwm If using mx,my,mz do we need degenerate n(xyz)grid=1 cases??
!ajwm Much slower using partial array?
!      fac=1./(2.*dx)
      df=0.


      if (nxgrid/=1) then
         df(1     ,:,:) =  df(1    ,:,:) &
                           + (  4.*f(2,:,:,iux) &
                              - 3.*f(1,:,:,iux) &
                              -    f(3,:,:,iux))/(2.*dx) 
         df(2:mx-1,:,:) =  df(2:mx-1,:,:) &
                           + ( f(3:mx,:,:,iux)-f(1:mx-2,:,:,iux) ) / (2.*dx) 
         df(mx    ,:,:) =  df(mx    ,:,:) &
                           + (  3.*f(mx  ,:,:,iux) &
                              - 4.*f(mx-1,:,:,iux) &
                              +    f(mx-2,:,:,iux))/(2.*dx)
      endif

      if (nygrid/=1) then
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
      end if

      if (nzgrid/=1) then
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
      end if
    endsubroutine shock_divu

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
      call multm2_mn(sij,sij2)
!      if (headtt) print*,'viscous heating: ',ivisc

      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + TT1 * &
           (2.*nu*sij2  & 
             + nu_shock * shock * divu**2)
!
      if(ip==0) print*,glnrho,rho1,cs2,f !(to keep compiler quiet)
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
      real, dimension (nx,3) :: glnrho, del2u, graddivu, fvisc, sglnrho,tmp
      real, dimension (nx,3) :: gshock
      real, dimension (nx) :: rho1,divu,shock

      intent (in) :: f, glnrho, rho1
      intent (out) :: df,shock,gshock

!
!  viscosity operator
!  rho1 is pre-calculated in equ
!

      shock=f(l1:l2,m,n,ishock)
      call grad(f,ishock,gshock)
      nutotal_max=nu
!
!  shock viscosity
!
      if (nu_shock /= 0.) then
         !
         !  viscous force: nu*(del2u+graddivu/3+2S.glnrho)
         !  -- the correct expression for nu=const
         !
         call del2v_etc(f,iuu,del2u,GRADDIV=graddivu)
         if(ldensity) then
            call multmv_mn(sij,glnrho,sglnrho)
            call multsv_mn(divu,glnrho,fvisc)
            tmp=fvisc +graddivu
            call multsv_mn(nu_shock*shock,tmp,fvisc)
            call multsv_add_mn(fvisc,nu_shock*divu,gshock,tmp)
            fvisc = tmp + 2*nu*sglnrho+nu*(del2u+1./3.*graddivu)
            nutotal_max=nutotal_max+nu_shock*maxval(shock)
            maxdiffus=amax1(maxdiffus,nutotal_max)
         else
            if(lfirstpoint) &
                 print*,"ldensity better be .true. for ivisc='nu-const'"
         endif
            
         df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)+fvisc
      else ! (nu_shock=0)
         if (nu /= 0.) then
            !
            !  viscous force: nu*(del2u+graddivu/3+2S.glnrho)
            !  -- the correct expression for nu=const
            !
            if (headtt.and.lroot) print*,'viscous force: nu*(del2u+graddivu/3+2S.glnrho)'
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
            if (headtt.and.lroot) print*,'no viscous force: (nu=0)'
         endif
      endif
!
      if (ldiagnos) then
        if (i_dtnu/=0) call max_mn_name(spread(nutotal_max,1,nx)/dxmin**2/cdtvDim,i_dtnu,l_dt=.true.)
        if (i_shockmax/=0) call max_mn_name(shock,i_shockmax)
      endif
!
      if(ip==0) print*,rho1 !(to keep compiler quiet)
    end subroutine calc_viscous_force
endmodule Viscosity
