! $Id: visc_shock.f90,v 1.14 2003-01-16 20:27:45 mee Exp $

!  This modules implements viscous heating and diffusion terms
!  here for shock viscosity nu_total = nu + nu_shock * dx * smooth(max5(-(div u)))) 

module Viscosity

  use Cparam
  use Cdata
  use Density

  implicit none

!  real :: nu=0.
  real :: nu_shock = 0.
  character (len=labellen) :: ivisc=''
  integer :: icalculated = -1
  real :: maxeffectivenu

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
      ishock = mvar + naux + 1
      naux = naux + 1 
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_viscosity: constant viscosity nvar = ', nvar
        print*, 'ishock = ', ishock
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: visc_shock.f90,v 1.14 2003-01-16 20:27:45 mee Exp $")


! Check we arn't registering too many auxilliary variables
      if (naux > maux) then
        if (lroot) write(0,*) 'naux = ', naux, ', maux = ', maux
        call stop_it('Register_viscosity: naux > maux')
      endif

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
      use IO

      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: tmp

      if (nu_shock /= 0.) then
         call shock_divu(f,f(:,:,:,ishock))
         f(:,:,:,ishock)=amax1(0., -f(:,:,:,ishock))

 
         call shock_max5(f(:,:,:,ishock),tmp)
         call shock_smooth(tmp,f(:,:,:,ishock))   

         f(:,:,:,ishock) = f(:,:,:,ishock) * dxmin**2 
      else
         f(:,:,:,ishock) = 0.
      end if

      icalculated = it

      maxeffectivenu = amax1(maxval(f(:,:,:,ishock))*nu_shock,nu)
!ajwm debug onyl line:-
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
!
      use Cdata
!
      real, dimension (mx,my,mz) :: f
      real, dimension (mx,my,mz) :: maxf
      
!
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
!  check for degeneracy
!
      if (nygrid/=1) then
         if (my.ge.5) then
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
         elseif (my.eq.4) then
            maxf(:,1,:)=amax1(f(:,1,:),f(:,2,:),f(:,3,:))
            maxf(:,2,:)=amax1(f(:,1,:),f(:,2,:),f(:,3,:),f(:,4,:))
            maxf(:,3,:)=maxf(:,2,:)
            maxf(:,4,:)=amax1(f(:,2,:),f(:,3,:),f(:,4,:))
         elseif (my.eq.3) then
            maxf(:,1,:)=amax1(f(:,1,:),f(:,2,:),f(:,3,:))
            maxf(:,2,:)=maxf(:,1,:)
            maxf(:,3,:)=maxf(:,1,:)
         elseif (my.eq.2) then
            maxf(:,1,:)=amax1(f(:,1,:),f(:,2,:))
            maxf(:,2,:)=maxf(:,1,:)
         else
            maxf=f
         endif
      else
         maxf=f
      endif
!
!  check for degeneracy
!
      if (nzgrid/=1) then
         if (mz.ge.5) then
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
         elseif (mz.eq.4) then
            maxf(:,:,1)=amax1(f(:,:,1),f(:,:,2),f(:,:,3))
            maxf(:,:,2)=amax1(f(:,:,1),f(:,:,2),f(:,:,3),f(:,:,4))
            maxf(:,:,3)=maxf(:,:,2)
            maxf(:,:,4)=amax1(f(:,:,2),f(:,:,3),f(:,:,4))
         elseif (mz.eq.3) then
            maxf(:,:,1)=amax1(f(:,:,1),f(:,:,2),f(:,:,3))
            maxf(:,:,2)=maxf(:,:,1)
            maxf(:,:,3)=maxf(:,:,1)
         elseif (mz.eq.2) then
            maxf(:,:,1)=amax1(f(:,:,1),f(:,:,2))
            maxf(:,:,2)=maxf(:,:,1)
         else
            maxf=f
         endif
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
      integer :: j,k

      ! FORGET corners and edges for the moment! Even faces
      ! shouldn't be needed but we'll calculate them anyway... 

      ! Smooth faces regions - lower x
      smoothf(1     ,2:my-1,2:mz-1) = f(1     ,2:my-1,2:mz-1) * 2.0 / 3.0
      smoothf(1     ,2:my-1,2:mz-1) = smoothf(1     ,2:my-1,2:mz-1) + &
                                      ( f(2     ,2:my-1,2:mz-1) + &
                                        f(1     ,1:my-2,2:mz-1) + &
                                        f(1     ,3:my  ,2:mz-1) + &
                                        f(1     ,2:my-1,1:mz-2) + &
                                        f(1     ,2:my-1,3:mz  )) / 3.0
      smoothf(1     ,2:my-1,2:mz-1) = smoothf(1     ,2:my-1,2:mz-1) + &
                                      ( f(2     ,1:my-2,2:mz-1) + &
                                        f(2     ,3:my  ,2:mz-1) + &
                                        f(2     ,2:my-1,1:mz-2) + &
                                        f(2     ,2:my-1,3:mz ) + &
                                        f(1     ,3:my  ,1:mz-2) + &
                                        f(1     ,3:my  ,3:mz  ) + &
                                        f(1     ,1:my-2,1:mz-2) + &
                                        f(1     ,1:my-2,3:mz  )) / 6.0
      smoothf(1     ,2:my-1,2:mz-1) = smoothf(1 ,2:my-1,2:mz-1) + &
                                      ( f(2     ,1:my-2,1:mz-2) + &
                                        f(2     ,1:my-2,3:mz  ) + &
                                        f(2     ,3:my  ,1:mz-2) + &
                                        f(2     ,3:my  ,3:mz  )) / 12.0

      ! Smooth faces regions - upper x
      smoothf(mx    ,2:my-1,2:mz-1) = f(mx    ,2:my-1,2:mz-1) * 2.0 / 3.0
      smoothf(mx    ,2:my-1,2:mz-1) = smoothf(mx    ,2:my-1,2:mz-1) + &
                                      ( f(mx-1  ,2:my-1,2:mz-1) + &
                                        f(mx    ,1:my-2,2:mz-1) + &
                                        f(mx    ,3:my  ,2:mz-1) + &
                                        f(mx    ,2:my-1,1:mz-2) + &
                                        f(mx    ,2:my-1,3:mz  )) / 3.0
      smoothf(mx    ,2:my-1,2:mz-1) = smoothf(mx    ,2:my-1,2:mz-1) + &
                                      ( f(mx-1  ,1:my-2,2:mz-1) + &
                                        f(mx-1  ,3:my  ,2:mz-1) + &
                                        f(mx-1  ,2:my-1,1:mz-2) + &
                                        f(mx-1  ,2:my-1,3:mz ) + &
                                        f(mx    ,3:my  ,1:mz-2) + &
                                        f(mx    ,3:my  ,3:mz  ) + &
                                        f(mx    ,1:my-2,1:mz-2) + &
                                        f(mx    ,1:my-2,3:mz  )) / 6.0
      smoothf(mx    ,2:my-1,2:mz-1) = smoothf(mx,2:my-1,2:mz-1) + &
                                      ( f(mx-1  ,1:my-2,1:mz-2) + &
                                        f(mx-1  ,1:my-2,3:mz  ) + &
                                        f(mx-1  ,3:my  ,1:mz-2) + &
                                        f(mx-1  ,3:my  ,3:mz  )) / 12.0
      ! Smooth faces regions - lower y
      smoothf(2:mx-1,1     ,2:mz-1) = f(2:mx-1,1     ,2:mz-1) * 2.0 / 3.0
      smoothf(2:mx-1,1     ,2:mz-1) = smoothf(2:mx-1,1     ,2:mz-1) + &
                                      ( f(2:mx-1,2     ,2:mz-1) + &
                                        f(1:mx-2,1     ,2:mz-1) + &
                                        f(3:mx  ,1     ,2:mz-1) + &
                                        f(2:mx-1,1     ,1:mz-2) + &
                                        f(2:mx-1,1     ,3:mz  )) / 3.0
      smoothf(2:mx-1,1     ,2:mz-1) = smoothf(2:mx-1,1     ,2:mz-1) + &
                                      ( f(1:mx-2,2     ,2:mz-1) + &
                                        f(3:mx  ,2     ,2:mz-1) + &
                                        f(2:mx-1,2     ,1:mz-2) + &
                                        f(2:mx-1,2     ,3:mz ) + &
                                        f(3:mx  ,1     ,1:mz-2) + &
                                        f(3:mx  ,1     ,3:mz  ) + &
                                        f(1:mx-2,1     ,1:mz-2) + &
                                        f(1:mx-2,1     ,3:mz  )) / 6.0
      smoothf(2:mx-1,1     ,2:mz-1) = smoothf(2:mx-1,1 ,2:mz-1) + &
                                      ( f(1:mx-2,2     ,1:mz-2) + &
                                        f(1:mx-2,2     ,3:mz  ) + &
                                        f(3:mx  ,2     ,1:mz-2) + &
                                        f(3:mx  ,2     ,3:mz  )) / 12.0

      ! Smooth faces regions - upper y
      smoothf(2:mx-1,my    ,2:mz-1) = f(2:mx-1,my    ,2:mz-1) * 2.0 / 3.0
      smoothf(2:mx-1,my    ,2:mz-1) = smoothf(2:mx-1,my    ,2:mz-1) + &
                                      ( f(2:mx-1,my-1  ,2:mz-1) + &
                                        f(1:mx-2,my    ,2:mz-1) + &
                                        f(3:mx  ,my    ,2:mz-1) + &
                                        f(2:mx-1,my    ,1:mz-2) + &
                                        f(2:mx-1,my    ,3:mz  )) / 3.0
      smoothf(2:mx-1,my    ,2:mz-1) = smoothf(2:mx-1,my    ,2:mz-1) + &
                                      ( f(1:mx-2,my-1  ,2:mz-1) + &
                                        f(3:mx  ,my-1  ,2:mz-1) + &
                                        f(2:mx-1,my-1  ,1:mz-2) + &
                                        f(2:mx-1,my-1  ,3:mz ) + &
                                        f(3:mx  ,my    ,1:mz-2) + &
                                        f(3:mx  ,my    ,3:mz  ) + &
                                        f(1:mx-2,my    ,1:mz-2) + &
                                        f(1:mx-2,my    ,3:mz  )) / 6.0
      smoothf(2:mx-1,my    ,2:mz-1) = smoothf(2:mx-1,my,2:mz-1) + &
                                      ( f(1:mx-2,my-1  ,1:mz-2) + &
                                        f(1:mx-2,my-1  ,3:mz  ) + &
                                        f(3:mx  ,my-1  ,1:mz-2) + &
                                        f(3:mx  ,my-1  ,3:mz  )) / 12.0
      ! Smooth faces regions - lower z
      smoothf(2:mx-1,2:my-1,1     ) = f(2:mx-1,2:my-1,1     ) * 2.0 / 3.0
      smoothf(2:mx-1,2:my-1,1     ) = smoothf(2:mx-1,2:my-1,1     ) + &
                                      ( f(2:mx-1,2:my-1,2     ) + &
                                        f(1:mx-2,2:my-1,1     ) + &
                                        f(3:mx  ,2:my-1,1     ) + &
                                        f(2:mx-1,1:my-2,1     ) + &
                                        f(2:mx-1,3:my  ,1     )) / 3.0
      smoothf(2:mx-1,2:my-1,1     ) = smoothf(2:mx-1,2:my-1,1     ) + &
                                      ( f(1:mx-2,2:my-1,2     ) + &
                                        f(3:mx  ,2:my-1,2     ) + &
                                        f(2:mx-1,1:my-2,2     ) + &
                                        f(2:mx-1,3:my  ,2     ) + &
                                        f(3:mx  ,1:my-2,1     ) + &
                                        f(3:mx  ,3:my  ,1     ) + &
                                        f(1:mx-2,1:my-2,1     ) + &
                                        f(1:mx-2,3:my  ,1     )) / 6.0
      smoothf(2:mx-1,2:my-1,1     ) = smoothf(2:mx-1,2:my-1,1 ) + &
                                      ( f(1:mx-2,1:my-2,2     ) + &
                                        f(1:mx-2,3:my  ,2     ) + &
                                        f(3:mx  ,1:my-2,2     ) + &
                                        f(3:mx  ,3:my  ,2     )) / 12.0

      ! Smooth faces regions - upper z
      smoothf(2:mx-1,2:my-1,mz    ) = f(2:mx-1,2:my-1,mz    ) * 2.0 / 3.0
      smoothf(2:mx-1,2:my-1,mz    ) = smoothf(2:mx-1,2:my-1,mz    ) + &
                                      ( f(2:mx-1,2:my-1,mz-1  ) + &
                                        f(1:mx-2,2:my-1,mz    ) + &
                                        f(3:mx  ,2:my-1,mz    ) + &
                                        f(2:mx-1,1:my-2,mz    ) + &
                                        f(2:mx-1,3:my  ,mz    )) / 3.0
      smoothf(2:mx-1,2:my-1,mz    ) = smoothf(2:mx-1,2:my-1,mz    ) + &
                                      ( f(1:mx-2,2:my-1,mz-1  ) + &
                                        f(3:mx  ,2:my-1,mz-1  ) + &
                                        f(2:mx-1,1:my-2,mz-1  ) + &
                                        f(2:mx-1,3:my  ,mz-1  ) + &
                                        f(3:mx  ,1:my-2,mz    ) + &
                                        f(3:mx  ,3:my  ,mz    ) + &
                                        f(1:mx-2,1:my-2,mz    ) + &
                                        f(1:mx-2,3:my  ,mz    )) / 6.0
      smoothf(2:mx-1,2:my-1,mz    ) = smoothf(2:mx-1,2:my-1,mz) + &
                                      ( f(1:mx-2,1:my-2,mz-1  ) + &
                                        f(1:mx-2,3:my  ,mz-1  ) + &
                                        f(3:mx  ,1:my-2,mz-1  ) + &
                                        f(3:mx  ,3:my  ,mz-1  )) / 12.0

      ! Smooth central region - pencil by pencil
      smoothf(2:mx-1,:,:) = f(2:mx-1,:,:) / 8.0
      do j=2,my-1
         do k=2,mz-1
           smoothf(2:mx-1,j,k) = smoothf(2:mx-1,j,k) + &
                                         ( f(1:mx-2,j,k) + &
                                           f(3:mx,j,k) + &
                                           f(2:mx-1,j-1,k) + &
                                           f(2:mx-1,j+1,k) + &
                                           f(2:mx-1,j,k-1) + &
                                           f(2:mx-1,j,k+1)) / 16.0
  
           smoothf(2:mx-1,j,k) = smoothf(2:mx-1,j,k) + &
                                         ( f(1:mx-2,j-1,k) + &
                                           f(1:mx-2,j+1,k) + &
                                           f(1:mx-2,j  ,k+1) + &
                                           f(1:mx-2,j  ,k-1) + &
                                           f(3:mx  ,j-1,k) + &
                                           f(3:mx  ,j+1,k) + &
                                           f(3:mx  ,j  ,k-1) + &
                                           f(3:mx  ,j  ,k+1) + &
                                           f(2:mx-1,j-1,k-1) + &
                                           f(2:mx-1,j+1,k-1) + &
                                           f(2:mx-1,j-1,k+1) + &
                                           f(2:mx-1,j+1,k+1)) / 32.0
           smoothf(2:mx-1,j,k) = smoothf(2:mx-1,j,k) + &
                                         ( f(1:mx-2,j-1,k-1) + &
                                           f(1:mx-2,j-1,k+1) + &
                                           f(1:mx-2,j+1,k-1) + &
                                           f(1:mx-2,j+1,k+1) + &
                                           f(3:mx  ,j-1,k-1) + &
                                           f(3:mx  ,j-1,k+1) + &
                                           f(3:mx  ,j+1,k-1) + &
                                           f(3:mx  ,j+1,k+1)) / 64.0

         end do
      end do
!!$!
!!$!  check for degeneracy
!!$!
!!$      if (nxgrid/=1) then
!!$         if (mx.ge.3) then
!!$            smoothfx(1     ,:,:) = 0.25 * (3.*f(1     ,:,:) +  &
!!$                                             f(2     ,:,:))
!!$
!!$            smoothfx(2:mx-1,:,:) = 0.25 * (   f(1:mx-2,:,:) + &
!!$                                          2.*f(2:mx-1,:,:) + &
!!$                                             f(3:mx  ,:,:))
!!$                                
!!$            smoothfx(mx    ,:,:) = 0.25 * (   f(mx-1  ,:,:) + &
!!$                                          3.*f(mx    ,:,:))
!!$         else
!!$            smoothfx=f
!!$         endif
!!$      else
!!$         smoothfx=f
!!$      endif
!!$!
!!$!  check for degeneracy
!!$!
!!$      if (nygrid/=1) then
!!$         if (my.ge.3) then
!!$            smoothfy(:,1     ,:) = 0.25 * (3.*f(:,1     ,:) +  &
!!$                                       f(:,2     ,:))
!!$
!!$            smoothfy(:,2:my-1,:) = 0.25 * (   f(:,1:my-2,:) + &
!!$                                    2.*f(:,2:my-1,:) + &
!!$                                       f(:,3:my  ,:))
!!$                                
!!$            smoothfy(:,my    ,:) = 0.25 * (   f(:,my-1  ,:) + &
!!$                                    3.*f(:,my    ,:))
!!$         else
!!$            f=smoothf
!!$         endif
!!$      else
!!$         f=smoothf
!!$      endif
!!$!
!!$!  check for degeneracy
!!$!
!!$      if (nzgrid/=1) then
!!$         if (mz.ge.3) then
!!$            smoothf(:,:,1     ) = 0.25 * (3.*f(:,:,1     ) +  &
!!$                                             f(:,:,2     ))
!!$
!!$            smoothf(:,:,2:mz-1) = 0.25 * (   f(:,:,1:mz-2) + &
!!$                                          2.*f(:,:,2:mz-1) + &
!!$                                             f(:,:,3:mz  ))
!!$                                
!!$            smoothf(:,:,mz    ) = 0.25 * (   f(:,:,mz-1  ) + &
!!$                                          3.*f(:,:,mz    ))
!!$         else
!!$            smoothf=f
!!$         endif
!!$      else
!!$         smoothf=f
!!$      end if

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
!
!ajwm If using mx,my,mz do we need degenerate n(xyz)grid=1 cases??
!ajwm Much slower using partial array?
!      fac=1./(2.*dx)
      df=0.

      if (nxgrid/=1) then
         df(1     ,:,:) =   df(:,1     ,:) &
                            + (  4.*f(2,:,:,iux) &
                               - 3.*f(1,:,:,iux) &
                               -    f(3,:,:,iux) &
                               )/(2.*dx) 
         df(2:mx-1,:,:) =   (f(3:mx,:,:,iux)-f(1:mx-2,:,:,iux))/(2.*dx) 
         df(mx    ,:,:) =   (  3.*f(mx  ,:,:,iux) &
                             - 4.*f(mx-1,:,:,iux) &
                             +    f(mx-2,:,:,iux) &
                             )/(2.*dx)
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
      real, dimension (nx) :: rho1,TT1,cs2
      real, dimension (nx) :: sij2, divu
      real, dimension (nx,3) :: glnrho

      if ( icalculated<it ) call calc_viscosity(f)
!
!  traceless strainmatrix squared
!
      call multm2_mn(sij,sij2)
!      if (headtt) print*,'viscous heating: ',ivisc

      df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) + TT1 * &
           (2.*nu*sij2  & 
           + nu_shock * f(l1:l2,m,n,ishock) * divu**2)

      maxheating=amax1(maxheating,df(l1:l2,m,n,ient))
!
      if(ip==0) print*,glnrho,rho1,cs2 !(to keep compiler quiet)
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

      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3) :: glnrho, del2u, graddivu, fvisc, sglnrho,tmp
      real, dimension (nx,3) :: gshock_characteristic
      real, dimension (nx) :: rho1, divu

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
!         if (headtt) print*,'viscous force: (nu_total) *(del2u+graddivu/3+2S.glnrho) + 2S.gnu_total)'
         call del2v_etc(f,iuu,del2u,GRADDIV=graddivu)
         if(ldensity) then
            call grad(f(:,:,:,ishock),gshock_characteristic)
            call multmv_mn(sij,glnrho,sglnrho)
            call multsv_mn(glnrho,divu,fvisc)
            tmp=fvisc +graddivu
            call multsv_mn(tmp,nu_shock*f(l1:l2,m,n,ishock),fvisc)
            call multsv_add_mn(fvisc,nu_shock * divu, gshock_characteristic,tmp)
            fvisc = tmp + 2*nu*sglnrho+nu*(del2u+1./3.*graddivu)
            maxdiffus=amax1(maxdiffus,maxeffectivenu)
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
!
      if(ip==0) print*,rho1 !(to keep compiler quiet)
    end subroutine calc_viscous_force
endmodule Viscosity
