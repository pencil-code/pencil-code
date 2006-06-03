! $Id: shock.f90,v 1.12 2006-06-03 20:51:12 ajohan Exp $

!  This modules implements viscous heating and diffusion terms
!  here for shock viscosity
!    nu_total = nu + nu_shock*dx^2*smooth(max5(-(div u)))) 
!
!  NOTE: this works and has been tested for periodic boundaries.
!  With the current version, if your shock fronts reach a non-periodic
!  boundary, unexpected things may happen, so you should monitor the
!  behavior on the boundaries in this case.
!

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lshock = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 1
! COMMUNICATED AUXILIARIES 1
!
! PENCILS PROVIDED shock,gshock
!
!***************************************************************

module Shock

  use Cparam
  use Cdata
  use Density
  use Messages

  implicit none

  include 'shock.h'

  logical :: lshock_first=.true.,lshock_max5=.false.
  logical :: lwith_extreme_div=.false.
  logical :: lmax_smooth=.false.
  logical :: lgauss_integral=.false.
  logical :: lcommunicate_uu=.true.
  real :: div_threshold=0., div_scaling=1.
  real, dimension (3,3,3) :: smooth_factor

  ! input parameters
  !namelist /viscosity_init_pars/ dummy

  ! run parameters
  namelist /shock_run_pars/ lshock_first, lshock_max5, div_threshold, div_scaling, &
                            lmax_smooth, lgauss_integral, lcommunicate_uu
 
  ! other variables (needs to be consistent with reset list below)
  integer :: idiag_shockmax=0
!
  interface shock_max3
    module procedure shock_max3_farray
    module procedure shock_max3_pencil
  endinterface
!
  interface shock_smooth
    module procedure shock_smooth_farray
    module procedure shock_smooth_pencil
  endinterface
!
  interface shock_divu
    module procedure shock_divu_farray
    module procedure shock_divu_pencil
  endinterface
!
  contains

!***********************************************************************
    subroutine register_shock()
!
!  19-nov-02/tony: coded
!  24-jan-05/tony: modified from visc_shock.f90
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_shock called twice')
      first = .false.
!
      ishock = mvar + naux_com + 1
      naux = naux + 1
      naux_com = naux_com + 1
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_shock: shock viscosity nvar = ', nvar
        print*, 'ishock = ', ishock
      endif
!
!  Put variable name in array
!
      varname(ishock) = 'shock'
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: shock.f90,v 1.12 2006-06-03 20:51:12 ajohan Exp $")
!
! Check we aren't registering too many auxiliary variables
!
      if (naux > maux) then
        if (lroot) write(0,*) 'naux = ', naux, ', maux= ', maux
        call stop_it('register_shock: naux > maux')
      endif
      if (naux_com > maux_com) then
        if (lroot) write(0,*) 'naux_com = ', naux_com, ', maux_com = ', maux_com
        call stop_it('register_shock: naux_com > maux_com')
      endif
!
!  Writing files for use with IDL
!
      if (naux+naux_com <  maux+maux_com) aux_var(aux_count)=',shock $'
      if (naux+naux_com  == maux+maux_com) aux_var(aux_count)=',shock'
      aux_count=aux_count+1
      if (lroot) write(15,*) 'shock = fltarr(mx,my,mz)*one'
!
    endsubroutine register_shock
!***********************************************************************
    subroutine initialize_shock(lstarting)
!
!  20-nov-02/tony: coded
!
       use CData

       logical, intent(in) :: lstarting 
!
!  Calculate factors for polynomial smoothing
!
      smooth_factor=1.
      if (nxgrid/=1) then
        smooth_factor(1,:,:)=smooth_factor(1,:,:)*0.25
        smooth_factor(2,:,:)=smooth_factor(2,:,:) *0.5
        smooth_factor(3,:,:)=smooth_factor(3,:,:)*0.25
      else
        smooth_factor(1,:,:)=0.
        smooth_factor(3,:,:)=0.
      endif
!
      if (nygrid/=1) then
        smooth_factor(:,1,:)=smooth_factor(:,1,:)*0.25
        smooth_factor(:,2,:)=smooth_factor(:,2,:) *0.5
        smooth_factor(:,3,:)=smooth_factor(:,3,:)*0.25
      else
        smooth_factor(:,1,:)=0.
        smooth_factor(:,3,:)=0.
      endif
!
      if (nzgrid/=1) then
        smooth_factor(:,:,1)=smooth_factor(:,:,1)*0.25
        smooth_factor(:,:,2)=smooth_factor(:,:,2) *0.5
        smooth_factor(:,:,3)=smooth_factor(:,:,3)*0.25
      else
        smooth_factor(:,:,1)=0.
        smooth_factor(:,:,3)=0.
      endif
      if (lroot) print*,"initialize_shock: prenormalised shock_factor sum=",sum(smooth_factor)
      smooth_factor=smooth_factor/sum(smooth_factor)
!
      if (div_threshold/=0.) lwith_extreme_div=.true.

      if (.not.lgauss_integral) lmax_smooth=.true.
!        if (headtt.and.lroot.and.(.not.lstarting)) &
!                  print*,'viscosity: nu=',nu,', nu_shock=',nu_shock
      if (NO_WARN) print*,lstarting
!
    endsubroutine initialize_shock
!***********************************************************************
    subroutine read_shock_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat).and.NO_WARN) print*,iostat
      if (NO_WARN) print*,unit
!
    endsubroutine read_shock_init_pars
!***********************************************************************
    subroutine write_shock_init_pars(unit)
      integer, intent(in) :: unit
!
      if (NO_WARN) print*,unit
!
    endsubroutine write_shock_init_pars
!***********************************************************************
    subroutine read_shock_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=shock_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=shock_run_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_shock_run_pars
!***********************************************************************
    subroutine write_shock_run_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=shock_run_pars)

    endsubroutine write_shock_run_pars
!*******************************************************************
    subroutine rprint_shock(lreset,lwrite)
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
        idiag_shockmax=0
        lwith_extreme_div=.false.
        lmax_smooth=.false.
        lgauss_integral=.false.
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if (lroot.and.ip<14) print*,'rprint_shock: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),&
            'shockmax',idiag_shockmax)
      enddo
!
!  write column where which ionization variable is stored
!
      if (present(lwrite)) then
        if (lwrite) then
          write(3,*) 'i_shockmax=',idiag_shockmax
          write(3,*) 'ishock=',ishock
        endif
      endif
!   
      if (NO_WARN) print*,lreset  !(to keep compiler quiet)
    endsubroutine rprint_shock
!!***********************************************************************
    subroutine pencil_criteria_shock()
!    
!  All pencils that the Viscosity module depends on are specified here.
!
!  20-11-04/anders: coded
!
!   dummy
!
      if (idiag_shockmax/=0) then
          lpenc_diagnos(i_shock)=.true.
      endif
    endsubroutine pencil_criteria_shock
!***********************************************************************
    subroutine pencil_interdep_shock(lpencil_in)
!
!  Interdependency among pencils from the Viscosity module is specified here.
!
!  20-11-04/anders: coded
!
      use Cdata
!
      logical, dimension (npencils) :: lpencil_in
!      
      if (NO_WARN) print*, lpencil_in !(keep compiler quiet)
!
    endsubroutine pencil_interdep_shock
!***********************************************************************
    subroutine calc_pencils_shock(f,p)
!     
!  Calculate Viscosity pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  20-11-04/anders: coded
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
! shock
      if (lpencil(i_shock)) p%shock=f(l1:l2,m,n,ishock)
! gshock
      if (lpencil(i_gshock)) call grad(f,ishock,p%gshock)
!
!print*,maxval(f(l1:l2,m1:m2,n1:n2,ishock))
      if (ldiagnos.and.idiag_shockmax/=0) call max_mn_name(p%shock,idiag_shockmax)
!
    endsubroutine calc_pencils_shock
!***********************************************************************
    subroutine calc_shock_profile_simple(f)
!
!  Calculate divu based shock profile to be used in viscosity and 
!  diffusion type terms.
!
!  23-nov-02/tony: coded
!
!for debug      use IO
     use Cdata
!
     real, dimension (mx,my,mz,mvar+maux) :: f
     real, dimension (mx,my,mz) :: tmp
!
     if ((.not.lgauss_integral).and.(.not.lcommunicate_uu).and.((.not.lshock_first).or.lfirst)) then
!
!  calculate shock viscosity only when nu_shock /=0
!
!  calculate (-divu)_+
!
       call shock_divu(f,f(:,:,:,ishock))
!
       if (lwith_extreme_div) then
         where ((f(:,:,:,ishock) .gt. 0.) .and. (f(:,:,:,ishock) .lt. div_threshold)) &
             f(:,:,:,ishock)=0.
         where (f(:,:,:,ishock) .gt. 0.) &
             f(:,:,:,ishock)=f(:,:,:,ishock)*div_scaling
         f(:,:,:,ishock)=abs(f(:,:,:,ishock))
       else
         f(:,:,:,ishock)=max(0., -f(:,:,:,ishock))
       endif
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

       call shock_smooth(tmp,f(:,:,:,ishock))
!
!  scale with dxmin**2
!
       f(:,:,:,ishock) = f(:,:,:,ishock) * dxmin**2 
     endif
!
!debug only line:-
! if (ip=0) call output(trim(directory_snap)//'/shockvisc.dat',f(:,:,:,ishock),1)
    endsubroutine calc_shock_profile_simple
!!***********************************************************************
    subroutine calc_shock_profile(f)
!
!  Calculate divu based shock profile to be used in viscosity and 
!  diffusion type terms.
!
!  12-apr-05/tony: coded
!
      use Cdata
      use Boundcond
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: tmp 
      real, dimension(mx) :: penc
      integer :: jj,kk
!
      if ((.not.lshock_first).or.lfirst) then
        if (lgauss_integral) then
!
! Calculate contributions to other processors
!
          call shock_calc_externalboundary(f)
!
! Initiate asyncronous communication of contributions to other processors 
!
          if (nxgrid/=1) call bcshock_per_x(f)
          call initiate_isendrcv_bdry(f,ishock,ishock)
!
! Calculate all local shock profile contributions
!
!        call shock_calc_internalboundary(f)
          call shock_calc_body(f)
!
! Finalize all shock profile communications
!
          call finalize_isendrcv_bdry(f,ishock,ishock)
          if (nygrid/=1) call bcshock_per_y(f)
          if (nzgrid/=1) call bcshock_per_z(f)
!FIX ME
!! CALCULATE shock in REAL boundary locations... and neighbours!
!!!
!        call scale_and_chop_internalboundary(f)
         !f(:,:,:,ishock) = tmp * dxmin**2 
        elseif (lcommunicate_uu) then
!  Communicate uu ghost zones
          call initiate_isendrcv_bdry(f,iux,iuz)
          f(:,:,:,ishock)=0.
!
!  Divu over internal region
!          
          do n=n1+1,n2-1; do m=m1+1,m2-1
            call shock_divu(f,iuu,penc) 
            f(:,m,n,ishock)=max(0.,-penc)
          enddo; enddo
!
!  Max3 over internal region
!          
          do n=n1+2,n2-2; do m=m1+2,m2-2
            call shock_max3(f,ishock,penc) 
            tmp(:,m,n)=penc
          enddo; enddo
!
!  Smooth over internal region
!
          do n=n1+3,n2-3; do m=m1+3,m2-3
            call shock_smooth(tmp,penc) 
            f(:,m,n,ishock)=penc
          enddo; enddo
!
!  Set periodic x-boundaries in smooth(max3(shock))
!          
          if (nxgrid/=1) then 
            call bc_per_x(f,'top',ishock)
            call bc_per_x(f,'bot',ishock)
          endif
!  End communication of uu ghost zones.
          call finalize_isendrcv_bdry(f,iux,iuz)
          call boundconds_y(f)
          call boundconds_z(f)
!
! Divu over external region
!
          do n=2,mz-1; do jj=1,3
            m=1+jj
            call shock_divu(f,iuu,penc) 
            f(:,m,n,ishock)=max(-penc,0.)
            m=my-jj
            call shock_divu(f,iuu,penc) 
            f(:,m,n,ishock)=max(-penc,0.)
          enddo; enddo
          do kk=1,3; do m=5,my-4
            n=1+kk
            call shock_divu(f,iuu,penc) 
            f(:,m,n,ishock)=max(-penc,0.)
            n=mz-kk
            call shock_divu(f,iuu,penc) 
            f(:,m,n,ishock)=max(-penc,0.)
          enddo; enddo
!
! Max over external region
!
          do n=3,mz-2; do jj=2,4
            m=1+jj
            call shock_max3(f,ishock,penc) 
            tmp(:,m,n)=penc
            m=my-jj
            call shock_max3(f,ishock,penc) 
            tmp(:,m,n)=penc
          enddo; enddo
          do kk=2,4; do m=6,my-5
            n=1+kk
            call shock_max3(f,ishock,penc) 
            tmp(:,m,n)=penc
            n=mz-kk
            call shock_max3(f,ishock,penc) 
            tmp(:,m,n)=penc
          enddo; enddo
!
! Smooth over external region
!
          do n=4,mz-3; do jj=3,5
            m=1+jj
            call shock_smooth(tmp,penc) 
            f(:,m,n,ishock)=penc
            m=my-jj
            call shock_smooth(tmp,penc) 
            f(:,m,n,ishock)=penc
          enddo; enddo
          do kk=3,5; do m=7,my-6
            n=1+kk
            call shock_smooth(tmp,penc) 
            f(:,m,n,ishock)=penc
            n=mz-kk
            call shock_smooth(tmp,penc) 
            f(:,m,n,ishock)=penc
          enddo; enddo
!
!  Non-MPI version:
!          
!          do n=2,mz-1; do m=2,my-1
!            call shock_divu(f,iuu,penc) 
!            f(:,m,n,ishock)=max(0.,-penc)
!          enddo; enddo
!          do n=3,mz-2; do m=3,my-2
!            call shock_max3(f,ishock,penc) 
!            tmp(:,m,n)=penc
!          enddo; enddo
!          do n=4,mz-3; do m=4,my-3
!            call shock_smooth(tmp,penc) 
!            f(:,m,n,ishock)=penc
!          enddo; enddo
!
!  Set global x-boundaries.
!
          call boundconds_x(f)
          f(:,:,:,ishock) = f(:,:,:,ishock) * dxmin**2 
        endif
      endif
!
!debug only line:-
! if (ip=0) call output(trim(directory_snap)//'/shockvisc.dat',f(:,:,:,ishock),1)
    endsubroutine calc_shock_profile
!***********************************************************************
!Utility routines - poss need moving elsewhere
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
            maxf(1     ,:,:) = max(f(1     ,:,:),  &
                                     f(2     ,:,:),  &
                                     f(3     ,:,:))
            maxf(2     ,:,:) = max(f(1     ,:,:), &
                                     f(2     ,:,:), &
                                     f(3     ,:,:), &
                                     f(4     ,:,:))
            maxf(3:mx-2,:,:) = max(f(1:mx-4,:,:), &
                                     f(2:mx-3,:,:), &
                                     f(3:mx-2,:,:), &
                                     f(4:mx-1,:,:), &
                                     f(5:mx  ,:,:))
            maxf(mx-1,:,:)   = max(f(mx-3  ,:,:), &
                                     f(mx-2  ,:,:), &
                                     f(mx-1  ,:,:), &
                                     f(mx    ,:,:))
            maxf(mx  ,:,:)   = max(f(mx-2  ,:,:), &
                                     f(mx-1  ,:,:), &
                                     f(mx    ,:,:))
         elseif (mx.eq.4) then
            maxf(1,:,:)=max(f(1,:,:),f(2,:,:),f(3,:,:))
            maxf(2,:,:)=max(f(1,:,:),f(2,:,:),f(3,:,:),f(4,:,:))
            maxf(3,:,:)=maxf(2,:,:)
            maxf(4,:,:)=max(f(2,:,:),f(3,:,:),f(4,:,:))
         elseif (mx.eq.3) then
            maxf(1,:,:)=max(f(1,:,:),f(2,:,:),f(3,:,:))
            maxf(2,:,:)=maxf(1,:,:)
            maxf(3,:,:)=maxf(1,:,:)
         elseif (mx.eq.2) then
            maxf(1,:,:)=max(f(1,:,:),f(2,:,:))
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
            tmp(:,1     ,:) = max(maxf(:,1     ,:),  &
                                  maxf(:,2     ,:),  &
                                  maxf(:,3     ,:))
            tmp(:,2     ,:) = max(maxf(:,1     ,:), &
                                  maxf(:,2     ,:), &
                                  maxf(:,3     ,:), &
                                  maxf(:,4     ,:))
            tmp(:,3:my-2,:) = max(maxf(:,1:my-4,:), &
                                  maxf(:,2:my-3,:), &
                                  maxf(:,3:my-2,:), &
                                  maxf(:,4:my-1,:), &
                                  maxf(:,5:my  ,:))
            tmp(:,my-1,:)   = max(maxf(:,my-3  ,:), &
                                  maxf(:,my-2  ,:), &
                                  maxf(:,my-1  ,:), &
                                  maxf(:,my    ,:))
            tmp(:,my  ,:)   = max(maxf(:,my-2  ,:), &
                                  maxf(:,my-1  ,:), &
                                  maxf(:,my    ,:))
         elseif (my.eq.4) then
            tmp(:,1,:)=max(maxf(:,1,:),maxf(:,2,:),maxf(:,3,:))
            tmp(:,2,:)=max(maxf(:,1,:),maxf(:,2,:),maxf(:,3,:),maxf(:,4,:))
            tmp(:,3,:)=tmp(:,2,:)
            tmp(:,4,:)=max(maxf(:,2,:),maxf(:,3,:),maxf(:,4,:))
         elseif (my.eq.3) then
            tmp(:,1,:)=max(maxf(:,1,:),maxf(:,2,:),maxf(:,3,:))
            tmp(:,2,:)=f(:,1,:)
            tmp(:,3,:)=f(:,1,:)
         elseif (my.eq.2) then
            tmp(:,1,:)=max(maxf(:,1,:),maxf(:,2,:))
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
            maxf(:,:,1     ) = max(tmp(:,:,1     ),  &
                                     tmp(:,:,2     ),  &
                                     tmp(:,:,3     ))
            maxf(:,:,2     ) = max(tmp(:,:,1     ), &
                                     tmp(:,:,2     ), &
                                     tmp(:,:,3     ), &
                                     tmp(:,:,4     ))
            maxf(:,:,3:mz-2) = max(tmp(:,:,1:mz-4), &
                                     tmp(:,:,2:mz-3), &
                                     tmp(:,:,3:mz-2), &
                                     tmp(:,:,4:mz-1), &
                                     tmp(:,:,5:mz  ))
            maxf(:,:,mz-1  ) = max(tmp(:,:,mz-3  ), &
                                     tmp(:,:,mz-2  ), &
                                     tmp(:,:,mz-1  ), &
                                     tmp(:,:,mz    ))
            maxf(:,:,mz    ) = max(tmp(:,:,mz-2  ), &
                                     tmp(:,:,mz-1  ), &
                                     tmp(:,:,mz    ))
         elseif (mz.eq.4) then
            maxf(:,:,1)=max(tmp(:,:,1),tmp(:,:,2),tmp(:,:,3))
            maxf(:,:,2)=max(tmp(:,:,1),tmp(:,:,2),tmp(:,:,3),tmp(:,:,4))
            maxf(:,:,3)=maxf(:,:,2)
            maxf(:,:,4)=max(tmp(:,:,2),tmp(:,:,3),tmp(:,:,4))
         elseif (mz.eq.3) then
            maxf(:,:,1)=max(tmp(:,:,1),tmp(:,:,2),tmp(:,:,3))
            maxf(:,:,2)=maxf(:,:,1)
            maxf(:,:,3)=maxf(:,:,1)
         elseif (mz.eq.2) then
            maxf(:,:,1)=max(tmp(:,:,1),tmp(:,:,2))
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
!Utility routines - poss need moving elsewhere
    subroutine shock_max3_farray(f,maxf)
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
            maxf(1     ,:,:) = max(f(1     ,:,:), &
                                     f(2     ,:,:))
            maxf(2:mx-1,:,:) = max(f(1:mx-2,:,:), &
                                     f(2:mx-1,:,:), &
                                     f(3:mx  ,:,:))
            maxf(  mx  ,:,:) = max(f(  mx-1,:,:), &
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
            tmp(:,1     ,:) = max(maxf(:,1     ,:),  &
                                    maxf(:,2     ,:))
            tmp(:,2:my-1,:) = max(maxf(:,1:my-2,:), &
                                    maxf(:,2:my-1,:), &
                                    maxf(:,3:my  ,:))
            tmp(:,  my  ,:) = max(maxf(:,  my-1,:), &
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
            maxf(:,:,1     ) = max(tmp(:,:,1     ),  &
                                     tmp(:,:,2     ))
            maxf(:,:,2:mz-1) = max(tmp(:,:,1:mz-2), &
                                     tmp(:,:,2:mz-1), &
                                     tmp(:,:,3:mz  ))
            maxf(:,:,mz    ) = max(tmp(:,:,  mz-1), &
                                     tmp(:,:,  mz  ))
         endif
      else
         maxf=tmp
      endif
!
    endsubroutine shock_max3_farray

!***********************************************************************
!Utility routines - poss need moving elsewhere
    subroutine shock_max3_pencil(f,j,maxf)
!
!  return array maxed with by 2 points either way
!  skipping 1 data point all round
!
!  27-apr-03/tony: adapted from shock_max3
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx) :: maxf
      integer :: j
      integer :: ii,jj,kk

      maxf=f(:,m,n,j)
      if ((nxgrid/=1).and.(nygrid/=1).and.(nzgrid/=1)) then
        do kk=-1,1,2
        do jj=-1,1,2
        do ii=-1,1,2
          maxf(3:mx-2)=max(maxf(3:mx-2),f(3+ii:mx-2+ii,m+jj,n+kk,j))
        enddo
        enddo
        enddo
      elseif ((nxgrid/=1).and.(nygrid/=1)) then
        do jj=-1,1,2
        do ii=-1,1,2
          maxf(3:mx-2)=max(maxf(3:mx-2),f(3+ii:mx-2+ii,m+jj,n1,j))
        enddo
        enddo
      elseif ((nxgrid/=1).and.(nzgrid/=1)) then
        do kk=-1,1,2
        do ii=-1,1,2
          maxf(3:mx-2)=max(maxf(3:mx-2),f(3+ii:mx-2+ii,m1,n+kk,j))
        enddo
        enddo
      elseif ((nygrid/=1).and.(nzgrid/=1)) then
        do kk=-1,1,2
        do jj=-1,1,2
          maxf(l1:l2)=max(maxf(l1:l2),f(l1:l2,m+jj,n+kk,j))
        enddo
        enddo
      elseif (nxgrid/=1) then 
        do ii=-1,1,2
          maxf(3:mx-2)=max(maxf(3:mx-2),f(3+ii:mx-2+ii,m1,n1,j))
        enddo
      elseif (nygrid/=1) then 
        do jj=-1,1,2
          maxf(l1:l2)=max(maxf(l1:l2),f(l1:l2,m+jj,n1,j))
        enddo
      elseif (nzgrid/=1) then 
        do kk=-1,1,2
          maxf(l1:l2)=max(maxf(l1:l2),f(l1:l2,m1,n+kk,j))
        enddo
      endif
!
    endsubroutine shock_max3_pencil
!***********************************************************************
    subroutine shock_max5_pencil(f,j,maxf)
!
!  return array maxed with by 2 points either way
!  skipping 1 data point all round
!
!  27-apr-03/tony: adapted from shock_max3
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx) :: maxf
      integer :: j
      integer :: ii,jj,kk

      maxf=f(:,m,n,j)
      if ((nxgrid/=1).and.(nygrid/=1).and.(nzgrid/=1)) then
        do kk=-2,2
          if (kk==0) cycle
        do jj=-2,2
          if (jj==0) cycle
        do ii=-2,2
          if (ii==0) cycle
          maxf(2+ii:mx-1+ii)=max(maxf(2+ii:mx-1+ii),f(2+ii:mx-1+ii,m+jj,n+kk,j))
        enddo
        enddo
        enddo
      elseif ((nxgrid/=1).and.(nygrid/=1)) then
        do jj=-2,2
          if (jj==0) cycle
        do ii=-2,2
          if (ii==0) cycle
          maxf(2+ii:mx-1+ii)=max(maxf(2+ii:mx-1+ii),f(2+ii:mx-1+ii,m+jj,n1,j))
        enddo
        enddo
      elseif ((nxgrid/=1).and.(nzgrid/=1)) then
        do kk=-2,2
          if (kk==0) cycle
        do ii=-2,2
          if (ii==0) cycle
          maxf(2+ii:mx-1+ii)=max(maxf(2+ii:mx-1+ii),f(2+ii:mx-1+ii,m1,n+kk,j))
        enddo
        enddo
      elseif ((nygrid/=1).and.(nzgrid/=1)) then
        do kk=-2,2
          if (kk==0) cycle
        do jj=-2,2
          if (jj==0) cycle
          maxf(l1:l2)=max(maxf(l1:l2),f(l1:l2,m+jj,n+kk,j))
        enddo
        enddo
      elseif (nxgrid/=1) then 
        do ii=-2,2
          if (ii==0) cycle
          maxf(2+ii:mx-1+ii)=max(maxf(2+ii:mx-1+ii),f(2+ii:mx-1+ii,m1,n1,j))
        enddo
      elseif (nygrid/=1) then 
        do jj=-2,2
          if (jj==0) cycle
          maxf(l1:l2)=max(maxf(l1:l2),f(l1:l2,m+jj,n1,j))
        enddo
      elseif (nzgrid/=1) then 
        do kk=-2,2
          if (kk==0) cycle
          maxf(l1:l2)=max(maxf(l1:l2),f(l1:l2,m1,n+kk,j))
        enddo
      endif
!
    endsubroutine shock_max5_pencil
!***********************************************************************
    subroutine shock_smooth_farray(f,smoothf)
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

    endsubroutine shock_smooth_farray
!***********************************************************************
    subroutine shock_smooth_pencil(f,smoothf)
!
!  return array smoothed with by 2 points either way
!  skipping 3 data point all round 
!  i.e. result valid ()
!
!  23-nov-02/tony: coded
!
      use Cdata, only: mx,my,mz,m,n
!
      real, dimension (mx,my,mz) :: f
      real, dimension (mx) :: smoothf
      integer :: ii,jj,kk
!
      smoothf=0.
      if ((nxgrid/=1).and.(nygrid/=1).and.(nzgrid/=1)) then
        do kk=-1,1
        do jj=-1,1
        do ii=-1,1
          smoothf(4:mx-3)=smoothf(4:mx-3)+smooth_factor(2+ii,2+jj,2+kk)*f(4+ii:mx-3+ii,m+jj,n+kk)
        enddo
        enddo
        enddo
      elseif ((nxgrid/=1).and.(nygrid/=1)) then
        do jj=-1,1
        do ii=-1,1
          smoothf(4:mx-3)=smoothf(4:mx-3)+smooth_factor(2+ii,2+jj,2)*f(4+ii:mx-3+ii,m+jj,n1)
        enddo
        enddo
      elseif ((nxgrid/=1).and.(nzgrid/=1)) then
        do kk=-1,1
        do ii=-1,1
          smoothf(4:mx-3)=smoothf(4:mx-3)+smooth_factor(2+ii,2,2+kk)*f(4+ii:mx-3+ii,m1,n+kk)
        enddo
        enddo
      elseif ((nygrid/=1).and.(nzgrid/=1)) then
        do kk=-1,1
        do jj=-1,1
          smoothf(l1:l2)=smoothf(l1:l2)+smooth_factor(2,2+jj,2+kk)*f(l1:l2,m+jj,n+kk)
        enddo
        enddo
      elseif (nxgrid/=1) then 
        do ii=-1,1
          smoothf(4:mx-3)=smoothf(4:mx-3)+smooth_factor(2+ii,2,2)*f(4+ii:mx-3+ii,m1,n1)
        enddo
      elseif (nygrid/=1) then 
        do jj=-1,1
          smoothf(l1:l2)=smoothf(l1:l2)+smooth_factor(2,2+jj,2)*f(l1:l2,m+jj,n1)
        enddo
      elseif (nzgrid/=1) then 
        do kk=-1,1
          smoothf(l1:l2)=smoothf(l1:l2)+smooth_factor(2,2,2+kk)*f(l1:l2,m1,n+kk)
        enddo
      endif
!
    endsubroutine shock_smooth_pencil
!***********************************************************************
    subroutine shock_divu_farray(f,df)
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
      real :: fac
!
      df=0.
!
      if (nxgrid/=1) then
         fac=1./(2.*dx)
         df(1     ,:,:) =  df(1    ,:,:) &
                           + (  4.*f(2,:,:,iux) &
                              - 3.*f(1,:,:,iux) &
                              -    f(3,:,:,iux))*fac
         df(2:mx-1,:,:) =  df(2:mx-1,:,:) &
                           + ( f(3:mx,:,:,iux)-f(1:mx-2,:,:,iux) ) * fac
         df(mx    ,:,:) =  df(mx    ,:,:) &
                           + (  3.*f(mx  ,:,:,iux) &
                              - 4.*f(mx-1,:,:,iux) &
                              +    f(mx-2,:,:,iux))*fac
      endif

      if (nygrid/=1) then
         fac=1./(2.*dy)
         df(:,1     ,:) = df(:,1     ,:) &
                          + (  4.*f(:,2,:,iuy) &
                             - 3.*f(:,1,:,iuy) &
                             -    f(:,3,:,iuy))*fac
         df(:,2:my-1,:) = df(:,2:my-1,:) &  
                          + (f(:,3:my,:,iuy)-f(:,1:my-2,:,iuy))*fac
         df(:,my    ,:) = df(:,my    ,:) & 
                          + (  3.*f(:,my  ,:,iuy) &
                             - 4.*f(:,my-1,:,iuy) &
                             +    f(:,my-2,:,iuy))*fac
      endif

      if (nzgrid/=1) then
         fac=1./(2.*dz)
         df(:,:,1     ) = df(:,:,1     ) &
                          + (  4.*f(:,:,2,iuz) &
                             - 3.*f(:,:,1,iuz) &
                             -    f(:,:,3,iuz))*fac
         df(:,:,2:mz-1) = df(:,:,2:mz-1) & 
                          + (f(:,:,3:mz,iuz)-f(:,:,1:mz-2,iuz))*fac
         df(:,:,mz    ) = df(:,:,mz    ) &   
                          + (  3.*f(:,:,mz  ,iuz) &
                             - 4.*f(:,:,mz-1,iuz) &
                             +    f(:,:,mz-2,iuz))*fac
      endif
!      
    endsubroutine shock_divu_farray
!***********************************************************************
    subroutine shock_divu_pencil(f,j,df)
!
!  calculate divergence of a vector U, get scalar
!  accurate to 2nd order, explicit,  centred an left and right biased 
!
!  23-nov-02/tony: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx) :: df
      real :: fac
      integer :: j
!
      df=0.
      if (nxgrid/=1) then
         fac=1./(2.*dx)
         df(2:mx-1) = df(2:mx-1)     &
               + (f(3:mx  ,m     ,n     ,j  ) &
               -  f(1:mx-2,m     ,n     ,j  ) ) &
               * fac
      endif

      if (nygrid/=1) then
         fac=1./(2.*dy)
         df = df  &  
               + (f(:     ,m+1   ,n     ,j+1)   &
               -  f(:     ,m-1   ,n     ,j+1) ) &
               * fac
      endif

      if (nzgrid/=1) then
         fac=1./(2.*dz)
         df = df       & 
               + (f(:     ,m     ,n+1   ,j+2)   &
               -  f(:     ,m     ,n-1   ,j+2) ) &
               * fac
      endif
    endsubroutine shock_divu_pencil
!***********************************************************************
    subroutine shock_smooth_cube_diamond7(f,df)
!
!  calculate divergence of a vector U smoothed over a 7x7x7 cube
!  by using avergaging by a volume integral and transforming to 
!  a surface flux integral using Gauss' Theorm
!
!  01-apr-05/tony: coded
!
      use Cdata
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: df
      real :: cube, diamond, fac
      integer :: i,j,k
!
      df=0.
!
      if ((nxgrid/=1).and.(nygrid/=1).and.(nzgrid/=1)) then
        do k=n1,n2
        do j=m1,m2
        do i=l1,l2
!          df(i,j,k)=sum(f(i-3:i+3, j-3:j+3, k+3    , iuz)) - sum(f(i-3:i+3, j-3:j+3, k-3    , iuz)) + &
!                    sum(f(i-3:i+3, j+3    , k-3:k+3, iuy)) - sum(f(i-3:i+3, j-3    , k-3:k+3, iuy)) + &
!                    sum(f(i+3    , j-3:j+3, k-3:k+3, iux)) - sum(f(i-3    , j-3:j+3, k-3:k+3, iux))
        enddo
        enddo
        enddo
        df=df*(dx*dy)
      elseif ((nxgrid/=1).and.(nzgrid/=1)) then
        do k=n1,n2
        do j=m1,m2
        do i=l1,l2
          fac = (-1./18.) / dx
          diamond = ( &
                       f(i  , j, k+3 , iuz)  &
                     + f(i+1, j, k+2 , iuz) + f(i+1, j, k+2 , iux)  &
                     + f(i+2, j, k+1 , iuz) + f(i+2, j, k+1 , iux)  &
                     + f(i+3, j, k   , iux) &
                     - f(i+2, j, k-1 , iuz) + f(i+2, j, k-1 , iux)  &
                     - f(i+1, j, k-2 , iuz) + f(i+1, j, k-2 , iux)  &
                     - f(i  , j, k-3 , iuz) & 
                     - f(i-1, j, k-2 , iuz) - f(i-1, j, k-2 , iux)  &
                     - f(i-2, j, k-1 , iuz) - f(i-2, j, k-1 , iux)  &
                     - f(i-3, j, k   , iux)  &
                     + f(i-2, j, k+1 , iuz) - f(i-2, j, k+1 , iux)  &
                     + f(i-1, j, k+2 , iuz) - f(i-1, j, k+2 , iux)  &
                    ) * fac

          fac = (1./16.) / dx
          cube    = (  sum(f(i-2:i+2, j, k-2    , iuz)) &
                     - sum(f(i-2:i+2, j, k+2    , iuz)) &
                     + sum(f(i-2    , j, k-2:k+2, iux)) &
                     - sum(f(i+2    , j, k-2:k+2, iux)) &
                    ) * fac

          if (lwith_extreme_div) then
            if ((diamond .lt. 0.) .and. (diamond .gt. -div_threshold)) then
              diamond=0.
            elseif (diamond.lt.-div_threshold) then
              diamond=diamond*div_scaling
            endif
            if ((cube .lt. 0.) .and. (cube .gt. -div_threshold)) then
              cube=0.
            elseif (cube.lt.-div_threshold) then
              cube=cube*div_scaling
            endif
            df(i,j,k)=0.5*(diamond + cube)
            !df(i,j,k)=sqrt(0.5*(diamond**2 + cube**2))
          else
            df(i,j,k)=0.5*(max(diamond,0.) + max(cube,0.))
            !df(i,j,k)=sqrt(0.5*(max(diamond,0.)**2 + max(cube,0.)**2))
          endif
        enddo
        enddo
        enddo
      else
        call stop_it('shock_smooth_cube_diamond7: CASE NOT IMPLEMENTED');
      endif
!
    endsubroutine shock_smooth_cube_diamond7
!***********************************************************************
    function scale_and_chop(value)
!
      real :: scale_and_chop
      real :: value
!
      if (lwith_extreme_div) then
        if ((value .lt. 0.) .and. (value .gt. -div_threshold)) then
          scale_and_chop=0.
        elseif (value.lt.-div_threshold) then
          scale_and_chop=-value*div_scaling
        else
          scale_and_chop=value
        endif
      else
        scale_and_chop=max(value,0.)
      endif
!
    endfunction scale_and_chop
!***********************************************************************
    subroutine scale_and_chop_internalboundary(f)
!
      use Cdata
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer :: i,j,k
!
      if ((nxgrid/=1).and.(nygrid/=1).and.(nzgrid/=1)) then
        call stop_it("scale_and_chop_internalboundary: 3D not implemented yet...")
      elseif ((nxgrid/=1).and.(nygrid/=1)) then
        do j=m1i+1,m2i-1
        do i=0,nghost-1
          f(l1+i,j,n1,ishock)=scale_and_chop(f(l1+i,j,n1,ishock))
          f(l2-i,j,n1,ishock)=scale_and_chop(f(l2-i,j,n1,ishock))
        enddo
        enddo
        do j=0,nghost-1
        do i=l1,l2
          f(i,m1+j,n1,ishock)=scale_and_chop(f(i,m1+j,n1,ishock))
          f(i,m2-j,n1,ishock)=scale_and_chop(f(i,m2-j,n1,ishock))
        enddo
        enddo
      elseif ((nxgrid/=1).and.(nzgrid/=1)) then
        do k=n1i+1,n2i-1
        do i=0,nghost-1
          f(l1+i,m1,k,ishock)=scale_and_chop(f(l1+i,m1,k,ishock))
          f(l2-i,m1,k,ishock)=scale_and_chop(f(l2-i,m1,k,ishock))
        enddo
        enddo
        do k=0,nghost-1
        do i=l1,l2
          f(i,m1,n1+k,ishock)=scale_and_chop(f(i,m1,n1+k,ishock))
          f(i,m1,n2-k,ishock)=scale_and_chop(f(i,m1,n2-k,ishock))
        enddo
        enddo
      elseif ((nygrid/=1).and.(nzgrid/=1)) then
        do k=0,nghost-1
        do j=m1i+1,m2i-1
          f(l1,j,n1+k,ishock)=scale_and_chop(f(l1,j,n1+k,ishock))
          f(l1,j,n2-k,ishock)=scale_and_chop(f(l1,j,n2-k,ishock))
        enddo
        enddo
        do k=n1,n2
        do j=0,nghost-1
          f(l1,m1+j,k,ishock)=scale_and_chop(f(l1,m1+j,k,ishock))
          f(l1,m2-j,k,ishock)=scale_and_chop(f(l1,m2-j,k,ishock))
        enddo
        enddo
      elseif (nxgrid/=1) then
        do i=0,nghost-1
          f(l1+i,m1,n1,ishock)=scale_and_chop(f(l1+i,m1,n1,ishock))
          f(l2-i,m1,n1,ishock)=scale_and_chop(f(l2-i,m1,n1,ishock))
        enddo
      elseif (nygrid/=1) then
        do j=0,nghost-1
          f(l1,m1+j,n1,ishock)=scale_and_chop(f(l1,m1+j,n1,ishock))
          f(l1,m2-j,n1,ishock)=scale_and_chop(f(l1,m2-j,n1,ishock))
        enddo
      elseif (nzgrid/=1) then
        do k=0,nghost-1
          f(n1,m1,n1+k,ishock)=scale_and_chop(f(l1,m1,n1+k,ishock))
          f(n1,m1,n2-k,ishock)=scale_and_chop(f(l1,m1,n2-k,ishock))
        enddo
      endif 
!
    endsubroutine scale_and_chop_internalboundary
!***********************************************************************
    subroutine shock_smooth_octagon7(f,df)
!
!  calculate divergence of a vector U smoothed over a 7x7x7 cube
!  by using avergaging by a volume integral and transforming to 
!  a surface flux integral using Gauss' Theorm
!
!  01-apr-05/tony: coded
!
      use Cdata
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: df
      real :: octagon, fac_diag, fac_straight
      integer :: i,j,k
!
      df=0.
!
      if ((nxgrid/=1).and.(nygrid/=1).and.(nzgrid/=1)) then
        do k=n1,n2
        do j=m1,m2
        do i=l1,l2
!          df(i,j,k)=sum(f(i-3:i+3, j-3:j+3, k+3    , iuz)) - sum(f(i-3:i+3, j-3:j+3, k-3    , iuz)) + &
!                    sum(f(i-3:i+3, j+3    , k-3:k+3, iuy)) - sum(f(i-3:i+3, j-3    , k-3:k+3, iuy)) + &
!                    sum(f(i+3    , j-3:j+3, k-3:k+3, iux)) - sum(f(i-3    , j-3:j+3, k-3:k+3, iux))
        enddo
        enddo
        enddo
        df=df*(dx*dy)
      elseif ((nxgrid/=1).and.(nzgrid/=1)) then
!
!  Top
!
        fac_diag     = -sqrt(dx**2+dz**2)/(25.*dx*dz)
        fac_straight = -(1./(25*dx))
        k=n2+1
        j=m1
        do i=l1,l2
          octagon = ( &
                     + f(i+2 , j, k+2 , iuz) + f(i+2 , j, k+2 , iux) &
                     + f(i+3 , j, k+1 , iuz) + f(i+3 , j, k+1 , iux) &
                     - f(i+1 , j, k-3 , iuz) + f(i+1 , j, k-3 , iux) &  ! Bottom right /
                     - f(i+2 , j, k-2 , iuz) + f(i+2 , j, k-2 , iux) &
                     - f(i+3 , j, k-1 , iuz) + f(i+3 , j, k-1 , iux) &
                     - f(i-1 , j, k-3 , iuz) - f(i-1 , j, k-3 , iux) &  ! Bottom left \
                     - f(i-2 , j, k-2 , iuz) - f(i-2 , j, k-2 , iux) &
                     - f(i-3 , j, k-1 , iuz) - f(i-3 , j, k-1 , iux) &
                     + f(i-2 , j, k+2 , iuz) - f(i-2 , j, k+2 , iux) &
                     + f(i-3 , j, k+1 , iuz) - f(i-3 , j, k+1 , iux) &
                    ) * fac_diag &
                  + (  &
                       f(i+3 , j, k-1 , iux) - f(i-3 , j, k-1 , iux) &  ! left and right |
                     + f(i+3 , j, k   , iux) - f(i-3 , j, k   , iux) &
                     + f(i+3 , j, k+1 , iux) - f(i-3 , j, k+1 , iux) &
                     + f(i-1 , j, k+2 , iuz) - f(i-1 , j, k-3 , iuz) &  ! top and bottom -
                     + f(i   , j, k+2 , iuz) - f(i   , j, k-3 , iuz) &
                     + f(i+1 , j, k+2 , iuz) - f(i+1 , j, k-3 , iuz) &
                     + f(i-2 , j, k+2 , iuz) &
                     + f(i+2 , j, k+2 , iuz) &  
                    ) * fac_straight 
!
          df(i,j,k)=scale_and_chop(octagon)
        enddo
!
!  Bottom
!
        fac_diag     = -sqrt(dx**2+dz**2)/(25.*dx*dz)
        fac_straight = -(1./(25*dx))
        k=n2+1
        j=m1
        do i=l1,l2
          octagon = ( &
                     + f(i+1 , j, k+3 , iuz) + f(i+1 , j, k+3 , iux) &
                     + f(i+2 , j, k+2 , iuz) + f(i+2 , j, k+2 , iux) &
                     + f(i+3 , j, k+1 , iuz) + f(i+3 , j, k+1 , iux) &
                     - f(i+2 , j, k-2 , iuz) + f(i+2 , j, k-2 , iux) &
                     - f(i+3 , j, k-1 , iuz) + f(i+3 , j, k-1 , iux) &
                     - f(i-2 , j, k-2 , iuz) - f(i-2 , j, k-2 , iux) &
                     - f(i-3 , j, k-1 , iuz) - f(i-3 , j, k-1 , iux) &
                     + f(i-2 , j, k+2 , iuz) - f(i-2 , j, k+2 , iux) &
                     + f(i-3 , j, k+1 , iuz) - f(i-3 , j, k+1 , iux) &
                    ) * fac_diag &
                  + (  &
                       f(i+3 , j, k-1 , iux) - f(i-3 , j, k-1 , iux) &  ! left and right |
                     + f(i+3 , j, k   , iux) - f(i-3 , j, k   , iux) &
                     + f(i+3 , j, k+1 , iux) - f(i-3 , j, k+1 , iux) &
                     + f(i-1 , j, k+3 , iuz) - f(i-1 , j, k-2 , iuz) &  ! top and bottom -
                     + f(i   , j, k+3 , iuz) - f(i   , j, k-2 , iuz) &
                     + f(i+1 , j, k+3 , iuz) - f(i+1 , j, k-2 , iuz) &
                     - f(i-2 , j, k-2 , iuz) &
                     - f(i+2 , j, k-2 , iuz) &
                    ) * fac_straight 
!
          df(i,j,k)=scale_and_chop(octagon)
        enddo
!
!  Bulk
!
        fac_diag     = -sqrt(dx**2+dz**2)/(28.*dx*dz)
        fac_straight = -(1./(28.*dx))
        j=m1
        do k=n1i,n2i
        do i=l1i,l2i
          octagon = ( &
                       f(i+1 , j, k+3 , iuz) + f(i+1 , j, k+3 , iux) &  ! Top right \
                     + f(i+2 , j, k+2 , iuz) + f(i+2 , j, k+2 , iux) &
                     + f(i+3 , j, k+1 , iuz) + f(i+3 , j, k+1 , iux) &
                     - f(i+1 , j, k-3 , iuz) + f(i+1 , j, k-3 , iux) &  ! Bottom right /
                     - f(i+2 , j, k-2 , iuz) + f(i+2 , j, k-2 , iux) &
                     - f(i+3 , j, k-1 , iuz) + f(i+3 , j, k-1 , iux) &
                     - f(i-1 , j, k-3 , iuz) - f(i-1 , j, k-3 , iux) &  ! Bottom left \
                     - f(i-2 , j, k-2 , iuz) - f(i-2 , j, k-2 , iux) &
                     - f(i-3 , j, k-1 , iuz) - f(i-3 , j, k-1 , iux) &
                     + f(i-1 , j, k+3 , iuz) - f(i-1 , j, k+3 , iux) &  ! Top left /
                     + f(i-2 , j, k+2 , iuz) - f(i-2 , j, k+2 , iux) &
                     + f(i-3 , j, k+1 , iuz) - f(i-3 , j, k+1 , iux) &
                    ) * fac_diag &
                  + (  &
                       f(i+3 , j, k-1 , iux) - f(i-3 , j, k-1 , iux) &  ! left and right |
                     + f(i+3 , j, k   , iux) - f(i-3 , j, k   , iux) &
                     + f(i+3 , j, k+1 , iux) - f(i-3 , j, k+1 , iux) &
                     + f(i-1 , j, k+3 , iuz) - f(i-1 , j, k-3 , iuz) &  ! top and bottom -
                     + f(i   , j, k+3 , iuz) - f(i   , j, k-3 , iuz) &
                     + f(i+1 , j, k+3 , iuz) - f(i+1 , j, k-3 , iuz) &
                    ) * fac_straight 
!
          df(i,j,k)=scale_and_chop(octagon)
        enddo
        enddo
      else
        call stop_it('shock_smooth_octagon7: CASE NOT IMPLEMENTED');
      endif
!
    endsubroutine shock_smooth_octagon7
!    
    include 'shock_profile.inc'
!
!***********************************************************************
    subroutine bcshock_per_x(f)
!
!  periodic boundary condition
!  11-nov-02/wolf: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      if (nprocx==1)  then 
        f(l1:l1i ,:,:,ishock) = f(l1:l1i,:,:,ishock) + f(l2+1:mx,:,:,ishock)
        f(l2+1:mx,:,:,ishock) = f(l1:l1i,:,:,ishock)
        f(l2i:l2 ,:,:,ishock) = f(l2i:l2,:,:,ishock) + f(1:l1-1 ,:,:,ishock) 
        f(1:l1-1 ,:,:,ishock) = f(l2i:l2,:,:,ishock)
      endif
!
    endsubroutine bcshock_per_x
!***********************************************************************
    subroutine bcshock_per_y(f)
!
!  periodic boundary condition
!  11-nov-02/wolf: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f

      if (nprocy==1) then
        f(:,m1:m1i,:,ishock)  = f(:,m1:m1i,:,ishock) + f(:,m2+1:my,:,ishock)
        f(:,m2+1:my,:,ishock) = f(:,m1:m1i,:,ishock)
        f(:,m2i:m2,:,ishock)  = f(:,m2i:m2,:,ishock) + f(:,1:m1-1 ,:,ishock)
        f(:,1:m1-1,:,ishock)  = f(:,m2i:m2,:,ishock)
      endif
!
    endsubroutine bcshock_per_y
!***********************************************************************
    subroutine bcshock_per_z(f)
!
!  periodic boundary condition
!  11-nov-02/wolf: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      if (nprocz==1)  then 
        f(:,:,n1:n1i ,ishock) = f(:,:,n1:n1i,ishock) + f(:,:,n2+1:mz,ishock)
        f(:,:,n2+1:mz,ishock) = f(:,:,n1:n1i,ishock)
        f(:,:,n2i:n2 ,ishock) = f(:,:,n2i:n2,ishock) + f(:,:,1:n1-1 ,ishock) 
        f(:,:,1:n1-1 ,ishock) = f(:,:,n2i:n2,ishock)
      endif
!
    endsubroutine bcshock_per_z
endmodule Shock
