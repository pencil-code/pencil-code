! $Id: cosmicray.f90,v 1.1 2003-10-09 16:31:57 mee Exp $

!  This modules solves the cosmic ray energy density advection difussion equation
!  it follows the description of Hanasz & Lesch (2002,2003) as used in their
!  ZEUS 3D implementation
!

module CosmicRay

  use Cparam
  use Cdata

  implicit none

  character (len=labellen) :: initecr='zero', initecr2='zero'
  character (len=40) :: tensor_pscalar_file

  real, parameter :: gammacr = 4./3., gammacr1=1./3.

  ! input parameters
  real :: amplecr=.1, widthecr=.5, ecr_min=0.
  real :: amplecr2=0.,kx_ecr=1.,ky_ecr=1.,kz_ecr=1.,radius_ecr=0.,epsilon_ecr=0.
!  real, dimension(3) :: gradC0=(/0.,0.,0./)

  namelist /cosmicray_init_pars/ &
       initlncc,initlncc2,ampllncc,ampllncc2,kx_lncc,ky_lncc,kz_lncc, &
       radius_lncc,epsilon_lncc,widthlncc,cc_min

  ! run parameters
  real :: cosmicray_diff=0.,tensor_cosmicray_diff=0.

  namelist /pscalar_run_pars/ &
       cosmicray_diff,tensor_pscalar_diff
       ! ,gradC0

  ! other variables (needs to be consistent with reset list below)
!  integer :: i_rhoccm=0,i_ccmax=0,i_lnccm=0,i_lnccmz=0
!  integer :: i_ucm=0,i_uudcm=0,i_Cz2m=0,i_Cz4m=0,i_Crmsm=0

  contains

!***********************************************************************
    subroutine register_cosmicray()
!
!  Initialise variables which should know that we solve for active 
!  scalar: iecr - the cosmic ray energy density; increase nvar accordingly
!
!  09-oct-03/tony: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_cosmicray called twice')
      first = .false.
!
      lcosmicray = .true.
      iecr = nvar+1            ! index to access lncc
      nvar = nvar+1             ! added 1 variable
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_cosmicray: nvar = ', nvar
        print*, 'register_cosmicray: iecr = ', iecr
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: cosmicray.f90,v 1.1 2003-10-09 16:31:57 mee Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('register_cosmicray: nvar > mvar')
      endif
!
!  Writing files for use with IDL
!
      if (maux == 0) then
         if (nvar < mvar) write(4,*) ',ecr $'
         if (nvar == mvar) write(4,*) ',ecr'
      else
         write(4,*) ',ecr $'
      endif
      write(5,*) 'ecr = fltarr(mx,my,mz)*one'
!
    endsubroutine register_cosmicray
!***********************************************************************
    subroutine initialize_cosmicray(f)
!
!  Perform any necessary post-parameter read initialization
!  Dummy routine
!
!  09-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f
! 
!  set to zero and then call the same initial condition
!  that was used in start.csh
!   
      if(ip==0) print*,'f=',f
    endsubroutine initialize_cosmicray
!***********************************************************************
    subroutine init_ecr(f,xx,yy,zz)
!
!  initialise cosmic ray energy density field; called from start.f90
!
!   09-oct-03/tony: coded
!
      use Cdata
      use Mpicomm
      use Density
      use Sub
      use Initcond
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz)      :: xx,yy,zz,prof
!
      select case(initecr)
        case('zero'); f(:,:,:,iecr)=0.
        case('gaussian-x'); call gaussian(amplecr,f,iecr,kx=kx_ecr)
        case('gaussian-y'); call gaussian(amplecr,f,iecr,ky=ky_ecr)
        case('gaussian-z'); call gaussian(amplecr,f,iecr,kz=kz_ecr)
        case('parabola-x'); call parabola(amplecr,f,iecr,kx=kx_ecr)
        case('parabola-y'); call parabola(amplecr,f,iecr,ky=ky_ecr)
        case('parabola-z'); call parabola(amplecr,f,iecr,kz=kz_ecr)
        case('gaussian-noise'); call gaunoise(amplecr,f,ilncc,iecr)
        case('wave-x'); call wave(amplecr,f,iecr,kx=kx_ecr)
        case('wave-y'); call wave(amplecr,f,iecr,ky=ky_ecr)
        case('wave-z'); call wave(amplecr,f,iecr,kz=kz_ecr)
        case('propto-ux'); call wave_uu(amplecr,f,iecr,kx=kx_ecr)
        case('propto-uy'); call wave_uu(amplecr,f,iecr,ky=ky_ecr)
        case('propto-uz'); call wave_uu(amplecr,f,iecr,kz=kz_ecr)
        case('tang-discont-z')
           print*,'init_lncc: widthlncc=',widthecr
        prof=.5*(1.+tanh(zz/widthecr))
        f(:,:,:,iecr)=-1.+2.*prof
        case('hor-tube'); call htube2(amplecr,f,iecr,iecr,xx,yy,zz,radius_ecr,epsilon_ecr)
        case default; call stop_it('init_lncc: bad initlncc='//trim(initecr))
      endselect
!
!  superimpose something else
!
      select case(initecr2)
        case('wave-x'); call wave(amplecr2,f,iecr,ky=5.)
      endselect
!
!  add floor value if cc_min is set
!
      if(cc_min/=0.) then
        lncc_min=alog(cc_min)
        if(lroot) print*,'set floor value for ecr; ecr_min=',ecr_min
        f(:,:,:,iecr)=amax1(ecr_min,f(:,:,:,iecr))
      endif
!
      if(ip==0) print*,xx,yy,zz !(prevent compiler warnings)
    endsubroutine init_ecr
!***********************************************************************
    subroutine decr_dt(f,df,uu,glnrho)
!
!  cosmic ray evolution
!  calculate decr/dt=
!-uu.glncc + pscaler_diff*[del2lncc + (glncc+glnrho).glncc]
!
!   09-oct-03/tony: coded
!
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: uu,glncc,glnrho
      real, dimension (nx) :: ecr,rho,uglncc,diff_op,del2lnecr
      integer :: j
!
      intent(in)  :: f,uu,glnrho
      intent(out) :: df
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'SOLVE decr_dt'
      if (headtt) call identify_bcs('ecr',iecr)

      call stop_it('decr_dt: NOT IMPLEMENTED YET')

!
!  gradient of 
!  allow for possibility to turn off passive scalar
!  without changing file size and recompiling everything.
!
!!      if (.not. nopscalar) then ! i.e. if (pscalar)
!!        call grad(f,ilncc,glncc)
!!        call dot_mn(uu,glncc,uglncc)
!
!  passive scalar equation
!
!!        if(lhydro) df(l1:l2,m,n,ilncc)=df(l1:l2,m,n,ilncc)-uglncc
!
!  diffusion operator
!
!!        if (pscalar_diff/=0.) then
!!          if(headtt) print*,'dlncc_dt: pscalar_diff=',pscalar_diff
!!          call dot_mn(glncc+glnrho,glncc,diff_op)
!!          call del2(f,ilncc,del2lncc)
!!          diff_op=diff_op+del2lncc
!!          df(l1:l2,m,n,ilncc)=df(l1:l2,m,n,ilncc)+pscalar_diff*diff_op
!!        endif
!
!  add advection of imposed constant gradient of lncc (called gradC0)
!  makes sense really only for periodic boundary conditions
!  This gradient can have arbitary direction.
!ajwm  temporarily removes... while modifying from pscalar
!        do j=1,3
!          if (gradC0(j)/=0.) then
!            df(l1:l2,m,n,ilncc)=df(l1:l2,m,n,ilncc)-gradC0(j)*uu(:,j)
!          endif
!        enddo
!
!  tensor diffusion (but keep the isotropic one)
!
!!        if (tensor_pscalar_diff/=0.) call tensor_diff(f,df,tensor_pscalar_diff,glncc)
!
!!      endif
!
!  diagnostics
!
!  output for double and triple correlators (assume z-gradient of cc)
!  <u_k u_j d_j c> = <u_k c uu.gradlncc>
!
!!      if (ldiagnos) then
!!        lncc=f(l1:l2,m,n,ilncc)
!!        cc=exp(lncc)
!!        rho=exp(f(l1:l2,m,n,ilnrho))
!!        if (i_rhoccm/=0) call sum_mn_name(rho*cc,i_rhoccm)
!!        if (i_ccmax/=0) call max_mn_name(cc,i_ccmax)
!!        if (i_lnccmz/=0) call xysum_mn_name_z(lncc,i_lnccmz)
!!        if (i_ucm/=0) call sum_mn_name(uu(:,3)*cc,i_ucm)
!!        if (i_uudcm/=0) call sum_mn_name(uu(:,3)*cc*uglncc,i_uudcm)
!!        if (i_Cz2m/=0) call sum_mn_name(rho*cc*z(n)**2,i_Cz2m)
!!        if (i_Cz4m/=0) call sum_mn_name(rho*cc*z(n)**4,i_Cz4m)
!!        if (i_Crmsm/=0) call sum_mn_name((rho*cc)**2,i_Crmsm,lsqrt=.true.)
!!      endif
!
    endsubroutine decr_dt
!***********************************************************************
    subroutine rprint_cosmicray(lreset)
!
!  reads and registers print parameters relevant for magnetic fields
!
!   6-jul-02/axel: coded
!
      use Sub
!
      integer :: iname,inamez
      logical :: lreset
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
!!        i_rhoccm=0; i_ccmax=0; i_lnccm=0; i_lnccmz=0
!!        i_ucm=0; i_uudcm=0; i_Cz2m=0; i_Cz4m=0; i_Crmsm=0
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
!!        call parse_name(iname,cname(iname),cform(iname),'rhoccm',i_rhoccm)
!!        call parse_name(iname,cname(iname),cform(iname),'ccmax',i_ccmax)
!!        call parse_name(iname,cname(iname),cform(iname),'lnccm',i_lnccm)
!!        call parse_name(iname,cname(iname),cform(iname),'ucm',i_ucm)
!!        call parse_name(iname,cname(iname),cform(iname),'uudcm',i_uudcm)
!!        call parse_name(iname,cname(iname),cform(iname),'Cz2m',i_Cz2m)
!!        call parse_name(iname,cname(iname),cform(iname),'Cz4m',i_Cz4m)
!!        call parse_name(iname,cname(iname),cform(iname),'Crmsm',i_Crmsm)
      enddo
!
!  check for those quantities for which we want xy-averages
!
      do inamez=1,nnamez
!        call parse_name(inamez,cnamez(inamez),cformz(inamez),'lnccmz',i_lnccmz)
      enddo
!
!  write column where which magnetic variable is stored
!
!!      write(3,*) 'i_rhoccm=',i_rhoccm
!!      write(3,*) 'i_ccmax=',i_ccmax
!!      write(3,*) 'i_lnccm=',i_lnccm
!!      write(3,*) 'i_ucm=',i_ucm
!!      write(3,*) 'i_uudcm=',i_uudcm
!!      write(3,*) 'i_lnccmz=',i_lnccmz
!!      write(3,*) 'i_Cz2m=',i_Cz2m
!!      write(3,*) 'i_Cz4m=',i_Cz4m
!!      write(3,*) 'i_Crmsm=',i_Crmsm
      write(3,*) 'iecr=',iecr
!
    endsubroutine rprint_cosmicray
!***********************************************************************
    subroutine tensor_diff(f,df,tensor_pscalar_diff,glncc)
!
!  reads file
!
!  11-jul-02/axel: coded
!
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, save, dimension (nx,ny,nz,3) :: bunit,hhh
      real, dimension (nx,3,3) :: g
      real, dimension (nx,3) :: glncc
      real, dimension (nx) :: tmp,scr
      real :: tensor_pscalar_diff
      integer :: iy,iz,i,j
      logical, save :: first=.true.
!
!  read H and Bunit arrays and keep them in memory
!
      if(first) then
        open(1,file=trim(directory)//'/bunit.dat',form='unformatted')
        print*,'read bunit.dat with dimension: ',nx,ny,nz,3
        read(1) bunit,hhh
        close(1)
        print*,'read bunit.dat; bunit(1,1,1,1)=',bunit(1,1,1,1)
      endif
!
!  tmp = (Bunit.G)^2 + H.G + Bi*Bj*Gij
!  for details, see tex/mhd/thcond/tensor_der.tex
!
      call dot_mn(bunit,glncc,scr)
      call dot_mn(hhh,glncc,tmp)
      tmp=tmp+scr**2
!
!  calculate Hessian matrix of lncc
!
      call g2ij(f,ilncc,g)
!
!  dot with bi*bj
!
      iy=m-m1+1
      iz=n-n1+1
      do j=1,3
      do i=1,3
        tmp=tmp+bunit(:,iy,iz,i)*bunit(:,iy,iz,j)*g(:,i,j)
      enddo
      enddo
!
!  and add result to the dlncc/dt equation
!
      df(l1:l2,m,n,ilncc)=df(l1:l2,m,n,ilncc)+tensor_pscalar_diff*tmp
!
      first=.false.
    endsubroutine tensor_diff
!***********************************************************************

endmodule Pscalar



