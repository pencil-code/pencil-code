! $Id: pscalar.f90,v 1.23 2003-05-13 17:24:04 pkapyla Exp $

!  This modules solves the passive scalar advection equation

module Pscalar

  use Cparam
  use Cdata

  implicit none

  integer :: ilncc=0
  character (len=labellen) :: initlncc='zero', initlncc2='zero'
  character (len=40) :: tensor_pscalar_file
  logical :: nopscalar=.false.

  ! input parameters
  real :: ampllncc=.1, widthlncc=.5, cc_min=0., lncc_min
  real :: ampllncc2=0.,kx_lncc=1.,ky_lncc=1.,kz_lncc=1.,radius_lncc=0.,epsilon_lncc=0.
  real :: gradC0=0.

  namelist /pscalar_init_pars/ &
       initlncc,initlncc2,ampllncc,ampllncc2,kx_lncc,ky_lncc,kz_lncc, &
       radius_lncc,epsilon_lncc,widthlncc,cc_min

  ! run parameters
  real :: pscalar_diff=0.,tensor_pscalar_diff=0.

  namelist /pscalar_run_pars/ &
       pscalar_diff,nopscalar,tensor_pscalar_diff,gradC0

  ! other variables (needs to be consistent with reset list below)
  integer :: i_rhoccm=0,i_ccmax=0,i_lnccm=0,i_lnccmz=0
  integer :: i_ucm=0,i_uudcm=0

  contains

!***********************************************************************
    subroutine register_pscalar()
!
!  Initialise variables which should know that we solve for passive
!  scalar: ilncc; increase nvar accordingly
!
!  6-jul-02/axel: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_lncc called twice')
      first = .false.
!
      lpscalar = .true.
      ilncc = nvar+1            ! index to access lncc
      nvar = nvar+1             ! added 1 variable
!
      if ((ip<=8) .and. lroot) then
        print*, 'Register_lncc:  nvar = ', nvar
        print*, 'ilncc = ', ilncc
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: pscalar.f90,v 1.23 2003-05-13 17:24:04 pkapyla Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('Register_lncc: nvar > mvar')
      endif
!
    endsubroutine register_pscalar
!***********************************************************************
    subroutine initialize_pscalar()
!
!  Perform any necessary post-parameter read initialization
! 
!
!  24-nov-02/tony: coded
!
      ! dummy
    endsubroutine initialize_pscalar
!***********************************************************************
    subroutine init_lncc(f,xx,yy,zz)
!
!  initialise passive scalar field; called from start.f90
!
!   6-jul-2001/axel: coded
!
      use Cdata
      use Mpicomm
      use Density
      use Sub
      use Initcond
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz)      :: xx,yy,zz,prof
!
      select case(initlncc)
        case('zero'); f(:,:,:,ilncc)=0.
        case('gaussian-x'); call gaussian(ampllncc,f,ilncc,kx=kx_lncc)
        case('gaussian-y'); call gaussian(ampllncc,f,ilncc,ky=ky_lncc)
        case('gaussian-z'); call gaussian(ampllncc,f,ilncc,kz=kz_lncc)
        case('parabola-x'); call parabola(ampllncc,f,ilncc,kx=kx_lncc)
        case('parabola-y'); call parabola(ampllncc,f,ilncc,ky=ky_lncc)
        case('parabola-z'); call parabola(ampllncc,f,ilncc,kz=kz_lncc)
        case('gaussian-noise'); call gaunoise(ampllncc,f,ilncc,ilncc)
        case('wave-x'); call wave(ampllncc,f,ilncc,kx=kx_lncc)
        case('wave-y'); call wave(ampllncc,f,ilncc,ky=ky_lncc)
        case('wave-z'); call wave(ampllncc,f,ilncc,kz=kz_lncc)
        case('propto-ux'); call wave_uu(ampllncc,f,ilncc,kx=kx_lncc)
        case('propto-uy'); call wave_uu(ampllncc,f,ilncc,ky=ky_lncc)
        case('propto-uz'); call wave_uu(ampllncc,f,ilncc,kz=kz_lncc)
        case('tang-discont-z')
           print*,'init_lncc: widthlncc=',widthlncc
        prof=.5*(1.+tanh(zz/widthlncc))
        f(:,:,:,ilncc)=-1.+2.*prof
        case('hor-tube'); call htube2(ampllncc,f,ilncc,ilncc,xx,yy,zz,radius_lncc,epsilon_lncc)
        case default; call stop_it('init_lncc: bad initlncc='//trim(initlncc))
      endselect
!
!  superimpose something else
!
      select case(initlncc2)
        case('wave-x'); call wave(ampllncc2,f,ilncc,ky=5.)
      endselect
!
!  add floor value if cc_min is set
!
      if(cc_min/=0.) then
        lncc_min=alog(cc_min)
        if(lroot) print*,'set floor value for cc; cc_min=',cc_min
        f(:,:,:,ilncc)=amax1(lncc_min,f(:,:,:,ilncc))
      endif
!
      if(ip==0) print*,xx,yy,zz !(prevent compiler warnings)
    endsubroutine init_lncc
!***********************************************************************
    subroutine dlncc_dt(f,df,uu,glnrho)
!
!  passive scalar evolution
!  calculate dc/dt=-uu.glncc + pscaler_diff*[del2lncc + (glncc+glnrho).glncc]
!
!   6-jul-02/axel: coded
!
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3) :: uu,glncc,glnrho
      real, dimension (nx) :: lncc,cc,rho,uglncc,diff_op,del2lncc
!
      intent(in)  :: f,uu,glnrho
      intent(out) :: df
!
!  identify module and boundary conditions
!
      if (nopscalar) then
        if (headtt.or.ldebug) print*,'not SOLVED: dlncc_dt'
      else
        if (headtt.or.ldebug) print*,'SOLVE dlncc_dt'
      endif
      if (headtt) call identify_bcs('cc',ilncc)
!
!  gradient of passive scalar
!  allow for possibility to turn off passive scalar
!  without changing file size and recompiling everything.
!
      if (.not. nopscalar) then ! i.e. if (pscalar)
        call grad(f,ilncc,glncc)
        call dot_mn(uu,glncc,uglncc)
!
!  passive scalar equation
!
        if(lhydro) df(l1:l2,m,n,ilncc)=df(l1:l2,m,n,ilncc)-uglncc
!
!  diffusion operator
!
        if (pscalar_diff/=0.) then
          if(headtt) print*,'dlncc_dt: pscalar_diff=',pscalar_diff
          call dot_mn(glncc+glnrho,glncc,diff_op)
          call del2(f,ilncc,del2lncc)
          diff_op=diff_op+del2lncc
          df(l1:l2,m,n,ilncc)=df(l1:l2,m,n,ilncc)+pscalar_diff*diff_op
        endif
!
!  add diffusion of imposed constant gradient of c
!  restrict ourselves (for the time being) to z-gradient only
!  makes sense really only for periodic boundary conditions
!
        if (gradC0/=0.) then
          cc=exp(lncc)
          df(l1:l2,m,n,ilncc)=df(l1:l2,m,n,ilncc)-gradC0*uu(:,3)/cc
        endif
!
!  tensor diffusion (but keep the isotropic one)
!
        if (tensor_pscalar_diff/=0.) call tensor_diff(f,df,tensor_pscalar_diff,glncc)
!
      endif
!
!  diagnostics
!
!  output for double and triple correlators (assume z-gradient of cc)
!  <u_k u_j d_j c> = <u_k c uu.gradlncc>
!
      if (ldiagnos) then
        lncc=f(l1:l2,m,n,ilncc)
        cc=exp(lncc)
        rho=exp(f(l1:l2,m,n,ilnrho))
        if (i_rhoccm/=0) call sum_mn_name(rho*cc,i_rhoccm)
        if (i_ccmax/=0) call max_mn_name(cc,i_ccmax)
        if (i_lnccmz/=0) call xysum_mn_name_z(lncc,i_lnccmz)
        if (i_ucm/=0) call sum_mn_name(uu(:,3)*cc,i_ucm)
        if (i_uudcm/=0) call sum_mn_name(uu(:,3)*cc*uglncc,i_uudcm)
      endif
!
    endsubroutine dlncc_dt
!***********************************************************************
    subroutine rprint_pscalar(lreset)
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
        i_rhoccm=0; i_ccmax=0; i_lnccm=0; i_lnccmz=0
        i_ucm=0; i_uudcm=0
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'rhoccm',i_rhoccm)
        call parse_name(iname,cname(iname),cform(iname),'ccmax',i_ccmax)
        call parse_name(iname,cname(iname),cform(iname),'lnccm',i_lnccm)
        call parse_name(iname,cname(iname),cform(iname),'ucm',i_ucm)
        call parse_name(iname,cname(iname),cform(iname),'uudcm',i_uudcm)
      enddo
!
!  check for those quantities for which we want xy-averages
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'lnccmz',i_lnccmz)
      enddo
!
!  write column where which magnetic variable is stored
!
      write(3,*) 'i_rhoccm=',i_rhoccm
      write(3,*) 'i_ccmax=',i_ccmax
      write(3,*) 'i_lnccm=',i_lnccm
      write(3,*) 'i_ucm=',i_ucm
      write(3,*) 'i_uudcm=',i_uudcm
      write(3,*) 'i_lnccmz=',i_lnccmz
      write(3,*) 'ilncc=',ilncc
!
    endsubroutine rprint_pscalar
!***********************************************************************
    subroutine calc_mpscalar
!
!  calculate mean magnetic field from xy- or z-averages
!
!  14-apr-03/axel: adaped from calc_mfield
!
      use Cdata
      use Sub
!
      logical,save :: first=.true.
      real :: lnccm
!
!  Magnetic energy in horizontally averaged field
!  The bxmz and bymz must have been calculated,
!  so they are present on the root processor.
!
      if (i_lnccm/=0) then
        if (i_lnccmz==0) then
          if(first) print*
          if(first) print*,"NOTE: to get lnccm, lnccmz must also be set in xyaver"
          if(first) print*,"      We proceed, but you'll get lnccm=0"
          lnccm=0.
        else
          lnccm=sqrt(sum(fnamez(:,:,i_lnccmz)**2)/(nz*nprocz))
        endif
        call save_name(lnccm,i_lnccm)
      endif
!
    endsubroutine calc_mpscalar
!***********************************************************************
    subroutine tensor_diff(f,df,tensor_pscalar_diff,glncc)
!
!  reads file
!
!  11-jul-02/axel: coded
!
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, save, dimension (nx,ny,nz,3) :: bunit,hhh
      real, dimension (nx,3,3) :: g
      real, dimension (nx,3) :: glncc
      real, dimension (nx) :: tmp,scr
      real :: tensor_pscalar_diff
      integer :: mm,nn,i,j
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
      mm=m-m1+1
      nn=n-n1+1
      do j=1,3
      do i=1,3
        tmp=tmp+bunit(:,mm,nn,i)*bunit(:,mm,nn,j)*g(:,i,j)
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



