! $Id: pscalar.f90,v 1.3 2002-07-11 00:24:33 brandenb Exp $

!  This modules solves the passive scalar advection equation

module Pscalar

  use Cparam
  use Cdata

  implicit none

  integer :: ilncc=0
  character(len=labellen) :: initlncc='wave-y',initlncc2='zero'
  logical :: nopscalar=.false.

  ! input parameters
  real :: ampllncc=.1, widthlncc=.5
  real :: ampllncc2=0.,kx_lncc=1.,ky_lncc=1.,kz_lncc=1.

  namelist /pscalar_init_pars/ &
       initlncc,initlncc2,ampllncc,ampllncc2,kx_lncc,ky_lncc,kz_lncc

  ! run parameters
  real :: pscalar_diff=0.

  namelist /pscalar_run_pars/ &
       pscalar_diff,nopscalar,tensor_pscalar_diff

  ! other variables (needs to be consistent with reset list below)
  integer :: i_ccm=0,i_ccmax=0

  contains

!***********************************************************************
    subroutine register_lncc()
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
           "$Id: pscalar.f90,v 1.3 2002-07-11 00:24:33 brandenb Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('Register_lncc: nvar > mvar')
      endif
!
    endsubroutine register_lncc
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
      real, dimension (mx,my,mz)      :: xx,yy,zz
!
      select case(initlncc)
        case('gaussian-noise'); call gaunoise(ampllncc,f,ilncc,ilncc)
        case('wave-x'); call wave(ampllncc,f,ilncc,kx=kx_lncc)
        case('wave-y'); call wave(ampllncc,f,ilncc,ky=ky_lncc)
        case('wave-z'); call wave(ampllncc,f,ilncc,kz=kz_lncc)
        case default; call stop_it('init_lncc: bad initlncc='//trim(initlncc))
      endselect
!
!  superimpose something else
!
      select case(initlncc2)
        case('wave-x'); call wave(ampllncc2,f,ilncc,ky=5.)
      endselect
!
      if(ip==0) print*,xx,yy,zz !(prevent compiler warnings)
    endsubroutine init_lncc
!***********************************************************************
    subroutine dlncc_dt(f,df,uu,glnrho)
!
!  passive scalar evolution
!  calculate dc/dt=-uu.glncc + DD*[del2lncc + (glncc+glnrho).glncc]
!
!   6-jul-02/axel: coded
!
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3) :: uu,glncc,glnrho
      real, dimension (nx) :: cc,uglncc,diff_op,del2lncc
!
      intent(in)  :: f,uu,glnrho
      intent(out) :: df
!
!  gradient of passive scalar
!  allow for possibility to turn off passive scalar
!  without changing file size etc.
!
      if (nopscalar) then
        if (headtt.or.ldebug) print*,'not SOLVED: dlncc_dt'
      else
        if (headtt.or.ldebug) print*,'SOLVE dlncc_dt'
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
      endif
!
!  tensor diffusion (keep the isotropic one)
!  dlnc/dt = ... + tensor_pscalar_diff*[Tij*Gi*Gj+Tij*Gij],
!  where Gj=dlnc/dxj and Gij=d2lnc/(dxi*dxj) and 
!
        if (tensor_pscalar_diff/=0.) then
          if(headtt) print*,'dlncc_dt: pscalar_diff=',pscalar_diff
          call dot_mn(glncc+glnrho,glncc,diff_op)
          call del2(f,ilncc,del2lncc)
          diff_op=diff_op+del2lncc
          df(l1:l2,m,n,ilncc)=df(l1:l2,m,n,ilncc)+pscalar_diff*diff_op
        endif
      endif
!
!  diagnostics
!
      if (ldiagnos) then
        cc=exp(f(l1:l2,m,n,ilncc))
        if (i_ccm/=0) call sum_mn_name(cc,i_ccm)
        if (i_ccmax/=0) call max_mn_name(cc,i_ccmax)
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
      integer :: iname
      logical :: lreset
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        i_ccm=0; i_ccmax=0
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'ccm',i_ccm)
        call parse_name(iname,cname(iname),cform(iname),'ccmax',i_ccmax)
      enddo
!
!  write column where which magnetic variable is stored
!
      write(3,*) 'i_ccm=',i_ccm
      write(3,*) 'i_ccmax=',i_ccmax
      write(3,*) 'ilncc=',ilncc
!
    endsubroutine rprint_pscalar
!***********************************************************************

endmodule Pscalar
