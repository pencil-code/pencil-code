module Magnetic

  use Cparam

  implicit none

  real :: grav=1.
  integer :: ient

  contains

!***********************************************************************
    subroutine register_aa()
!
!  Initialise variables which should know that we solve for the vector
!  potential: iaa, etc; increase nvar accordingly
!
!  6-nov-01/wolf: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_aa called twice')
      first = .false.
!
      lmagnetic = .true.
!      iaa = nvar+1              ! indices to access aa
      iax = iaa
      iay = iaa+1
      iaz = iaa+2
      nvar = nvar+3             ! added 3 variables
!
      if ((ip<=8) .and. lroot) then
        print*, 'Register_aa:  nvar = ', nvar
        print*, 'iaa,iax,iay,iaz = ', iaa,iax,iay,iaz
      endif
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('Register_aa: nvar > mvar')
      endif
!
    endsubroutine register_aa
!***********************************************************************
    subroutine init_aa(f,init,ampl,xx,yy,zz)
!
!  initialise magnetic field; called from start.f90
!  7-nov-2001/wolf: coded
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: tmp,r,p,xx,yy,zz
      real :: ampl
      integer :: init
!
      f(:,:,:,iax:iaz) = 0.
!
    endsubroutine init_aa
!***********************************************************************
    subroutine daa_dt(f,df,uu,uij,divu,glnrho,gpprho,cs2)
!
!  calculate right hand side of entropy equation
!
!  7-nov-01/wolf: [to be] coded
!
!  WARNING: this is just a verbatim copy of dssdt; needs to be adapted to
! aa
!
!      use Mpicomm
      use Cdata
      use Slices
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3,3) :: uij,sij
      real, dimension (nx,3) :: uu,gss,glnrho,gpprho
      real, dimension (nx) :: ugss,thdiff,del2ss,divu,sij2,cs2,ss,lnrho,TT1
      real :: gamma1
      integer :: i,j
!
      call grad(f,ient,gss)
      call del2(f,ient,del2ss)
!
!  sound speed squared
!
      gamma1=gamma-1.
      ss=f(l1:l2,m,n,ient)
      lnrho=f(l1:l2,m,n,ilnrho)
      cs2=cs20*exp(gamma1*lnrho+gamma*ss)
!
!  pressure gradient term
!
      !gpprho=cs20*glnrho  !(in isothermal case)
      do j=1,3
        gpprho(:,j)=cs2*(glnrho(:,j)+gss(:,j))
      enddo
!
!  advection term
!
      ugss=uu(:,1)*gss(:,1)+uu(:,2)*gss(:,2)+uu(:,3)*gss(:,3)
!
!  calculate rate of strain tensor
!
      do j=1,3
        do i=1,3
          sij(:,i,j)=.5*(uij(:,i,j)+uij(:,j,i))
        enddo
        sij(:,j,j)=sij(:,j,j)-.333333*divu
      enddo
!
      sij2=0.
      do j=1,3
      do i=1,3
        sij2=sij2+sij(:,i,j)**2
      enddo
      enddo
!
!  this is a term for numerical purposes only
!
      thdiff=nu*del2ss
!
      TT1=gamma1/cs2
      df(l1:l2,m,n,ient)=df(l1:l2,m,n,ient)+TT1*(-ugss+2.*nu*sij2)+thdiff
!
!  add gravity
!
      !if (headt) print*,'add gravity'
      df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)-grav
!
    endsubroutine daa_dt
!***********************************************************************

endmodule Magnetic
