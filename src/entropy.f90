module Entropy

  use Cparam

  implicit none

! integer, parameter :: mvar=nvar+1
  integer :: iss=5
  real :: grav=1.

  contains

!***********************************************************************
    subroutine dssdt(f,df,uu,uij,divu,glnrho,gpprho,cs2)
!
!  calculate right hand side of entropy equation
!
!  17-sep-01/axel: coded
!
      use Mpicomm
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
      call grad(f,iss,gss)
      call del2(f,iss,del2ss)
!
!  sound speed squared
!
      gamma1=gamma-1.
      ss=f(l1:l2,m,n,iss)
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
      df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss)+TT1*(-ugss+2.*nu*sij2)+thdiff
!
!  add gravity
!
      !if (headt) print*,'add gravity'
      df(l1:l2,m,n,iu+2)=df(l1:l2,m,n,iu+2)-grav
!
    endsubroutine dssdt
!***********************************************************************

endmodule Entropy
