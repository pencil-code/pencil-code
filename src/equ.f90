module Equ

  implicit none

  contains

!***********************************************************************
    subroutine cread
!
!  read input parameters
!  14-sep-01/axel: inserted from run.f90
!
      use Cdata
!
      open(1,file='run.in',form='formatted')
      read(1,*) nt,it1,dt,isave,itorder
      read(1,*) dsnap,dvid,dforce
      read(1,*) ip,ix,iy,iz
      read(1,*) cs,nu,ivisc
      read(1,*) iforce,force,relhel
      read(1,*) ibc
      read(1,*) form1
      close(1)
      cs20=cs**2 !(goes into cdata module)
!  
!  make sure ix,iy,iz are not outside the boundaries
!
      ix=min(ix,mx); iy=min(iy,my); iz=min(iz,mz)
      ix=max(ix,0); iy=max(iy,0); iz=max(iz,0)
!
    endsubroutine cread
!***********************************************************************
    subroutine cprint
!
!  print input parameters
!  14-sep-01/axel: inserted from run.f90
!
      use Cdata
!
      if (lroot) then
        print*,'nt,it1,dt,isave,itorder=',nt,it1,dt,isave,itorder
        print*,'dsnap,dvid,dforce',dsnap,dvid,dforce
        print*,'ip,ix,iy,iz',ip,ix,iy,iz
        print*, cs,nu,ivisc
        print*, iforce,force,relhel
        print*, ibc
        print*, form1
        print*, 'cs20=',cs20
      endif
!
    endsubroutine cprint
!***********************************************************************
    subroutine calc_UUmax
!
!  calculate diagnostic quantities
!   2-sep-01/axel: coded
!
      use Mpicomm
      use Cdata
      use Sub
!
      real, dimension(1) :: fmax_tmp,fmax
!
!  communicate over all processors
!  the result is present only on the root processor
!  need to take sqare root; reassemble using old names
!
      fmax_tmp(1)=UUmax
      call mpireduce_max(fmax_tmp,fmax,1)
      if(lroot) UUmax=sqrt(fmax(1))
!
    endsubroutine calc_uumax
!***********************************************************************
    subroutine diagnostic
!
!  calculate diagnostic quantities
!   2-sep-01/axel: coded
!
      use Mpicomm
      use Cdata
      use Sub
!
      real, dimension (mreduce) :: fmax_tmp,fsum_tmp,fmax,fsum
!
      fmax_tmp(1)=umax
      fmax_tmp(2)=omax
      fmax_tmp(3)=oumax
      fmax_tmp(4)=divumax
      fmax_tmp(5)=rmax
!
      fsum_tmp(1)=urms
      fsum_tmp(2)=orms
      fsum_tmp(3)=ourms
      fsum_tmp(4)=divurms
      fsum_tmp(5)=rrms
!
!  communicate over all processors
!
      call mpireduce_max(fmax_tmp,fmax,mreduce)
      call mpireduce_sum(fsum_tmp,fsum,mreduce)
!
!  the result is present only on the root processor
!
      if(lroot) then
!
!   need to take sqare root
!
        fsum=fsum/(nw*ncpus)
!
!  reassemble using old names
!
        umax=sqrt(fmax(1))
        omax=sqrt(fmax(2))
        oumax=fmax(3)
        divumax=sqrt(fmax(4))
        rmax=fmax(5)
!
        urms=sqrt(fsum(1))
        orms=sqrt(fsum(2))
        ourms=fsum(3)
        divurms=sqrt(fsum(4))
        rrms=fsum(5)
!
      endif
!
    endsubroutine diagnostic
!***********************************************************************
    subroutine prints
!
!  calculate and print some useful quantities during the run
!  29-sep-97/axel: coded
!
      use Cdata
      use Sub
!
!  print max and rms values of u and divu (calculated in pde)
!  do this only if we are on the root processor
!
      if(lroot) then
        write(6,form1) it,t,urms,umax,orms,omax,ourms,oumax,rrms,rmax,divurms,divumax
        write(20,form1) it,t,urms,umax,orms,omax,ourms,oumax,rrms,rmax,divurms,divumax
        open(3,file='check')
        write(3,'(a)')'it,t,urms,umax,orms,omax,ourms,oumax,rrms,rmax,divurms,divumax'
        write(3,form1) it,t,urms,umax,orms,omax,ourms,oumax,rrms,rmax,divurms,divumax
        close(3)
      endif
!
    endsubroutine prints
!***********************************************************************
    subroutine wvid(chdir)
!
!  write into video file
!  20-oct-97/axel: coded
!
      use Cdata
      use Slices
      use Sub
!
      real, save :: tvid
      integer, save :: ifirst,nvid
!
      character ch*4,file*9,chdir*(*)
      logical lvid
!
!  Output vid-data in 'tvid' time intervals
!
      file='tvid.dat'
      if (ifirst==0) then
        open(41,file=chdir//'/divu.dat',form='unformatted')
        open(43,file=chdir//'/ux.dat',form='unformatted')
        open(44,file=chdir//'/uz.dat',form='unformatted')
        call out1 (file,tvid,nvid,dvid,t)
        ifirst=1
      else
        call out2 (file,tvid,nvid,dvid,t,lvid,ch,.false.)
        if (lvid) then
          write(41) divu_slice(:,:),t
          write(43) uu_slice(:,:,1),t
          write(44) uu_slice(:,:,3),t
        endif
      endif
    endsubroutine wvid
!***********************************************************************
    subroutine pde(f,df)
!
!  Hydro equation
!  du/dt = -u.gradu -cs2*gradlam +nu*del2(u) ;(if graddivu term ignored)
!  dlnrho/dt = -u.grad(lnrho) - divu
!
!  10-sep-01/axel: coded
!
      use Mpicomm
      use Cdata
      use Slices
      use Sub
      use Entropy
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3,3) :: uij
      real, dimension (nx,3) :: uu,del2u,glnrho,ugu,oo,graddivu,fvisc,gpprho
      real, dimension (nx) :: divu,uglnrho,u2,o2,ou,divu2,rho,rho1,nurho1,cs2
      logical :: headtt,lfirstpoint
      integer :: i,j
!
!  print statements when they are first executed
!
      headtt = headt .and. lfirst .and. lroot
      if (headtt) print*,'$Id: equ.f90,v 1.2 2001-11-06 20:08:03 dobler Exp $'
!
!  initiate communication
!
      call initiate_isendrcv_bdry(f)
!
!  do loop over y and z
!  set indices and check whether communication must now be completed
!
      lfirstpoint=.true.
      do imn=1,ny*nz
        n=nn(imn)
        m=mm(imn)
        if (necessary(imn)) call finalise_isendrcv_bdry(f)
!
!  do all the neccessary derivatives here
!
        call gij(f,iuu,uij)
        call grad(f,ilnrho,glnrho)
!
!  viscosity operator
!
        if (ivisc==1) then
          if (headtt) print*,'full viscous force'
          rho1=exp(-f(l1:l2,m,n,ilnrho))
          nurho1=nu*rho1
          call del2v_etc(f,iuu,del2u,graddiv=graddivu)
          do i=1,3
            fvisc(:,i)=nurho1*(del2u(:,i)+1./3*graddivu(:,i))
          enddo
        else
          if (headtt) print*,'reduced viscous force'
          call del2v(f,iuu,del2u)
          fvisc=nu*del2u
        endif
!
!  abbreviations
!
        uu=f(l1:l2,m,n,iux:iuz)
!
!  auxiliary terms
!
        u2=uu(:,1)**2+uu(:,2)**2+uu(:,3)**2
        divu=uij(:,1,1)+uij(:,2,2)+uij(:,3,3)
        uglnrho=uu(:,1)*glnrho(:,1)+uu(:,2)*glnrho(:,2)+uu(:,3)*glnrho(:,3)
        do i=1,3
          ugu(:,i)=uu(:,1)*uij(:,i,1)+uu(:,2)*uij(:,i,2)+uu(:,3)*uij(:,i,3)
        enddo
!
!  entropy equation
!
        call dss_dt(f,df,uu,uij,divu,glnrho,gpprho,cs2)
!
!  momentum equation (forcing is now done in timestep)
!
        df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)-ugu-gpprho+fvisc
!
!  continuity equation
!
        df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho)-uglnrho-divu
!
!  write slices for animation
!
        if (n.eq.iz) then
          divu_slice(:,m-m1+1)=divu
          do j=1,3
            uu_slice(:,m-m1+1,j)=uu(:,j)
          enddo
        endif
!
!  In max_mn maximum values of u^2 (etc) are determined sucessively
!  In rms_mn sum of all u^2 (etc) is accumulated
!  Calculate maximum advection speed for timestep; needs to be done at every step
!
        if (lfirst.and.ldt) call max_mn(u2+cs2,UUmax)
!
!  Calculate maxima and rms values for diagnostic purposes
!
        if (lfirst.and.lout) then
          oo(:,1)=uij(:,3,2)-uij(:,2,3)
          oo(:,2)=uij(:,1,3)-uij(:,3,1)
          oo(:,3)=uij(:,2,1)-uij(:,1,2)
          o2=oo(:,1)**2+oo(:,2)**2+oo(:,3)**2
          ou=oo(:,1)*uu(:,1)+oo(:,2)*uu(:,2)+oo(:,3)*uu(:,3)
          divu2=divu**2
          rho=exp(f(l1:l2,m,n,ilnrho))
          call max_mn(u2,umax)
          call rms_mn(u2,urms)
          call max_mn(o2,omax)
          call rms_mn(o2,orms)
          call max_mn(ou,oumax)
          call rms_mn(ou,ourms)
          call max_mn(divu2,divumax)
          call rms_mn(divu2,divurms)
          call max_mn(rho,rmax)
          call rms_mn(rho,rrms)
        endif
!
!  end of loops over m and n
!
        headtt=.false.
        lfirstpoint=.false.
      enddo
!
!  diagnostic quantities
!  calculate maximum speed for time step
!
      if (lfirst.and.ldt) call calc_uumax
      if (lfirst.and.lout) call diagnostic
!
    endsubroutine pde
!***********************************************************************

endmodule Equ
