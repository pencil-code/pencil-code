module Equ

  implicit none

  contains

!***********************************************************************
    subroutine cread(print)
!
!  read input parameters
!  14-sep-01/axel: inserted from run.f90
!
      use Cdata
!
      logical, optional :: print
!
      open(1,file='run.in',form='formatted')
      read(1,*) nt,it1,dt,isave,itorder
      read(1,*) dsnap,dvid,dforce
      read(1,*) tinit,tdamp,dampu
      read(1,*) dampuext,rdamp,wdamp
      read(1,*) ip,ix,iy,iz
      read(1,*) cs0,nu,ivisc
      read(1,*) chi0,chi2
      read(1,*) cdiffrho
      read(1,*) gravz
      read(1,*) cheat,wheat,cool,wcool
      read(1,*) iforce,force,relhel
      read(1,*) bcx
      read(1,*) bcy
      read(1,*) bcz
      read(1,*) form1
      close(1)
      cs20=cs0**2 !(goes into cdata module)

      if (present(print) .and. print) then
        call cprint()
      endif
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
        print*, 'nt,it1,dt,isave,itorder=', nt,it1,dt,isave,itorder
        print*, 'dsnap,dvid,dforce=', dsnap,dvid,dforce
        print*, 'tinit,tdamp,dampu=', tinit,tdamp,dampu
        print*, 'dampuext,rdamp,wdamp=', dampuext,rdamp,wdamp
        print*, 'ip,ix,iy,iz=', ip,ix,iy,iz
        print*, 'cs0,nu,ivisc=', cs0,nu,ivisc
        print*, 'chi,chi2=', chi0,chi2
        print*, 'cdiffrho=', cdiffrho
        print*, 'gravz=', gravz
        print*, 'cheat,wheat,cool,wcool=', cheat,wheat,cool,wcool
        print*, 'iforce,force,relhel=', iforce,force,relhel
        print*, 'bcx=', bcx
        print*, 'bcy=', bcy
        print*, 'bcz=', bcz
        print*, 'form1=', form1
        print*, 'cs20=', cs20
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
      fmax_tmp(1)=u2max
      fmax_tmp(2)=o2max
      fmax_tmp(3)=oumax
      fmax_tmp(4)=divu2max
      fmax_tmp(5)=rmax
!
      fsum_tmp(1)=urms
      fsum_tmp(2)=orms
      fsum_tmp(3)=ourms
      fsum_tmp(4)=divurms
      fsum_tmp(5)=rmean
      fsum_tmp(6)=rrms
!
!  communicate over all processors
!
      call mpireduce_max(fmax_tmp,fmax,mreduce)
      call mpireduce_sum(fsum_tmp,fsum,mreduce)
!
!  the result is present only on the root processor
!
      if(lroot) then
        fsum=fsum/(nw*ncpus)
!
!  reassemble using old names
!  need to take sqare root
!
        umax=sqrt(fmax(1))
        omax=sqrt(fmax(2))
        oumax=fmax(3)
        divumax=sqrt(fmax(4))
        rmax=fmax(5)
!
        urms=sqrt(fsum(1))
        orms=sqrt(fsum(2))
        ourms=sqrt(fsum(3))
        divurms=sqrt(fsum(4))
        rmean=sqrt(fsum(5))
        rrms=sqrt(fsum(6))
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
        write( 6,form1) it-1,t_diag,urms,umax,orms,omax,ourms,oumax, &
             rmean,rrms,rmax,divurms,divumax
        write(20,form1) it-1,t_diag,urms,umax,orms,omax,ourms,oumax, &
             rmean,rrms,rmax,divurms,divumax
        open(3,file='check')
        write(3,'(a)')'it,t_diag,urms,umax,orms,omax,ourms,oumax,rmean,rrms,rmax,divurms,divumax'
        write( 3,form1) it-1,t_diag,urms,umax,orms,omax,ourms,oumax, &
             rmean,rrms,rmax,divurms,divumax
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
      character (LEN=4) :: ch
      character (LEN=9) :: file
      character (LEN=*) :: chdir
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
      use Global
      use Gravity
      use Entropy
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3,3) :: uij
      real, dimension (nx,3) :: uu,del2u,glnrho,ugu,oo,graddivu,fvisc,gpprho
      real, dimension (nx) :: divu,uglnrho,u2,o2,ou,divu2
      real, dimension(nx) :: rho,rho1,nurho1,cs2,del2lam
      real, dimension(nx) :: r,pdamp
      real :: diffrho
      integer :: i,j
!
!  print statements when they are first executed
!
      headtt = headt .and. lfirst .and. lroot
      if (headtt) call cvs_id( &
           "$RCSfile: equ.f90,v $", &
           "$Revision: 1.13 $", &
           "$Date: 2002-01-23 19:56:13 $")
!
!  initiate communication
!
      call initiate_isendrcv_bdry(f)
      diffrho = cdiffrho*dxmin*cs0
!
!  do loop over y and z
!  set indices and check whether communication must now be completed
!
      lfirstpoint=.true.        ! true for very first m-n loop
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
        if (lentropy) call dss_dt(f,df,uu,uij,divu,glnrho,gpprho,cs2)
!
!  momentum equation (forcing is now done in timestep)
!
        df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)-ugu-gpprho+fvisc
!
!  damping terms:
!
!  1. damp motion during time interval 0<t<tdamp.
!  damping coefficient is dampu (if >0) or |dampu|/dt (if dampu <0)
!
        if ((dampu .ne. 0.) .and. (t < tdamp)) then
          ! damp motion provided t<tdamp
          if (dampu > 0) then
            df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) &
                                    - dampu*f(l1:l2,m,n,iux:iuz)
          else
            if (dt > 0) then    ! dt known and good
              df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) &
                                      + dampu/dt*f(l1:l2,m,n,iux:iuz)
            endif
          endif
        endif
!
!  2. damp motions for r>1
!
        if (lgravr) then
!        r = rr(l1:l2,m,n)
          pdamp = 0.5*(1+tanh((r-rdamp)/wdamp)) ! damping profile
          do i=iux,iuz
            df(l1:l2,m,n,i) = df(l1:l2,m,n,i) - dampuext*pdamp*f(l1:l2,m,n,i)
          enddo
        endif
!
!  add gravity
!
        if (lgrav) call duu_dt_grav(f,df)
!
!  continuity equation
!
        df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho)-uglnrho-divu
!
!  mass diffusion
!
        if (cdiffrho /= 0.) then
          call del2(f,ilnrho,del2lam)
          df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + diffrho*del2lam
        endif
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
!  Calculate maximum advection speed for timestep; needs to be done at
! every step
!
        if (lfirst.and.ldt) call max_mn(u2+cs2,UUmax)
!
!  Calculate maxima and rms values for diagnostic purposes
!
        if (lfirst.and.lout) then
          t_diag = t            ! diagnostic quantities are for this time
          oo(:,1)=uij(:,3,2)-uij(:,2,3)
          oo(:,2)=uij(:,1,3)-uij(:,3,1)
          oo(:,3)=uij(:,2,1)-uij(:,1,2)
          o2=oo(:,1)**2+oo(:,2)**2+oo(:,3)**2
          ou=oo(:,1)*uu(:,1)+oo(:,2)*uu(:,2)+oo(:,3)*uu(:,3)
          divu2=divu**2
          rho=exp(f(l1:l2,m,n,ilnrho))
          call max_mn (u2,u2max)
          call rms2_mn(u2,urms)
          call max_mn (o2,o2max)
          call rms2_mn(o2,orms)
          call max_mn (ou,oumax)
          call rms_mn (ou,ourms)
          call max_mn (divu2,divu2max)
          call rms2_mn(divu2,divurms)
          call mean_mn (rho,rmean)
          call max_mn (rho,rmax)
          call rms_mn (rho,rrms)
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
