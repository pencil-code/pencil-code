! $Id: equ.f90,v 1.50 2002-06-01 09:36:38 brandenb Exp $

module Equ

  use Cdata

  implicit none

  contains

!***********************************************************************
      subroutine calc_UUmax
!
!  This routine is used for calculating the maximum effective advection
!  velocity in the domain for determining dt at each timestep 
!  
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
      endsubroutine calc_UUmax
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
      integer :: iname,imax_count,isum_count,nmax_count,nsum_count
      real, dimension (mname) :: fmax_tmp,fsum_tmp,fmax,fsum
!
!  go through all print names, and sort into communicators
!  corresponding to their type
!
      imax_count=0
      isum_count=0
      do iname=1,nname
        if(itype_name(iname)<0) then
          imax_count=imax_count+1
          fmax_tmp(imax_count)=fname(iname)
        elseif(itype_name(iname)>0) then
          isum_count=isum_count+1
          fsum_tmp(isum_count)=fname(iname)
        endif
      enddo
      nmax_count=imax_count
      nsum_count=isum_count
!
!  communicate over all processors
!
      call mpireduce_max(fmax_tmp,fmax,nmax_count)
      call mpireduce_sum(fsum_tmp,fsum,nsum_count)
!
!  the result is present only on the root processor
!
      if(lroot) then
        fsum=fsum/(nw*ncpus)
!
!  sort back into original array
!  need to take sqare root if |itype|=2
!  THIS COULD BE SIMPLIFIED!!
!
      imax_count=0
      isum_count=0
      do iname=1,nname
        if(itype_name(iname)<0) then
          imax_count=imax_count+1
          if(itype_name(iname)==-1) fname(iname)=fmax(imax_count)
          if(itype_name(iname)==-2) fname(iname)=sqrt(fmax(imax_count))
        elseif(itype_name(iname)>0) then
          isum_count=isum_count+1
          if(itype_name(iname)==+1) fname(iname)=fsum(isum_count)
          if(itype_name(iname)==+2) fname(iname)=sqrt(fsum(isum_count))
        endif
      enddo
      !nmax_count=imax_count
      !nsum_count=isum_count
!
      endif
!
    endsubroutine diagnostic
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
      character (len=4) :: ch
      character (len=9) :: file
      character (len=*) :: chdir
      logical lvid
!
!  Output vid-data in 'tvid' time intervals
!
      file='tvid.dat'
      if (ifirst==0) then
        open(41,file=chdir//'/divu.xy',form='unformatted')
        open(43,file=chdir//'/ux.xy',form='unformatted')
        open(44,file=chdir//'/uz.xy',form='unformatted')
        open(45,file=chdir//'/uz.xz',form='unformatted')
        open(46,file=chdir//'/lnrho.xz',form='unformatted')
        call out1 (file,tvid,nvid,dvid,t)
        ifirst=1
      else
        call out2 (file,tvid,nvid,dvid,t,lvid,ch,.false.)
        if (lvid) then
          write(41) divu_xy(:,:),t
          write(43) uu_xy(:,:,1),t
          write(44) uu_xy(:,:,3),t
          write(45) uu_xz(:,:,3),t
          write(46) lnrho_xz(:,:),t
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
      use Cdata
      use Mpicomm
      use Slices
      use Sub
      use Global
      use Hydro
      use Gravity
      use Entropy
      use Magnetic
      use IO
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3,3) :: uij
      real, dimension (nx,3) :: uu,del2u,glnrho,ugu,oo,graddivu,fvisc,gpprho
      real, dimension (nx) :: divu,lnrho,uglnrho,u2,o2,ou
      real, dimension(nx) :: rho1,nu_var,chi,diff,del2lam
      real, dimension(nx) :: pdamp
      real :: diffrho,fac
      integer :: i,j
!
!  print statements when they are first executed
!
      headtt = headt .and. lfirst .and. lroot

      if (headtt) call cvs_id( &
           "$RCSfile: equ.f90,v $", &
           "$Revision: 1.50 $", &
           "$Date: 2002-06-01 09:36:38 $")
!
!  initialize counter for calculating and communicating print results
!
      ldiagnos= lfirst .and. lout
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
!  abbreviations
!
        uu=f(l1:l2,m,n,iux:iuz)
        lnrho=f(l1:l2,m,n,ilnrho)
!
!  do all the neccessary derivatives here
!
        call gij(f,iuu,uij)
        call grad(f,ilnrho,glnrho)
!
!  coordinates are needed all the time
!
        x_mn = x(l1:l2)
        y_mn = spread(y(m),1,nx)
        z_mn = spread(z(n),1,nx)
        r_mn = sqrt(x_mn**2+y_mn**2+z_mn**2)
!
!  rho1 (=1/rho) is needed for viscous term, heat conduction, and Lorentz force
!
        rho1=exp(-lnrho)
!
!  viscosity operator
!
        if (ivisc==1) then
          if (headtt) print*,'full viscous force'
          nu_var=nu*rho0*rho1   ! spatially varying nu
          call del2v_etc(f,iuu,del2u,GRADDIV=graddivu)
          do i=1,3
            fvisc(:,i)=nu_var*(del2u(:,i)+1./3*graddivu(:,i))
          enddo
        else
          if (headtt) print*,'reduced viscous force'
          call del2v(f,iuu,del2u)
          fvisc=nu*del2u
        endif
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
!  pure hydro part of eq. of motion (forcing is now done in timestep)
!
        df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) - ugu + fvisc
!
!  entropy equation
!  needs to be called every time even with noentropy,
!  because it is here that we set cs2, TT1, and gpprho
!
!       if (lentropy .or. headtt) then
          call dss_dt(f,df,uu,uij,divu,rho1,glnrho,gpprho,cs2,TT1,chi)
!       endif
        df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) - gpprho
        !!if (lentropy) df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) - gpprho
!
!  thermal part of eq. of motion (pressure force)
!
!
!  magnetic part
!
        if (lmagnetic) call daa_dt(f,df,uu,rho1,TT1)
!
!  damping terms (artificial, but sometimes useful):
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
!  2. damp motions for r_mn>1
!
        if (lgravr) then
!          pdamp = 0.5*(1+tanh((r_mn-rdamp)/wdamp)) ! damping profile
          pdamp = step(r_mn,rdamp,wdamp) ! damping profile
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
          divu_xy(:,m-m1+1)=divu
          lnrho_xy(:,m-m1+1)=lnrho
          do j=1,3
            uu_xy(:,m-m1+1,j)=uu(:,j)
          enddo
        endif
!
        if (m.eq.iy) then
          lnrho_xz(:,n-n1+1)=lnrho
          do j=1,3
            uu_xz(:,n-n1+1,j)=uu(:,j)
          enddo
        endif
!
!  In max_mn maximum values of u^2 (etc) are determined sucessively
!  va2 is set in magnetic (or nomagnetic)
!  In rms_mn sum of all u^2 (etc) is accumulated
!  Calculate maximum advection speed for timestep; needs to be done at
!  the first substep of each time step
!
        if (lfirst.and.ldt) then
          fac=cdt/(cdtv*dxmin)
          diff = nu  !!(for the time being)
          if (lentropy)  diff = max(diff, chi)
          call max_mn(sqrt(u2+cs2+va2)+fac*diff,UUmax)
        endif
!
!  Calculate maxima and rms values for diagnostic purposes
!  (The corresponding things for magnetic fields etc happen inside magnetic etc)
!  The length of the timestep is not known here (--> moved to prints.f90)
!
        if (ldiagnos) then
          tdiagnos = t !(diagnostics are for THIS time)
          oo(:,1)=uij(:,3,2)-uij(:,2,3)
          oo(:,2)=uij(:,1,3)-uij(:,3,1)
          oo(:,3)=uij(:,2,1)-uij(:,1,2)
          if (i_oum/=0) then
            call dot_mn(oo,uu,ou)
            call sum_mn_name(ou,i_oum)
          endif
          if (i_o2m/=0) then
            call dot2_mn(oo,o2)
            call sum_mn_name(o2,i_o2m)
          endif
          if (i_u2m/=0) call sum_mn_name(u2,i_u2m)
          if (i_um2/=0) call max_mn_name(u2,i_um2)
!         rho=exp(f(l1:l2,m,n,ilnrho))
!         call max_mn (u2,u2max)
!         call rms2_mn(u2,urms)
!         call max_mn (o2,o2max)
!         call rms2_mn(o2,orms)
!         call max_mn (ou,oumax)
!         call rms_mn (ou,ourms)
!         call mean_mn(rho,rmean)
!         call max_mn (rho,rmax)
!         call rms_mn (rho,rrms)
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
      if (lfirst.and.ldt) call calc_UUmax
      if (ldiagnos) then
        call diagnostic
        !call diagnostic_old !!(should still be intact)
      endif
!
    endsubroutine pde
!***********************************************************************

endmodule Equ
