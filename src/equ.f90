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
      use Sub, only: wparam2,parse_bc
      use Forcing
      use Entropy
      use Magnetic
!
      logical, optional :: print
!
      open(1,file='run.in',form='formatted')
      read(1,*) nt,it1,dt,isave,itorder
      read(1,*) dsnap,dvid,dforce,dtmin
      read(1,*) tinit,tdamp,dampu
      read(1,*) dampuext,rdamp,wdamp
      read(1,*) ip,ix,iy,iz
      read(1,*) cs0,nu,ivisc,cdtv
      read(1,*) cdiffrho
      read(1,*) gravz
      read(1,*) bcx
      read(1,*) bcy
      read(1,*) bcz
      read(1,*) form1
      !
      if (lforcing ) read(1,NML=forcing_run_pars )
      if (lentropy ) read(1,NML=entropy_run_pars )
      if (lmagnetic) read(1,NML=magnetic_run_pars)
      close(1)
      cs20=cs0**2 !(goes into cdata module)
      ss0 = (alog(cs20) - gamma1*alog(rho0)-alog(gamma))/gamma   !!AB: this looks like it belongs to entropy
!
!  Write data to file for IDL
!
      call wparam2()
!
!  Give online feedback if called with the PRINT optional argument
!  Note: Some compiler's [like Compaq's] code crashes with the more
!  compact `if (present(print) .and. print)' 
!
      if (present(print)) then
        if (print) then
          call cprint()
        endif
      endif
!  
!  make sure ix,iy,iz are not outside the boundaries
!
      ix=min(ix,mx); iy=min(iy,my); iz=min(iz,mz)
      ix=max(ix,0); iy=max(iy,0); iz=max(iz,0)
!
!  parse boundary conditions; compound conditions of the form `a:s' allow
!  to have different variables at the lower and upper boundaries
!
      call parse_bc(bcx,bcx1,bcx2)
      call parse_bc(bcy,bcy1,bcy2)
      call parse_bc(bcz,bcz1,bcz2)
      if (lroot) then
        print*, 'bcx1,bcx2= ', bcx1," : ",bcx2
        print*, 'bcy1,bcy2= ', bcy1," : ",bcy2
        print*, 'bcz1,bcz2= ', bcz1," : ",bcz2
      endif
!
!  timestep
!
      ldt=dt.lt.0.
      if (ldt) cdt=abs(dt)
!
    endsubroutine cread
!***********************************************************************
    subroutine cprint
!
!  print input parameters
!  14-sep-01/axel: inserted from run.f90
!
      use Cdata
      use Forcing
      use Entropy
      use Magnetic
!
      if (lroot) then
        print*, 'nt,it1,dt,isave,itorder=', nt,it1,dt,isave,itorder
        print*, 'dsnap,dvid,dforce,dtmin=', dsnap,dvid,dforce,dtmin
        print*, 'tinit,tdamp,dampu=', tinit,tdamp,dampu
        print*, 'dampuext,rdamp,wdamp=', dampuext,rdamp,wdamp
        print*, 'ip,ix,iy,iz=', ip,ix,iy,iz
        print*, 'cs0,nu,ivisc,cdtv=', cs0,nu,ivisc,cdtv
        print*, 'cdiffrho=', cdiffrho
        print*, 'gravz=', gravz
        print*, 'bcx=', bcx
        print*, 'bcy=', bcy
        print*, 'bcz=', bcz
        print*, 'form1=', form1
        print*, 'cs20=', cs20
        !
        if (lforcing ) write(*,NML=forcing_run_pars )
        if (lentropy ) write(*,NML=entropy_run_pars )
        if (lmagnetic) write(*,NML=magnetic_run_pars)
      endif
!
    endsubroutine cprint
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
      real, dimension (nx) :: divu,uglnrho,u2,o2,ou
      real, dimension(nx) :: rho1,nu_var,chi,diff,del2lam
      real, dimension(nx) :: pdamp
      real :: diffrho
      integer :: i,j
!
!  print statements when they are first executed
!
      headtt = headt .and. lfirst .and. lroot

      if (headtt) call cvs_id( &
           "$RCSfile: equ.f90,v $", &
           "$Revision: 1.47 $", &
           "$Date: 2002-05-30 16:31:35 $")
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
        rho1=exp(-f(l1:l2,m,n,ilnrho))
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
          divu_slice(:,m-m1+1)=divu
          do j=1,3
            uu_slice(:,m-m1+1,j)=uu(:,j)
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
          diff = nu
          if (lentropy)  diff = max(diff, chi)
          call max_mn(u2+cs2+va2+cdtv*diff/dxmin,UUmax)
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
