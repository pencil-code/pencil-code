! $Id: equ.f90,v 1.56 2002-06-08 11:17:15 brandenb Exp $

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
    subroutine zaverages
!
!  calculate z-averages
!   6-jun-02/axel: coded
!
      use Mpicomm
      use Cdata
      use Sub
!
      real, dimension (nz,mnamez) :: fsumz
!
!  communicate over all processors
!
      call mpireduce_sum(fnamez,fsumz,nnamez*nz)
!
!  the result is present only on the root processor
!
      if(lroot) then
        fnamez=fsumz/(nx*ny*ncpus)
      endif
!
    endsubroutine zaverages
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
        open(42,file=chdir//'/ux.xy',form='unformatted')
        open(43,file=chdir//'/uz.xy',form='unformatted')
        open(44,file=chdir//'/ux.xz',form='unformatted')
        open(45,file=chdir//'/uz.xz',form='unformatted')
        open(46,file=chdir//'/lnrho.xz',form='unformatted')
        open(47,file=chdir//'/ss.xz',form='unformatted')
        call out1 (file,tvid,nvid,dvid,t)
        ifirst=1
      else
        call out2 (file,tvid,nvid,dvid,t,lvid,ch,.false.)
        if (lvid) then
          write(41) divu_xy(:,:),t
          write(42) uu_xy(:,:,1),t
          write(43) uu_xy(:,:,3),t
          write(44) uu_xz(:,:,1),t
          write(45) uu_xz(:,:,3),t
          write(46) lnrho_xz(:,:),t
          write(47) ss_xz(:,:),t
        endif
      endif
    endsubroutine wvid
!***********************************************************************
    subroutine pde(f,df)
!
!  call the different evolution equations (now all in their own modules)
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
      real, dimension (nx,3,3) :: uij,sij
      real, dimension (nx,3) :: uu,glnrho,oo,gpprho
      real, dimension (nx) :: lnrho,divu,u2,o2,ou,rho,ee
      real, dimension(nx) :: rho1
      real :: fac
      integer :: j
!
!  print statements when they are first executed
!
      headtt = headt .and. lfirst .and. lroot

      if (headtt) call cvs_id( &
           "$RCSfile: equ.f90,v $", &
           "$Revision: 1.56 $", &
           "$Date: 2002-06-08 11:17:15 $")
!
!  initialize counter for calculating and communicating print results
!
      ldiagnos= lfirst .and. lout
!
!  initiate communication
!
      call initiate_isendrcv_bdry(f)
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
!  coordinates are needed all the time
!  (but not for isotropic turbulence!)
!
        x_mn = x(l1:l2)
        y_mn = spread(y(m),1,nx)
        z_mn = spread(z(n),1,nx)
        r_mn = sqrt(x_mn**2+y_mn**2+z_mn**2)
!
!  for each pencil, accummulate through the different routines
!  maximum diffusion and maximum advection (keep as nx-array)
!
        maxdiffus=0.
        maxadvec2=0.
!
!  hydro and density parts
!
        if (lhydro)   call duu_dt   (f,df,uu,divu,sij,uij,u2)
        if (ldensity) call dlnrho_dt(f,df,uu,divu,sij,lnrho,glnrho,rho1)
!
!  entropy equation: needs to be called EVERY time EVEN with noentropy,
!  because it is here that we set cs2, TT1, and gpprho
!
        call dss_dt(f,df,uu,sij,lnrho,glnrho,gpprho,cs2,TT1)
        df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)-gpprho
!
!  magnetic part
!
        if (lmagnetic) call daa_dt(f,df,uu,rho1,TT1)
!
!  add gravity
!
        if (lgrav) call duu_dt_grav(f,df)
!
!  write slices for animation
!
        if (n.eq.iz) then
          divu_xy(:,m-m1+1)=divu
          if (ldensity) lnrho_xy(:,m-m1+1)=f(l1:l2,m,n,ilnrho)
          do j=1,3
            uu_xy(:,m-m1+1,j)=f(l1:l2,m,n,iuu+j-1)
          enddo
        endif
!
        if (m.eq.iy) then
          if (ldensity) lnrho_xz(:,n-n1+1)=f(l1:l2,m,n,ilnrho)
          if (lentropy) ss_xz(:,n-n1+1)=f(l1:l2,m,n,ient)
          do j=1,3
            uu_xz(:,n-n1+1,j)=f(l1:l2,m,n,iuu+j-1)
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
          call max_mn(sqrt(maxadvec2)+fac*maxdiffus,UUmax)
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
!
!  calculate density diagnostics: mean density
!
          ee=cs2/(gamma*gamma1)
          rho=exp(f(l1:l2,m,n,ilnrho))
          if (i_eth/=0) call max_mn_name(rho*ee,i_eth)
          if (i_ekin/=0) call max_mn_name(.5*rho*u2,i_ekin)
          if (i_rhom/=0) call sum_mn_name(rho,i_rhom)
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
        call zaverages
      endif
!
    endsubroutine pde
!***********************************************************************

endmodule Equ
