! $Id: equ.f90,v 1.79 2002-07-06 20:29:17 brandenb Exp $

module Equ

  use Cdata

  implicit none

  contains

!***********************************************************************
      subroutine collect_UUmax
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
!  the result is then present only on the root processor
!  reassemble using old names
!
      fmax_tmp(1)=UUmax
      call mpireduce_max(fmax_tmp,fmax,1)
      if(lroot) UUmax=fmax(1)
!
      endsubroutine collect_UUmax
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
!  (in current version, don't need itype=2 anymore)
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
    subroutine xyaverages
!
!  calculate xy-averages (as functions of z)
!   6-jun-02/axel: coded
!
      use Mpicomm
      use Cdata
      use Sub
!
      real, dimension (nz,nprocz,mnamez) :: fsumz
!
!  communicate over all processors
!  the result is only present on the root processor
!
      if(nnamez>0) then
        call mpireduce_sum(fnamez,fsumz,nnamez*nz*nprocz)
        if(lroot) fnamez=fsumz/(nx*ny*nprocy)
      endif
!
    endsubroutine xyaverages
!***********************************************************************
    subroutine zaverages
!
!  calculate z-averages (as functions of x and y)
!  19-jun-02/axel: coded
!
      use Mpicomm
      use Cdata
      use Sub
!
      real, dimension (nx,ny,nprocy,mnamexy) :: fsumxy
!
!  communicate over all processors
!  the result is only present on the root processor
!
      if (nnamexy>0) then
        call mpireduce_sum(fnamexy,fsumxy,nnamexy*nx*ny*nprocy)
        if(lroot) fnamexy=fsumxy/(nz*nprocz)
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
      use Pscalar
      use Boundcond
      use IO
      use Shear
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3) :: uu,glnrho
      real, dimension (nx) :: lnrho,divu,u2,rho,ee=0.,rho1
      real :: fac
      integer :: j
!
!  print statements when they are first executed
!
      headtt = headt .and. lfirst .and. lroot

      if (headtt.or.ldebug) print*,'ENTER: pde'
      if (headtt) call cvs_id( &
           "$Id: equ.f90,v 1.79 2002-07-06 20:29:17 brandenb Exp $")
!
!  initialize counter for calculating and communicating print results
!
      ldiagnos=lfirst.and.lout
      if (ldiagnos) tdiagnos=t !(diagnostics are for THIS time)
!
!  initiate communication and do boundary conditions
!  need to call boundconds, because it also deals with x-boundaries!
!
      if (ldebug) print*,'bef. initiate_isendrcv_bdry'
      call initiate_isendrcv_bdry(f)
      call boundconds(f)
!
!  do loop over y and z
!  set indices and check whether communication must now be completed
!
      lfirstpoint=.true.        ! true for very first m-n loop
      do imn=1,ny*nz
        n=nn(imn)
        m=mm(imn)
        if (necessary(imn)) then 
           call finalise_isendrcv_bdry(f)
        endif
!
!  coordinates are needed all the time
!  (but not for isotropic turbulence!)
!  There are many other ciscumstances where this is not needed.
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
!  calculate inverse density: currenly only needed
!  in magnetic (provided lhydro=.true.) and for the
!  viscosity (currently in ldensity)
!  Better check that this is still true (after modifications)
!  Otherwise, better always calculate rho1.
!
!        if(ivisc==1.or.(lmagnetic.and.lhydro)) rho1=exp(-f(l1:l2,m,n,ilnrho))
!
!  WD: Also needed with heat conduction, so we better calculate it in all
!  cases. Could alternatively have a switch lrho1known and check for it,
!  or initialise to 1e35.
        if (ldensity) rho1=exp(-f(l1:l2,m,n,ilnrho))
!
!  hydro, density, and entropy evolution
!  They all are needed for setting some variables even
!  if their evolution is turned off.
!
        call duu_dt   (f,df,uu,glnrho,divu,rho1,u2)
        call dlnrho_dt(f,df,uu,glnrho,divu,lnrho)
        call dss_dt   (f,df,uu,glnrho,rho1,lnrho,cs2,TT1)
        call dlncc_dt (f,df,uu,glnrho)
!
!  Add gravity, if present
!  Shouldn't we call this one in hydro itself?
!  WD: there is some virtue in calling all of the dXX_dt in equ.f90
!  AB: but it is not really a new dXX_dt, because XX=uu.
!
        if (lhydro) then
          if(lgrav) call duu_dt_grav(f,df)
        endif
!
!  Magnetic field evolution
!
        if (lmagnetic) call daa_dt(f,df,uu,rho1,TT1)
!
!  Add shear if precent
!
        if (lshear) call shearing(f,df)
!
!  write slices for animation
!  this needs to be put somewhere else for compactness
!
        if (lhydro) then
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
        endif  !(end from lhydro)
!
!  In max_mn maximum values of u^2 (etc) are determined sucessively
!  va2 is set in magnetic (or nomagnetic)
!  In rms_mn sum of all u^2 (etc) is accumulated
!  Calculate maximum advection speed for timestep; needs to be done at
!  the first substep of each time step
!  Note that we are (currently) accumulating the maximum value,
!  not the maximum squared!
!
        if (lfirst.and.ldt) then
          fac=cdt/(cdtv*dxmin)
          call max_mn(sqrt(maxadvec2)+(fac*maxdiffus),UUmax)
        endif
!
!  calculate density diagnostics: mean density
!  Note that p/rho = gamma1*e = cs2/gamma, so e = cs2/(gamma1*gamma).
!
        if (ldiagnos.and.lhydro.and.ldensity) then
          if (gamma1/=0.) ee=cs2/(gamma*gamma1)
          rho=exp(f(l1:l2,m,n,ilnrho))
          if (i_eth/=0)  call sum_mn_name(rho*ee,i_eth)
          if (i_ekin/=0) call sum_mn_name(.5*rho*u2,i_ekin)
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
!  collect from different processors UUmax for the time step
!
      if (lfirst.and.ldt) call collect_UUmax
      if (ldiagnos) then
        call diagnostic
        call xyaverages
        call zaverages
      endif
!
    endsubroutine pde
!***********************************************************************

endmodule Equ
