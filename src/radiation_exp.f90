! $Id: radiation_exp.f90,v 1.40 2003-07-02 15:19:43 theine Exp $

!!!  NOTE: this routine will perhaps be renamed to radiation_feautrier
!!!  or it may be combined with radiation_ray.

module Radiation

!  Radiation (solves transfer equation along rays)
!  The direction of the ray is given by the vector (lrad,mrad,nrad),
!  and the parameters radx0,rady0,radz0 gives the maximum number of
!  steps of the direction vector in the corresponding direction.

  use Cparam
!
  implicit none
!
  !integer, parameter :: radx0=3,rady0=3,radz0=3
  integer, parameter :: radx0=1,rady0=1,radz0=1
  real, dimension(radx0,my,mz) :: Irad_yz,Irad0_yz,tau_yz
  real, dimension(mx,rady0,mz) :: Irad_zx,Irad0_zx,tau_zx
  real, dimension(mx,my,radz0) :: Irad_xy,Irad0_xy,tau_xy
  integer, dimension(radx0,my,mz,3) :: pos_yz
  integer, dimension(mx,rady0,mz,3) :: pos_zx
  integer, dimension(mx,my,radz0,3) :: pos_xy
  real, dimension (mx,my,mz) :: Srad,kaprho
  integer :: directions
  real :: frac
!
!  default values for one pair of vertical rays
!
  integer :: radx=0,rady=0,radz=1,rad2max=1
!
  logical :: nocooling=.false.,test_radiation=.false.,output_Qrad=.false.
  logical :: lkappa_es=.false.
!
!  definition of dummy variables for FLD routine
!
  real :: DFF_new=0.  !(dum)
  integer :: i_frms=0,i_fmax=0,i_Erad_rms=0,i_Erad_max=0
  integer :: i_Egas_rms=0,i_Egas_max=0

  namelist /radiation_init_pars/ &
       radx,rady,radz,rad2max,output_Qrad,test_radiation,lkappa_es

  namelist /radiation_run_pars/ &
       radx,rady,radz,rad2max,output_Qrad,test_radiation,lkappa_es,nocooling

  contains

!***********************************************************************
    subroutine register_radiation()
!
!  initialise radiation flags
!
!  24-mar-03/axel+tobi: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if(.not. first) call stop_it('register_radiation called twice')
      first = .false.
!
      lradiation=.true.
      lradiation_ray=.true.
!
!  set indices for auxiliary variables
!
      iQrad = mvar + naux +1; naux = naux + 1
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_radiation: radiation naux = ', naux
        print*, 'iQrad = ', iQrad
      endif
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: radiation_exp.f90,v 1.40 2003-07-02 15:19:43 theine Exp $")
!
!  Check that we aren't registering too many auxilary variables
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'naux = ', naux, ', maux = ', maux
        call stop_it('register_radiation: naux > maux')
      endif
!
!  Writing files for use with IDL
!
      if (naux < maux) aux_var(aux_count)=',Qrad $'
      if (naux == maux) aux_var(aux_count)=',Qrad'
      aux_count=aux_count+1
      write(5,*) 'Qrad = fltarr(mx,my,mz)*one'
!
    endsubroutine register_radiation
!***********************************************************************
    subroutine initialize_radiation()
!
!  Calculate number of directions of rays
!  Do this in the beginning of each run
!
!  16-jun-03/axel+tobi: coded
!
  integer :: lrad,mrad,nrad,rad2
!
!  check that the number of rays does not exceed maximum
!
      if(radx>radx0) stop "radx0 is too small"
      if(rady>rady0) stop "rady0 is too small"
      if(radz>radz0) stop "radz0 is too small"
!
!  count
!
      directions=0
      do nrad=-radz,radz
      do mrad=-rady,rady
      do lrad=-radx,radx
        rad2=lrad**2+mrad**2+nrad**2
        if(rad2>0 .and. rad2<=rad2max) then 
          directions=directions+1
        endif
      enddo
      enddo
      enddo
      print*,'initialize_radiation: directions=',directions
!
    endsubroutine initialize_radiation
!***********************************************************************
    subroutine radcalc(f)
!
!  calculate source function and opacity
!
!  24-mar-03/axel+tobi: coded
!
      use Cdata
      use Ionization
!
      real, dimension(mx,my,mz,mvar+maux), intent(in) :: f
!
!  test
!
      if(test_radiation) then
        if(lroot) print*,'radcalc: put Srad=kaprho=1 (as a test)'
        Srad=1.
        kaprho=1.
        return
      endif
!
!  At the moment we don't calculate ghost zones (ok for vertical arrays)  
!
      do n=n0,n3
      do m=m0,m3
         call sourcefunction(f,Srad)
!
!  opacity: if lkappa_es then take electron scattering opacity only;
!  otherwise use Hminus opacity (but may need to add kappa_es as well).
!
         if (lkappa_es) then
            kaprho(l0:l3,m,n)=kappa_es*exp(f(l0:l3,m,n,ilnrho))
         else
            call opacity(f,kaprho)
         endif
      enddo
      enddo
!
    endsubroutine radcalc
!***********************************************************************
    subroutine radtransfer(f)
!
!  Integration radioation transfer equation along rays
!
!  This routine is called before the communication part
!  (certainly needs to be given a better name)
!  All rays start with zero intensity
!
!  16-jun-03/axel+tobi: coded
!
      use Cdata
      use Sub
!
      real, dimension(mx,my,mz,mvar+maux) :: f
      integer :: lrad,mrad,nrad,rad2
!
!  identifier
!
      if(lroot.and.headt) print*,'radtransfer'
!
!  calculate source function and opacity
!
      call radcalc(f)
!
!  Accumulate the result for Qrad=(J-S),
!  First initialize Qrad=-S. 
!
      f(:,:,:,iQrad)=-Srad
!
!  calculate weights
!
      frac=1./directions
!
!  loop over rays
!
      do nrad=-radz,radz
      do mrad=-rady,rady
      do lrad=-radx,radx
        rad2=lrad**2+mrad**2+nrad**2
        if (rad2>0 .and. rad2<=rad2max) then 
        call radtransfer_intrinsic  (f,lrad,mrad,nrad)
        call radtransfer_communicate(  lrad,mrad,nrad)
        call radtransfer_revision   (f,lrad,mrad,nrad)
        endif
      enddo
      enddo
      enddo
!
    endsubroutine radtransfer
!***********************************************************************
    subroutine radtransfer_intrinsic(f,lrad,mrad,nrad)
!
!  Integration radiation transfer equation along rays
!
!  This routine is called before the communication part
!  All rays start with zero intensity
!
!  16-jun-03/axel+tobi: coded
!
      use Cdata
      use Sub
!
      real, dimension(mx,my,mz,mvar+maux) :: f
      integer, intent(in) :: lrad,mrad,nrad
      real, dimension(mx,my,mz) :: tau,Irad
!
!  identifier
!
      if(lroot.and.headt) print*,'radtransfer_intrinsic'
!
!  calculate intrinsicinsic intensity and optical depth
!
      call intensity_intrinsic(lrad,mrad,nrad,tau,Irad)
!
!  add contribution to the heating rate Q
!
      f(:,:,:,iQrad)=f(:,:,:,iQrad)+frac*Irad
!
!  store boundary values for communication
!
     if (lrad<0) then
        tau_yz (:,:,:)=tau (l1:l1+radx0-1,:,:)
        Irad_yz(:,:,:)=Irad(l1:l1+radx0-1,:,:)
     endif
     if (lrad>0) then
        tau_yz (:,:,:)=tau (l2-radx0+1:l2,:,:)
        Irad_yz(:,:,:)=Irad(l2-radx0+1:l2,:,:)
     endif
     if (mrad<0) then
        tau_zx (:,:,:)=tau (:,m1:m1+rady0-1,:)
        Irad_zx(:,:,:)=Irad(:,m1:m1+rady0-1,:)
     endif
     if (mrad>0) then
        tau_zx (:,:,:)=tau (:,m2-rady0+1:m2,:)
        Irad_zx(:,:,:)=Irad(:,m2-rady0+1:m2,:)
     endif
     if (nrad<0) then
        tau_xy (:,:,:)=tau (:,:,n1:n1+radz0-1)
        Irad_xy(:,:,:)=Irad(:,:,n1:n1+radz0-1)
     endif
     if (nrad>0) then
        tau_xy (:,:,:)=tau (:,:,n2-radz0+1:n2)
        Irad_xy(:,:,:)=Irad(:,:,n2-radz0+1:n2)
     endif
!
    endsubroutine radtransfer_intrinsic
!***********************************************************************
    subroutine intensity_intrinsic(lrad,mrad,nrad,tau,Irad)
!
!  Integration radiation transfer equation along all rays
!
!  This routine is called before the communication part
!  (certainly needs to be given a better name)
!  All rays start with zero intensity
!
!  16-jun-03/axel+tobi: coded
!
      use Cdata
!
      integer, intent(in) :: lrad,mrad,nrad
      real, dimension(mx,my,mz), intent(out) :: tau,Irad
      integer :: lstart,lstop,lsgn
      integer :: mstart,mstop,msgn
      integer :: nstart,nstop,nsgn
      real :: dlength,dtau,emdtau
      integer :: l
      logical, save :: first=.true.
!
!  identifier
!
      if(first) then
        print*,'intensity_intrinsic'
        first=.false.
      endif
!
!  calculate start and stop values
!
      if(lrad>=0) then; lstart=l1; lstop=l2; else; lstart=l2; lstop=l1; endif
      if(mrad>=0) then; mstart=m1; mstop=m2; else; mstart=m2; mstop=m1; endif
      if(nrad>=0) then; nstart=n1; nstop=n2; else; nstart=n2; nstop=n1; endif
!
!  make sure the loop is executed at least once, even when
!  lrad,mrad,nrad=0.
!
      if(lrad>=0) then; lsgn=1; else; lsgn=-1; endif
      if(mrad>=0) then; msgn=1; else; msgn=-1; endif
      if(nrad>=0) then; nsgn=1; else; nsgn=-1; endif
!
!  line elements
!
      dlength=sqrt((dx*lrad)**2+(dy*mrad)**2+(dz*nrad)**2)
!
!  set optical depth and intensity initially to zero
!
      tau=0.
      Irad=0.
!
!  loop
!
      do n=nstart,nstop,nsgn
      do m=mstart,mstop,msgn
      do l=lstart,lstop,lsgn 
          dtau=.5*(kaprho(l-lrad,m-mrad,n-nrad)+kaprho(l,m,n))*dlength
          tau(l,m,n)=tau(l-lrad,m-mrad,n-nrad)+dtau
          emdtau=exp(-dtau)
          Irad(l,m,n)=Irad(l-lrad,m-mrad,n-nrad)*emdtau &
                      +(1.-emdtau)*Srad(l-lrad,m-mrad,n-nrad) &
                      +(emdtau-1+dtau)*(Srad(l,m,n) &
                                       -Srad(l-lrad,m-mrad,n-nrad))/dtau
      enddo
      enddo
      enddo
!
    endsubroutine intensity_intrinsic
!***********************************************************************
    subroutine radtransfer_communicate(lrad,mrad,nrad)
!
!  Integration radioation transfer equation along rays
!
!  This routine is called after the communication part
!  The true boundary intensities I0 are now known and
!    the correction term I0*exp(-tau) is added
!  16-jun-03/axel+tobi: coded
!
      use Cdata
      use Sub
!
      integer, intent(in) :: lrad,mrad,nrad
      real, dimension(mx,my,mz) :: Irad0
      integer :: tag_xy=301,tag_yz=302,tag_zx=303
!
!  identifier
!
      if(lroot.and.headt) print*,'radtransfer_communicate'
!
!  set ghost zones, data from preceeding processor
!
      call radcomm_yz_recv(lrad,radx0,Irad0_yz,tag_yz)
      if (lrad>0) Irad0(l1-radx0:l1-1,:,:)=Irad0_yz(:,:,:)
      if (lrad<0) Irad0(l2+1:l2+radx0,:,:)=Irad0_yz(:,:,:)
      call radcomm_zx_recv(mrad,rady0,Irad0_zx,tag_zx)
      if (mrad>0) Irad0(:,m1-rady0:m1-1,:)=Irad0_zx(:,:,:)
      if (mrad<0) Irad0(:,m2+1:m2+rady0,:)=Irad0_zx(:,:,:)
      call radcomm_xy_recv(nrad,radz0,Irad0_xy,tag_xy)
      if (nrad>0) Irad0(:,:,n1-radz0:n1-1)=Irad0_xy(:,:,:)
      if (nrad<0) Irad0(:,:,n2+1:n2+radz0)=Irad0_xy(:,:,:)
      !
      !  propagate boundary values
      !
      call intensity_communicate(lrad,mrad,nrad,Irad0)
      !
      !  set boundary values for following processor
      !
      if (lrad<0) Irad0_yz(:,:,:)=Irad0(l1:l1+radx0-1,:,:)
      if (lrad>0) Irad0_yz(:,:,:)=Irad0(l2-radx0+1:l2,:,:)
      call radcomm_yz_send(lrad,radx0,Irad_yz,tag_yz)
      if (mrad<0) Irad0_zx(:,:,:)=Irad0(:,m1:m1+rady0-1,:)
      if (mrad>0) Irad0_zx(:,:,:)=Irad0(:,m2-rady0+1:m2,:)
      call radcomm_zx_send(mrad,rady0,Irad_zx,tag_zx)
      if (nrad<0) Irad0_xy(:,:,:)=Irad0(:,:,n1:n1+radz0-1)
      if (nrad>0) Irad0_xy(:,:,:)=Irad0(:,:,n2-radz0+1:n2)
      call radcomm_xy_send(nrad,radz0,Irad_xy,tag_xy)
!
    endsubroutine radtransfer_communicate
!***********************************************************************
    subroutine intensity_communicate(lrad,mrad,nrad,Irad0)
!
!  Integration radiation transfer equation along all rays
!
!  16-jun-03/axel+tobi: coded
!
      use Cdata
!
      integer, intent(in) :: lrad,mrad,nrad
      real, dimension(mx,my,mz), intent(inout) :: Irad0
      integer :: lstart,lstop,lsgn
      integer :: mstart,mstop,msgn
      integer :: nstart,nstop,nsgn
      integer :: l
      logical, save :: first=.true.
!
!  identifier
!
      if(first) then
        print*,'intensity_communicate'
        first=.false.
      endif
!
!  calculate start and stop values
!
      if(lrad>=0) then; lstart=l1; lstop=l2; else; lstart=l2; lstop=l1; endif
      if(mrad>=0) then; mstart=m1; mstop=m2; else; mstart=m2; mstop=m1; endif
      if(nrad>=0) then; nstart=n1; nstop=n2; else; nstart=n2; nstop=n1; endif
!
!  make sure the loop is executed at least once, even when
!  lrad,mrad,nrad=0.
!
      if(lrad>=0) then; lsgn=1; else; lsgn=-1; endif
      if(mrad>=0) then; msgn=1; else; msgn=-1; endif
      if(nrad>=0) then; nsgn=1; else; nsgn=-1; endif
!
!  loop
!
      do n=nstart,nstop,nsgn
      do m=mstart,mstop,msgn
      do l=lstart,lstop,lsgn
          Irad0(l,m,n)=Irad0(l-lrad,m-mrad,n-nrad)
      enddo
      enddo
      enddo
!
    endsubroutine intensity_communicate
!***********************************************************************
    subroutine radtransfer_revision(f,lrad,mrad,nrad)
!
!  Integration radioation transfer equation along rays
!
!  This routine is called after the communication part
!  The true boundary intensities I0 are now known and
!    the correction term I0*exp(-tau) is added
!  16-jun-03/axel+tobi: coded
!
      use Cdata
      use Sub
!
      integer, intent(in) :: lrad,mrad,nrad
      real, dimension(mx,my,mz,mvar+maux) :: f
      real, dimension(mx,my,mz) :: Irad0,Irad
!
!  identifier
!
      if(lroot.and.headt) print*,'radtransfer_revision'
!
!  set ghost zones, data from next processor (or opposite boundary)
!
      if(lrad>0) Irad0(l1-radx0:l1-1,:,:)=Irad0_yz(:,:,:)
      if(lrad<0) Irad0(l2+1:l2+radx0,:,:)=Irad0_yz(:,:,:)
      if(mrad>0) Irad0(:,m1-rady0:m1-1,:)=Irad0_zx(:,:,:)
      if(mrad<0) Irad0(:,m2+1:m2+rady0,:)=Irad0_zx(:,:,:)
      if(nrad>0) Irad0(:,:,n1-radz0:n1-1)=Irad0_xy(:,:,:)
      if(nrad<0) Irad0(:,:,n2+1:n2+radz0)=Irad0_xy(:,:,:)
!
!  do the ray, and add corresponding contribution to Q
!
      call intensity_revision(lrad,mrad,nrad,Irad0,Irad)
      f(:,:,:,iQrad)=f(:,:,:,iQrad)+frac*Irad
!
    endsubroutine radtransfer_revision
!***********************************************************************
    subroutine intensity_revision(lrad,mrad,nrad,Irad0,Irad)
!
!  Integration radiation transfer equation along all rays
!
!  16-jun-03/axel+tobi: coded
!
      use Cdata
!
      integer, intent(in) :: lrad,mrad,nrad
      real, dimension(mx,my,mz), intent(inout) :: Irad0
      real, dimension(mx,my,mz), intent(out) :: Irad
      real, dimension(mx,my,mz) :: tau
      integer :: lstart,lstop,lsgn
      integer :: mstart,mstop,msgn
      integer :: nstart,nstop,nsgn
      real :: dlength,dtau
      integer :: l
      logical, save :: first=.true.
!
!  identifier
!
      if(first) then
        print*,'intensity_revision'
        first=.false.
      endif
!
!  calculate start and stop values
!
      if(lrad>=0) then; lstart=l1; lstop=l2; else; lstart=l2; lstop=l1; endif
      if(mrad>=0) then; mstart=m1; mstop=m2; else; mstart=m2; mstop=m1; endif
      if(nrad>=0) then; nstart=n1; nstop=n2; else; nstart=n2; nstop=n1; endif
!
!  make sure the loop is executed at least once, even when
!  lrad,mrad,nrad=0.
!
      if(lrad>=0) then; lsgn=1; else; lsgn=-1; endif
      if(mrad>=0) then; msgn=1; else; msgn=-1; endif
      if(nrad>=0) then; nsgn=1; else; nsgn=-1; endif
!
!  line elements
!
      dlength=sqrt((dx*lrad)**2+(dy*mrad)**2+(dz*nrad)**2)
!
!  set optical depth initially to zero
!
      tau=0.
!
!  loop
!
      do n=nstart,nstop,nsgn
      do m=mstart,mstop,msgn
      do l=lstart,lstop,lsgn
          dtau=.5*(kaprho(l-lrad,m-mrad,n-nrad)+kaprho(l,m,n))*dlength
          tau(l,m,n)=tau(l-lrad,m-mrad,n-nrad)+dtau
          Irad0(l,m,n)=Irad0(l-lrad,m-mrad,n-nrad)
          Irad(l,m,n)=Irad0(l,m,n)*exp(-tau(l,m,n))
      enddo
      enddo
      enddo
!
    endsubroutine intensity_revision
!***********************************************************************
    subroutine radiative_cooling(f,df)
!
!  calculate source function
!
!  25-mar-03/axel+tobi: coded
!
      use Cdata
      use Ionization
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real :: formfactor=0.5
!
!  Add radiative cooling
!
      do n=n1,n2
      do m=m1,m2
         if(.not. nocooling) then
            df(l1:l2,m,n,ient)=df(l1:l2,m,n,ient) &
                              +4.*pi*kaprho(l1:l2,m,n) &
                               *f(l1:l2,m,n,iQrad) &
                               /f(l1:l2,m,n,iTT)*formfactor &
                               *exp(-f(l1:l2,m,n,ilnrho))
         endif
      enddo
      enddo
!
    endsubroutine radiative_cooling
!***********************************************************************
    subroutine output_radiation(lun)
!
!  Optional output of derived quantities along with VAR-file
!  Called from wsnap
!
!   5-apr-03/axel: coded
!
      use Cdata
      use Ionization
!
      integer, intent(in) :: lun
!
!  identifier
!
      !if(lroot.and.headt) print*,'output_radiation',Qrad(4,4,4)
      !if(output_Qrad) write(lun) Qrad,Srad,kaprho,TT
!
    endsubroutine output_radiation
!***********************************************************************
    subroutine init_rad(f,xx,yy,zz)
!
!  Dummy routine for Flux Limited Diffusion routine
!  initialise radiation; called from start.f90
!
!  15-jul-2002/nils: dummy routine
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz)      :: xx,yy,zz
!
      if(ip==0) print*,f,xx,yy,zz !(keep compiler quiet)
    endsubroutine init_rad
!***********************************************************************
   subroutine de_dt(f,df,rho1,divu,uu,uij,TT1,gamma)
!
!  Dummy routine for Flux Limited Diffusion routine
!
!  15-jul-2002/nils: dummy routine
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: uu
      real, dimension (nx) :: rho1,TT1
      real, dimension (nx,3,3) :: uij
      real, dimension (nx) :: divu
      real :: gamma
!
      if(ip==0) print*,f,df,rho1,divu,uu,uij,TT1,gamma !(keep compiler quiet)
    endsubroutine de_dt
!*******************************************************************
    subroutine rprint_radiation(lreset)
!
!  Dummy routine for Flux Limited Diffusion routine
!  reads and registers print parameters relevant for radiative part
!
!  16-jul-02/nils: adapted from rprint_hydro
!
      use Cdata
      use Sub
!  
      logical :: lreset
!
!  write column where which radiative variable is stored
!
      write(3,*) 'i_frms=',i_frms
      write(3,*) 'i_fmax=',i_fmax
      write(3,*) 'i_Erad_rms=',i_Erad_rms
      write(3,*) 'i_Erad_max=',i_Erad_max
      write(3,*) 'i_Egas_rms=',i_Egas_rms
      write(3,*) 'i_Egas_max=',i_Egas_max
      write(3,*) 'nname=',nname
      write(3,*) 'ie=',ie
      write(3,*) 'ifx=',ifx
      write(3,*) 'ify=',ify
      write(3,*) 'ifz=',ifz
      write(3,*) 'iQrad=',iQrad
!   
      if(ip==0) print*,lreset  !(to keep compiler quiet)
    endsubroutine rprint_radiation
!***********************************************************************
    subroutine  bc_ee_inflow_x(f,topbot)
!
!  Dummy routine for Flux Limited Diffusion routine
!
!  8-aug-02/nils: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      if (ip==1) print*,topbot,f(1,1,1,1)  !(to keep compiler quiet)
!
    end subroutine bc_ee_inflow_x
!***********************************************************************
    subroutine  bc_ee_outflow_x(f,topbot)
!
!  Dummy routine for Flux Limited Diffusion routine
!
!  8-aug-02/nils: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      if (ip==1) print*,topbot,f(1,1,1,1)  !(to keep compiler quiet)
!
    end subroutine bc_ee_outflow_x
!***********************************************************************

endmodule Radiation
