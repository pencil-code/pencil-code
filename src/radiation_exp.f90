! $Id: radiation_exp.f90,v 1.17 2003-06-19 21:28:17 brandenb Exp $

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
  integer, parameter :: radx0=3,rady0=3,radz0=3
  real, dimension(radx0,my,mz,-radx0:radx0,-rady0:rady0,-radz0:radz0) :: Irad_yz
  real, dimension(mx,rady0,mz,-radx0:radx0,-rady0:rady0,-radz0:radz0) :: Irad_zx
  real, dimension(mx,my,radz0,-radx0:radx0,-rady0:rady0,-radz0:radz0) :: Irad_xy
  real, dimension (mx,my,mz) :: Srad,kaprho
  integer :: directions
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
           "$Id: radiation_exp.f90,v 1.17 2003-06-19 21:28:17 brandenb Exp $")
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
!  side boundaries : initially I=0
!
      do nrad=-radz,radz
      do mrad=-rady,rady
      do lrad=-radx,radx
        Irad_yz(:,:,:,lrad,mrad,nrad)=0
        Irad_zx(:,:,:,lrad,mrad,nrad)=0
      enddo
      enddo
      enddo
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
      real, dimension(nx) :: lnrho,ss,yH,TT,kappa
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
      do n=n1,n2
      do m=m1,m2
         lnrho=f(l1:l2,m,n,ilnrho)
         ss=f(l1:l2,m,n,ient)
!
!  Use the ionization module to calculate temperature.
!  If noionization is used, yH and TT do not exist as 3-D arrays,
!  but we still want Srad and kaprho to be full arrays.
!
         call ionset(f,ss,lnrho,yH,TT)
!
!  opacity: if lkappa_es then take electron scattering opacity only;
!  otherwise use Hminus opacity (but may need to add kappa_es as well).
!
         if(lkappa_es) then
           kappa=kappa_es
         else
           kappa=.25*exp(lnrho-lnrho_ion_)*(TT_ion_/TT)**1.5 &
                 *exp(TT_ion_/TT)*yH*(1.-yH)*kappa0
         endif
!
!  save the 3-D arrays Srad and kaprho (=kappa*rho) because they are used
!  many times, depending on the number of rays in the radiation module.
!
         Srad(l1:l2,m,n)=sigmaSB*TT**4/pi
         kaprho(l1:l2,m,n)=kappa*exp(lnrho)
      enddo
      enddo
!
    endsubroutine radcalc
!***********************************************************************
    function mean_intensity(f)
!
!  Solves the transfer equation using Feautrier's method
!  At the moment for vertical rays only
!
!  12-may-03/tobi: coded
!
      use Cdata
!
      real, dimension(mx,my,mz,mvar+maux) :: f
      real, dimension(mx,my,mz) :: mean_intensity,Iup,Idown
      real :: dtau,emdtau
      integer :: lr,mr,nr
!
!  loop fastest over first index --> could do this for l-vectors
!
      do mr=m1,m2
      do lr=l1,l2
        Iup(lr,mr,n1)=Srad(lr,mr,n1)
        do nr=n1,n2
          dtau=.5*(kaprho(lr,mr,nr)+kaprho(lr,mr,nr+1))*dz
          emdtau=exp(-dtau)
          Iup(lr,mr,nr+1)=Iup(lr,mr,nr)*emdtau &
                          +(1.-emdtau)*Srad(lr,mr,nr) &
                          +(emdtau-1+dtau) &
                           *(Srad(lr,mr,nr+1)-Srad(lr,mr,nr))/dtau
        enddo
        Idown(lr,mr,n2)=0
        do nr=n2,n1,-1
          dtau=.5*(kaprho(lr,mr,nr)+kaprho(lr,mr,nr-1))*dz
          emdtau=exp(-dtau)
          Idown(lr,mr,nr-1)=Idown(lr,mr,nr)*emdtau &
                            +(1.-emdtau)*Srad(lr,mr,nr) &
                            +(emdtau-1+dtau) &
                             *(Srad(lr,mr,nr-1)-Srad(lr,mr,nr))/dtau
        enddo
      enddo
      enddo
      mean_intensity=.5*(Iup+Idown)
    endfunction mean_intensity
!***********************************************************************
    subroutine radtransfer_old(f)
!
!  Integration radiation transfer equation along rays
!
!  24-mar-03/axel+tobi: coded
!
      use Cdata
      use Sub
!
      real, dimension(mx,my,mz,mvar+maux) :: f
!
!  identifier
!
      if(lroot.and.headt) print*,'radtransfer'
!
      call radcalc(f)
      f(:,:,:,iQrad)=-Srad+mean_intensity(f)
!
    endsubroutine radtransfer_old
!***********************************************************************
    subroutine radtransfer(f)
!
!  Integration radioation transfer equation along rays
!
!  16-jun-03/axel+tobi: coded
!
      use Cdata
      use Sub
!
      real, dimension(mx,my,mz,mvar+maux) :: f
      real, dimension(mx,my,mz) :: Irad
      real :: frac
      integer :: lrad,mrad,nrad,rad2,i
      integer :: counter=0
!
!  identifier
!
      if(lroot.and.headt) print*,'radtransfer'
!
!  calculate source function and opacity
!
      call radcalc(f)
!
!  set boundary values
!  bottom boundary (rays point upwards): I=S
!
      do nrad=+1,+radz
      do mrad=-rady,rady
      do lrad=-radx,radx
        Irad_xy(:,:,:,lrad,mrad,nrad)=Srad(:,:,n1-radz0:n1-1)
      enddo
      enddo
      enddo
!
!  top boundary (rays point downwards): I=0
!
      do nrad=-radz,-1
      do mrad=-rady,rady
      do lrad=-radx,radx
        !in principle ok
        !Irad_xy(:,:,:,lrad,mrad,nrad)=0.
        ! to make sure we reproduce old results
        do i=1,nghost
          Irad_xy(:,:,i,lrad,mrad,nrad)=1.-exp(dz*i)
        enddo
      enddo
      enddo
      enddo
!
!  Accumulate the result for Qrad=(J-S),
!  First initialize Qrad=-S. 
!
      f(:,:,:,iQrad)=-Srad
!
!  loop over rays
!
      frac=1./directions
!
      do nrad=-radz,radz
      do mrad=-rady,rady
      do lrad=-radx,radx
        rad2=lrad**2+mrad**2+nrad**2
        if(rad2>0 .and. rad2<=rad2max) then 
          !
          !  set ghost zones, data from next processor (or opposite boundary)
          !
          if(lrad>0) Irad(l1-radx0:l1-1,:,:)=Irad_yz(:,:,:,lrad,mrad,nrad)
          if(lrad<0) Irad(l2+1:l2+radx0,:,:)=Irad_yz(:,:,:,lrad,mrad,nrad)
          if(mrad>0) Irad(:,m1-rady0:m1-1,:)=Irad_zx(:,:,:,lrad,mrad,nrad)
          if(mrad<0) Irad(:,m2+1:m2+rady0,:)=Irad_zx(:,:,:,lrad,mrad,nrad)
          if(nrad>0) Irad(:,:,n1-radz0:n1-1)=Irad_xy(:,:,:,lrad,mrad,nrad)
          if(nrad<0) Irad(:,:,n2+1:n2+radz0)=Irad_xy(:,:,:,lrad,mrad,nrad)
          !
          !  do the ray, and add corresponding contribution to Q
          !
          call transfer(lrad,mrad,nrad,Irad)
          f(:,:,:,iQrad)=f(:,:,:,iQrad)+frac*Irad
          !
          !  safe boundary values for next processor (or opposite boundary)
          !
          if(lrad<0) Irad_yz(:,:,:,lrad,mrad,nrad)=Irad(l1-radx0:l1-1,:,:)
          if(lrad>0) Irad_yz(:,:,:,lrad,mrad,nrad)=Irad(l2+1:l2+radx0,:,:)
          if(mrad<0) Irad_zx(:,:,:,lrad,mrad,nrad)=Irad(:,m1-rady0:m1-1,:)
          if(mrad>0) Irad_zx(:,:,:,lrad,mrad,nrad)=Irad(:,m2+1:m2+rady0,:)
          if(nrad<0) Irad_xy(:,:,:,lrad,mrad,nrad)=Irad(:,:,n1-radz0:n1-1)
          if(nrad>0) Irad_xy(:,:,:,lrad,mrad,nrad)=Irad(:,:,n2+1:n2+radz0)
        endif
      enddo
      enddo
      enddo
!
      !print*,'Number of directions in this run:',counter
!
    endsubroutine radtransfer
!***********************************************************************
    subroutine transfer(lrad,mrad,nrad,Irad)
!
!  Integration radiation transfer equation along all rays
!
!  16-jun-03/axel+tobi: coded
!
      use Cdata
!
      integer :: lrad,mrad,nrad
      real, dimension(mx,my,mz) :: Irad
      integer :: lstart,lstop,lrad1
      integer :: mstart,mstop,mrad1
      integer :: nstart,nstop,nrad1
      real :: dlength,dtau,emdtau
      integer :: l
      logical, save :: first=.true.
!
!  identifier
!
      if(first) then
        print*,'transfer'
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
      if(lrad==0) then; lrad1=1; else; lrad1=lrad; endif
      if(mrad==0) then; mrad1=1; else; mrad1=mrad; endif
      if(nrad==0) then; nrad1=1; else; nrad1=nrad; endif
!
!  line elements
!
      dlength=sqrt((dx*lrad)**2+(dy*mrad)**2+(dz*nrad)**2)
!
!
!  loop
!
      do n=nstart,nstop,nrad1
      do m=mstart,mstop,mrad1
      do l=lstart,lstop,lrad1
          dtau=.5*(kaprho(l-lrad,m-mrad,n-nrad)+kaprho(l,m,n))*dlength
          emdtau=exp(-dtau)
          Irad(l,m,n)=Irad(l-lrad,m-mrad,n-nrad)*emdtau &
                      +(1.-emdtau)*Srad(l-lrad,m-mrad,n-nrad) &
                      +(emdtau-1+dtau)*(Srad(l,m,n) &
                                       -Srad(l-lrad,m-mrad,n-nrad))/dtau
      enddo
      enddo
      enddo
!
    endsubroutine transfer
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
