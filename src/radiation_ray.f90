! $Id: radiation_ray.f90,v 1.18 2003-06-13 11:56:11 nilshau Exp $

module Radiation

!  Radiation (solves transfer equation along rays)
!  The direction of the ray is given by the vector (lrad,mrad,nrad),
!  and the parameters radx0,rady0,radz0 gives the maximum number of
!  steps of the direction vector in the corresponding direction.

  use Cparam

  implicit none

  logical :: nocooling=.false.,output_Qrad=.false.

  integer :: directions
  integer, parameter :: radx0=3,rady0=3,radz0=3
  real, dimension (radx0,my,mz,-radx0:radx0,-rady0:rady0,-radz0:radz0) :: Intensity_yz
  real, dimension (mx,rady0,mz,-radx0:radx0,-rady0:rady0,-radz0:radz0) :: Intensity_zx
  real, dimension (mx,my,radz0,-radx0:radx0,-rady0:rady0,-radz0:radz0) :: Intensity_xy
!
!  default values for one pair of vertical rays
!
  integer :: radx=0,rady=0,radz=1,rad2max=1
!
!  definition of dummy variables for FLD routine
!
  real :: DFF_new=0.  !(dum)
  integer :: i_frms=0,i_fmax=0,i_Erad_rms=0,i_Erad_max=0
  integer :: i_Egas_rms=0,i_Egas_max=0

  namelist /radiation_init_pars/ &
       radx,rady,radz,rad2max,output_Qrad

  namelist /radiation_run_pars/ &
       radx,rady,radz,rad2max,output_Qrad,nocooling

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
      lradiation = .true.
      lradiation_ray = .true.
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: radiation_ray.f90,v 1.18 2003-06-13 11:56:11 nilshau Exp $")
!
    endsubroutine register_radiation
!***********************************************************************
    subroutine initialize_radiation()
!
!  Calculate number of directions of rays
!  Do this in the beginning of each run
!
!  25-mar-03/axel+tobi: coded
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
    subroutine source_function(f)
!
!  calculate source function
!
!  24-mar-03/axel+tobi: coded
!
      use Cdata
      use Ionization
!
      real, dimension(mx,my,mz,mvar+maux) :: f
      real, dimension(nx) :: lnrho,ss,yH,kappa_,cs2,TT1,cp1tilde
!
!  Use the ionization module to calculate temperature
!  At the moment we don't calculate ghost zones (ok for vertical arrays)  
!  Need Source function even in the ghost zones (for lower bc)
!
print*,'ss_border=',f(l1,m1,1:7,ient)
print*,'lnrho_border=',f(l1,m1,1:7,ilnrho)
!!!!  do n=n1,n2
      do n=1,mz
      do m=m1,m2
         lnrho=f(l1:l2,m,n,ilnrho)
         ss=f(l1:l2,m,n,ient)
!         yH=yyH(l1:l2,m,n) yh is beeing calculated in ioncalc
         call ioncalc(lnrho,ss,yH,kappa=kappa_)
         call thermodynamics(lnrho,ss,cs2,TT1,cp1tilde)
         f(l1:l2,m,n,iSrad)=sigmaSB/(pi*TT1**4)
         f(l1:l2,m,n,iTT)=1./TT1
         f(l1:l2,m,n,ikappa)=kappa_
      enddo
      enddo
print*,'Srad_border=',f(l1,m1,1:7,iSrad)
!
    endsubroutine source_function
!***********************************************************************
    subroutine radtransfer(f)
!
!  Integration radioation transfer equation along rays
!
!  24-mar-03/axel+tobi: coded
!
      use Cdata
      use Sub
!
      real :: frac
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: Intensity
      integer :: lrad,mrad,nrad,rad2
!
!  identifier
!
      if(lroot.and.headt) print*,'radtransfer'
!
      call source_function(f)
!
!  set boundary values
!  bottom boundary (rays point upwards): I=S
!
      do nrad=+1,+radz
      do mrad=-rady,rady
      do lrad=-radx,radx
        Intensity_xy(:,:,:,lrad,mrad,nrad)=f(:,:,n1-radz0:n1-1,iSrad)
      enddo
      enddo
      enddo
!
!  top boundary (rays point downwards): I=0
!
      do nrad=-radz,-1
      do mrad=-rady,rady
      do lrad=-radx,radx
        Intensity_xy(:,:,:,lrad,mrad,nrad)=0.
      enddo
      enddo
      enddo
!
!  Accumulate the result for Qrad=(J-S),
!  First initialize Qrad=-S. 
!
      f(:,:,:,iQrad)=-f(:,:,:,iSrad)
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
          if(lrad>0) Intensity(l1-radx0:l1-1,:,:)=Intensity_yz(:,:,:,lrad,mrad,nrad)
          if(lrad<0) Intensity(l2+1:l2+radx0,:,:)=Intensity_yz(:,:,:,lrad,mrad,nrad)
          if(mrad>0) Intensity(:,m1-rady0:m1-1,:)=Intensity_zx(:,:,:,lrad,mrad,nrad)
          if(mrad<0) Intensity(:,m2+1:m2+rady0,:)=Intensity_zx(:,:,:,lrad,mrad,nrad)
          if(nrad>0) Intensity(:,:,n1-radz0:n1-1)=Intensity_xy(:,:,:,lrad,mrad,nrad)
          if(nrad<0) Intensity(:,:,n2+1:n2+radz0)=Intensity_xy(:,:,:,lrad,mrad,nrad)
          !
          !  do the ray, and add corresponding contribution to Q
          !
          call transfer(f,Intensity,lrad,mrad,nrad)
          f(:,:,:,iQrad)=f(:,:,:,iQrad)+frac*Intensity
write(28) Intensity,nrad
          !
          !  safe boundary values for next processor (or opposite boundary)
          !
          if(lrad<0) Intensity_yz(:,:,:,lrad,mrad,nrad)=Intensity(l1-radx0:l1-1,:,:)
          if(lrad>0) Intensity_yz(:,:,:,lrad,mrad,nrad)=Intensity(l2+1:l2+radx0,:,:)
          if(mrad<0) Intensity_zx(:,:,:,lrad,mrad,nrad)=Intensity(:,m1-rady0:m1-1,:)
          if(mrad>0) Intensity_zx(:,:,:,lrad,mrad,nrad)=Intensity(:,m2+1:m2+rady0,:)
          if(nrad<0) Intensity_xy(:,:,:,lrad,mrad,nrad)=Intensity(:,:,n1-radz0:n1-1)
          if(nrad>0) Intensity_xy(:,:,:,lrad,mrad,nrad)=Intensity(:,:,n2+1:n2+radz0)
        endif
      enddo
      enddo
      enddo
!
    endsubroutine radtransfer
!***********************************************************************
    subroutine transfer(f,Intensity,lrad,mrad,nrad)
!
!  Integration radiation transfer equation along all rays
!
!  24-mar-03/axel+tobi: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: Intensity
      real :: dlength,lnrhom,fnew,fold,dtau05
      integer :: l,lrad,mrad,nrad
      integer :: lstart,lstop,lrad1
      integer :: mstart,mstop,mrad1
      integer :: nstart,nstop,nrad1
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
!  dI/dlength = -kappa*rho*(I-S)
!
!  in discretized form, regardless of the signs of lrad,mrad,nrad,
!  [I(l,m,n)-I(l-lrad,m-mrad,n-nrad)]/dlength = -kappa*rho*(I-S)
!
!  Define dtau=kappa*rho*dlength (always positive!), so
!  I(l,m,n)*(1+.5*dtau) = I(l-lrad,m-mrad,n-nrad)*(1-.5*dtau) + dtau*Smean
!
!  loop
!
      do n=nstart,nstop,nrad1
      do m=mstart,mstop,mrad1
      do l=lstart,lstop,lrad1
        lnrhom=0.5*(f(l,m,n,ilnrho)+f(l-lrad,m-mrad,n-nrad,ilnrho))
        dtau05=0.5*f(l,m,n,ikappa)*exp(lnrhom)*dlength
        fnew=1.+dtau05
        fold=1.-dtau05
        Intensity(l,m,n)=(fold*Intensity(l-lrad,m-mrad,n-nrad) &
           +dtau05*(f(l,m,n,iSrad)+f(l-lrad,m-mrad,n-nrad,iSrad)))/fnew
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
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
!
!  Add radiative cooling
!
      do n=n1,n2
      do m=m1,m2
         if(.not. nocooling) then
            df(l1:l2,m,n,ient)=df(l1:l2,m,n,ient) &
                              +4.*pi*f(l1:l2,m,n,ikappa) &
                               *f(l1:l2,m,n,iQrad) &
                               /f(l1:l2,m,n,iTT)
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
!
      integer, intent(in) :: lun
!
!  identifier
!
!      if(lroot.and.headt) print*,'output_radiation',Qrad(4,4,4)
!      if(output_Qrad) write(lun) f(:,:,:,iQrad),f(:,:,:,iSrad), &
!           f(:,:,:,ikappa),f(:,:,:,iTT)
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
      real, dimension (mx,my,mz,mvar) :: f
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
      real, dimension (mx,my,mz,mvar) :: f,df
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
      write(3,*) 'iSrad=',iSrad
      write(3,*) 'ikappa=',ikappa
      write(3,*) 'iTT=',iTT
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
