! $Id: radiation_ray.f90,v 1.2 2003-03-25 05:42:57 brandenb Exp $

module Radiation

!  Radiation (solves transfer equation along rays)

  use Cparam

  implicit none

  real, dimension (mx,my,mz) :: Qrad,Source
  real :: arad=0.,kappa=0.
  integer :: radx=0,rady=0,radz=1,rad2max=1
!
!  definition of dummy variables for FLD routine
!
  real :: DFF_new=0.  !(dum)
  integer :: i_frms=0,i_fmax=0,i_Erad_rms=0,i_Erad_max=0
  integer :: i_Egas_rms=0,i_Egas_max=0

  namelist /radiation_init_pars/ &
       kappa,arad,radx,rady,radz,rad2max

  namelist /radiation_run_pars/ &
       kappa,arad,radx,rady,radz,rad2max

  contains

!***********************************************************************
    subroutine register_radiation()
!
!  initialise radiation flags
!
! 24-mar-03/axel+tobi: coded
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
           "$Id: radiation_ray.f90,v 1.2 2003-03-25 05:42:57 brandenb Exp $")
!
    endsubroutine register_radiation
!***********************************************************************
    subroutine initialize_radiation()
!
!  Set up some variables for radiation; do nothing
!  Possible examples: reading in of tables etc
!
    endsubroutine initialize_radiation
!***********************************************************************
    subroutine source_function(f)
!
!  calculate source function
!
! 24-mar-03/axel+tobi: coded
!
      use Cdata
      use Ionization
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (nx) :: lnrho,ss,cs2,TT,TT1,cp1tilde
!
      intent(in)  :: f
!
!  Use the thermodynamics module to calculate temperature
!  At the moment we don't calculate ghost zones (ok for vertical arrays)  
!
      do n=1,mz
      do m=1,my
        lnrho=f(l1:l2,m,n,ilnrho)
        ss=f(l1:l2,m,n,ient)
        call thermodynamics(lnrho,ss,cs2,TT1,cp1tilde,Temperature=TT)
        Source(l1:l2,m,n)=arad*TT**4
        !Source(l1:l2,m,n)=arad*1./TT1**4
      enddo
      enddo
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
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: Intensity
      real :: frac
      integer :: lrad,mrad,nrad,rad2
      integer :: directions=2  !(at the moment)
!
!  identifier
!
      if(lroot.and.headt) print*,'radtransfer'
!
      call source_function(f)
!
!  set boundary values, I=0 at the top, I=S at the bottom
!
      Intensity(:,:,n2+1:mz)=0.
      Intensity(:,:,1:n1-1)=Source(:,:,1:n1-1)
!
!  Accumulate the result for Qrad=(J-S),
!  First initialize Qrad=-S. 
!
      Qrad=-Source
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
          call transfer(f,Intensity,lrad,mrad,nrad)
          Qrad=Qrad+frac*Intensity
        endif
      enddo
      enddo
      enddo
write(28) Qrad
!
    endsubroutine radtransfer
!***********************************************************************
    subroutine transfer(f,Intensity,lrad,mrad,nrad)
!
!  Integration radioation transfer equation along all rays
!
!  24-mar-03/axel+tobi: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: Intensity
      real :: dlength,lnrhom,fnew,fold,dtau05
      integer :: l,lrad,mrad,nrad,rad2
      integer :: lstart,lstop,lrad1
      integer :: mstart,mstop,mrad1
      integer :: nstart,nstop,nrad1
!
!  identifier
!
      if(lroot.and.headt) print*,'transfer'
!
!  calculate start and stop values
!
      if(lrad>0) then; lstart=l1; lstop=l2; else; lstart=l2; lstop=l1; endif
      if(mrad>0) then; mstart=m1; mstop=m2; else; mstart=m2; mstop=m1; endif
      if(nrad>0) then; nstart=n1; nstop=n2; else; nstart=n2; nstop=n1; endif
!
!  avoid problems with lrad,mrad,nrad=0.
!
      if(lrad==0) then; lstart=0; lstop=0; lrad1=1; else; lrad1=lrad; endif
      if(mrad==0) then; mstart=0; mstop=0; mrad1=1; else; mrad1=mrad; endif
      if(nrad==0) then; nstart=0; nstop=0; nrad1=1; else; nrad1=nrad; endif
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
        dtau05=0.5*kappa*exp(lnrhom)*dlength
        fnew=1.+dtau05
        fold=1.-dtau05
        Intensity(l,m,n)=(fold*Intensity(l-lrad,m-mrad,n-nrad) &
           +dtau05*(Source(l,m,n)+Source(l-lrad,m-mrad,n-nrad)))
      enddo
      enddo
      enddo
!
    endsubroutine transfer
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
      real, dimension (mx,my,mz,mvar) :: f
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
      real, dimension (mx,my,mz,mvar) :: f
!
      if (ip==1) print*,topbot,f(1,1,1,1)  !(to keep compiler quiet)
!
    end subroutine bc_ee_outflow_x
!***********************************************************************

endmodule Radiation
