! $Id: feautrier.f90,v 1.23 2003-06-13 09:28:58 nilshau Exp $

!!!  NOTE: this routine will perhaps be renamed to radation_feautrier
!!!  or it may be combined with radiation_ray.

module Radiation

!  Radiation (solves transfer equation along rays)
!  The direction of the ray is given by the vector (lrad,mrad,nrad),
!  and the parameters radx0,rady0,radz0 gives the maximum number of
!  steps of the direction vector in the corresponding direction.

  use Cparam

  implicit none

  logical :: nocooling=.false.,output_Qrad=.false.
!
!  default values for one pair of vertical rays
!
  integer :: radx=0,rady=0,radz=1,rad2max=1
!
! init parameteres
!
  character (len=labellen) :: initrad='equil',pertee='none'
  real :: amplee=0
  real :: ampl_pert=0
!
!  definition of dummy variables for FLD routine
!
  real :: DFF_new=0.  !(dum)
  integer :: i_frms=0,i_fmax=0,i_Erad_rms=0,i_Erad_max=0
  integer :: i_Egas_rms=0,i_Egas_max=0

  namelist /radiation_init_pars/ &
       radx,rady,radz,rad2max,output_Qrad,initrad,amplee,pertee,ampl_pert

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
      iQrad = mvar + naux + 1
      iSrad = mvar + naux + 2
      ikappa = mvar + naux + 3
      iTT = mvar + naux + 4
      naux = naux + 4 
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: feautrier.f90,v 1.23 2003-06-13 09:28:58 nilshau Exp $")
!
! Check we arn't registering too many auxilliary variables
!
      if (naux > maux) then
        if (lroot) write(0,*) 'naux = ', naux, ', maux = ', maux
        call stop_it('register_viscosityfeautrier: naux > maux')
      endif
!
    endsubroutine register_radiation
!***********************************************************************
    subroutine initialize_radiation()
!
!  nothing to be done here
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
      real, dimension(mx,my,mz,mvar + maux) :: f
      real, dimension(nx) :: lnrho,ss,yH,kappa_,cs2,TT1,cp1tilde
!
!  Use the ionization module to calculate temperature
!  At the moment we don't calculate ghost zones (ok for vertical arrays)  
!
      do n=n1,n2
      do m=m1,m2
         lnrho=f(l1:l2,m,n,ilnrho)
         ss=f(l1:l2,m,n,ient)
         call ioncalc(lnrho,ss,yH,kappa=kappa_)
         call thermodynamics(lnrho,ss,cs2,TT1,cp1tilde)
         f(l1:l2,m,n,iSrad)=sigmaSB/(pi*TT1**4)
         f(l1:l2,m,n,ikappa)=kappa_
         f(l1:l2,m,n,iTT)=1./TT1
      enddo
      enddo
!
    endsubroutine source_function
!***********************************************************************
    function feautrier(f)
!
!  Solves the transfer equation using Feautrier's method
!  At the moment for vertical rays only
!
!  01-apr/tobi: coded
!
      use Cdata
      use General
      use Ionization
!
      real, dimension(mx,my,mz,mvar + maux) :: f
      real, dimension(mx,my,mz) :: feautrier
      real, dimension(nz) :: kaprho,tau,Srad_,Prad_
      real, dimension(nz) :: a,b,c
      integer :: lrad,mrad,nrad
      logical :: err=.false.
!
!  loop fastest over first index --> could do this for l-vectors
!
      do mrad=m1,m2
      do lrad=l1,l2
         kaprho=f(lrad,mrad,n1:n2,ikappa)*exp(f(lrad,mrad,n1:n2,ilnrho))
         tau=spline_integral(z,kaprho)
         Srad_=f(lrad,mrad,n1:n2,iSrad)
         !print*,'kappa=',kappa(lrad,mrad,n1:n2)
         !print*,'tau=',tau
         !print*,'Srad=',Srad_
!
!  top boundary: P'=P, together with P1-P0=dtau*P'+.5*dtau^2*P" and P"=P-S
!  lower boundary: P=S
!
         b(1)=1.+2./(tau(2)-tau(1))+2./(tau(2)-tau(1))**2
         c(1)=-2./(tau(2)-tau(1))**2
         a(nz)=0.
         b(nz)=1.
!
!  interior, P"=P-S
!
         a(2:nz-1)=  -2./(tau(3:nz)-tau(1:nz-2))/(tau(2:nz-1)-tau(1:nz-2))
         b(2:nz-1)=1.+2./(tau(3:nz)-tau(1:nz-2))/(tau(2:nz-1)-tau(1:nz-2)) &
                     +2./(tau(3:nz)-tau(1:nz-2))/(tau(3:nz)-tau(2:nz-1))
         c(2:nz-1)=  -2./(tau(3:nz)-tau(1:nz-2))/(tau(3:nz)-tau(2:nz-1))
!
!  solve tri-diagonal matrix, and give detailed error output if problems
!
         call tridag(a,b,c,Srad_,Prad_,err=err)
         if (err) then
            print*,'lnrho=',f(lrad,mrad,n1:n1+5,ilnrho),'...', &
                            f(lrad,mrad,n2-5:n2,ilnrho)
            print*,'ss=',f(lrad,mrad,n1:n1+5,ient),'...', &
                         f(lrad,mrad,n2-5:n2,ient)
            print*,'tau=',tau(1:6),'...',tau(n2-n1-5:n2-n1)
            print*,'kappa=',f(lrad,mrad,n1:n1+5,ikappa),'...', &
                            f(lrad,mrad,n2-5:n2,ikappa)
            stop
         endif
!
         feautrier(lrad,mrad,n1:n2)=Prad_
!         print*,'Prad',Prad_
      enddo
      enddo
    endfunction feautrier
!***********************************************************************
    function feautrier_double(f)
!
!  Solves the transfer equation using Feautrier's method
!  At the moment for vertical rays only
!
!  01-apr/tobi: coded
!  11-apr/axel: turned some matrix stuff into double precision
!
      use Cdata
      use General
      use Ionization
!
      real, dimension(mx,my,mz,mvar + maux) :: f
      real, dimension(mx,my,mz) :: feautrier_double
      double precision, dimension(nz) :: kaprho,tau,Srad_,Prad_
      double precision, dimension(nz) :: a,b,c
      integer :: lrad,mrad,nrad
      logical :: err=.false.
!
!  loop fastest over first index --> could do this for l-vectors
!
      do mrad=m1,m2
      do lrad=l1,l2
         kaprho=f(lrad,mrad,n1:n2,ikappa)*exp(f(lrad,mrad,n1:n2,ilnrho))
         tau=spline_integral_double(z,kaprho)
         Srad_=f(lrad,mrad,n1:n2,iSrad)
         !print*,'kappa=',kappa(lrad,mrad,n1:n2)
         !print*,'tau=',tau
         !print*,'Srad=',Srad_
!
!  top boundary: P'=P, together with P1-P0=dtau*P'+.5*dtau^2*P" and P"=P-S
!  lower boundary: P=S
!
         b(1)=1.+2./(tau(2)-tau(1))+2./(tau(2)-tau(1))**2
         c(1)=-2./(tau(2)-tau(1))**2
         a(nz)=0.
         b(nz)=1.
!
!  interior, P"=P-S
!
         a(2:nz-1)=  -2./(tau(3:nz)-tau(1:nz-2))/(tau(2:nz-1)-tau(1:nz-2))
         b(2:nz-1)=1.+2./(tau(3:nz)-tau(1:nz-2))/(tau(2:nz-1)-tau(1:nz-2)) &
                     +2./(tau(3:nz)-tau(1:nz-2))/(tau(3:nz)-tau(2:nz-1))
         c(2:nz-1)=  -2./(tau(3:nz)-tau(1:nz-2))/(tau(3:nz)-tau(2:nz-1))
!
!  solve tri-diagonal matrix, and give detailed error output if problems
!
         call tridag_double(a,b,c,Srad_,Prad_,err=err)
         if (err) then
            print*,'lnrho=',f(lrad,mrad,n1:n1+5,ilnrho),'...', &
                            f(lrad,mrad,n2-5:n2,ilnrho)
            print*,'ss=',f(lrad,mrad,n1:n1+5,ient),'...', &
                         f(lrad,mrad,n2-5:n2,ient)
            print*,'tau=',tau(1:6),'...',tau(n2-n1-5:n2-n1)
            print*,'kappa=',f(lrad,mrad,n1:n1+5,ikappa),'...', &
                            f(lrad,mrad,n2-5:n2,ikappa)
            stop
         endif
!
         feautrier_double(lrad,mrad,n1:n2)=Prad_
!         print*,'Prad',Prad_
      enddo
      enddo

!print*,'her2'
!print*,maxval(kaprho),maxval(tau),maxval(Srad_)
!print*,minval(kaprho),minval(tau),minval(Srad_)
!print*,maxval(Prad_),minval(Prad_)
    endfunction feautrier_double
!***********************************************************************
    subroutine radtransfer(f)
!
!  Integration radiation transfer equation along rays
!
!  24-mar-03/axel+tobi: coded
!
      use Cdata
      use Sub
      use Ionization
!
      real, dimension(mx,my,mz,mvar + maux) :: f
!
!  identifier
!
      if(lroot.and.headt) print*,'radtransfer'
!
      call source_function(f)
      f(:,:,:,iQrad)=-f(:,:,:,iSrad)
      f(:,:,:,iQrad)=f(:,:,:,iQrad)+feautrier_double(f)
!print*,'her',maxval(feautrier_double(f)),maxval(Qrad(l1:l2,m1:m2,n1:n2)),maxval(-Srad)
!
    endsubroutine radtransfer
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
      real, dimension (mx,my,mz,mvar + maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real :: formfactor=0.5
!
!  Add radiative cooling
!
      do n=n1,n2
      do m=m1,m2
         if(.not. nocooling) then
            df(l1:l2,m,n,ient)=df(l1:l2,m,n,ient) &
                              +4.*pi*f(l1:l2,m,n,ikappa) &
                               *f(l1:l2,m,n,iQrad) &
                               /f(l1:l2,m,n,iTT)*formfactor
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
!!$
!!$print*,'output_Qrad,lun=',output_Qrad,lun,size(Qrad)
!!$
!!$      if(lroot.and.headt) print*,'output_radiation',Qrad(4,4,4)
!!$      if(output_Qrad) write(2) 'help'
!!$      if(output_Qrad) write(90) 'help'
!!$!
    endsubroutine output_radiation
!***********************************************************************
    subroutine init_rad(f,xx,yy,zz)
!
!  initialise radiation; called from start.f90
!  We have an init parameter (initrad) to stear radiation i.c. independently.
!
!   15-jul-2002/nils: coded
!
      use Cdata
      use Mpicomm
      use Sub
      use Initcond
!
      real, dimension (mx,my,mz,mvar + maux) :: f
      real, dimension (mx,my,mz)      :: xx,yy,zz
      real :: nr1,nr2
      integer :: l12
!
      select case(initrad)

      case('zero', '0') 
         f(:,:,:,iQrad     ) = 1.
      case('gaussian-noise','1'); call gaunoise(amplee,f,iQrad)
      case('equil','2'); call init_equil(f)
      case ('cos', '3')
         f(:,:,:,iQrad) = -amplee*(cos(sqrt(3.)*0.5*xx)*(xx-Lx/2)*(xx+Lx/2)-1)
      case ('step', '4')
         l12=(l1+l2)/2
         f(1    :l12,:,:,iQrad) = 1.
         f(l12+1: mx,:,:,iQrad) = 2.
      case ('substep', '5')
         l12=(l1+l2)/2
         nr1=1.
         nr2=2.
         f(1    :l12-2,:,:,iQrad) = nr1
         f(l12-1      ,:,:,iQrad) = ((nr1+nr2)/2+nr1)/2
         f(l12+0      ,:,:,iQrad) = (nr1+nr2)/2
         f(l12+1      ,:,:,iQrad) = ((nr1+nr2)/2+nr2)/2
         f(l12+2: mx  ,:,:,iQrad) = nr2
      case ('lamb', '6')
         f(:,:,:,iQrad) = 2+(sin(2*pi*xx)*sin(2*pi*zz))
      case default
        !
        !  Catch unknown values
        !
        if (lroot) print*, 'No such such value for initrad: ', trim(initrad)
        call stop_it("")
      endselect
!
!  Pertubations
!
      select case(pertee)
         
      case('none', '0') 
      case('left','1')
         l12=(l1+l2)/2
         f(l1:l12,m1:m2,n1:n2,iQrad) = ampl_pert*f(l1:l12,m1:m2,n1:n2,iQrad)
      case('whole','2')
         f(:,m1:m2,n1:n2,iQrad) = ampl_pert*f(:,m1:m2,n1:n2,iQrad)
      case('ent','3') 
         !
         !  For perturbing the entropy after haveing found the 
         !  equilibrium between radiation and entropy.
         !
         f(:,m1:m2,n1:n2,ient) = ampl_pert
         f(:,m1:m2,n1:n2,ilnrho) = ampl_pert
      case default
         !
         !  Catch unknown values
         !
         if (lroot) print*, 'No such such value for pertee: ', trim(pertee)
         call stop_it("")
         
      endselect
!
      if(ip==0) print*,yy !(keep compiler quiet)
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
    subroutine init_equil(f)
!
!  Routine for calculating equilibrium solution of radiation
!
!  18-jul-02/nils: coded
!
      use Cdata
      use Ionization
!
      real, dimension(mx,my,mz,mvar + maux) :: f
      real, dimension(nx) :: lnrho,ss,source,TT1,cs2,cp1tilde
!
!  Use the ionization module to calculate temperature
!  At the moment we don't calculate ghost zones (ok for vertical arrays)  
!
      do n=n1,n2
      do m=m1,m2
         lnrho=f(l1:l2,m,n,ilnrho)
         ss=f(l1:l2,m,n,ient)
         call thermodynamics(lnrho,ss,cs2,TT1,cp1tilde)
         source=sigmaSB/(pi*TT1**4)
         f(l1:l2,m,n,iQrad) = source
      enddo
      enddo
!
    end subroutine init_equil
!***********************************************************************
endmodule Radiation
