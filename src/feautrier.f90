! $Id: feautrier.f90,v 1.41 2003-10-24 13:17:31 dobler Exp $

!!!  NOTE: this routine will perhaps be renamed to radiation_feautrier
!!!  or it may be combined with radiation_ray.

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 1
!
!***************************************************************

module Radiation

!  Radiation (solves transfer equation along rays)
!  The direction of the ray is given by the vector (lrad,mrad,nrad),
!  and the parameters radx0,rady0,radz0 gives the maximum number of
!  steps of the direction vector in the corresponding direction.

  use Cparam

  implicit none

  real, dimension (mx,my,mz) :: Srad,kaprho
  logical :: nocooling=.false.,test_radiation=.false.,output_Qrad=.false.
  logical :: lkappa_es=.false.
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
      if(.not. first) call stop_it('register_radiation: called twice')
      first = .false.
!
      lradiation = .true.
      lradiation_ray = .true.
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
           "$Id: feautrier.f90,v 1.41 2003-10-24 13:17:31 dobler Exp $")
!
! Check we aren't registering too many auxiliary variables
!
      if (naux > maux) then
        if (lroot) write(0,*) 'naux = ', naux, ', maux = ', maux
        call stop_it('register_radiation: naux > maux')
      endif
!
!  Writing files for use with IDL
!
      if (naux < maux) aux_var(aux_count)=',Qrad $'
      if (naux == maux) aux_var(aux_count)=',Qrad'
      aux_count=aux_count+1
      if (lroot) write(15,*) 'Qrad = fltarr(mx,my,mz)*one'
!
    endsubroutine register_radiation
!***********************************************************************
    subroutine initialize_radiation()
!
!  nothing to be done here
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
      real, dimension(nx) :: lnrho
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
!  Use the ionization module to calculate temperature
!  At the moment we don't calculate ghost zones (ok for vertical arrays)  
!
      do n=n1,n2
      do m=m1,m2
         Srad(l1:l2,m,n)=sourcefunction(f)
!
!  opacity: if lkappa_es then take electron scattering opacity only;
!  otherwise use Hminus opacity (but may need to add kappa_es as well).
!
         if (lkappa_es) then
            lnrho=f(l1:l2,m,n,ilnrho)
            kaprho(l1:l2,m,n)=kappa_es*exp(lnrho)
         else
            kaprho(l1:l2,m,n)=opacity(f)
         endif
      enddo
      enddo
!
    endsubroutine radcalc
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
      real, dimension(mx,my,mz,mvar+maux) :: f
      real, dimension(mx,my,mz) :: feautrier
      real, dimension(nz) :: kaprho_,tau,Srad_,Prad_
      real, dimension(nz) :: a,b,c
      integer :: lrad,mrad,nrad
      logical :: err=.false.
!
!  loop fastest over first index --> could do this for l-vectors
!
      do mrad=m1,m2
      do lrad=l1,l2
         kaprho_=kaprho(lrad,mrad,n1:n2)
         tau=spline_integral(z,kaprho_)
         Srad_=Srad(lrad,mrad,n1:n2)
!
!  top boundary: P'=P, together with P1-P0=dtau*P'+.5*dtau^2*P" and P"=P-S
!
         b(1)=1.+2./(tau(2)-tau(1))+2./(tau(2)-tau(1))**2
         c(1)=-2./(tau(2)-tau(1))**2
!
!  lower boundary: P=S
!
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
            print*,'feautrier: lnrho=',f(lrad,mrad,n1:n1+5,ilnrho),'...', &
                            f(lrad,mrad,n2-5:n2,ilnrho)
            print*,'feautrier: ss=',f(lrad,mrad,n1:n1+5,iss),'...', &
                         f(lrad,mrad,n2-5:n2,iss)
            print*,'feautrier: tau=',tau(1:6),'...',tau(n2-n1-5:n2-n1)
            stop
         endif
!
         feautrier(lrad,mrad,n1:n2)=Prad_
!         print*,'feautrier: Prad',Prad_
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
      real, dimension(mx,my,mz,mvar+maux) :: f
      real, dimension(mx,my,mz) :: feautrier_double
      double precision, dimension(nz) :: kaprho_,tau,Srad_,Prad_
      double precision, dimension(nz) :: a,b,c
      integer :: lrad,mrad,nrad
      logical :: err=.false.
!
!  loop fastest over first index --> could do this for l-vectors
!
      do mrad=m1,m2
      do lrad=l1,l2
         do nrad=n1,n2
           kaprho_(nrad)=kaprho(lrad,mrad,n1+n2-nrad)
           Srad_(nrad)=Srad(lrad,mrad,n1+n2-nrad)
         enddo
         tau=spline_integral_double(z,-kaprho_)
!
!  top boundary: P'=P, together with P1-P0=dtau*P'+.5*dtau^2*P" and P"=P-S
!
         b(1)=1.+2./(tau(2)-tau(1))+2./(tau(2)-tau(1))**2
         c(1)=-2./(tau(2)-tau(1))**2
!
         !b(1)=1.+3./(tau(2)-tau(1))
         !c(1)=-4./(tau(2)-tau(1))
!
!  lower boundary: P=S
!
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
            print*,'feautrier_double: lnrho=',f(lrad,mrad,n1:n1+5,ilnrho),'...', &
                            f(lrad,mrad,n2-5:n2,ilnrho)
            print*,'feautrier_double: ss=',f(lrad,mrad,n1:n1+5,iss),'...', &
                         f(lrad,mrad,n2-5:n2,iss)
            print*,'feautrier_double: tau=',tau(1:6),'...',tau(n2-n1-5:n2-n1)
            print*,'feautrier_double: kappa=',f(lrad,mrad,n1:n1+5,ikappa),'...', &
                            f(lrad,mrad,n2-5:n2,ikappa)
            stop
         endif
!
         do nrad=n1,n2
            feautrier_double(lrad,mrad,nrad)=Prad_(n1+n2-nrad)
         enddo
      enddo
      enddo
!
    endfunction feautrier_double
!***********************************************************************
    function transx(chi,S)
!
!  Solve the transfer equation, given optical depth and source function.
!
!  The answer is q=p-s, where p is "Feautriers P".
!
!  "Steins trick" is used! storing the sum of the three elements in each
!  row of the tridiagonal equation system, instead of the diagonal ele-
!  ment.  This prevents the deterioration of the numerical precision for
!  small optical depths.
!
!  This is a second order version.  For simulations, with rapidly
!  varying absorption coefficients and source functions, this is to be
!  preferred over spline and Hermitean versions because it is positive
!  definite, in the sense that a positive source function is guaranteed
!  to result in a positive average intensity p.  Also, the flux
!  divergence is exactly equal to q, for the conventional definition
!  of the flux.
!
!  Operation count: 5d+5m+11a = 21 flops
!
!  Timings:
!            Alliant: 0.58 * 21 / 1.25 = 9.7 Mfl @ 31*31*31
!
!  Update history:
!
!  28-oct-87/aake: vector-concurrent version
!  03-nov-87/aake: signs moved around to simplify stores, comments added
!  05-nov-87/aake: first point is mmmz=1 again, was mmmz=2
!  05-nov-87/aake: tested if/then/else timing on Alliant
!  08-jan-88/aake: added surface intensity calculation
!
!***********************************************************************
!
      use Cdata
!
      real, dimension(nz), intent(in) :: chi,S
      real, dimension(nz) :: transx
      real, dimension(nz) :: sp1,sp2,sp3,Q
      real :: dinv,dtau2
!
!  k=nz
!
      dtau2=.5*(chi(nz-1)+chi(nz))*dz
      sp2(nz)=dtau2*(1.+.5*dtau2)
      sp3(nz)=-1.
      Q(nz)=S(nz-1)-(1.+dtau2)*S(nz)
!
!  2<k<nz [3d+2m+6a]
!
      do n=nz-1,2,-1
        dinv=4./(chi(n-1)+2.*chi(n)+chi(n+1))/dz
        sp1(n)=-2.*dinv/(chi(n)+chi(n+1))/dz
        sp2(n)=1.
        sp3(n)=-2.*dinv/(chi(n-1)+chi(n))/dz
        Q(n)=sp1(n)*(S(n)-S(n+1))+sp3(n)*(S(n)-S(n-1))
      enddo
!
!  k=1
!
      sp2(1)=1.
      Q(1)=0.
!
!  eliminate subdiagonal, save factors in sp1 [1d+2m+4a]
!
      do n=nz,3,-1
        sp1(n)=-sp1(n-1)/(sp2(n)-sp3(n))
        Q(n-1)=Q(n-1)+sp1(n)*Q(n)
        sp2(n-1)=sp2(n-1)+sp1(n)*sp2(n)
        sp2(n)=sp2(n)-sp3(n)
      enddo
      sp2(2)=sp2(2)-sp3(2)
!
!  backsubstitute [1d+1m+1a]
!
      do n=2,nz
        Q(n)=(Q(n)-sp3(n)*Q(n-1))/sp2(n)
      enddo
!
!  return value
!
      transx=Q
    endfunction transx
!***********************************************************************
    subroutine radtransfer1(f)
!
!  Integration radiation transfer equation along rays
!
!  24-mar-03/axel+tobi: coded
!
      use Cdata
!
      real, dimension(mx,my,mz,mvar+maux) :: f
      integer :: l
!
!  identifier
!
      if(lroot.and.headt) print*,'radtransfer1: ENTER'
!
      call radcalc(f)
      do l=l1,l2
      do m=m1,m2
        f(l,m,n1:n2,iQrad)=transx(kaprho(l,m,n1:n2),Srad(l,m,n1:n2))
      enddo
      enddo
!
    endsubroutine radtransfer1
!***********************************************************************
    subroutine radtransfer2(f)
!
!  dummy routine
!
      use Cdata
      use Sub
!
      real, dimension(mx,my,mz,mvar+maux) :: f
!
      if(ip==0) print*,f !(keep compiler quiet)
    endsubroutine radtransfer2
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
            df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss) &
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
!!$print*,'output_radiation: lun=',output_Qrad,lun,size(Qrad)
!!$
!!$      if(lroot.and.headt) print*,'output_radiation: Qrad(4,4,4) =',Qrad(4,4,4)
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
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz)      :: xx,yy,zz
      real :: nr1,nr2
      integer :: l12
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
    subroutine rprint_radiation(lreset,lwrite)
!
!  Writes iQrad, etc, to index.pro file
!  Also: dummy routine for Flux Limited Diffusion routine
!  reads and registers print parameters relevant for radiative part
!
!  16-jul-02/nils: adapted from rprint_hydro
!  14-jun-03/axel: moved iTT to rprint_ionization
!
      use Cdata
      use Sub
!  
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  write column where which radiative variable is stored
!
      if (lwr) then
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
      endif
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
      if (ip==0) print*,topbot,f(1,1,1,1)  !(to keep compiler quiet)
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
      if (ip==0) print*,topbot,f(1,1,1,1)  !(to keep compiler quiet)
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
      real, dimension(mx,my,mz,mvar+maux) :: f
      real, dimension(nx) :: lnrho,ss,yH,TT
!
!  Use the ionization module to calculate temperature
!  At the moment we don't calculate ghost zones (ok for vertical arrays)  
!
      if(lionization) call ioncalc(f)
      do n=n1,n2
      do m=m1,m2
         lnrho=f(l1:l2,m,n,ilnrho)
         ss=f(l1:l2,m,n,iss)
         call ionset(f,ss,lnrho,yH,TT)
         f(l1:l2,m,n,iQrad)=sigmaSB*TT**4/pi
      enddo
      enddo
!
    end subroutine init_equil
!***********************************************************************
endmodule Radiation
