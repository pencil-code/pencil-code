! $Id: reflect.f90,v 1.1 2011-04-04 12:03:22 piyali Exp $
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!***************************************************************
!
!
module InitialCondition
!
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include '../initial_condition.h'
!
  real :: r1=1.0d36,r2=1.0d36,cs2top=impossible
  real :: rescale_ff=1.
  logical :: lreflectua=.false.,lperturb=.false.
  real, dimension (mx,my,mz,3),save :: aa_save
  real, dimension (mx,my,mz),save :: ss_save,rho_save

  namelist /initial_condition_pars/ &
           r1,r2,lreflectua,cs2top,rescale_ff,lperturb
!
  contains
!***********************************************************************
    subroutine register_initial_condition()
!
!  Register variables associated with this module; likely none.
!
!  07-may-09/wlad: coded
!
      if (lroot) call svn_id( &
           "$Id: reflect.f90,v 1.1 2011-04-04 12:03:22 piyali Exp $")
!
    endsubroutine register_initial_condition
!***********************************************************************
    subroutine initialize_initial_condition(f)

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  Initialize any module variables which are parameter dependent.
!
!  07-may-09/wlad: coded
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine initial_condition_uu(f)
!
!  Initialize the velocity field.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (3) :: tmpvec
      real, dimension (mx,my,mz,mfarray) :: fcopy
!
!  SAMPLE IMPLEMENTATION
!
          fcopy(:,:,:,iuu:iuu+2)=rescale_ff*f(:,:,:,iuu:iuu+2)
! copy the magnetic vector potential too. 
          if(lperturb) then
            aa_save(:,:,:,1:3)=rescale_ff*f(:,:,:,iaa:iaa+2)
            rho_save(:,:,:)=rescale_ff*f(:,:,:,ilnrho)
            ss_save(:,:,:)=rescale_ff*f(:,:,:,iss)
            f(:,:,:,iaa:iaa+2) = 0.
            f(:,:,:,ilnrho) = 0.
            f(:,:,:,iss) = 0.
          endif 
          if (lreflectua) then
           do iz=1,mz; do iy=1, my/2; do ix=1, mx
             tmpvec = f(ix,my-iy+1,iz,iux:iuz)
             f(ix,my-iy+1,iz,iux)= f(ix,iy,iz,iux)
             f(ix,iy,iz,iux)=tmpvec(1)
             f(ix,my-iy+1,iz,iuy)= -f(ix,iy,iz,iuy)
             f(ix,iy,iz,iuy)=-tmpvec(2)
             f(ix,my-iy+1,iz,iuz)= f(ix,iy,iz,iuz)
             f(ix,iy,iz,iuz)=tmpvec(3)
            enddo; enddo; enddo
          endif

      f(:,:,:,iuu:iuu+2)=r1*rescale_ff*f(:,:,:,iuu:iuu+2)+r2*fcopy(:,:,:,iuu:iuu+2)
!
    endsubroutine initial_condition_uu
!***********************************************************************
    subroutine initial_condition_lnrho(f)
!
!  Initialize logarithmic density. init_lnrho will take care of
!  converting it to linear density if you use ldensity_nolog.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx,my,mz,mfarray) :: fcopy
      real :: tmpscl

      if(lperturb) then 
          fcopy(:,:,:,ilnrho)=rho_save(:,:,:)
          if (lreflectua) then
           do iz=1,mz; do iy=1, my/2; do ix=1, mx
             tmpscl = rho_save(ix,my-iy+1,iz)
             rho_save(ix,my-iy+1,iz)= rho_save(ix,iy,iz)
             rho_save(ix,iy,iz)=tmpscl
            enddo; enddo; enddo
          endif
        f(:,:,:,ilnrho)= f(:,:,:,ilnrho)+ & 
           r1*rho_save(:,:,:)+r2*fcopy(:,:,:,ilnrho)
     else
          fcopy(:,:,:,ilnrho)=rescale_ff*f(:,:,:,ilnrho)
          if (lreflectua) then
           do iz=1,mz; do iy=1, my/2; do ix=1, mx
             tmpscl = f(ix,my-iy+1,iz,ilnrho)
             f(ix,my-iy+1,iz,ilnrho)= f(ix,iy,iz,ilnrho)
             f(ix,iy,iz,ilnrho)=tmpscl
            enddo; enddo; enddo
          endif
         f(:,:,:,ilnrho)=log(r1*exp(rescale_ff*f(:,:,:,ilnrho))+r2*exp(fcopy(:,:,:,ilnrho)))
     endif
!
    endsubroutine initial_condition_lnrho
!***********************************************************************
    subroutine initial_condition_ss(f)
!
!  Initialize entropy.
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx,my,mz,mfarray) :: fcopy
      real :: tmpscl

      if(lperturb) then 
          fcopy(:,:,:,iss)=ss_save(:,:,:)
          if (lreflectua) then
           do iz=1,mz; do iy=1, my/2; do ix=1, mx
             tmpscl = ss_save(ix,my-iy+1,iz)
             ss_save(ix,my-iy+1,iz)= ss_save(ix,iy,iz)
             ss_save(ix,iy,iz)=tmpscl
            enddo; enddo; enddo
          endif
        f(:,:,:,iss)= f(:,:,:,iss)+ & 
           r1*ss_save(:,:,:)+r2*fcopy(:,:,:,iss)
      else  
        fcopy(:,:,:,iss)=rescale_ff*f(:,:,:,iss)
          if (lreflectua) then
           do iz=1,mz; do iy=1, my/2; do ix=1, mx
             tmpscl = f(ix,my-iy+1,iz,iss)
             f(ix,my-iy+1,iz,iss)= f(ix,iy,iz,iss)
             f(ix,iy,iz,iss)=tmpscl
            enddo; enddo; enddo
          endif
        f(:,:,:,iss)=r1*rescale_ff*f(:,:,:,iss)+r2*fcopy(:,:,:,iss)
      endif
      if (lroot) print *, fcopy(64,64,133,iss),f(64,64,133,iss)

    endsubroutine initial_condition_ss
!********************************************************************
    subroutine initial_condition_aa(f)
!
! Initialize magnetic potential
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx,my,mz,mfarray) :: fcopy
      real, dimension (3) :: tmpvec

      if(lperturb) then 
          fcopy(:,:,:,iaa:iaa+2)=aa_save(:,:,:,1:3)
          if (lreflectua) then
           do iz=1,mz; do iy=1, my/2; do ix=1, mx
             tmpvec = aa_save(ix,my-iy+1,iz,1:3)
             aa_save(ix,my-iy+1,iz,1)= aa_save(ix,iy,iz,1)
             aa_save(ix,iy,iz,1)=tmpvec(1)
             aa_save(ix,my-iy+1,iz,2)= -aa_save(ix,iy,iz,2)
             aa_save(ix,iy,iz,2)=-tmpvec(2)
             aa_save(ix,my-iy+1,iz,3)= aa_save(ix,iy,iz,3)
             aa_save(ix,iy,iz,3)=tmpvec(3)
            enddo; enddo; enddo
          endif
          f(:,:,:,iaa:iaa+2)=f(:,:,:,iaa:iaa+2)+ & 
                     r1*aa_save(:,:,:,1:3)+r2*fcopy(:,:,:,iaa:iaa+2)
       else
          fcopy(:,:,:,iaa:iaa+2)=rescale_ff*f(:,:,:,iaa:iaa+2)
          if (lreflectua) then
           do iz=1,mz; do iy=1, my/2; do ix=1, mx
             tmpvec = f(ix,my-iy+1,iz,iax:iaz)
             f(ix,my-iy+1,iz,iax)= f(ix,iy,iz,iax)
             f(ix,iy,iz,iax)=tmpvec(1)
             f(ix,my-iy+1,iz,iay)= -f(ix,iy,iz,iay)
             f(ix,iy,iz,iay)=-tmpvec(2)
             f(ix,my-iy+1,iz,iaz)= f(ix,iy,iz,iaz)
             f(ix,iy,iz,iaz)=tmpvec(3)
            enddo; enddo; enddo
         endif
            f(:,:,:,iaa:iaa+2)=r1*rescale_ff*f(:,:,:,iaa:iaa+2)+r2*fcopy(:,:,:,iaa:iaa+2)
         endif
!
    endsubroutine initial_condition_aa
!********************************************************************
    subroutine read_initial_condition_pars(unit,iostat)
!
!  07-may-09/wlad: coded
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=initial_condition_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=initial_condition_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_initial_condition_pars
!***********************************************************************
    subroutine write_initial_condition_pars(unit)
!     
      integer, intent(in) :: unit
!
      write(unit,NML=initial_condition_pars)
!
    endsubroutine write_initial_condition_pars
!***********************************************************************
!
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from noinitial_condition.f90 for any    **
!**  InitialCondition routines not implemented in this file        **
!**                                                                **
    include '../initial_condition_dummies.inc'
!********************************************************************
endmodule InitialCondition
