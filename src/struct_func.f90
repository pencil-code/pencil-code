! $Id: struct_func.f90,v 1.12 2003-01-21 12:42:28 nilshau Exp $
!
!  Calculates 2-point structure functions and/or PDFs
!  and saves them during the run.
!
!  For the time being, the structure functions (or PDFs) are
!  called from power, so the output frequency is set by dspec.
!
!  The save files are under data/proc# under the names
!  sfz1_sum_ or sfz1_sum_transp_ .
!
!-----------------------------------------------------------------------
!   23-dec-02/nils: adapted from postproc/src/struct_func_mpi.f90
!

module struct_func
  !
  implicit none
  !
  contains

!***********************************************************************
    subroutine structure(f,ivec,b_vec,variabl)
!
!  The following parameters may need to be readjusted:
!  qmax should be set to the largest moment to be calculated
!  n_pdf gives the number of bins of the PDF
!
!   23-dec-02/nils: adapted from postproc/src/struct_func_mpi.f90
!   28-dec-02/axel: need also n_pdf in normalization
!
  use Cdata
  use Sub
  use General
  use Mpicomm
  !
  implicit none
  !
  integer, parameter :: qmax=8, imax=lb_nxgrid
  integer, parameter :: n_pdf=101
  real, dimension (mx,my,mz,mvar) :: f
  real, dimension (nx,ny,nz) :: vect,b_vec
  real, dimension (imax,qmax,3) :: sf,sf_sum
  real, dimension (ny,nz) :: dvect1,dvect2
  real, dimension(n_pdf,imax,3) :: p_du,p_du_sum
  real, dimension(n_pdf) :: x_du
  integer, dimension (ny,nz) :: i_du1,i_du2
  integer :: l,ll,j,q,direction,ll1,ll2
  integer :: i,ivec,lb_ll,separation
  real :: pdf_max,pdf_min,normalization,dx_du
  character (len=4) :: var
  character (len=*) :: variabl
  character (len=20):: filetowrite
  logical :: llsf=.false., llpdf=.false.
  ! 
  !
  !
  if (iproc==root) print*,'Doing structure functions'
  !
  if (variabl .eq. 'u') then
     vect(:,:,:)=f(l1:l2,m1:m2,n1:n2,iuu+ivec-1)
     filetowrite='/sfu-'
     sf=0.
     llsf=.true.
     llpdf=.false.
  elseif (variabl .eq. 'b') then
     vect(:,:,:)=b_vec(:,:,:)
     filetowrite='/sfb-'
     sf=0.
     llsf=.true.
     llpdf=.false.
  elseif (variabl .eq. 'z1') then
     vect(:,:,:)=f(l1:l2,m1:m2,n1:n2,iuu+ivec-1)+b_vec(:,:,:)
     filetowrite='/sfz1-'
     sf=0.
     llsf=.true.
     llpdf=.false.
  elseif (variabl .eq. 'z2') then
     vect(:,:,:)=f(l1:l2,m1:m2,n1:n2,iuu+ivec-1)-b_vec(:,:,:)
     filetowrite='/sfz2-'
     sf=0.
     llsf=.true.
     llpdf=.false.
  end if
  !
  !  Setting some variables depending on wether we want to
  !  calculate pdf or structure functions.
  !
  if (variabl .eq. 'pdf') then 
     vect(:,:,:)=f(l1:l2,m1:m2,n1:n2,iuu+ivec-1)
     pdf_max= 1.  !(for the time being; assumes |u|<1)
     pdf_min=-pdf_max
     dx_du=(pdf_max-pdf_min)/n_pdf
     do l=1,n_pdf
        x_du(l)=(l-.5)*dx_du+pdf_min
     enddo
     p_du=0.
     llpdf=.true.
     llsf=.false.
  endif
  !
  !  Beginning the loops
  !
  do direction=1,nr_directions
     do l=1,nx
        if ((iproc==root) .and. (lpostproc)) print*,'l=',l
        do lb_ll=1,lb_nxgrid
           !separation=min(mod(ll-l+nx,nx),mod(l-ll+nx,nx))
           separation=2**(lb_ll-1)
           ll1=mod(l+separation-1,nx)+1
           ll2=mod(l-separation+nx-1,nx)+1
           dvect1=vect(l,:,:)-vect(ll1,:,:)
           dvect2=vect(l,:,:)-vect(ll2,:,:)
           if (llpdf) then !if pdf=.true.
              i_du1=1+int((dvect1-pdf_min)*n_pdf/(pdf_max-pdf_min))
              i_du1=min(max(i_du1,1),n_pdf)  !(make sure its inside array bdries)
              i_du2=1+int((dvect2-pdf_min)*n_pdf/(pdf_max-pdf_min))
              i_du2=min(max(i_du2,1),n_pdf)  !(make sure its inside array bdries)
              !
              !  Calculating pdf
              !
              do m=1,ny
                 do n=1,nz
                    p_du(i_du1(m,n),lb_ll,direction) &
                         =p_du(i_du1(m,n),lb_ll,direction)+1
                    p_du(i_du2(m,n),lb_ll,direction) &
                         =p_du(i_du2(m,n),lb_ll,direction)+1
                 enddo
              enddo
           endif
           !
           if (llsf) then
              !
              !  Calculates sf
              !
              do q=1,qmax   
                 sf(lb_ll,q,direction) &
                      =sf(lb_ll,q,direction) &
                      +sum(abs(dvect1(:,:))**q)+sum(abs(dvect2(:,:))**q)
              enddo
           endif
        enddo
     enddo
     if (nr_directions .gt. 1) then
        if (direction .eq. 1) then
           !Doing transpose of y direction
           call transp(vect(:,:,:),'y')
        endif
        if (direction .eq. 2) then
           !Doing transpose of z direction
           call transp(vect(:,:,:),'z')
        endif
     endif
  enddo
  !
  !  Collecting all data on root processor and normalizing pdf and sf
  !
  if(llpdf) then
     call mpireduce_sum(p_du,p_du_sum,n_pdf*imax*3)  !Is this safe???
     do i=1,imax
        do direction=1,nr_directions
           normalization=1./(n_pdf*dx_du*sum(p_du_sum(:,i,direction)))
           p_du_sum(:,i,direction)=normalization*p_du_sum(:,i,direction)
        enddo
     enddo
  endif
  !
  if(llsf) then
     call mpireduce_sum(sf,sf_sum,imax*qmax*3)  !Is this safe???
     sf_sum=sf_sum/(nw*ncpus*2)
     !sf_sum(imax,:,:)=2*sf_sum(imax,:,:)
  endif
  !
  !  Writing output file
  !
  if (iproc==root) then
     call chn(ivec,var)
     if(llpdf) then
        if (ip<10) print*,'Writing pdf of variable',var &
             ,'to ',trim(datadir)//'/pdf-'//trim(var)//'.dat'
        open(1,file=trim(datadir)//'/pdf-'//trim(var) &
             //'.dat',position='append')
        write(1,*) t,n_pdf
        write(1,'(1p,8e10.2)') p_du_sum(:,:,:)
        write(1,'(1p,8e10.2)') x_du
        close(1)
     endif
     !
     if(llsf) then
        if (ip<10) print*,'Writing structure functions of variable',var &
             ,'to ',trim(datadir)//trim(filetowrite)//trim(var)//'.dat'
        open(1,file=trim(datadir)//trim(filetowrite)//trim(var) &
             //'.dat',position='append')
        write(1,*) t,qmax
        write(1,'(1p,8e10.2)') sf_sum(:,:,:)
        close(1)
     endif
  endif
  !
end subroutine structure
!***********************************************************************
end module struct_func
