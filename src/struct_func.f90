! $Id: struct_func.f90,v 1.5 2003-01-07 12:59:46 nilshau Exp $
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
    subroutine structure(f)
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
  integer, parameter :: qmax=8, imax=nx/2
  integer, parameter :: n_pdf=101
  real, dimension (mx,my,mz,mvar) :: f
  real, dimension (nx,ny,nz,3) :: u_vec,b_vec,z1_vec,z2_vec
  real, dimension (nx) :: bb 
  real, dimension (imax,3,qmax,3) :: sf,sf_sum,sfz1,sfz1_sum,sfz2,sfz2_sum
  real, dimension (imax) :: totall 
  real, dimension (ny,nz,3) :: du,dz1,dz2
  real, dimension(n_pdf,imax,3,3) :: p_du,p_du_sum
  real, dimension(n_pdf) :: x_du
  integer, dimension (ny,nz,3) :: i_du
  integer :: l,ll,j,q,direction
  integer :: separation,i,ivec,im,in
  real :: pdf_max,pdf_min,normalization,dx_du
  character (len=4) :: var
  character (len=20):: filetowriteu,filetowritez1,filetowritez2
  ! 
  !
  !
  if (iproc==root) print*,'Doing structure functions'
  do ivec=1,3
     !
     u_vec(:,:,:,ivec)=f(l1:l2,m1:m2,n1:n2,iuu+ivec-1)
     do n=n1,n2
        do m=m1,m2
           call curli(f,iaa,bb,ivec)
           im=m-nghost
           in=n-nghost
           b_vec(:,im,in,ivec)=bb
        enddo
     enddo
     z1_vec(:,:,:,ivec)=u_vec(:,:,:,ivec)+b_vec(:,:,:,ivec)
     z2_vec(:,:,:,ivec)=u_vec(:,:,:,ivec)-b_vec(:,:,:,ivec)
  enddo
  !
  !  Setting some variables depending on wether we want to
  !  calculate pdf or structure functions.
  !
  if (lpdf) then 
     pdf_max= 1.  !(for the time being; assumes |u|<1)
     pdf_min=-pdf_max
     dx_du=(pdf_max-pdf_min)/n_pdf
     do l=1,n_pdf
        x_du(l)=(l-.5)*dx_du+pdf_min
     enddo
     p_du=0.
  endif
  !
  !  Initialize structure function
  !
  if (lsf) then 
     sf=0.
     totall=0
     if (nr_directions .eq. 1) filetowriteu='/sf_sum_'
     if (nr_directions .ne. 1) filetowriteu='/sf_sum_transp_'           
     if (nr_directions .eq. 1) filetowritez1='/sfz1_sum_'
     if (nr_directions .ne. 1) filetowritez1='/sfz1_sum_transp_'           
     if (nr_directions .eq. 1) filetowritez2='/sfz2_sum_'
     if (nr_directions .ne. 1) filetowritez2='/sfz2_sum_transp_'           
  endif
  !
  !  Beginning the loops
  !
  do direction=1,nr_directions
     do l=1,nx
        !if (iproc==root) print*,'l=',l
        do ll=l+1,nx
           separation=min(mod(ll-l+nx,nx),mod(l-ll+nx,nx))
!           du=f(l,m1:m2,n1:n2,1:3)-f(ll,m1:m2,n1:n2,1:3)
           du=u_vec(l,:,:,:)-u_vec(ll,:,:,:)
           dz1=z1_vec(l,:,:,:)-z1_vec(ll,:,:,:)
           dz2=z2_vec(l,:,:,:)-z2_vec(ll,:,:,:)
           if (lpdf) then !if pdf=.true.
              i_du=1+int((du-pdf_min)*n_pdf/(pdf_max-pdf_min))
              i_du=min(max(i_du,1),n_pdf)  !(make sure its inside array bdries)
              !
              !  Calculating pdf
              !
              do j=1,3
                 do m=1,ny
                    do n=1,nz
                       p_du(i_du(m,n,j),separation,j,direction) &
                            =p_du(i_du(m,n,j),separation,j,direction)+1
                    enddo
                 enddo
              enddo
           endif
           !
           if (lsf) then
              !
              !  Calculates sf
              !
              totall(separation)=totall(separation)+1
              do j=1,3
                 do q=1,qmax                          
                    sf(separation,j,q,direction) &
                         =sf(separation,j,q,direction) &
                         +sum(abs(du(:,:,j))**q)
                    sfz1(separation,j,q,direction) &
                         =sfz1(separation,j,q,direction) &
                         +sum(abs(dz1(:,:,j))**q)
                    sfz2(separation,j,q,direction) &
                         =sfz2(separation,j,q,direction) &
                         +sum(abs(dz2(:,:,j))**q)
                 enddo
              enddo
           endif
        enddo
     enddo
     if (nr_directions .gt. 1) then
        if (direction .eq. 1) then
           !Doing transpose of y direction
           call transp(u_vec(:,:,:,1),'y')
           call transp(u_vec(:,:,:,2),'y')
           call transp(u_vec(:,:,:,3),'y')
           call transp(z1_vec(:,:,:,1),'y')
           call transp(z1_vec(:,:,:,2),'y')
           call transp(z1_vec(:,:,:,3),'y')
           call transp(z2_vec(:,:,:,1),'y')
           call transp(z2_vec(:,:,:,2),'y')
           call transp(z2_vec(:,:,:,3),'y')
        endif
        if (direction .eq. 2) then
           !Doing transpose of z direction
           call transp(u_vec(:,:,:,1),'z')
           call transp(u_vec(:,:,:,2),'z')
           call transp(u_vec(:,:,:,3),'z')
           call transp(z1_vec(:,:,:,1),'z')
           call transp(z1_vec(:,:,:,2),'z')
           call transp(z1_vec(:,:,:,3),'z')
           call transp(z2_vec(:,:,:,1),'z')
           call transp(z2_vec(:,:,:,2),'z')
           call transp(z2_vec(:,:,:,3),'z')
        endif
     endif
  enddo
  !
  !  Collecting all data on root processor and normalizing pdf and sf
  !
  if(lpdf) then
     call mpireduce_sum(p_du,p_du_sum,n_pdf*imax*3*3)  !Is this safe???
     do i=1,imax
        do j=1,3
           do direction=1,nr_directions
              normalization=1./(n_pdf*dx_du*sum(p_du_sum(:,i,j,direction)))
              p_du_sum(:,i,j,direction)=normalization*p_du_sum(:,i,j,direction)
           enddo
        enddo
     enddo
  endif
  !
  if(lsf) then
     call mpireduce_sum(sf,sf_sum,imax*3*qmax*3)  !Is this safe???
     sf_sum=sf_sum/(nw*ncpus)
     sf_sum(imax,:,:,:)=2*sf_sum(imax,:,:,:)
     call mpireduce_sum(sfz1,sfz1_sum,imax*3*qmax*3)  !Is this safe???
     sfz1_sum=sfz1_sum/(nw*ncpus)
     sfz1_sum(imax,:,:,:)=2*sfz1_sum(imax,:,:,:)
     call mpireduce_sum(sfz2,sfz2_sum,imax*3*qmax*3)  !Is this safe???
     sfz2_sum=sfz2_sum/(nw*ncpus)
     sfz2_sum(imax,:,:,:)=2*sfz2_sum(imax,:,:,:)
  endif
  !
  !  Writing output file
  !
  if (iproc==root) then
     do j=1,3              
        call chn(j,var)
        if(lpdf) then
           if (ip<10) print*,'Writing pdf of variable',var &
                ,'to ',trim(datadir)//'/pdf_'//trim(var)//'.dat'
           open(1,file=trim(datadir)//'/pdf_'//trim(var) &
                //'.dat',position='append')
           write(1,*) t,n_pdf
           write(1,'(1p,8e10.2)') p_du_sum(:,:,j,:)
           write(1,'(1p,8e10.2)') x_du
           close(1)
        endif
        !
        if(lsf) then
           if (ip<10) print*,'Writing structure functions of variable',var &
                ,'to ',trim(datadir)//trim(filetowriteu)//trim(var)//'.dat'
           open(1,file=trim(datadir)//trim(filetowriteu)//trim(var) &
                //'.dat',position='append')
           write(1,*) t,qmax
           write(1,'(1p,8e10.2)') sf_sum(:,j,:,:)
           close(1)
!
           if (ip<10) print*,'Writing structure functions of variable',var &
                ,'to ',trim(datadir)//trim(filetowritez1)//trim(var)//'.dat'
           open(1,file=trim(datadir)//trim(filetowritez1)//trim(var) &
                //'.dat',position='append')
           write(1,*) t,qmax
           write(1,'(1p,8e10.2)') sfz1_sum(:,j,:,:)
           close(1)
!
           if (ip<10) print*,'Writing structure functions of variable',var &
                ,'to ',trim(datadir)//trim(filetowritez2)//trim(var)//'.dat'
           open(1,file=trim(datadir)//trim(filetowritez2)//trim(var) &
                //'.dat',position='append')
           write(1,*) t,qmax
           write(1,'(1p,8e10.2)') sfz2_sum(:,j,:,:)
           close(1)
        endif
     enddo
  endif
  !
end subroutine structure
!***********************************************************************
end module struct_func
