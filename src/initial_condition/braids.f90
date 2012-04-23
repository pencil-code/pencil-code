!  Initial condition (density, magnetic field, velocity) 
!  for a particular configuration of a braided magnetic field.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!***************************************************************
module InitialCondition
!
  use Cdata
  use Cparam
  use Messages
  use Streamlines
  use Fixed_point
  use Sub, only: keep_compiler_quiet  
  use Boundcond ! for the core boundary communication
  use mpicomm
!
  implicit none
!
  include '../initial_condition.h'
!
! ampl = amplitude of the magnetic field
! width_tube = width of the flux tube
! braid_margin = margin of outer most strands to the borders
! l_sigma = length of the twist region
! steepnes = steepnes of the braiding
! B_bkg = strength of the background field in z-direction
! word = sequence of the braid group
! prof = the amplitude profile across the tube
!
! n_blobs = number of blobs for the blob configuration
! xc, yc, zc = position of the blobs
! blob_sgn = sign of the twist in the blob
! l_blob = length in z-direction of the blob
! blob_scale = scaling factor for the Gaussian
!
  real :: ampl = 1.0, width_tube = 0.3, braid_margin = -1., l_sigma = 0.3, steepnes = 1.0  
  real :: B_bkg = 0.0
  character (len=labellen) :: prof='gaussian'
  character (len=labellen) :: word = "AA"
! variables for the blob configurations  
  integer :: n_blobs = 0
  real, dimension (6) :: xc, yc, zc, blob_sgn, l_blob
  real :: blob_scale = 1.
!
  namelist /initial_condition_pars/ &
    ampl,width_tube,braid_margin,l_sigma,steepnes,B_bkg,word,prof, &
!   blob variables
    n_blobs,xc,yc,zc,blob_sgn,l_blob,blob_scale
!
  contains
!***********************************************************************
  subroutine register_initial_condition()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  01-july-11/simon: coded
!
!  Identify CVS/SVN version information.
!
    if (lroot) call svn_id( &
        "$Id: braids.f90,v 1.9 2011-08-02 16:43:18 iomsn Exp $")
!
  endsubroutine register_initial_condition
!***********************************************************************
  subroutine initial_condition_uu(f)
!
!  Initialize the velocity field.
!
!  07-may-09/wlad: coded
!
    real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
    call keep_compiler_quiet(f)
!
  endsubroutine initial_condition_uu
!***********************************************************************
  subroutine initial_condition_lnrho(f)
!
!  Initialize logarithmic density. init_lnrho 
!  will take care of converting it to linear 
!  density if you use ldensity_nolog
!
!  07-may-09/wlad: coded
!
    real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
    call keep_compiler_quiet(f)
!
  endsubroutine initial_condition_lnrho
!***********************************************************************
  subroutine initial_condition_aa(f)
!
!  Initialize the magnetic vector potential.
!
!  01-july-11/simon: coded
!
!  Braided magnetic flux tube starting from the lower xy-plane and
!  ending at the top plane.
!
!  Created 2011-07-01 by Simon Candelaresi (Iomsn)
!
    use Mpicomm, only: stop_it
    use Poisson
    use Sub
!    
    real, dimension (mx,my,mz,mfarray) :: f
!        
    real :: distance_tubes ! distance between the flux tubes
    real :: l_straight ! length of the regions with a straight magnetic field
    real :: a(5) ! the coefficients for the rotation polynomial
!    
    real :: phi ! the rotation angle for the braid
    real :: phi_offset ! auxiliary variable for the rotation
    real :: z_start ! auxiliary variable for the creation of the braid
    integer :: rotation_sign ! sign of the rotation in the braid
    real :: rotation_center(2) ! the center of the rotation for the braid
!    
    integer :: word_len, idx, idx_strand
    character (len=len_trim(word)) :: wordn ! the trimmed word
!
!   keeps track of the current position of the strand
    integer, dimension(9) :: strand_position
!   position of the other strand in the braid
    integer :: other_position, ierr
    integer :: braid_num ! auxiliary variable for the braid position
!
!   step length for some curve parameters
    real :: delta_tube_param, delta_circle_param, delta_circle_radius
    real :: tube_param, circle_param, circle_radius
!   vectors for the creation of the flux tubes
    real, dimension(3) :: tube_pos, circle_pos, tangent, normal
!   working variables for the strands
    integer :: n_strands, ascii_code
    real, dimension(3) :: new_bb
!
!   The next variables are used for the uncurling.
    integer :: l, j, ju, k
    real, dimension (nx,ny,nz,3) :: jj, tmpJ  ! This is phi for poisson.f90
!
!   In case field line tracing is applied, use this array.
!   The last dimension is used for the following:
!   2 for the initial seed position (x0,y0)
!   3 for the current position in 3d space
!   1 for the total length
!   1 for the integrated quantity
!     real, dimension (nx*ny,7) :: tracers
    real, pointer, dimension (:,:) :: tracers
!   the traced field
    real, pointer, dimension (:,:,:,:) :: vv
!   filename for the tracer and fixed point output
    character(len=1024) :: filename, str_tmp
!
!
    if (trace_field == 'bb' .and. ipz == 0) then
!     allocate memory for the traced field
      allocate(vv(nx,ny,nz,3))
    endif
!
!   check the word
    wordn = word
    word_len = len(wordn)
    if (verify(wordn, 'ABCDEFGHabcdefgh') >= word_len) then
      call fatal_error("braids", "invalid word. Note that the highest braid is 'h'")
    endif
!
!   determine the number of strands
    n_strands = 0
    do idx = 1, word_len
      ascii_code = ichar(wordn(idx:idx))
      if ((ascii_code >= 65) .and. (ascii_code <= 72)) then
        if (n_strands <= ascii_code - 64) then
          n_strands = ascii_code - 64
        endif
      elseif ((ascii_code >= 97) .and. (ascii_code <= 104)) then
        if (n_strands <= ascii_code - 96) then
          n_strands = ascii_code - 96
        endif
      endif
    enddo
    n_strands = n_strands + 1
!
!   compute the distance between the braids
    l_straight = (Lz - word_len*l_sigma) / (word_len+1)
    if (l_straight <= 0.) then
      call fatal_error("braids", "distance between braids is negative, check l_sigma")
    endif
!
    if (braid_margin < 0) braid_margin = 0.5*Lx/n_strands
    if (braid_margin > Lx/2.) &
        call fatal_error("braids", "braid margin is larger then Lx/2")
!
    distance_tubes = (Lx-2*braid_margin)/(n_strands-1)
!
!   clear the magnetic field to zero
!     f(l1:l2,m1:m2,n1:n2,iax:iaz) = 0.
!     f(1:l2+nghost,1:m2+nghost,1:n2+nghost,iax:iaz) = 0.
    f(:,:,:,iax:iaz) = 0.
!
!
!   set the coefficients for the rotation polynomial
    a(1) = steepnes
    a(2) = 0
    a(3) = -(a(1)*l_sigma + 5*l_sigma*a(1) - 15*pi) / (3*l_sigma*(l_sigma/2.)**2)
    a(4) = 0      
    a(5) = (l_sigma*a(1) - 3*pi) / (2*(l_sigma/2.)**5)
!
!   Calculate the minimum step size of the curve parameters
!   to avoid discretization issues, like mesh points without magnetic field
    delta_tube_param = min(dx, dy, dz)
!   correct for the braid steepnes
    delta_tube_param = delta_tube_param * l_sigma / (steepnes * pi * distance_tubes * 8)
    delta_circle_radius = delta_tube_param*8
    delta_circle_param = delta_circle_radius/(width_tube/2.)
!
!   loop over all strands
    if (n_blobs == 0) then
      do idx_strand = 1,n_strands
!       compute the initial position of this strand
        tube_pos(1) = x(l1)-ipx*nx*dx + (idx_strand-1)*distance_tubes + braid_margin
        tube_pos(2) = Ly/2. + y(m1)-ipy*ny*dy
        tube_pos(3) = z(n1)-ipz*nz*dz
        idx = 1
!
!       reset the strand_position vector
        strand_position = (/1,2,3,4,5,6,7,8,9/)
!
!       loop over all braids
        do
!         create straight lines
          do
            if (tube_pos(3) > (l_straight*idx + l_sigma*(idx-1) + z(n1)-ipz*nz*dz)) exit
            tangent = (/0,0,1/)
!
            circle_radius = 0.
!           loop which changes the circle's radius
            do
              if (circle_radius > width_tube/2.) exit
              circle_param = 0.
!             loop which goes around the circle
              do
                if (circle_param > 2.*pi) exit
                circle_pos(1) = tube_pos(1) + circle_radius*cos(circle_param)
                circle_pos(2) = tube_pos(2) + circle_radius*sin(circle_param)
                circle_pos(3) = tube_pos(3)
!
!               Find the corresponding mesh point to this position.
                l = nint((circle_pos(1) - x(1))/dx) + 1
                m = nint((circle_pos(2) - y(1))/dy) + 1
                n = nint((circle_pos(3) - z(1))/dz) + 1
!
                if ((l > mx .or. m > my .or. n > mz .or. l < 1 .or. m < 1 .or. n < 1) .eqv. .false.) then
!                 Write the magnetic field B.
!                 Note that B is written in the f-array where A is stored.
!                 This is corrected further in the code.
                    if (prof == 'gaussian') then
                      new_bb = tangent*ampl*(exp(-(2*circle_radius/width_tube)**2)-exp(-1.)) / (1-exp(-1.))
                    else if (prof == 'constant') then
                      new_bb = tangent*ampl
                    endif
!                   Avoid issues in spots with high curvature.
                    if ((f(l,m,n,iax)**2 + f(l,m,n,iay)**2 + f(l,m,n,iaz)**2) < &
                        (new_bb(1)**2 + new_bb(2)**2 + new_bb(3)**2)) &
                        f(l,m,n,iax:iaz) = new_bb
                  endif
                circle_param = circle_param + delta_circle_param
              enddo
              circle_radius = circle_radius + delta_circle_radius
            enddo
            tube_pos(3) = tube_pos(3) + delta_tube_param*8
          enddo
          if (idx > word_len) exit
!
!         create the braid
!         determine the number of the braid
          if (ichar(word(idx:idx)) > ichar('H')) then
            braid_num = ichar(word(idx:idx)) - ichar('a') + 1
          else
            braid_num = ichar(word(idx:idx)) - ichar('A') + 1
          endif
!         check if this braid affects this strand
          if ((braid_num == strand_position(idx_strand)) .or. &
              (braid_num == strand_position(idx_strand) - 1)) then
!
!           switch the position of the strands
            if (braid_num == strand_position(idx_strand)) then
              strand_position(idx_strand) = strand_position(idx_strand) + 1
              strand_position(idx_strand+1) = strand_position(idx_strand+1) - 1
!
              rotation_center(1) = tube_pos(1) + distance_tubes/2.
              rotation_center(2) = tube_pos(2)
!
              phi_offset = 0.
              if (ichar(word(idx:idx)) > ichar('H')) then
                rotation_sign = -1
              else
                rotation_sign = 1
              endif
            else
              strand_position(idx_strand) = strand_position(idx_strand) - 1
              strand_position(idx_strand-1) = strand_position(idx_strand-1) + 1
!
              rotation_center(1) = tube_pos(1) - distance_tubes/2.
              rotation_center(2) = tube_pos(2)
!
              phi_offset = pi
              if (ichar(word(idx:idx)) > ichar('H')) then
                rotation_sign = -1
              else
                rotation_sign = 1
              endif
            endif
!
!           do the rotation
            tube_param = 0.
            z_start = tube_pos(3)
            do while (tube_param < (l_sigma))
              phi = ((a(5)*(tube_param-l_sigma/2.)**5 + a(3)*(tube_param-l_sigma/2.)**3 &
                  + a(1)*(tube_param-l_sigma/2.) - pi)/2. + phi_offset) * rotation_sign
              tube_pos(1) = rotation_center(1) + distance_tubes/2.*cos(phi)
              tube_pos(2) = rotation_center(2) + distance_tubes/2.*sin(phi)
              tube_pos(3) = tube_param + z_start
!
              tangent(1) = -distance_tubes/2. * sin(phi) * &
                  (5*a(5)*(tube_param-l_sigma/2.)**4 + &
                  3*a(3)*(tube_param-l_sigma/2.)**2 + a(1)) * rotation_sign
              tangent(2) = distance_tubes/2. * cos(phi) * &
                  (5*a(5)*(tube_param-l_sigma/2.)**4 + &
                  3*a(3)*(tube_param-l_sigma/2.)**2 + a(1)) * rotation_sign
              tangent(3) = 2.
              tangent = tangent / sqrt(tangent(1)**2+tangent(2)**2+tangent(3)**2)
!
!             Find vector which is orthonormal to tangent vector.
              if (abs(tangent(1)) <= 0.5) then
                normal(1) = tangent(1)**2 - 1.0
                normal(2) = tangent(2)*tangent(1)
                normal(3) = tangent(3)*tangent(1)
              elseif (abs(tangent(2)) <= 0.5) then
                normal(1) = tangent(1)*tangent(2)
                normal(2) = tangent(2)**2 - 1.0
                normal(3) = tangent(3)*tangent(2)
              else
                normal(1) = tangent(1)*tangent(3)
                normal(2) = tangent(2)*tangent(3)
                normal(3) = tangent(3)**2 - 1.0
              endif
              normal = normal / sqrt(normal(1)**2+normal(2)**2+normal(3)**2)
!
              circle_radius = 0.
!             loop which changes the circle's radius
              do
                if (circle_radius > width_tube/2.) exit
!
                circle_param = 0.
!               loop which goes around the circle
                do
                  if (circle_param > 2.*pi) exit
                  circle_pos(1) = tube_pos(1) + circle_radius * &
                      ((tangent(1)*tangent(1)*(1-cos(circle_param))+cos(circle_param))*normal(1) + &
                      (tangent(1)*tangent(2)*(1-cos(circle_param))-tangent(3)*sin(circle_param))*normal(2) + &
                      (tangent(1)*tangent(3)*(1-cos(circle_param))+tangent(2)*sin(circle_param))*normal(3))
                  circle_pos(2) = tube_pos(2) + circle_radius * &
                      ((tangent(1)*tangent(2)*(1-cos(circle_param))+tangent(3)*sin(circle_param))*normal(1) + &
                      (tangent(2)*tangent(2)*(1-cos(circle_param))+cos(circle_param))*normal(2) + &
                      (tangent(2)*tangent(3)*(1-cos(circle_param))-tangent(1)*sin(circle_param))*normal(3))
                  circle_pos(3) = tube_pos(3) + circle_radius * &
                      ((tangent(1)*tangent(3)*(1-cos(circle_param))-tangent(2)*sin(circle_param))*normal(1) + &
                      (tangent(2)*tangent(3)*(1-cos(circle_param))+tangent(1)*sin(circle_param))*normal(2) + &
                      (tangent(3)*tangent(3)*(1-cos(circle_param))+cos(circle_param))*normal(3))
!
!                 Find the corresponding mesh point to this position.
                  l = nint((circle_pos(1) - x(1))/dx) + 1
                  m = nint((circle_pos(2) - y(1))/dy) + 1
                  n = nint((circle_pos(3) - z(1))/dz) + 1
!
                  if ((l > mx .or. m > my .or. n > mz .or. l < 1 .or. m < 1 .or. n < 1) .eqv. .false.) then
!                   Write the magnetic field B.
!                   Note that B is written in the f-array where A is stored.
!                   This is corrected further in the code.
                    if (prof == 'gaussian') then
                      new_bb = tangent*ampl*(exp(-(2*circle_radius/width_tube)**2)-exp(-1.)) / (1-exp(-1.))
                    else if (prof == 'constant') then
                      new_bb = tangent*ampl
                    endif
!                   Avoid issues in spots with high curvature.
                    if ((f(l,m,n,iax)**2 + f(l,m,n,iay)**2 + f(l,m,n,iaz)**2) < &
                        (new_bb(1)**2 + new_bb(2)**2 + new_bb(3)**2)) &
                        f(l,m,n,iax:iaz) = new_bb
                  endif
                  circle_param = circle_param + delta_circle_param
                enddo
                circle_radius = circle_radius + delta_circle_radius
              enddo
              tube_param = tube_param + delta_tube_param
            enddo
          else
          endif
          idx = idx + 1
        enddo
      enddo
!
!     Transform the magnetic field into a vector potential
!
!     Compute curl(B) = J for the Poisson solver
      do m=m1,m2
        do n=n1,n2
          call curl(f,iaa,jj(:,m-nghost,n-nghost,:))
        enddo
      enddo
      tmpJ = -jj
!
!     Use the Poisson solver to solve \nabla^2 A = -J for A
      do j=1,3
        call inverse_laplacian(f,tmpJ(:,:,:,j))
      enddo
!
!     Overwrite the f-array with the correct vector potential A
      do j=1,3
        ju=iaa-1+j
        f(l1:l2,m1:m2,n1:n2,ju) = tmpJ(:,:,:,j)
      enddo
!
!     Add a background field to the braid
      do l=1,mx
        do m=1,my
          f(l,m,:,iax) = f(l,m,:,iax) - y(m)*B_bkg/2.
          f(l,m,:,iay) = f(l,m,:,iay) + x(l)*B_bkg/2.
        enddo
      enddo
!
!     communicate the core boundaries for taking the curl
      call MPI_BARRIER(MPI_comm_world, ierr)
      call boundconds_x(f)
      call initiate_isendrcv_bdry(f)
      call finalize_isendrcv_bdry(f)
      call boundconds_y(f)
      call boundconds_z(f)
      call MPI_BARRIER(MPI_comm_world, ierr)
!
!     convert the magnetic vector potential into the magnetic field
      if (trace_field == 'bb' .and. ipz == 0) then
        do m=m1,m2
          do n=n1,n2
            call curl(f,iaa,vv(:,m-nghost,n-nghost,:))
          enddo
        enddo
      endif
    endif
!
!   In case the blob configuration is wished create it.
!
    if (n_blobs > 0) then
      do j=1,n_blobs
!       create the field in the blobs
        do l=1,mx
          do m=1,my
            do n=1,mz
              f(l,m,n,iaz) = f(l,m,n,iaz) + ampl * blob_sgn(j) * &
                  exp(-((x(l)-xc(j))**2 + (y(m)-yc(j))**2)/blob_scale**2 - ((z(n)-zc(j))/l_blob(j))**2)
            enddo
          enddo
        enddo
      enddo      
!
!     Add a background field to the braid
      do l=1,mx
        do m=1,my
          f(l,m,:,iax) = f(l,m,:,iax) - y(m)*B_bkg/2.
          f(l,m,:,iay) = f(l,m,:,iay) + x(l)*B_bkg/2.
        enddo
      enddo
!
!     convert the magnetic vector potential into the magnetic field
      if (trace_field == 'bb' .and. ipz == 0) then
        do m=m1,m2
          do n=n1,n2
            call curl(f,iaa,vv(:,m-nghost,n-nghost,:))
          enddo
        enddo
      endif
    endif
!
!   Trace the specified field lines
!
    if (trace_field == 'bb' .and. ipz == 0) then
!
!     allocate the memory for the tracers
      allocate(tracers(nx*ny*trace_sub**2,7))
!     create the initial seeds at z(1+nghost)-ipz*nz*dz+dz
      do j=1,nx*trace_sub
        do k=1,ny*trace_sub
          tracers(j+(k-1)*(nx*trace_sub),1) = x(1+nghost) + (dx/trace_sub)*(j-1)
          tracers(j+(k-1)*(nx*trace_sub),2) = y(1+nghost) + (dy/trace_sub)*(k-1)
          tracers(j+(k-1)*(nx*trace_sub),3) = tracers(j+(k-1)*(nx*trace_sub),1)
          tracers(j+(k-1)*(nx*trace_sub),4) = tracers(j+(k-1)*(nx*trace_sub),2)
          tracers(j+(k-1)*(nx*trace_sub),5) = z(n1)-ipz*nz*dz+dz
          tracers(j+(k-1)*(nx*trace_sub),6) = 0.
          tracers(j+(k-1)*(nx*trace_sub),7) = 0.
        enddo
      enddo
!
!     find the tracers
!
      call trace_streamlines(f,tracers,nx*ny*trace_sub**2,vv)
!     write into output file
      write(str_tmp, "(I10.1,A)") iproc, '/tracers.dat'
      write(filename, *) 'data/proc', adjustl(trim(str_tmp))
      open(unit = 1, file = adjustl(trim(filename)), form = "unformatted")
      write(1) 0., tracers(:,:)
      close(1)
!
!     find the fixed points
!
      if (int_q == 'curlyA') then
        call get_fixed_points(f,tracers,vv)
        write(str_tmp, "(I10.1,A)") iproc, '/fixed_points.dat'
        write(filename, *) 'data/proc', adjustl(trim(str_tmp))
        open(unit = 1, file = adjustl(trim(filename)), form = "unformatted")
        write(1) 0.
        write(1) float(fidx)
        do l=1,fidx
          write(1) fixed_points(l,:)
        enddo
        close(1)
      endif
!
!     free allocated memory
      deallocate(tracers)
      deallocate(vv)
    endif
!
  endsubroutine initial_condition_aa
!***********************************************************************
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
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from noinitial_condition.f90 for any    **
!**  InitialCondition routines not implemented in this file        **
!**                                                                **
  include '../initial_condition_dummies.inc'
!********************************************************************
endmodule InitialCondition
