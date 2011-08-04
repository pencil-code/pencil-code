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
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include '../initial_condition.h'
!
! ampl = amplitude of the magnetic field
! width_tube = width of the flux tube
! l_sigma = length of the twist region
! word = sequence of the braid group
! prof = the amplitude profile across the tube

  real :: ampl = 1.0, width_tube = 0.3, l_sigma = 0.3, steepnes = 1.0
  character (len=labellen) :: prof='gaussian'
  character (len=labellen) :: word = "AA"
!
  namelist /initial_condition_pars/ &
      ampl,width_tube,l_sigma,steepnes,word,prof
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
!  NB: Currently this works only for one CPU. For multi CPU usage
!  initialize on one CPU and do the remeshing.
!  NB2: Works best with the same resolution (box_size_xyz/mesh_points_xyz) in all directions.
!  This has to do with the curl operator I guess.
!
!  Created 2011-07-01 by Simon Candelaresi (Iomsn)
!
        use Mpicomm, only: stop_it
        use Poisson
        use Sub
        
        real, dimension (mx,my,mz,mfarray) :: f
        
        real :: distance_tubes ! distance between the flux tubes
        real :: l_straight ! length of the regions with a straight magnetic field
        real :: a(5) ! the coefficients for the rotation polynomial
        real :: z0 ! auxiliary variable for the creation of the magnetic field
        
        real :: phi ! = phi, the rotation angle for the braid
        real :: phi_offset ! auxiliary variable for the rotation
        integer :: rotation_sign ! sign of the rotation in the braid
        real :: rotation_center(2) ! the center of the rotation for the braid
        
        integer :: word_len, idx, idx_strand
        real :: x_size, y_size, z_size ! the dimensions of the box
        character (len=len_trim(word)) :: wordn ! the trimmed word
        
        ! keeps track of the current position of the strand
        integer, dimension(9) :: strand_position
        ! position of the other strand in the braid
        integer :: other_position
        integer :: braid_num ! auxiliary variable for the braid position
        
        ! step length for some curve parameters
        real :: delta_tube_param, delta_circle_param, delta_circle_radius
        real :: tube_param, circle_param, circle_radius
        ! vectors for the creation of the flux tubes
        real, dimension(3) :: tube_pos, circle_pos, tangent, normal
        
        ! The next 3 variables are used for the uncurling.
        integer :: l, j, ju, n_strands, ascii_code
        real, dimension (nx,ny,nz,3) :: jj, tmpJ  ! This is phi for poisson.f90


        !  check the word
        
        wordn = word
        word_len = len(wordn)
        if (verify(wordn, 'ABCDEFGHabcdefgh') .ge. word_len) then
            write(*,*) "error: invalid word. Note that the highest braid is 'h'"
        endif
            
        !  determine the number of strands
        n_strands = 0
        do idx = 1, word_len
            ascii_code = ichar(wordn(idx:idx))
            if ((ascii_code .ge. 65) .and. (ascii_code .le. 72)) then
                if (n_strands .le. ascii_code - 64) then
                    n_strands = ascii_code - 64
                endif
            elseif ((ascii_code .ge. 97) .and. (ascii_code .le. 104)) then
                if (n_strands .le. ascii_code - 96) then
                    n_strands = ascii_code - 96
                endif
            endif
        enddo
        n_strands = n_strands + 1
 
        !  determine the size of the box
        x_size = x(l2) - x(l1)
        y_size = y(m2) - y(m1)
        z_size = z(n2) - z(n1)
      
        !  compute the distance between the braids
        l_straight = (z_size - word_len*l_sigma) / (word_len+1)  
        if (l_straight .le. 0.) then
            write(*,*) "error: distance between braids is negative, check l_sigma"
        endif
        
        distance_tubes = x_size / (n_strands+1)
      

        !  clear the magnetic vector potential in the farray

        f(l1:l2,m1:m2,n1:n2,iax:iaz) = 0.0      

        !  set the coefficients for the rotation polynomial
        a(1) = steepnes
        a(2) = 0
        a(3) = -(a(1)*l_sigma + 5*l_sigma*a(1) - 15*pi) / (3*l_sigma*(l_sigma/2.)**2)
        a(4) = 0      
        a(5) = (l_sigma*a(1) - 3*pi) / (2*(l_sigma/2.)**5)

        !  Calculate the minimum step size of the curve parameters 
        !  to avoid discretization issues, like mesh points without magnetic field
        delta_tube_param = min(x(l1+1)-x(l1), y(m1+1)-y(m1), z(n1+1)-z(n1))
        ! correct for the braid steepnes
        delta_tube_param = delta_tube_param * l_sigma / (steepnes * pi*distance_tubes)
        delta_circle_radius = delta_tube_param
        delta_circle_param = delta_circle_radius/(width_tube/2.)

        !  loop over all strands      
        do idx_strand = 1,n_strands
            !  compute the initial position if this strand
            tube_pos(1) = x(l1) + idx_strand * x_size / (n_strands + 1)
            tube_pos(2) = (y(m2) + y(m1))/2.
            tube_pos(3) = z(n1)
            idx = 1
            
            ! reset the strand_position vector
            strand_position = (/1,2,3,4,5,6,7,8,9/)
            
            ! create a field in the ghost zones for numerical stability
            tangent = (/0,0,1/)
            !  loop which changes the circle's radius
            circle_radius = 0.
            do
                if (circle_radius .gt. width_tube/2.) exit
                !  loop which goes around the circle
                circle_param = 0.
                do
                    if (circle_param .gt. 2.*pi) exit
                    circle_pos(1) = tube_pos(1) + circle_radius*cos(circle_param)
                    circle_pos(2) = tube_pos(2) + circle_radius*sin(circle_param)

                    !  Find the corresponding mesh point to this position.
                    l = nint((circle_pos(1) - x(l1))/x_size * (l2-l1)) + l1
                    m = nint((circle_pos(2) - y(m1))/y_size * (m2-m1)) + m1
                    n = nint((circle_pos(3) - z(n1))/z_size * (n2-n1)) + n1

                    !  Write the magnetic field B.
                    !  Note that B is written in the f-array where A is stored. This is
                    !  corrected further in the code.
                    do j = 1,n1
                        if (prof == 'gaussian') then
                            f(j,m,n,iax:iaz) = tangent*ampl*exp(-(2*circle_radius/width_tube)**2)
                        else if (prof == 'constant') then
                            f(j,m,n,iax:iaz) = tangent*ampl
                        else
                            write(*,*) "error: invalid magnetic field profile"
                        endif
                    enddo
                    circle_param = circle_param + delta_circle_param
                enddo
                circle_radius = circle_radius + delta_circle_radius
            enddo            
            
            !  loop over all braids
            do
                !  create straight lines
                do
                    if (tube_pos(3) .gt. (l_straight*idx + l_sigma*(idx-1) + z(n1))) exit
                    tangent = (/0,0,1/)

                    !  loop which changes the circle's radius
                    circle_radius = 0.
                    do
                        if (circle_radius .gt. width_tube/2.) exit
                        !  loop which goes around the circle
                        circle_param = 0.
                        do
                            if (circle_param .gt. 2.*pi) exit
                            circle_pos(1) = tube_pos(1) + circle_radius*cos(circle_param)
                            circle_pos(2) = tube_pos(2) + circle_radius*sin(circle_param)
                            circle_pos(3) = tube_pos(3)

                            !  Find the corresponding mesh point to this position.
                            l = nint((circle_pos(1) - x(l1))/x_size * (l2-l1)) + l1
                            m = nint((circle_pos(2) - y(m1))/y_size * (m2-m1)) + m1
                            n = nint((circle_pos(3) - z(n1))/z_size * (n2-n1)) + n1

                            !  Write the magnetic field B.
                            !  Note that B is written in the f-array where A is stored. This is
                            !  corrected further in the code.
                            if (prof == 'gaussian') then
                                f(l,m,n,iax:iaz) = tangent*ampl*exp(-(2*circle_radius/width_tube)**2)
                            else if (prof == 'constant') then
                                f(l,m,n,iax:iaz) = tangent*ampl
                            endif
                            circle_param = circle_param + delta_circle_param
                        enddo
                        circle_radius = circle_radius + delta_circle_radius
                    enddo            
                    tube_pos(3) = tube_pos(3) + delta_tube_param
!                     tube_pos(3) = tube_pos(3) + z(l1+1) - z(l1)
                enddo
                if (idx .gt. word_len) exit

                !  create the braid
                ! determine the number of the braid
                if (ichar(word(idx:idx)) .gt. ichar('H')) then
                    braid_num = ichar(word(idx:idx)) - ichar('a') + 1
                else
                    braid_num = ichar(word(idx:idx)) - ichar('A') + 1
                endif                
                ! check if this braid affects this strand
                if ((braid_num == strand_position(idx_strand)) .or. &
                    (braid_num == strand_position(idx_strand) - 1)) then
                    ! create the braid
                    
                    ! switch the position of the strands
                    if (braid_num == strand_position(idx_strand)) then
                        strand_position(idx_strand) = strand_position(idx_strand) + 1
                        strand_position(idx_strand+1) = strand_position(idx_strand+1) - 1
                        
                        rotation_center(1) = tube_pos(1) + distance_tubes/2.
                        rotation_center(2) = tube_pos(2)
                        
                        phi_offset = 0.
                        if (ichar(word(idx:idx)) .gt. ichar('H')) then
                            rotation_sign = -1
                        else
                            rotation_sign = 1
                        endif                    
                    else
                        strand_position(idx_strand) = strand_position(idx_strand) - 1
                        strand_position(idx_strand-1) = strand_position(idx_strand-1) + 1

                        rotation_center(1) = tube_pos(1) - distance_tubes/2.
                        rotation_center(2) = tube_pos(2)
                        
                        phi_offset = pi
                        if (ichar(word(idx:idx)) .gt. ichar('H')) then
                            rotation_sign = -1
                        else
                            rotation_sign = 1
                        endif                    
                    endif
                
                    ! do the rotation
                    tube_param = 0.
                    z0 = tube_pos(3)
                    do while (tube_param .lt. (l_sigma))
                        phi = ((a(5)*(tube_param-l_sigma/2.)**5 + a(3)*(tube_param-l_sigma/2.)**3 &
                            + a(1)*(tube_param-l_sigma/2.) - pi)/2. + phi_offset) * rotation_sign
                        tube_pos(1) = rotation_center(1) + distance_tubes/2.*cos(phi)
                        tube_pos(2) = rotation_center(2) + distance_tubes/2.*sin(phi)
                        tube_pos(3) = tube_param + z0
                        
                        tangent(1) = -distance_tubes/2. * sin(phi) * &
                                    (5*a(5)*(tube_param-l_sigma/2.)**4 + &
                                    3*a(3)*(tube_param-l_sigma/2.)**2 + a(1)) * rotation_sign
                        tangent(2) = distance_tubes/2. * cos(phi) * &
                                    (5*a(5)*(tube_param-l_sigma/2.)**4 + &
                                    3*a(3)*(tube_param-l_sigma/2.)**2 + a(1)) * rotation_sign
                        tangent(3) = 2.
                        tangent = tangent / sqrt(tangent(1)**2+tangent(2)**2+tangent(3)**2)
                        
                        !  Find vector which is orthonormal to tangent vector.
                        if (abs(tangent(1)) .le. 0.5) then
                            normal(1) = tangent(1)**2 - 1.0
                            normal(2) = tangent(2)*tangent(1)
                            normal(3) = tangent(3)*tangent(1)
                        elseif (abs(tangent(2)) .le. 0.5) then
                            normal(1) = tangent(1)*tangent(2)
                            normal(2) = tangent(2)**2 - 1.0
                            normal(3) = tangent(3)*tangent(2)
                        else
                            normal(1) = tangent(1)*tangent(3)
                            normal(2) = tangent(2)*tangent(3)
                            normal(3) = tangent(3)**2 - 1.0
                        endif
                        normal = normal / sqrt(normal(1)**2+normal(2)**2+normal(3)**2)

                        !  loop which changes the circle's radius
                        circle_radius = 0.
                        do
                            if (circle_radius .gt. width_tube/2.) exit

                            !  loop which goes around the circle
                            circle_param = 0.
                            do
                                if (circle_param .gt. 2.*pi) exit
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

                                !  Find the corresponding mesh point to this position.

                                l = nint((circle_pos(1) - x(l1))/x_size * (l2-l1)) + l1
                                m = nint((circle_pos(2) - y(m1))/y_size * (m2-m1)) + m1
                                n = nint((circle_pos(3) - z(n1))/z_size * (n2-n1)) + n1

                                !  Write the magnetic field B.
                                !  Note that B is written in the f-array where A is stored. This is
                                !  corrected further in the code.
                        
                                if (prof == 'gaussian') then
                                    f(l,m,n,iax:iaz) = tangent*ampl*exp(-(2*circle_radius/width_tube)**2)
                                else if (prof == 'constant') then
                                    f(l,m,n,iax:iaz) = tangent
                                endif
                                circle_param = circle_param + delta_circle_param
                            enddo
                            circle_radius = circle_radius + delta_circle_radius
                        enddo            

                        tube_param = tube_param + delta_tube_param                                
                    enddo
                else
                    ! create a straight field
                    tube_pos(3) = tube_pos(3)
                endif
                idx = idx + 1
            enddo
            
            ! create a field in the ghost zones for numerical stability
            tangent = (/0,0,1/)
            !  loop which changes the circle's radius
            circle_radius = 0.
            do
                if (circle_radius .gt. width_tube/2.) exit
                !  loop which goes around the circle
                circle_param = 0.
                do
                    if (circle_param .gt. 2.*pi) exit
                    circle_pos(1) = tube_pos(1) + circle_radius*cos(circle_param)
                    circle_pos(2) = tube_pos(2) + circle_radius*sin(circle_param)

                    !  Find the corresponding mesh point to this position.
                    l = nint((circle_pos(1) - x(l1))/x_size * (l2-l1)) + l1
                    m = nint((circle_pos(2) - y(m1))/y_size * (m2-m1)) + m1
                    n = nint((circle_pos(3) - z(n1))/z_size * (n2-n1)) + n1

                    !  Write the magnetic field B.
                    !  Note that B is written in the f-array where A is stored. This is
                    !  corrected further in the code.
                    do j = 1,n1
                        if (prof == 'gaussian') then
                            f(j+n2,m,n,iax:iaz) = tangent*ampl*exp(-(2*circle_radius/width_tube)**2)
                        else if (prof == 'constant') then
                            f(j+n2,m,n,iax:iaz) = tangent*ampl
                        endif
                    enddo
                    circle_param = circle_param + delta_circle_param
                enddo
                circle_radius = circle_radius + delta_circle_radius
            enddo                        
        enddo

!
!  Transform the magnetic field into a vector potential
!
!  Compute curl(B) = J for the Poisson solver

      do m=m1,m2
         do n=n1,n2
            call curl(f,iaa,jj(:,m-nghost,n-nghost,:))
         enddo
      enddo
      tmpJ = -jj
!  Use the Poisson solver to solve \nabla^2 A = -J for A
      do j=1,3
        call inverse_laplacian(f,tmpJ(:,:,:,j))
      enddo
      
!  Overwrite the f-array with the correct vector potential A
      do j=1,3
          ju=iaa-1+j
          f(l1:l2,m1:m2,n1:n2,ju) = tmpJ(:,:,:,j)
      enddo
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
