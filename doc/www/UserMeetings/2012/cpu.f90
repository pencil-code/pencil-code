! This is a pseudo code that describes the basic functionality of the
! nm loop for comparison between the GPU solution.
!
! f, df, and p are all "CPU" memory
!
! p is the 1D pencil with length nx
!=====================================================================
subroutine time_step(f,df,p)

  ...
  call pde(f,df,p)
  ...

endsubroutine time_step

!=====================================================================
subroutine pde(f,df,p)

  intent(inout)  :: f
  intent(out)    :: df, p

  mn_loop: do n and m

    p%uu=f(l1:l2,m,n,iux:iuz)
    ...

    p%ugu=...
    ...

    df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)-p%ugu
    ...

  enddo mn_loop

endsubroutine pde
