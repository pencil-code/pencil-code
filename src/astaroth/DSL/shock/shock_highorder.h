
      call smoothing_kernel(smooth_factor,lgaussian_smooth)
      
    	dt_div_pow = dtfactor**(shock_div_pow-1)
    
call div(f,iuu,penc)

	if (lconvergence_only) then
          f(l1:l2,m,n,ishock) = max(0.0,-penc)
        else
          if (lconvergence_bias) then
            f(l1:l2,m,n,ishock) = max(0.0,-penc) + max(0.0,con_bias*penc)
          else
            f(l1:l2,m,n,ishock) = abs(penc)
          endif
        endif
        if (shock_div_pow /= 1.) f(l1:l2,m,n,ishock)=dt_div_pow*f(l1:l2,m,n,ishock)**shock_div_pow

call smooth_shock(f,ishock,tmp)

	f(l1:l2,m1:m2,n1:n2,ishock) = tmp(l1:l2,m1:m2,n1:n2)*dxmin**2


      if (lshear) then
        call boundconds_y(f,ishock,ishock)
        call initiate_isendrcv_bdry(f,ishock,ishock)
        call finalize_isendrcv_bdry(f,ishock,ishock)
      endif
    call copy_addr(ishock_max   ,p_par(1))  ! int
    call copy_addr(div_threshold,p_par(2))
    call copy_addr(shock_linear ,p_par(3))
    call copy_addr(shock_div_pow,p_par(4))
    call copy_addr(dt_div_pow   ,p_par(5))
    call copy_addr(con_bias     ,p_par(6))
