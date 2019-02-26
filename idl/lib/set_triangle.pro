    pro set_triangle,i1,i2,j1,j2,values,f,iyy,yz

      isvec = (size(values))[0] eq 3
      if isvec then values=transform_thph_yy_other(yz, values)

      istep = sign(i2-i1) & jstep = sign(j2-j1)
      ind=0

      for i=i1,i2,istep do $
        for j=j1,j2,jstep do $
          if jstep*(j-j2)+istep*(i-i1) le 0 then begin
            if isvec then $
              f[*,i,j,*,iyy]=values[*,ind,*] $
            else $
              f[*,i,j,iyy]=values[*,ind]
            ind++
          endif
    end
