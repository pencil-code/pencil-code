  pro merge_yin_yang, m1, m2, n1, n2, y, z, dy, dz, yz, inds, yghosts=yghosts, zghosts=zghosts
;
;  Merges Yin and Yang grids: Yin unchanged while Yang is clipped and transformed in Yin basis.
;  Output is merged coordinate array yz[2,*] and index vector inds into 1D coordinate arrays of Yang
;  for selection of unclipped points.
;
    yin2yang_coors, y[m1:m2], z[n1:n2], yzyang
    inds = where(yzyang(0,*) lt min(y[m1:m2])-(.0)*dy or yzyang(0,*) gt max(y[m1:m2])+(.0)*dy or $
                 yzyang(1,*) lt min(z[n1:n2])-(.0)*dz or yzyang(1,*) gt max(z[n1:n2])+(.0)*dz)
    yzyang=yzyang(*,inds)

    ny=m2-m1+1 & nz=n2-n1+1            ;long(n_elements(y)) & nz=n_elements(z)
    yz=fltarr(2,ny*nz)*dy              ; *dy to enforce correct numerical precision

    ind=0L
    for i=0,ny-1 do begin
      yz(0,ind:ind+nz-1) = y(i+m1)
      yz(1,ind:ind+nz-1) = z[n1:n2]
      ind+=nz
    endfor
    
    yz=[[yz],[yzyang]]   

    if arg_present(yghosts) then begin
      yin2yang_coors, y[0:m1-1], z, yghosts
      yin2yang_coors, y[m2+1:*], z, yzyang
      yghosts=[[yghosts],[yzyang]]
    endif

    if arg_present(zghosts) then begin
      yin2yang_coors, y, z[0:n1-1], zghosts
      yin2yang_coors, y, z[n2+1:*], yzyang
      zghosts=[[zghosts],[yzyang]]
    endif

  end
