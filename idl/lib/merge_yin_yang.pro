  pro merge_yin_yang, y, z, dy, dz, yz, inds
;
;  Merges Yin and Yang grids: Yin unchanged while Yang is clipped and transformed in Yin basis.
;  Output is merged coordinate array yz[2,*] and index vector inds into 1D coordinate arrays of Yang
;  for selection of unclipped points.
;
    yin2yang_coors, y, z, yzyang
    inds = where(yzyang(0,*) lt min(y)-dy or yzyang(0,*) gt max(y)+dy or $
                 yzyang(1,*) lt min(z)-dz or yzyang(1,*) gt max(z)+dz)
    yzyang=yzyang(*,inds)

    ny=n_elements(y) & nz=n_elements(z)
    yz=fltarr(2,ny*nz)

    ind=0
    for i=0,ny-1 do begin
      yz(0,ind:ind+nz-1) = y(i)
      yz(1,ind:ind+nz-1) = z
      ind+=nz
    endfor
    
    yz=[[yz],[yzyang]]   

  end
