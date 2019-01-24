
    function yy_coors_transform, sth, cth, sph, cph
;
;  Rotate by Pi about z axis, then by Pi/2 about x axis.
;  No distinction between Yin and Yang as transformation matrix is self-inverse.
;
      xprime = -cph*sth
      yprime = -cth
      zprime = -sph*sth

      sprime = sqrt(xprime^2 + yprime^2)

      yz=fltarr(2)*sth
      yz(0) = atan(sprime,zprime)
      yz(1) = atan(yprime,xprime)
      if yz(1) lt 0. then yz(1) += 2.*!pi

      return, yz

    end
;*******************************************************************************
    pro yin2yang_coors_tri,i1,i2,j1,j2,y,z,yz
 
      yz=fltarr(2,(abs(i2-i1)+1)*(abs(j2-j1)+2)/2)*y[0]

      istep = sign(i2-i1) & jstep = sign(j2-j1)

      ind=0
      for i=i1,i2,istep do begin
        for j=j1,j2,jstep do begin
          if ( jstep*(j-j2)+istep*(i-i1) le 0 ) then begin
            yz(*,ind)=yy_coors_transform(sin(y(i)),cos(y(i)),sin(z(j)),cos(z(j)))
            ind++
          endif
        endfor
      endfor

    end
;*******************************************************************************
  pro yin2yang_coors, y, z, yz

      yz=fltarr(2,n_elements(y)*n_elements(z))*y[0]    ; 'y[0] to enforce correct numerical precision

      ind=0L
      for i=0,n_elements(y)-1 do begin

        sth=sin(y(i)) & cth=cos(y(i))

        for j=0,n_elements(z)-1 do begin

          yz(*,ind)=yy_coors_transform(sth,cth,sin(z(j)),cos(z(j)))
          ind++
        
        endfor
      endfor

    end
