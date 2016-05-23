  pro yin2yang_coors, y, z, yz

      yz=fltarr(2,n_elements(y)*n_elements(z))

      ind=0
      for i=0,n_elements(y)-1 do begin

        sth=sin(y(i)) & cth=cos(y(i))

        for j=0,n_elements(z)-1 do begin
;
;  Rotate by Pi about z axis, then by Pi/2 about x axis.
;  No distinction between Yin and Yang as transformation matrix is self-inverse.
;
          xprime = -cos(z(j))*sth
          yprime = -cth
          zprime = -sin(z(j))*sth

          sprime = sqrt(xprime^2 + yprime^2)

          yz(0,ind) = atan(sprime,zprime)
          yz(1,ind) = atan(yprime,xprime)
          if yz(1,ind) lt 0. then yz(1,ind) += 2.*!pi
          ind +=1
        
        endfor
      endfor

    end

