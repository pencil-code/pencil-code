   function trans_coords_yy, th, ph, ith1, ith2, iph1, iph2
;
; Performs coordinate transformation between Yin and Yang bases for
; strip defined by theta and phi values th(ith1:ith2) and ph(iph1:iph2), respectively.
; Transformed coordinates are returned as 3D array thphprime with the first dimension
; containing the pairs (theta',phi').
;
      thphprime=fltarr(2,ith2-ith1+1,iph2-iph1+1)

      for i=ith1,ith2 do begin

        sth=sin(th(i)) & cth=cos(th(i))

        for j=iph1,iph2 do begin
;
;  Rotate by Pi about z axis, then by Pi/2 about x axis.
;  No distinction between Yin and Yang as transformation matrix is self-inverse.
;
          xprime = -cos(ph(j))*sth
          yprime = -cth
          zprime = -sin(ph(j))*sth

          sprime = sqrt(xprime^2 + yprime^2)

          itp = i-ith1 & jtp = j-iph1

          thphprime(0,itp,jtp) = atan(sprime,zprime)    ; theta'
          thphprime(1,itp,jtp) = atan(yprime,xprime)    ; phi'
          if thphprime(1,itp,jtp) lt 0. then thphprime(1,itp,jtp) += 2.*!pi

        endfor
      endfor

      return, thphprime

    end

