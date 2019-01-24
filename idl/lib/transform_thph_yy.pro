    function transform_thph_yy, y, z, arr
;
;  Transforms theta and phi components of a vector arr defined with the Yang(Yin) grid basis
;  to the Yin(Yang) grid basis using theta and phi grid coordinates y and z of the Yang(Yin) grid.
;
      result=arr

      for i=0,n_elements(z)-1 do begin

        sinp=sin(z(i))
        for j=0,n_elements(y)-1 do begin

          sisisq=sqrt(1.-(sin(y(j))*sinp)^2)
          if sisisq  eq 0. then begin             ; i.e. at pole of other grid -> theta and phi components indefined
            a=0. & b=0.
          endif else begin
            sinth1=1./sisisq
            a=cos(z(i))*sinth1 & b=sinp*cos(y(j))*sinth1
          endelse

          result[*,j,i,1] = b*arr[*,j,i,1] + a*arr[*,j,i,2]
          result[*,j,i,2] =-a*arr[*,j,i,1] + b*arr[*,j,i,2]

        endfor
      endfor

      return, result

    end
