    function transform_thph_yy_other, yz, arr
;
;  Transforms theta and phi components of a vector arr defined with the Yang(Yin) grid basis
;  to the Yin(Yang) grid basis using arbitrary theta and phi coordinates yz of the Yang/(Yin) grid.
;     
      result=arr

      for i=0,(size(yz))[2]-1 do begin

          sinp=sin(yz[1,i])
          sisisq=sqrt(1.-(sin(yz[0,i])*sinp)^2)
          if sisisq  eq 0. then begin             ; i.e. at pole of other grid -> theta and phi components indefined
            a=0. & b=0.
          endif else begin
            sinth1=1./sisisq
            a=cos(yz[1,i])*sinth1 & b=sinp*cos(yz[0,i])*sinth1
          endelse

          result[*,i,1] = b*arr[*,i,1] + a*arr[*,i,2]
          result[*,i,2] =-a*arr[*,i,1] + b*arr[*,i,2]

      endfor

      return, result

    end

