    function transform_thph_yy, y, z, arr
;
;  Transforms theta and phi components of a vector arr defined with the Yang grid basis
;  to the Yin grid basis using theta and phi coordinates of the Yang grid.
;
      result=arr

      for i=0,n_elements(z)-1 do begin
        for j=0,n_elements(y)-1 do begin

          sinth1=1./sqrt(cos(y(j))^2+sin(y(j))^2*cos(z(i))^2)
          a=-cos(z(i))*sinth1 & b=sin(z(i))*cos(y(j))*sinth1

          result[*,j,i,1] = b*arr[*,j,i,1] - a*arr[*,j,i,2]
          result[*,j,i,2] = a*arr[*,j,i,1] + b*arr[*,j,i,2]

        endfor
      endfor

      return, result

    end
