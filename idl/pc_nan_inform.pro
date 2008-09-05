; $Id$
pro pc_nan_inform, object=object

  vars = tag_names(object)

  for kk=0,n_elements(vars)-1 do begin

    string = 'var = object.'+vars[kk]
    res=execute(string)

    if (min(finite(var)) ne 1) then begin
      print, 'NaN or inf found in ', vars[kk]
    endif

  endfor

end
