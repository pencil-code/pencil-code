; $Id: pc_nan_inform.pro,v 1.1 2004-06-11 13:30:16 ajohan Exp $
pro pc_nan_inform, object=object

  vars = tag_names(object)

  for kk=0,n_elements(vars)-1 do begin

    string = 'var = object.'+vars[kk]
    res=execute(string)

    if (mean(finite(var)) ne 1.) then begin
      print, 'NaN or inf found in ', vars[kk]
    endif

  endfor

end
