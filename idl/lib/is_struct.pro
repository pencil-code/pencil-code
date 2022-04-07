function is_struct, var

; TRUE, if var is structure, FALSE otherwise,
; if var undefined, FALSE is realised as -2

  if is_defined(var) then begin

    strsiz = size(var)

    return, strsiz(n_elements(strsiz)-2) eq 8

  endif else $

    return, -2

end

