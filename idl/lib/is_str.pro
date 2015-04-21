function is_str, var

; TRUE, if var is string, FALSE otherwise,
; if var undefined, FALSE is realised as -2

  if is_defined(var) then begin

    strsiz = size(var)

    return, strsiz(n_elements(strsiz)-2) eq 7

  endif else $

    return, -2

end

