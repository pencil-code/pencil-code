function is_in, feld, val, replace=replace, abbrev=abbrev

; 10/03/27: MR coded
;
; detects whether value val is contained in array feld;
; if so, the return value is the index of the first position in feld, at which
; val occurs, if not -1;
; if feld undefined, return value = -1.
; if feld is of type string, switch /abbrev causes
; abbreviations of val to be detected
; if replace is set, val is everywhere in feld replaced by
; replace ersetzt;
; in case of ambiguities due to abbrev only those hits are replaced, which
; are completely identical with the first hit

  len = n_elements(feld)

  strabbr = keyword_set(abbrev) and is_str(feld)
  if len gt 0 then begin

    if strabbr then begin

      feldt = strtrim(feld,2)
      minlen = strlen(val) < strlen(feldt)

      for i=0,len-1 do $
        if strmid( val, 0, minlen(i) ) eq strmid( feldt(i), 0, minlen(i) ) then $
          inds = add2vec( i, inds )

      if not is_defined( inds ) then inds = -1

    endif else $

      inds = where( val eq feld )

    if (inds(0) ne -1) and is_defined(replace) then $
      if strabbr then $
        feld(inds( where( feld(inds) eq feld(inds(0)) ))) = replace $
      else $
        feld(inds) = replace

    return, inds(0)

  endif else $

    return, -1

end

