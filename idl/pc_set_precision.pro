; Set the precision parameters in the pc_precision common block.
;
; 14-Sep-2015/Bourdin.KIS: redesigned completely

pro pc_set_precision, precision=new, dim=dim, datadir=datadir, QUIET=QUIET

COMPILE_OPT IDL2, HIDDEN

	common pc_precision, zero, one, precision, data_type, data_bytes, type_idl

        cd, current=wdir & wdir=strtrim(wdir,2)

        if (is_defined(precision)) then begin
          if (strtrim(getenv('IDL_LAST_WORKDIR'),2) eq wdir) then begin
;
;  Working directory was not changed.
;
            if is_defined(dim) then begin
              if (precision eq strupcase(dim.precision)) then return
            endif
          endif else begin
;
;  Working directory was changed: invalidate dim object.
;
            undefine, dim
          endelse
        endif

        setenv, 'IDL_LAST_WORKDIR='+wdir

	; *** WORK HERE: ***
	; [PABourdin, 07-Oct-2015]
	; Remove the following block after all scripts have been shifted to new config finder:

        if is_defined(new) then begin
          precision=strupcase(strtrim(new,2))
          if precision ne 'D' and precision ne 'S' then begin
            print, 'pc_set_precision: Error - new precision '+strtrim(new,2)+' is neither "S" nor "D"!!!'
            stop
          endif
        endif else begin
	  if not is_defined(dim) then pc_read_dim, obj=dim, datadir=datadir, /quiet
          precision = strupcase(strmid(strtrim(dim.precision, 2), 0, 1))
        endelse

	if (precision eq 'D') then begin
		; double precision
		zero = 0.d0
		one = 1.d0
		data_type = 'double'
		data_bytes = 8
		type_idl = 5
	end else begin
		; single precision
		precision = 'S'
		zero = 0.e0
		one = 1.e0
		data_type = 'real'
		data_bytes = 4
		type_idl = 4
	end
end
