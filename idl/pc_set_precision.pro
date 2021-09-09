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
            if arg_present(new) then new=precision
            return
          endif else begin
;
;  Working directory was changed: invalidate precision and dim object.
;
            undefine, precision
            undefine, dim
          endelse
        endif else begin
;
;  Very first call of pc_set_precision during a session.
;
          undefine, dim
        endelse

        setenv, 'IDL_LAST_WORKDIR='+wdir

	; *** WORK HERE: ***
	; [PABourdin, 07-Oct-2015]
	; Remove the following block after all scripts have been shifted to new config finder:

	if (not is_defined(new)) then begin
		pc_read_dim, obj=dim, datadir=datadir, /quiet
		new = strupcase (dim.precision)
	end

	if (not is_defined(new)) then message, "ERROR: precision is a mandatory parameter"

	precision = strupcase (strmid (strtrim (new, 2), 0, 1))
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
        if arg_present(new) then new=precision
end
