;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   input_table.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;; Author:  wd (dobler@uni-sw.gwdg.de)
;;; Date:    9-Feb-2001
;;; Version: 0.1
;;;
;;; Usage:
;;;   IDL> data = input_table(FILE)
;;;   IDL> data = input_table(FILE, OPTIONS)
;;;
;;; Options:
;;;   COMMENT_CHAR, DOUBLE, VERBOSE
;;;
;;; Description:
;;;   Reads tabulated ASCII-formatted data in a rectangular layout
;;;   from FILE.
;;;   Similar to read_table, but reformulated as a function.
;;;
;;; - COMMENT_CHAR=cchar specifies the character (or string) initiating
;;;   comments (`#' by default). For efficiency reasons, this
;;;   string must not be preceeded by white space in the data file,
;;;   i.e. the line must begin directly with it.
;;; - /DOUBLE tells INPUT_TABLE to read double precision values; the
;;;   default are FLOATs
;;; - /VERBOSE tells INPUT_TABLE to waffle a lot.
;;;
;;;   The routine relies on a quadrangular layout of the data,
;;;   i.e. there must be the same number of columns in all lines
;;;   (except the comment lines).
;;;

function input_table, filename, $
                COMMENT_CHAR=cchar, DOUBLE=double, VERBOSE=verb
  ;on_error,2
  if (n_params() ne 1) then begin
    message, 'input_table: wrong number of arguments'
  endif
  if (n_elements(cchar)  eq 0) then cchar  = '#'
  if (n_elements(double) eq 0) then double = 0
  if (n_elements(verb)   eq 0) then verb   = 0

  plural = ['s', '']            ; English plural suffix

  ;;; Obtain number of columns and an optimistic guess to number of
  ;;; lines. This information is needed to pre-allocate the array --
  ;;; extending arrays is inefficient.
  ;;
  ;; Call wc to get the number of fields in a line.
  ;; Assumes there are < 100 comment lines..
  spawn, 'head -100 ' + filename + '| grep -v "^' + cchar + '"' $
      + ' | head -1 | wc -w', ans
  N_cols = long(ans(0))
  ;; Determine the number of lines (either with wc or with egrep)
  ;spawn, 'wc -l ' + filename, ans
  ;spawn, 'grep -v "^' + cchar + '" ' + filename + '| wc -l ', ans
  spawn, 'wc -l ' + filename, ans
  ;spawn,'egrep -c ^ ' + filename, ans
  N_lines = long(ans(0))
  
  if (double) then begin
    row  = dblarr(N_cols)
    data = dblarr(N_cols, N_lines)
  endif else begin
    row  = fltarr(N_cols)
    data = fltarr(N_cols, N_lines)
  endelse
  
  if (verb) then begin
    print, FORMAT= $
        '("Found ", I0, " column", A, " and ",I0 , " line", A, " in file ", A)', $
        N_cols, plural(N_cols eq 1), $
        N_lines, plural(N_lines eq 1), $
        filename
  
    ;; Check whether there are any data at all
    if ((N_cols lt 1) or (N_lines lt 1)) then begin
      message, 'No data! Stopping program execution.'
    endif
  endif
  
  ;; Open the file
  if (verb) then print, format='(A,$)', 'Opening file for input..  '
  openr, in_file, filename, /GET_LUN, ERROR=err
  if (err ne 0) then begin
    free_lun, in_file
    message, !ERR_STRING
  endif
  
  ;; Read the data entries
  first = 1
  line = ''
  iline = 0L                    ; line number
  idat = 0L                     ; datum number (iline + # comment lines)
  clen = strlen(cchar)
  success = 0
  on_ioerror, read_err            ; Catch read errors
  while (not eof(in_file)) do begin
    readf, in_file, line
    is_comm = (strmid(line,0,clen) eq cchar)
    if (not is_comm) then begin ; If this is not a comment line
      reads, line, row
      data[*,idat] = row
      idat = idat + 1
    endif
    iline = iline + 1
  endwhile
  success = 1                     ; If we get here, everything went well
  
  read_err: if (success ne 1) then begin
    print
    print, $
        'input_table: Problems reading data from file ', $
        filename, ', line ', strtrim(iline,2)
  endif
  
  ;; Clean up
  data = data[*,0:idat-1]
  
  close, in_file
  free_lun, in_file
  if (verb) then begin
    print, '..Done reading.'
    print, 'Factually read:'
    help, data
  endif
  
  return, data

end

; End of file 'input_table.pro'
