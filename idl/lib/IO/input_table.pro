;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   input_table.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;; Author:  wd (dobler@uni-sw.gwdg.de)
;;; $Date: 2004-03-12 12:36:09 $
;;; $Revision: 1.7 $
;;;
;;; 21/08/2003 - ajwm (A.J.Mee@ncl.ac.uk) 
;;;   Added STOP_AT and resume with FILEPOSITION behaviour to handle
;;;   'RELOAD' cases where print.in changes in the Pencil-Code
;;;
;;; Usage:
;;;   IDL> data = input_table(FILE)
;;;   IDL> data = input_table(FILE, COMMENT_CHAR='%', /DOUBLE)
;;;   IDL> data = input_table(FILE, STOP_AT='(\?\?\?|\*\*\*)', FILEPOS=nan_pos)
;;;
;;; Options:
;;;   COMMENT_CHAR, DOUBLE, VERBOSE, STOP_AT, FILEPOSITION
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
;;; - STOP_AT may be set on entry to a regular expression to test against
;;;   each line of input.  When a line matches, INPUT_TABLE stops
;;;   reading data, returning the current file position in
;;;   FILEPOSITION and the content of the matching line in STOP_AT
;;; - FILEPOSITION contains the position in the file from which to start
;;;   processing the data and upon exit contains the file position at
;;;   which processing was terminated using STOP_AT or -1 if the end of
;;;   the file was reached.
;;;
;;;   The routine relies on a quadrangular layout of the data,
;;;   i.e. there must be the same number of columns in all lines
;;;   (except the comment lines).
;;;

function input_table, filename, $
                STOP_AT=stop_AT,FILEPOSITION=fileposition, $
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
  ;spawn, 'head -100 ' + filename + '| grep -v "^' + cchar + '"' $
  ;    + ' | head -1 | wc -w', ans
  ;N_cols = long(ans(0))
  ;
  ; Initialisation now done on first loop
  N_cols = 0

  ;; Determine the number of lines (either with wc or with egrep)
  ;;spawn, 'wc -l ' + filename, ans
  ;;spawn, 'grep -v "^' + cchar + '" ' + filename + '| wc -l ', ans
  ;;spawn,'egrep -c ^ ' + filename, ans

  ;; Favour the new IDL intrinsic if available
  if (!version.release ge 5.6) then begin
    N_lines = file_lines(filename) 
  endif else begin
    spawn, 'wc -l ' + filename, ans
    N_lines = long(ans(0))
  endelse

  ;; Open the file
  if (verb) then print, format='(A,$)', 'Opening file for input..  '
  openr, in_file, filename, /GET_LUN, ERROR=err
  if (err ne 0) then begin
    free_lun, in_file
    message, !err_string
  endif
  
  ;; Read the data entries
  line = ''
  iline = 0L                    ; line number
  idat = 0L                     ; datum number (iline - # comment lines)
  clen = strlen(cchar)

  slen = -1
  
  if (keyword_set(stop_at)) then begin
    slen = strlen(stop_at)
  endif 

  success = 0
  on_ioerror, read_err            ; Catch read errors

  if (keyword_set(fileposition)) then begin
    point_lun, in_file, fileposition
    iline = fileposition
  endif
  fileposition=-1 
  found_stop=0
  while ((not eof(in_file)) and (idat lt N_lines) $
                            and (not found_stop)) do begin
    readf, in_file, line
    is_comm = (strmid(line,0,clen) eq cchar)

    if (not is_comm) then begin
      ;; If this is first data line, determine number of columns
      if (N_cols eq 0) then begin
        N_cols = n_elements(strsplit(line,'[\ ]',/REGEX,/EXTRACT))

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
      endif
    
      if (slen gt 0) then begin
        if (stregex(line,STOP_AT,/BOOLEAN)) then begin
          point_lun,-in_file,fileposition ; Save file position
          stop_at=line                    ; Return the line
          found_stop=-1                   ; Exit the loop
          if (verb) then begin
            print, 'Found stop regexp at'
            print, '  position ', strtrim(fileposition,2)
            print, '  line no. ', strtrim(iline,2), ' (starting from 0)'
            print, '  line = <'+line+'>'
          endif
        endif
      endif

      ;; Read one line and store in data[:,:]
      if (not found_stop) then begin
        is_empty = (strlen(line) eq 0)
        if (not (is_empty)) then begin ; a data line
          reads, line, row
          data[*,idat] = row
          idat = idat + 1
        endif
      endif
    endif

    iline = iline + 1
  endwhile
  success = 1                   ; If we get here, everything went well
  
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
