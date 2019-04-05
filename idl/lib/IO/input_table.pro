;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   input_table.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;; Author:  wd (dobler@uni-sw.gwdg.de)
;;; $Date: 2007-06-12 13:19:08 $
;;; $Revision: 1.12 $
;;;
;;; 21/08/2003 - ajwm (A.J.Mee@ncl.ac.uk) 
;;;   Added STOP_AT and resume with FILEPOSITION behaviour to handle
;;;   'RELOAD' cases where print.in changes in the Pencil-Code
;;; 03/10/2013 - MR:
;;;   Added detection of complex quantities in output, new output parameter
;;;   for their positions in line
;;; 10/01/2017 - MR:
;;;   Added keyword parameter SEPMINUS.
;;;   Allowed STOP_AT to be an array of regular expressions.
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
;;; Keywords:

;;;   COMMENT_CHAR -- specify the character (or string) initiating
;;;                   comments (`#' by default). For efficiency
;;;                   reasons, this string must not be preceeded by
;;;                   white space in the data file, i.e. the line must
;;;                   begin directly with it.
;;;   DOUBLE       -- if set, read double precision values; the
;;;                   default are FLOATs
;;;   NOFILE_OK    -- return one NaN if file does not exist
;;;   VERBOSE      -- waffle a lot
;;;   STOP_AT      -- regular expression to test against each line of
;;;                   input. When a line matches, INPUT_TABLE stops
;;;                   reading data, returning the current file
;;;                   position in FILEPOSITION and the content of the
;;;                   matching line in STOP_AT
;;;   FILEPOSITION -- start reading data at this position. Upon exit
;;;                   contains the position where reading was
;;;                   terminated due to STOP_AT, or -1 if the end of
;;;                   the file was reached.
;;;  INDS_COMPL    -- returns list of positions of complex quantities in line 
;;;                   (detected by opening bracket); = [-1] if none are present
;;;  SEPMINUS      -- minus sign of a numeral allowed to separate columns (e.g. .123E+00-.456E+00 would be read as two quantities)
;;;
;;; Note:
;;;   INPUT_TABLE relies on a quadrangular layout of the data, i.e.
;;;   there must be the same number of columns in all lines (except
;;;   the comment lines).
;;;

function input_table, filename, $
                      STOP_AT=stop_at, FILEPOSITION=fileposition, $
                      COMMENT_CHAR=cchar, DOUBLE=double, $
                      NOFILE_OK=nofile_ok, VERBOSE=verb, $
                      HELP=help, INDS_COMPL=inds_compl, SEPMINUS=sepminus
  ;on_error,2
  default, sepminus, 1

  if (keyword_set(help)) then extract_help, 'input_table'

  if (n_params() ne 1) then begin
    message, 'input_table: wrong number of arguments'
  endif
  if (n_elements(cchar)  eq 0) then cchar  = '#'
  if (n_elements(cchar2) eq 0) then cchar2 = '%'
  if (n_elements(double) eq 0) then double = 0
  if (n_elements(verb)   eq 0) then verb   = 0

  if (keyword_set(nofile_ok)) then begin
    if (not file_test(filename)) then begin
      message, /INFO, 'No such file: ' + filename
      if (double) then return, !values.d_nan else return, !values.f_nan
    endif
  endif

  plural = ['s', '']            ; English plural suffix
  N_cols = 0                    ; Initial value

  ;; Determine the number of lines (either with wc or with egrep)
  ;; Favour the new IDL intrinsic if available
  if ((!version.release ge 5.6) and (running_gdl() eq 0)) then begin
    N_lines = file_lines(filename) 
  endif else begin
    spawn, 'wc -l ' + filename, ans, /SH
    ans_lines = n_elements(ans)
    if (ans_lines gt 1) then begin
      message, /INFO, "`wc -l' returned more than 1 line"
      message, /INFO, "Is your shell clean?"
      message, /INFO, $
          "`bash -c /bin/true' and `csh -c /bin/true' should produce no output"
      print, 'Trying to proceed anyway..'
    endif
    N_lines = long(ans(ans_lines-1))
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

  success = 0
  on_ioerror, read_err            ; Catch read errors

  if (keyword_set(fileposition)) then point_lun, in_file, fileposition
  found_stop=0

  while ((not eof(in_file)) and (idat lt N_lines) $
                            and (not found_stop)) do begin
    readf, in_file, line

    if (keyword_set(fileposition)) then $
      if (fileposition ne -1) then begin
        iline = round(fileposition/strlen(line))
        fileposition=-1 
      endif

    if keyword_set(stop_at) then begin

      stopit=0
      for i=0,n_elements(stop_at)-1 do $
        if strlen(stop_at(i)) gt 0 then stopit += stregex(line,stop_at(i),/BOOLEAN)

      if (stopit gt 0) then begin

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

    is_comm = (strmid(line,0,clen) eq cchar or strmid(line,0,clen) eq cchar2)

    if (not is_comm) then begin

      if keyword_set(sepminus) then begin
        inds_fields=strsplit(line,'[0-9]-',length=len,/REGEX,/PRESERVE_NULL)
        for i=n_elements(inds_fields)-1,1,-1 do $
          line = strmid(line,0,inds_fields(i)-1)+' '+strmid(line,inds_fields(i)-1)
      endif

      ;; If this is first data line, determine number of columns
      if (N_cols eq 0) then begin

        ;; split the line, accept whitespace and (|,|) as separators (for complex output)
        inds_fields = strsplit(line,'[(,)]| ',/REGEX,count=N_cols)      
 
        ;; obtain number of complex quantities
        inds_bracks = strsplit(line,'[(] *',/REGEX,count=N_compl)
        N_compl=N_compl-1
        ;; reduce by one as strsplit returns at least one substring
        ;; yet untreated case: already first quantity is complex
        inds_compl = [-1]
        if N_compl gt 0 then begin
          ;; determine positions of complex quantities
          for i=1,N_compl do $
            inds_compl = [inds_compl,where( inds_fields eq inds_bracks(i) )]
          inds_compl = inds_compl[1:*]
        endif 

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
    

      ;; Read one line and store in data[:,:]
      if (not found_stop) then begin
        is_empty = (strlen(line) eq 0)
        if (not (is_empty)) then begin ; a data line
          if N_compl gt 0 then $
            line = strjoin(strsplit(line,' |[(,)]',/REGEX,/EXTRACT),' ')
          if (not sepminus) then begin

;  If the minus sign is not supposed to be a separator, the situation "+ or - in between two decimal digits"
;  is most likely due to Fortran output of numbers in exponential format with 3-digit exponents *without* eEdD.
;  So it is inserted here. To avoid zeros, the sub must be called with /double.

            ipos=0
            while ipos ge 0 do begin
              ipos=stregex(line,'[0-9][+-][0-9]')
              if ipos ge 0 then $
                line=strmid(line,0,ipos)+'E'+strmid(line,ipos+1)
            endwhile
          endif 
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
