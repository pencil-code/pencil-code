; $Id$

;;;;;;;;;;;;;;;
;;;  r.pro  ;;;
;;;;;;;;;;;;;;;

;;; Read the data produced on one processor
;;; You should have run `start.pro' once before.
;;; $Id$

function param2
COMPILE_OPT HIDDEN 
; Dummy to keep IDL from complaining. The real param2() routine will be
; compiled below
end

;
;  read data
;
if ((n_elements(started) le 0) or (n_elements(read_all) gt 0)) then begin
  message, "You need to run start.pro first: use `.rnew start'"
endif
undefine, read_all
;
if (n_elements(datatopdir) eq 0) then datatopdir=pc_get_datadir()
default, varfile, 'var.dat'
;
;  Read startup parameters
;
@readstartpars

;
;  Read data
;
varcontent=pc_varcontent(QUIET=quiet)
totalvars=(size(varcontent))[1]
; Prepare for read
res=''
content=''
for iv=0L,totalvars-1L do begin
  res     = res + ',' + varcontent[iv].idlvar
  content = content + ', ' + varcontent[iv].variable
  ; Initialise variable
  if (varcontent[iv].variable eq 'UNKNOWN') then $
           message, 'Unknown variable at position ' + str(iv)  $
                    + ' needs declaring in pc_varcontent.pro', /INFO   
  cmd = varcontent[iv].idlvar + '='+varcontent[iv].idlinit
  if (execute(cmd) ne 1) then $
      message, 'Error initialising ' + varcontent[iv].variable $
                                     +' - '+ varcontent[iv].idlvar, /INFO
  ; For vector quantities skip the required number of elements
  iv=iv+varcontent[iv].skip
end

content = strmid(content,2)
if (quiet le 2) then print,'File '+varfile+' contains: ', content

close,1
openr,1, datadir+'/'+varfile, /F77
  if (execute('readu,1'+res) ne 1) then $
           message, 'Error reading: ' + 'readu,1'+res
 
;
if (lshear) then begin
  readu,1, t, x, y, z, dx, dy, dz, deltay
end else begin
  readu,1, t, x, y, z, dx, dy, dz
end
;
xx = spread(x, [1,2], [my,mz])
yy = spread(y, [0,2], [mx,mz])
zz = spread(z, [0,1], [mx,my])
rr = sqrt(xx^2+yy^2+zz^2)
;
;  Summarize data
;
;
@varcontent_stats

close,1
;
END

; End of file
