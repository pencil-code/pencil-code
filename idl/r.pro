; $Id: r.pro,v 1.53 2003-06-24 11:25:04 dobler Exp $

;;;;;;;;;;;;;;;
;;;  r.pro  ;;;
;;;;;;;;;;;;;;;

;;; Read the data produced on one processor
;;; You should have run `start.pro' once before.
;;; $Id: r.pro,v 1.53 2003-06-24 11:25:04 dobler Exp $

;
;  read data
;
if ((n_elements(started) le 0) or (n_elements(read_all) gt 0)) then begin
  message, "You need to run start.pro first: use `.rnew start'"
endif
undefine, read_all
;
default, datadir, 'data'
default, varfile, 'var.dat'
;

;
;  Read startup parameters
;
@readstartpars

;
;  Read data
;
@varcontent

; Prepare for read
res=''
content=''
for iv=1,totalvars do begin
  res     = res + ',' + varcontent[iv].idlvar
  content = content + ', ' + varcontent[iv].variable
  ; Initialise variable
  if (varcontent[iv].variable eq 'UNKNOWN') then $
           message, 'Unknown variable at position ' + str(iv)  $
                    + ' needs declaring in varcontent.pro', /INFO   
  cmd = varcontent[iv].idlvar + '='+varcontent[iv].idlinit
  if (execute(cmd) ne 1) then $
      message, 'Error initialising ' + varcontent[iv].variable $
                                     +' - '+ varcontent[iv].idlvar, /INFO
  ; For vector quantities skip the required number of elements
  iv=iv+varcontent[iv].skip
end

dummy=0.

content = strmid(content,2)
if (quiet le 2) then print,'File contains: '+content

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
