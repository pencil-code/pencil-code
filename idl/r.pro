; $Id: r.pro,v 1.50 2003-06-19 22:36:16 mee Exp $

;;;;;;;;;;;;;;;
;;;  r.pro  ;;;
;;;;;;;;;;;;;;;

;;; Read the data produced on one processor
;;; You should have run `start.pro' once before.
;;; $Id: r.pro,v 1.50 2003-06-19 22:36:16 mee Exp $

function param2
; Dummy to keep IDL from complaining. The real param() routine will be
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
for i=1,totalvars do begin
  res     = res + ',' + varcontent[i].idlvar
  content = content + ', ' + varcontent[i].variable
  ; Initialise variable
  if (varcontent[i].variable eq 'UNKNOWN') then $
           message, 'Unknown variable at position ' + str(i)  $
                                    + ' needs declaring in varcontent.pro', /INFO   
  if (execute(varcontent[i].idlvar+'='+varcontent[i].idlinit,0) ne 1) then $
           message, 'Error initialising ' + varcontent[i].variable $
                                    +' - '+ varcontent[i].idlvar, /INFO
;If it's a vector quantity skip the required number of elements
  i=i+varcontent[i].skip
end

dummy=0.

content = strmid(content,2)
print,'File contains: '+content

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
;  Summarise data
;
;
@varcontent_stats

close,1
;
END

; End of file
