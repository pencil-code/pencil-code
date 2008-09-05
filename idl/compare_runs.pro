;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   compare_runs.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   10-Oct-2002
;;;  $Id$
;;;
;;;  Description:
;;;   Compare serial and parallel runs with reference data. called by
;;;   run-ser-par.pro

pro compare_runs, var1, var2, lmn12, REFERENCE=ref

  ; extract l1, l2, etc
  l1=lmn12.l1
  l2=lmn12.l2
  m1=lmn12.m1
  m2=lmn12.m2
  n1=lmn12.n1
  n2=lmn12.n2

  comp_vars = ['par']

  if (keyword_set(ref)) then begin
    comp_vars = [comp_vars, 'ser_orig','par_orig']
    restore, 'vars-ser_orig.sav'
    restore, 'vars-par_orig.sav'
  endif

  for i=0,n_elements(comp_vars)-1 do begin
    print, 'Comparing serial data with ', comp_vars[i]
    case (comp_vars[i]) of
      'par':      var = var2
      'par_orig': var = varp
      'ser_orig': var = vars
    endcase
    duu = var1.uu   -var.uu
    dll = var1.lnrho-var.lnrho
    dss = var1.ss   -var.ss

    print, '          ', '  max |df|   ', '  max |df|   ', '  max |df|   '
    print, '          ', '  global     ', '  no corners ', '  inner      '

    mvar = 5
    nghost = 3
    for i=0,mvar-1 do begin
      case i of
        0: begin & dvar=duu[*,*,*,i] & label='ux   ' & end
        1: begin & dvar=duu[*,*,*,i] & label='uy   ' & end
        2: begin & dvar=duu[*,*,*,i] & label='uz   ' & end
        3: begin & dvar=dll[*,*,*  ] & label='lnrho' & end
        4: begin & dvar=dss[*,*,*  ] & label='ss   ' & end
      endcase
      dv_glob   = max(abs(dvar))
      dv_inner  = max(abs(dvar[l1:l2,m1:m2,n1:n2]))
      dvar[0:l1-1,0:m1-1,0:n1-1] = 0.
      dvar[0:l1-1,0:m1-1,n2+1:*] = 0.
      dvar[0:l1-1,m2+1:*,0:n1-1] = 0.
      dvar[0:l1-1,m2+1:*,n2+1:*] = 0.
      dvar[l2+1:*,0:m1-1,0:n1-1] = 0.
      dvar[l2+1:*,0:m1-1,n2+1:*] = 0.
      dvar[l2+1:*,m2+1:*,0:n1-1] = 0.
      dvar[l2+1:*,m2+1:*,n2+1:*] = 0.
      dv_nocorn = max(abs(dvar))

      ;; Now tell us where the differences are located
      loc = [' ', ' ',  ' ',  '  ',  '  ',  '  ',  '   ']
      ;; faces:
      if ((where(dvar[0:l1-1,m1:m2,n1:n2] ne 0 or $
                 dvar[l2+1:*,m1:m2,n1:n2] ne 0))[0] ge 0) then loc[0] = 'x'
      if ((where(dvar[l1:l2,0:m1-1,n1:n2] ne 0 or $
                 dvar[l1:l2,m2+1:*,n1:n2] ne 0))[0] ge 0) then loc[1] = 'y'
      if ((where(dvar[l1:l2,m1:m2,0:n1-1] ne 0 or $
                 dvar[l1:l2,m1:m2,n2+1:*] ne 0))[0] ge 0) then loc[2] = 'z'
      ;; edges:
      if ((where(dvar[0:l1-1,0:m1-1,n1:n2] ne 0 or $
                 dvar[0:l1-1,m2+1:*,n1:n2] ne 0 or $
                 dvar[l2+1:*,0:m1-1,n1:n2] ne 0 or $
                 dvar[l2+1:*,m2+1:*,n1:n2] ne 0))[0] ge 0) then loc[3] = 'xy'
      if ((where(dvar[l1:l2,0:m1-1,0:n1-1] ne 0 or $
                 dvar[l1:l2,0:m1-1,n2+1:*] ne 0 or $
                 dvar[l1:l2,m2+1:*,0:n1-1] ne 0 or $
                 dvar[l1:l2,m2+1:*,n2+1:*] ne 0))[0] ge 0) then loc[4] = 'yz'
      if ((where(dvar[0:l1-1,m1:m2,0:n1-1] ne 0 or $
                 dvar[0:l1-1,m1:m2,n2+1:*] ne 0 or $
                 dvar[l2+1:*,m1:m2,0:n1-1] ne 0 or $
                 dvar[l2+1:*,m1:m2,n2+1:*] ne 0))[0] ge 0) then loc[5] = 'zx'
      ;; corners
      if ((where(dvar[0:l1-1,0:m1-1,0:n1-1] ne 0 or $
                 dvar[0:l1-1,0:m1-1,n2+1:*] ne 0 or $
                 dvar[0:l1-1,m2+1:*,0:n1-1] ne 0 or $
                 dvar[0:l1-1,m2+1:*,n2+1:*] ne 0 or $
                 dvar[l2+1:*,0:m1-1,0:n1-1] ne 0 or $
                 dvar[l2+1:*,0:m1-1,n2+1:*] ne 0 or $
                 dvar[l2+1:*,m2+1:*,0:n1-1] ne 0 or $
                 dvar[l2+1:*,m2+1:*,n2+1:*] ne 0))[0] ge 0) then loc[6] = 'xyz'

      print, label+': ', dv_glob, dv_nocorn, dv_inner, '  ', loc
    endfor

endfor




end
; End of file
