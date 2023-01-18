pro pc_read_saffman,tt,iir,choice=choice
;
;  Simple routine for reading Saffman integrals.
;  Here with hardwired file name and file size.
;
;  Choices:
;
;  data/Iv_bc_saffman_mag.dat
;
;  w(1) = 8*rr**3*(j0**2)      !  box counting, cubic, Eq.(2.7)
;  data/Iv_bcc_saffman_mag.dat
;
;  w(2) = 48*pi*rr**3*(j1**2)  !  box counting, spheric, Eq.(2.8)
;  data/Iv_bcs_saffman_mag.dat
;
;  w(3) = 8*rr**3*j0           !  spectral, cubic, Eq.(2.12)
;  data/Iv_spc_saffman_mag.dat
;
;  w(4) = 8*pi*rr**3*j1        !  spectral, spheric, Eq.(2.13)
;  data/Iv_sps_saffman_mag.dat
;
if keyword_set(choice) then begin
  CASE choice of
    1: file='data/Iv_bcc_saffman_mag.dat'
    2: file='data/Iv_bcs_saffman_mag.dat'
    3: file='data/Iv_spc_saffman_mag.dat'
    4: file='data/Iv_sps_saffman_mag.dat'
    ELSE: print,'invalid choice'
  ENDCASE
endif else begin
  file='data/Iv_bcs_saffman_mag.dat'
endelse
;
;  determine resolution
;
pc_read_dim,obj=dim
n=dim.nx/2
print,'AXEL: n=',n
;
ir=fltarr(n)
openr,1,file
;
;  Initialize counter, read, and append, until end of file (eof):
;
it=0L
while not eof(1) do begin
  readf,1,t
  readf,1,ir
  print,t
  if it eq 0 then begin
    tt=t
    iir=ir
  endif else begin
    tt=[tt,t]
    iir=[iir,ir]
  endelse
  it=it+1
endwhile
nt=it
;
;  Reform
;
iir=reform(iir,n,nt)
close,1
END
