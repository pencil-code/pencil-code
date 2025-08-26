;
; $Id$
;+
;  Read in all snapshots stored in video slice.
;
;  Typical calling sequence:
;    rvid_plane_multiplot,'uz',times=times,planes=planes
;
;-
pro pc_read_videoslices,field,quiet=quiet, $
    extension=extension,datadir=datadir,oldfile=oldfile,debug=debug, $
    proc=proc,planes=planes,times=times,xx=xx,yy=yy, lplot=lplot
;
common pc_precision, zero, one, precision, data_type, data_bytes, type_idl
;
if (keyword_set(help)) then begin
  doc_library, 'rvid_plane'
  return
endif

default,quiet,0
default,lplot,0
default,single,0
default,extension,'xz'
;
default,swap_endian,0
default,field,'lnrho'
datadir = pc_get_datadir(datadir)
default,stride,0
default,dimfile,'dim.dat'
default,varfile,'var.dat'
default,title,'rvid_plane'
;
sample = ~keyword_set(interp)
;
if not check_slices_par(field, stride, arg_present(proc) ? datadir+'/proc'+str(proc) : datadir, s) then return
;
if (not any (tag_names (s) eq strupcase (strtrim (extension,2)+'read'))) then begin
  print, "rvid_plane: ERROR: slice '"+extension+"' is missing!"
  return
endif
;
tiny=1e-30 ; a small number
;
;  Read the dimensions from "dim.dat".
;
pc_read_dim, obj=dim, datadir=datadir, proc=proc, /quiet
pc_set_precision, dim=dim, /quiet
nx=dim.nx
ny=dim.ny
nz=dim.nz
mx=dim.mx
my=dim.my
mz=dim.mz
nghostx=dim.nghostx
nghosty=dim.nghosty
nghostz=dim.nghostz
nprocx=dim.nprocx
nprocy=dim.nprocy
nprocz=dim.nprocz
ncpus=nprocx*nprocy*nprocz
;
;  Read grid data.
;
pc_read_grid, obj=grid, proc=proc, swap_endian=swap_endian, /quiet, single=single
xx=grid.x(dim.l1:dim.l2)
yy=grid.y(dim.m1:dim.m2)
zz=grid.z(dim.n1:dim.n2)
;
;  Set reasonable extension for 2-D runs.
;
if ( (nx ne 1) and (ny ne 1) and (nz eq 1) ) then extension='xy'
if ( (nx ne 1) and (ny eq 1) and (nz ne 1) ) then extension='xz'
if ( (nx eq 1) and (ny ne 1) and (nz ne 1) ) then extension='yz'
;
if (is_defined(proc)) then begin
  file_slice=datadir+'/proc'+str(proc)+'/slice_'+field+'.'+extension
endif else begin
  file_slice=datadir+'/slice_'+field+'.'+extension
endelse
;
if not file_test(file_slice) then begin
  print, 'Slice file "'+file_slice+'" does not exist!!!'
  pos=strpos(file_slice,'.'+extension)
  compfile=strmid(file_slice,0,pos)+'1'+'.'+extension
  if file_test(compfile) then $
    print, 'Field name "'+field+'" refers to a vectorial quantity -> select component!!!'
  return
endif
if (not quiet) then print,'Will read "'+file_slice+'".'
;
;
;
t=zero
islice=0
;
if (extension eq 'xy') then plane=make_array(nx,ny, type=type_idl)
if (extension eq 'xy2') then plane=make_array(nx,ny, type=type_idl)
if (extension eq 'xy3') then plane=make_array(nx,ny, type=type_idl)
if (extension eq 'xz') then plane=make_array(nx,nz, type=type_idl)
if (extension eq 'yz') then plane=make_array(ny,nz, type=type_idl)
size_plane=size(plane)
;
slice_xpos=zero
slice_ypos=zero
slice_zpos=zero
slice_z2pos=zero
;
window, xsize=700, ysize=700, title=title
;
; Count number of times that are stored in the slices
;
nt=0
print, 'Reading "'+file_slice+'".'
openr, lun, file_slice, /f77, /get_lun, swap_endian=swap_endian
while (not eof(lun)) do begin
   if (keyword_set(oldfile)) then begin ; For files without position
      readu, lun, plane, t
   endif else begin
      readu, lun, plane, t, slice_z2pos
   endelse
   nt=nt+1
end
close, lun
free_lun, lun
;
print, 'Reading "'+file_slice+'".'
openr, lun, file_slice, /f77, /get_lun, swap_endian=swap_endian

if (extension eq 'xy') then  planes=make_array(nx,ny,nt, type=type_idl)
if (extension eq 'xy2') then planes=make_array(nx,ny,nt, type=type_idl)
if (extension eq 'xy3') then planes=make_array(nx,ny,nt, type=type_idl)
if (extension eq 'xz') then  planes=make_array(nx,nz,nt, type=type_idl)
if (extension eq 'yz') then  planes=make_array(ny,nz,nt, type=type_idl)
times=fltarr(nt)

it=0
while (not eof(lun)) do begin
   if (keyword_set(oldfile)) then begin ; For files without position
      readu, lun, plane, t
   end else begin
      readu, lun, plane, t, slice_z2pos
   end

   planes(*,*,it)=plane
   times(it)=t
   if (lplot) then begin
      contour,planes(*,*,it),xx,yy,nlev=50,/fill
   endif
   
   if (keyword_set(debug)) then begin
      print, t, min([plane,xy,xz,yz]), max([plane,xy,xz,yz])
   endif
   it=it+1 
end
close, lun
free_lun, lun
;
END
