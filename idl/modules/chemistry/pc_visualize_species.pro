;
; Visualize the species
;
; 2010.05.23/nils: coded
;
device,decompose=0
loadct,5
!x.style=1
!x.title='x [cm]'
!p.charsize=2
!p.multi=[0,1,1]
;
; Read in some data
;
pc_read_dim,obj=dim
;
; Read in data about species
;
@./data/pc_constants
print,'The available species are:',specname
;
; Set defaults
;
default,specie,'H2'
default,zoom,1
default,tmin,0.0
print,'specie=',specie
;
; Choose which variable to plot
;
index=-1
for i=0, n_elements(specname)-1 do begin
    if specname[i] eq specie then begin
        index=i
    endif
end
if (index eq -1) then begin
    print,'No such specie:', specie
    print,'The available species are:',specname
    stop
end
if (index eq 0) then begin
    variable='chemspec'
endif else begin
    variable='chemspec'+str(index+1)
end
;
; Show results
;
!y.title='Y [-]'
!p.title=specie
if (dim.nygrid eq 1 and dim.nzgrid eq 1) then begin
    rvid_line,variable,/nocontour,/xgrid,/global_scaling,tmin=tmin,zoom=zoom
endif else begin
    rvid_plane,variable,/xgrid,/global_scaling,tmin=tmin,zoom=zoom
end

END
