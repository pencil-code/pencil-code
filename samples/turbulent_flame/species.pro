;
; Visualize the species
;
device,decompose=0
loadct,5
!p.charsize=2
!p.multi=[0,1,1]
;
; Read in data about species
;
@./data/pc_constants
print,'The available species are:',specname
;
; Set defaults
;
default,specie,'H2'
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
!p.title=specie
rvid_plane,variable,/nocontour,/xgrid,/global_scaling



END
