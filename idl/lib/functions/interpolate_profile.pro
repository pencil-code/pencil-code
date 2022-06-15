function interpolate_profile,data,data_z,n_data,z

;linear interpolation of data
	num_below = 0
	num_over = 0
	n_z = n_elements(z)
	out_profile = replicate(!Values.D_NaN, n_z)
	for j = 0, n_z-1 do begin
		if (z[j] lt data_z[0] ) then begin
			;extrapolate linarily below bottom
			num_below = num_below + 1
			out_profile[j] = data[0] + (data[1]-data[0])/(data_z[1]-data_z[0]) * (z[j]-data_z[0])
		end else if (z[j] gt data_z[n_data-1]) then begin
			;extrapolate linarily over top
			num_over = num_over + 1
			out_profile[j] = data[n_data-1] + (data[n_data-1]-data[n_data-2])/(data_z[n_data-1]-data_z[n_data-2]) * (z[j]-data_z[n_data-1])
		end else if (z[j] eq data_z[n_data-1]) then begin
			out_profile[j] = data[n_data-1]
		end 
	for i = 0, n_data-2 do begin
		if ((z[j] ge data_z[i]) and (z[j] lt data_z[i+1])) then begin
			;y = m*(x-x1) + y1
			out_profile[j] = (data[i+1]-data[i]) / (data_z[i+1]-data_z[i]) * (z[j]-data_z[i]) + data[i]
			break 
		endif
	endfor
	endfor

	if (num_below gt 0) then begin
		message,"extrapolated "+str(num_below)+" grid points below bottom", /info
	endif
	if (num_over gt 0) then begin
		message,"extrapolated "+str(num_over)+" grid points over top", /info
	endif

return, out_profile
 
end
 

