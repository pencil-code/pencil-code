; General settings:
device, decompose=0
loadct, 5

; Process videofiles
default, videofiles_processed, 0
if (not videofiles_processed) then begin
	print, "Reading all videofiles."
	spawn, 'src/read_all_videofiles.x'
	videofiles_processed = 1
end

; Visualization
print, "Showing 'u_y'..."
rvid_plane, 'uu2', title='u_y', min=-7000, max=7000, zoom=4, /quiet
wait, 1
print, "Showing 'u_x'..."
rvid_plane, 'uu1', title='u_x', min=-7000, max=7000, zoom=4, /quiet
wait, 1
print, "Showing 'u_z'..."
rvid_plane, 'uu3', title='u_z', min=-7000, max=7000, zoom=4, /quiet
print, "Done."

END
