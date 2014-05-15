; General settings:
default, wait, 0.1
default, zoom, 4

device, decompose=0
loadct, 5

; Process videofiles
spawn, 'pc_videofiles'

; Visualization
print, "Showing 'u_y'..."
rvid_plane, 'uu2', title='u_y', min=-7000, max=7000, zoom=zoom, wait=wait, /quiet
wait, 1
print, "Showing 'u_x'..."
rvid_plane, 'uu1', title='u_x', min=-7000, max=7000, zoom=zoom, wait=wait, /quiet
wait, 1
print, "Showing 'u_z'..."
rvid_plane, 'uu3', title='u_z', min=-7000, max=7000, zoom=zoom, wait=wait, /quiet
print, "Done."

END
