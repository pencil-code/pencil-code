;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   write_bb_for_dx.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ucalgary.ca)
;;;  Date:   11-Jun-2007
;;;
;;;  Description:
;;;    Write magntic field with extrapolated potential field for
;;;    visualizaton with OpenDX.
;;;  Usage:
;;;    .r start
;;;    ; varfile='VAR5'    ; or such, if desired
;;;    .r rall
;;;    .r write_bb_for_dx

pc_potentialfield_exterior_z, aa, $
                              AAA=aa_expanded, ZZZ=z_expanded
bb_expanded = curl(aa_expanded)
dx_export, bb_expanded, X=x, Y=y, Z=z_expanded, T=t, $
    BASENAME='bb', LABELS='BB', GHOSTS=3



end
; End of file write_bb_for_dx.pro
