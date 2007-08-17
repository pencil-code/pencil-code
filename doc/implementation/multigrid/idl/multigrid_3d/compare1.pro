;;;;;;;;;;;;;;;;;;;;
;;;   toto.pro   ;;;
;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ucalgary.ca)
;;;  Date:   29-May-2007
;;;
;;;  Description:
;;;    Compare Laplacian(f) to density for multigrid

.r multigrid
rr = sqrt( (grid.xx-cent[0])^2 + (grid.yy-cent[1])^2 + (grid.zz-cent[2])^2)

plot_binned, rr, -la, PSYM=4, $
             XR=[-1,1]*0.1, $
             YR=[10,600], YSTYLE=1, /YLOG
plot_binned, -rr, -la, PSYM=4, /OVER

plot_binned, rr, -g, PSYM=1, COLOR=120, /OVER

; End of file toto.pro
