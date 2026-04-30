;
;  plot_breakup.pro
;
;  Visualise the falling-droplet breakup run defined by start.in / run.in.
;
;  The current setup releases a single 20 cm water drop at rest under gravity
;  in still air, with KH and RT breakup enabled (particles_breakup.f90,
;  Aguerre & Nigro 2019).  The expected outcome was a steady shower of
;  rain-sized droplets (radius ~ 0.05-0.3 cm).  This script reads the
;  Pencil Code time-series diagnostics and shows what the run actually
;  produces:
;
;    * panel 1 - droplet radius a_p versus time (log scale).  A grey band
;                marks the typical raindrop range (radius 0.05-0.3 cm) and
;                a dashed line marks the maximum stable raindrop radius
;                (~0.4 cm, Pruppacher & Klett).  The trace plunges
;                straight through the rain band into the cloud-droplet
;                regime.
;    * panel 2 - vertical particle velocity, showing the collapse of the
;                terminal velocity once the radius has been driven below
;                rain-droplet size.
;    * panel 3 - gas-phase Weber number We_g = rho_g U_r^2 (2 a_p) / sigma
;                with the canonical breakup-regime thresholds overlaid
;                (We=12 bag, We=50 bag-stamen, We=350 catastrophic).
;                The trace shows that the run sits in the bag-breakup
;                regime for a long time without anything happening,
;                then suddenly jumps into the RT-active range.
;    * panel 4 - number of physical droplets represented by the swarm,
;                npar_found * npswarmm.  RT breakup keeps a single
;                parcel and inflates npswarm, so this is the right
;                count to look at.
;
;  All quantities are plotted in CGS - the units used by start.in/run.in.
;
;  Usage (from the run directory):
;
;    idl
;    .r plot_breakup
;
;  A PNG copy is also written to plot_breakup.png.
;
;-----------------------------------------------------------------------------

pc_read_ts, obj=ts

; Physical constants of the setup (CGS), copied from start.in / run.in.
rho_g    = 1.225e-3       ; air density            [g/cm^3]
sigma    = 72.0           ; water surface tension  [dyn/cm]
ap_rain_lo = 0.05         ; raindrop radius lower edge [cm]  (~ 1 mm diameter)
ap_rain_hi = 0.30         ; raindrop radius upper edge [cm]  (~ 6 mm diameter)
ap_rain_max = 0.40        ; max stable raindrop radius [cm]  (~ 8 mm diameter)

t        = ts.t
ap_max   = ts.apmax
ap_min   = ts.apmin
ap_mean  = ts.apm
vpz      = ts.vpzm
np_phys  = ts.npar_found * ts.npswarmm

; Weber number based on the *gas* phase, mean radius and mean parcel speed
; (in this setup vpzm = vpzmin = vpzmax because there is one active parcel).
ur       = abs(vpz) > 1.0e-30
we_g     = rho_g * ur^2 * (2.0*ap_mean) / sigma  > 1.0e-30

; ---------------------------------------------------------------------------
; Plot to screen.
; ---------------------------------------------------------------------------
set_plot, 'x'
device, decomposed=0
loadct, 39, /silent
!p.multi      = [0, 1, 4]
!p.charsize   = 2.4
!p.charthick  = 1.6
!p.thick      = 1.6
!x.thick      = 1.6
!y.thick      = 1.6
window, 0, xsize=900, ysize=900, title='Falling water drop with KH/RT breakup'

; Panel 1 - radius vs time, log scale, with raindrop reference band.
ymin = min([ap_min, ap_rain_lo/3.0]) > 1.0e-5
ymax = max([ap_max, ap_rain_max*2.0])
plot, t, ap_max, /ylog, yrange=[ymin, ymax], ystyle=1, $
      xtitle='t  [s]', ytitle='a_p  [cm]', $
      title='droplet radius (target: rain-drop band, grey)'
; rain-drop band
polyfill, [t[0], t[n_elements(t)-1], t[n_elements(t)-1], t[0]], $
          [ap_rain_lo, ap_rain_lo, ap_rain_hi, ap_rain_hi], $
          color=220, /data, noclip=0
oplot, t, ap_max,  thick=2, color=250
oplot, t, ap_mean, thick=2, color=100
oplot, t, ap_min,  thick=2, color=60
oplot, t, replicate(ap_rain_max, n_elements(t)), linestyle=2, color=0
xyouts, t[1], ap_rain_max*1.15, '!17 max stable raindrop !6', charsize=1.6
xyouts, t[1], sqrt(ap_rain_lo*ap_rain_hi), '!17 raindrop range !6', charsize=1.6

; Panel 2 - vertical velocity.
plot, t, vpz, xtitle='t  [s]', ytitle='v_z  [cm/s]', $
      title='vertical particle velocity (terminal v collapses with a_p)'
oplot, t, replicate(0.0, n_elements(t)), linestyle=1

; Panel 3 - Weber number with breakup-regime thresholds.
ymax_we = max([we_g, 1000.0])
plot, t, we_g, /ylog, yrange=[1.0e-2, ymax_we], ystyle=1, $
      xtitle='t  [s]', ytitle='We_g', $
      title='gas Weber number (bag 12 / bag-stamen 50 / catastrophic 350)'
oplot, t, replicate(12.0,  n_elements(t)), linestyle=2, color=80
oplot, t, replicate(50.0,  n_elements(t)), linestyle=2, color=160
oplot, t, replicate(350.0, n_elements(t)), linestyle=2, color=250
xyouts, t[1], 13.0,  'bag', charsize=1.6
xyouts, t[1], 55.0,  'bag-stamen', charsize=1.6
xyouts, t[1], 380.0, 'catastrophic', charsize=1.6

; Panel 4 - number of physical droplets represented (log scale).
ymin_n = (min(np_phys) > 1.0)*0.5
ymax_n = max(np_phys)*2.0
plot, t, np_phys, /ylog, yrange=[ymin_n, ymax_n], ystyle=1, $
      xtitle='t  [s]', ytitle='N_phys = npar_found * npswarm', $
      title='physical droplet count carried by the swarm (RT inflates npswarm)'

; ---------------------------------------------------------------------------
; Also dump a PNG so the figure can be inspected outside an IDL session.
; ---------------------------------------------------------------------------
img = tvrd(true=1)
write_png, 'plot_breakup.png', img

; Restore defaults.
!p.multi     = 0
!p.charsize  = 1.0
!p.charthick = 1.0
!p.thick     = 1.0
!x.thick     = 1.0
!y.thick     = 1.0

; ---------------------------------------------------------------------------
; A short text summary of the cascade endpoint, printed to the IDL console.
; ---------------------------------------------------------------------------
print, ''
print, '=== particles-breakup run summary (CGS) ============================='
print, format='("initial radius      a_p(t=0)   = ",e10.3," cm")', ap_max[0]
print, format='("final radius        a_p(t_end) = ",e10.3," cm")', ap_max[n_elements(t)-1]
print, format='("typical raindrop    a_p        = ",f6.3,"  cm  (1-6 mm dia)")', $
       0.5*(ap_rain_lo+ap_rain_hi)
print, format='("final v_z                       = ",e10.3," cm/s")', vpz[n_elements(t)-1]
print, format='("final We_g                      = ",e10.3)', we_g[n_elements(t)-1]
print, format='("final N_phys                    = ",e10.3)', np_phys[n_elements(t)-1]
print, '====================================================================='
print, 'wrote plot_breakup.png'

END
