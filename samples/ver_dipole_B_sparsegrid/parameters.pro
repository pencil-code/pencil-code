t1=2000.
Om0=3.10652
nu=5e-5
eta=5e-5
chit=5e-5
nplot=120
run='Om3_Pr100_Pm100'
GM=3.
nt=2667000L
nt0=1000L
ntu=1000L
;Unit of time is years
tscale=0.62/2.7e-6/(365.*86400.)
;Unit of mgf is kG
Bscale=2.7e-6*7.e8/0.62*sqrt(4.*!pi*1.e-7*220.)*1.e4/1000.
; Unit of velocity is in km/s
Uscale=2.7e-6*7.e8/0.62/1000.
;
Omscale=1./(tscale*365.*86400.)
;For latitude strips
nt1=0 & nt2=127 & nt3=128 & nt4=255
;nt1=70 & nt2=80 & nt3=128 & nt4=255
;For radial strips
nr1=0 & nr2=127
wdsmth=147
av_ulim=1.1
av_llim=0.9
