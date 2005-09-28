;;;;;;;;;;;;;;;;;;;;;;;;
;;;   boxbtest.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ucalgary.ca)
;;;  Date:   28-Sep-2005
;;;
;;;  Description:
;;;    Test boxbotex_scl (Anders reports to get different layout in IDL 6.0,
;;;    6.1 and 6.2)

;; Construct some data (in interval [0,1]) to plot
n=50
seed=!pi/2
s1=my_smooth(randomu(seed,n,n)-0.5,15)
m=minmax(s1)
s1=(s1-m[0])/(m[1]-m[0])
s2=s1^2
s3=s2^2
s4=s3^2

wdwset, 0              ; need to have window open ere calling boxbote_scl
boxbotex_scl,s1,s2,s3,s4,2,2,ip=3,amin=0.,amax=1.


end
; End of file boxbtest.pro
