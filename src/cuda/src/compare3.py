import numpy as np
import pylab as plt
import read_ts as rts 

#testing
tsfilex = "src/data/ts.ac"
tsfiley = "../3dgausy2/src/data/ts.ac"
tsfilez = "../3dgausz2/src/data/ts.ac"
tsx = rts.tsdata(tsfilex)
tsy = rts.tsdata(tsfiley)
tsz = rts.tsdata(tsfilez)

plt.figure()
#plt.plot(ts.t, ts.urms, 'k--', label='urms')
plt.plot(tsx.t, tsx.urms, 'r-', label='urms x')
plt.plot(tsy.t, tsy.urms, 'g-', label='urms y')
plt.plot(tsz.t, tsz.urms, 'b-', label='urms z')
#plt.plot(ts.t, ts.uxrms, 'r-', label='uxrms')
#plt.plot(ts.t, ts.uyrms, 'g-', label='uyrms')
#plt.plot(ts.t, ts.uzrms, 'b-', label='uzrms')
plt.xlabel("t")
plt.ylabel("u")
legend = plt.legend(loc='upper right', shadow=True)
legend.get_frame().set_facecolor('#00FFCC')

plt.figure()
plt.plot(tsx.t, tsx.uxrms, 'r-', label='uxrms x')
plt.plot(tsy.t, tsy.uxrms, 'g-', label='uxrms y')
plt.plot(tsz.t, tsz.uxrms, 'b-', label='uxrms z')
plt.xlabel("t")
plt.ylabel("u")
legend = plt.legend(loc='upper right', shadow=True)
legend.get_frame().set_facecolor('#00FFCC')

plt.figure()
plt.plot(tsx.t, tsx.uyrms, 'r-', label='uyrms x')
plt.plot(tsy.t, tsy.uyrms, 'g-', label='uyrms y')
plt.plot(tsz.t, tsz.uyrms, 'b-', label='uyrms z')
plt.xlabel("t")
plt.ylabel("u")
legend = plt.legend(loc='upper right', shadow=True)
legend.get_frame().set_facecolor('#00FFCC')

plt.figure()
plt.plot(tsx.t, tsx.uzrms, 'r-', label='uzrms x')
plt.plot(tsy.t, tsy.uzrms, 'g-', label='uzrms y')
plt.plot(tsz.t, tsz.uzrms, 'b-', label='uzrms z')
plt.xlabel("t")
plt.ylabel("u")
legend = plt.legend(loc='upper right', shadow=True)
legend.get_frame().set_facecolor('#00FFCC')

plt.figure()
plt.plot(tsx.t, np.abs(tsx.uyrms-tsx.uzrms), 'r-', label='uyrms-uzrms x')
plt.plot(tsy.t, np.abs(tsy.uxrms-tsy.uzrms), 'g-', label='uxrms-uzrms y')
plt.plot(tsz.t, np.abs(tsz.uxrms-tsz.uyrms), 'b-', label='uxrms-uyrms z')
plt.xlabel("t")
plt.ylabel("u")
legend = plt.legend(loc='upper right', shadow=True)
legend.get_frame().set_facecolor('#00FFCC')

plt.figure()
plt.plot(tsx.t, np.abs(tsx.uymax-tsx.uzmax), 'r-', label='uymax-uzmax x')
plt.plot(tsy.t, np.abs(tsy.uxmax-tsy.uzmax), 'g-', label='uxmax-uzmax y')
plt.plot(tsz.t, np.abs(tsz.uxmax-tsz.uymax), 'b-', label='uxmax-uymax z')
plt.xlabel("t")
plt.ylabel("u")
legend = plt.legend(loc='upper right', shadow=True)
legend.get_frame().set_facecolor('#00FFCC')


plt.show()


