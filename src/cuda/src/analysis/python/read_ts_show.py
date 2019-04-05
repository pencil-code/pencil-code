import numpy as np
import pylab as plt
import read_ts as rts 

#testing
tsfile = "data/ts.ac"
ts = rts.tsdata(tsfile)

#Some test plots
plt.figure()
plt.plot(ts.t, ts.dt, 'k-.', label='dt')
plt.xlabel("t")
plt.ylabel("dt")
legend = plt.legend(loc='upper right', shadow=True)
# Put a nicer background color on the legend.
legend.get_frame().set_facecolor('#00FFCC')

plt.figure()
plt.plot(ts.t, ts.urms, 'k-.', label='urms')
plt.plot(ts.t, ts.umax, 'k:', label='umax')
plt.xlabel("t")
plt.ylabel("u")
legend = plt.legend(loc='upper right', shadow=True)
# Put a nicer background color on the legend.
legend.get_frame().set_facecolor('#00FFCC')

plt.figure()
plt.plot(ts.t, ts.urms, 'k-.', label='urms')
plt.plot(ts.t, ts.uxrms, 'r-', label='uxrms')
plt.plot(ts.t, ts.uyrms, 'g--', label='uyrms')
plt.plot(ts.t, ts.uzrms, 'b:', label='uzrms')
plt.xlabel("t")
plt.ylabel("u")
legend = plt.legend(loc='upper right', shadow=True)
legend.get_frame().set_facecolor('#00FFCC')

plt.figure()
plt.plot(ts.t, ts.umax, 'k-.', label='umax')
plt.plot(ts.t, ts.uxmax, 'r-', label='uxmax')
plt.plot(ts.t, ts.uymax, 'g--', label='uymax')
plt.plot(ts.t, ts.uzmax, 'b:', label='uzmax')
plt.xlabel("t")
plt.ylabel("u")
legend = plt.legend(loc='upper right', shadow=True)
legend.get_frame().set_facecolor('#00FFCC')

plt.figure()
plt.plot(ts.t, ts.umin, 'k-.', label='umin')
plt.plot(ts.t, ts.uxmin, 'r-', label='uxmin')
plt.plot(ts.t, ts.uymin, 'g--', label='uymin')
plt.plot(ts.t, ts.uzmin, 'b:', label='uzmin')
plt.xlabel("t")
plt.ylabel("u")
legend = plt.legend(loc='upper right', shadow=True)
legend.get_frame().set_facecolor('#00FFCC')

plt.figure()
plt.plot(ts.t, ts.umin, 'k-.', label='umin')
plt.xlabel("t")
plt.ylabel("u")
legend = plt.legend(loc='upper right', shadow=True)
legend.get_frame().set_facecolor('#00FFCC')

plt.figure()
plt.plot(ts.t, ts.rhorms, 'r-', label='rhorms')
plt.plot(ts.t, ts.rhomax, 'g--', label='rhomax')
plt.plot(ts.t, ts.rhomin, 'b:', label='rhomin')
plt.xlabel("t")
plt.ylabel("rho")
legend = plt.legend(loc='upper right', shadow=True)
legend.get_frame().set_facecolor('#00FFCC')

plt.figure()
plt.plot(ts.t, ts.uxmax-ts.uymax, 'r-',  label='uxmax-uymax')
plt.plot(ts.t, ts.uxmax-ts.uzmax, 'g--', label='uxmax-uzmax')
plt.plot(ts.t, ts.uymax-ts.uzmax, 'b:',  label='uymax-uzmax')
plt.xlabel("t")
plt.ylabel("delta u")
legend = plt.legend(loc='upper right', shadow=True)
legend.get_frame().set_facecolor('#00FFCC')

plt.figure()
plt.plot(ts.t, ts.uxrms-ts.uyrms, 'r-',  label='uxrms-uyrms')
plt.plot(ts.t, ts.uxrms-ts.uzrms, 'g--', label='uxrms-uzrms')
plt.plot(ts.t, ts.uyrms-ts.uzrms, 'b:',  label='uyrms-uzrms')
plt.xlabel("t")
plt.ylabel("delta u")
legend = plt.legend(loc='upper right', shadow=True)
legend.get_frame().set_facecolor('#00FFCC')

plt.figure()
plt.plot(ts.t, np.abs(ts.uxmax), 'r-', label='abs(uxmax)')
plt.plot(ts.t, np.abs(ts.uymax), 'g-', label='abs(uymax)')
plt.plot(ts.t, np.abs(ts.uzmax), 'b-', label='abs(uzmax)')
plt.plot(ts.t, np.abs(ts.uxmin), 'r--', label='abs(uxmin)')
plt.plot(ts.t, np.abs(ts.uymin), 'g--', label='abs(uymin)')
plt.plot(ts.t, np.abs(ts.uzmin), 'b--', label='abs(uzmin)')
plt.xlabel("t")
plt.ylabel("u")
legend = plt.legend(loc='upper right', shadow=True)
legend.get_frame().set_facecolor('#00FFCC')

plt.figure()
plt.plot(ts.t, np.abs(ts.uxmax)-np.abs(ts.uxmin), 'r-',  label='uxmax-uxmin')
plt.plot(ts.t, np.abs(ts.uymax)-np.abs(ts.uymin), 'g-',  label='uymax-uymin')
plt.plot(ts.t, np.abs(ts.uzmax)-np.abs(ts.uzmin), 'b-',  label='uzmax-uzmin')
plt.xlabel("t")
plt.ylabel("delta u")
legend = plt.legend(loc='upper right', shadow=True)
legend.get_frame().set_facecolor('#00FFCC')


plt.show()


