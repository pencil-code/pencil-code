import numpy as np
import pylab as plt
import read_k_used as rkk

#testing
kkfile = "src/data/kk_used.ac"
kk = rkk.kkdata(kkfile)

#print kk.kk_vec_x 
#print kk.kk_vec_y 
#print kk.kk_vec_z 
#print kk.phi 
#print kk.forcing_kk_part_x 
#print kk.forcing_kk_part_y 
#print kk.forcing_kk_part_z 

plt.figure()
plt.subplot(1,3,1)
plt.plot(kk.kk_vec_x, 'r-', label='kk_vec_x')
plt.plot(kk.kk_vec_y, 'g-', label='kk_vec_y')
plt.plot(kk.kk_vec_z, 'b-', label='kk_vec_z')
legend = plt.legend(loc='upper right', shadow=True)
legend.get_frame().set_facecolor('#00FFCC')

plt.subplot(1,3,2)
plt.plot(kk.phi, 'k-', label='phi')
legend = plt.legend(loc='upper right', shadow=True)
legend.get_frame().set_facecolor('#00FFCC')

plt.subplot(1,3,3)
plt.plot(kk.forcing_kk_part_x, 'r-', label='forcing_kk_part_x')
plt.plot(kk.forcing_kk_part_y, 'g-', label='forcing_kk_part_y')
plt.plot(kk.forcing_kk_part_z, 'b-', label='forcing_kk_part_z')
legend = plt.legend(loc='upper right', shadow=True)
legend.get_frame().set_facecolor('#00FFCC')



plt.figure()
plt.subplot(1,3,1)
n, bins, patches = plt.hist(kk.kk_vec_x, 5, normed=0, facecolor='red', alpha=0.75)
plt.title("kk_vec_x")

plt.subplot(1,3,2)
n, bins, patches = plt.hist(kk.kk_vec_y, 5, normed=0, facecolor='green', alpha=0.75)
plt.title("kk_vec_y")

plt.subplot(1,3,3)
n, bins, patches = plt.hist(kk.kk_vec_z, 5, normed=0, facecolor='blue', alpha=0.75)
plt.title("kk_vec_z")



plt.figure()
n, bins, patches = plt.hist(kk.phi, 40, normed=0, facecolor='green', alpha=0.75)
plt.title("phi")



plt.figure()
plt.subplot(1,3,1)
n, bins, patches = plt.hist(kk.forcing_kk_part_x, 40, normed=0, facecolor='red', alpha=0.75)
plt.ylim([0, 60])
plt.title("forcing_kk_part_x")

plt.subplot(1,3,2)
n, bins, patches = plt.hist(kk.forcing_kk_part_y, 50, normed=0, facecolor='green', alpha=0.75)
plt.ylim([0, 60])
plt.title("forcing_kk_part_y")

plt.subplot(1,3,3)
n, bins, patches = plt.hist(kk.forcing_kk_part_z, 50, normed=0, facecolor='blue', alpha=0.75)
plt.ylim([0, 60])
plt.title("forcing_kk_part_z")

plt.show()


