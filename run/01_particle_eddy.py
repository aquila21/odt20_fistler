import matplotlib.pyplot as plt
import numpy as np

from matplotlib import rcParams
rcParams['text.usetex']=True
rcParams['axes.labelsize'] =18
rcParams['xtick.labelsize'] = 14
rcParams['ytick.labelsize'] = 14
rcParams['axes.linewidth'] =1.4
params = {'text.latex.preamble' : [r'\usepackage{sfmath}']}
rcParams.update(params)

filename1 = 'log'
data1 = np.genfromtxt(filename1, comments='#')
time = data1[:,1]
eddyPos = data1[:,6]
eddyLength = data1[:,5]
eddyTau = 0.08/data1[:,10]
nEddies = time.size
diam = 0.007

plt.close()
plt.rc('text', usetex=True)
plt.rc('font', family='sans-serif')

plt.errorbar(eddyPos[:],time[:],xerr=eddyLength*0.5,fmt=",",color='grey', label='eddy', capsize=0, elinewidth=0.3)

nP = 1
for i in range(0,nP):
    if i < 10:
        filename2 = '../data/data_particle/particle_No_000'+str(i)+'.dat'
    elif i >= 10 and i< 100:
        filename2 = '../data/data_particle/particle_No_00'+str(i)+'.dat'
    elif i >= 100 and i < 1000:
        filename2 = '../data/data_particle/particle_No_0'+str(i)+'.dat'
    else:
        print("Modification required!!")
    data2 = np.genfromtxt(filename2, comments='#')
    partTime1 = data2[:,0]
    partPos1 = data2[:,1]
    plt.plot(partPos1[:],partTime1[:], color = 'k')
    #plt.plot(partPos1[:],partTime1[:], label='part '+str(i))


plt.xlim([-7,7])
plt.ylim([0,25])
plt.ylabel(r' z/D ')
plt.xlabel(r' r/D ')
#plt.title(r'\large Particle trajectory and eddies in jet test case')
#plt.legend(loc=1)
#plt.show()
plt.savefig('particle_traj.pdf')
plt.show()
