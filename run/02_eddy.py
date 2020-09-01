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
eddyTime1 = data1[:,1]
eddyPos1 = data1[:,6]
eddyLength1 = data1[:,5]
diam = 0.007

plt.rc('text', usetex=True)
plt.rc('font', family='sans-serif')
plt.errorbar(eddyPos1[:],eddyTime1[:],xerr=eddyLength1*0.5, fmt=',',color='black', label='Eddy events', capsize=0, elinewidth=1.0)
plt.xlabel(r' r ')
plt.ylabel(r' t ')
plt.legend(numpoints=1)
#plt.xlim([-7, 7])
plt.savefig('eddies.pdf')
plt.show()
