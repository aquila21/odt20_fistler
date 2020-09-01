import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from HSTpost import HSTpost

rcParams['text.usetex']=True
rcParams['axes.labelsize'] =18
rcParams['xtick.labelsize'] = 14
rcParams['ytick.labelsize'] = 14
rcParams['axes.linewidth'] =1.4
params = {'text.latex.preamble' : [r'\usepackage{sfmath}']}
rcParams.update(params)
plt.rc('text', usetex=True)
plt.rc('font', family ='sans-serif')

nRlz = 2

drc0 = '../../rC_HST/'
psF0 = HSTpost(drc0)
psF0.createGasPhase()
psF0.gas.postProc(nRlz)
psF0.gas.energSpect(1001,nRlz)

k0 = psF0.gas.k
ep0 = psF0.gas.d
nu0 = psF0.gas.visc

print('#### Input ####')
print('Density: '+str(psF0.gas.rho))
print('Viscosity: '+str(psF0.gas.visc))
print('Shear Rate: '+str(psF0.Srate))

print('#### Output ####')
print('Turbulent kinetic energy: '+str(psF0.gas.k))
print('Dissipation rate: '+str(psF0.gas.d))
print('Reynolds number: '+str(2*k0*np.sqrt(5/(3*nu0*ep0))))
print('<u2>/k: '+str(psF0.gas.u2/psF0.gas.k))
print('<v2>/k: '+str(psF0.gas.v2/psF0.gas.k))
print('<w2>/k: '+str(psF0.gas.w2/psF0.gas.k))

