import numpy as np
import matplotlib.pyplot as plt
from HSTpost import HSTpost

nRlz = 2

# Single-phase simulation
drc0 = '../../rC_HST/'
psF0 = HSTpost(drc0)
psF0.createGasPhase()
psF0.gas.postProc(nRlz)
psF0.gas.energSpect(1001,nRlz)

k0 = psF0.gas.k
ep0 = psF0.gas.d
nu0 = psF0.gas.visc

print('#### Single-phase HST ####')
print('#### Input ####')
print('Density: '+str(psF0.gas.rho))
print('Viscosity: '+str(psF0.gas.visc))
print('Shear Rate: '+str(psF0.Srate))

print('\n#### Output ####')
print('Turbulent kinetic energy: '+str(psF0.gas.k))
print('Dissipation rate: '+str(psF0.gas.d))
print('Reynolds number: '+str(2*k0*np.sqrt(5/(3*nu0*ep0))))
print('<u2>/k: '+str(psF0.gas.u2/psF0.gas.k))
print('<v2>/k: '+str(psF0.gas.v2/psF0.gas.k))
print('<w2>/k: '+str(psF0.gas.w2/psF0.gas.k))

tau_eta = np.power(psF0.gas.visc/psF0.gas.d, 0.5)

# Particle-laden simulation
drc = '../../rC_part_HST/'
psF = HSTpost(drc)
psF.createGasPhase()
psF.gas.postProc(nRlz)
psF.gas.energSpect(1001,nRlz)

psF.createParticlePhase()
psF.part.postProc(nRlz)

tau_p = np.power(psF.part.diam,2)*psF.part.rho/18/psF.gas.visc/psF.gas.rho
phi = (np.pi/6*np.power(psF.part.diam,3))*psF.part.Npart*psF.part.rho/psF.LStat/psF.gas.rho

print('\n#### Particle-laden HST ####')
print('#### Input ####')
print('No. of particles: '+str(psF.part.Npart))
print('Density of particle: '+str(psF.part.rho))
print('Particle diameter: ' +str(psF.part.diam))

print('\n#### Based on single-phase values ####')
print('Stokes no.: '+str(tau_p/tau_eta))
print('Mass loading: '+str(phi))

print('\n#### Output ####')
print('Dissipation rate: '+str(psF.gas.d))
print('TKE of gas phase: '+str(psF.gas.k))
print('TKE of particles: ' +str(psF.part.kp))


plt.figure(0)
plt.loglog(psF0.gas.frq[1:], (psF0.gas.E[1:,0]+psF0.gas.E[1:,1]+psF0.gas.E[1:,2]), label=r'Single-phase HST', linewidth=2.0)
plt.loglog(psF.gas.frq[1:], (psF.gas.E[1:,0]+psF.gas.E[1:,1]+psF.gas.E[1:,2]), label=r'Particle-phase HST', linewidth=2.0)
plt.xlabel(r'$\kappa$')
plt.ylabel(r'$E(\kappa)$')
plt.legend()
plt.show()




