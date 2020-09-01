import numpy as np

St = 5.0
ml = 0.8

k = 0.471
epsilon = 0.096
densRatio = 800
visc = 0.002 
L = 14

eta = np.power(np.power(visc,3)/epsilon, 0.25)
tauEta = np.sqrt(visc/epsilon)

d = np.sqrt(St*tauEta*visc/densRatio*18)
N = ml*L*6/(densRatio*np.power(d,3)*np.pi)

nP = np.power(N,2.0/3.0)
Nodt = ml*6*L/(densRatio*nP*np.pi*np.power(d,3))

print([d,nP,Nodt]) 

m = densRatio*np.power(d,3)*np.pi*nP*Nodt/14/6
Stoke = np.power(d,2)*densRatio/18/visc
print([m,St])

