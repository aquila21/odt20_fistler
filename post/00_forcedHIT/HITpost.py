import yaml
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import rcParams

class HITpost:
    def __init__(self,drc):
        self.drc  = drc
        stream    = open(drc+'CPU_1/input/shearFlow/odt_input.yaml')
        self.docs      = yaml.load(stream)
        stream.close()
        self.tEnd = self.docs['params']['tEnd']
        self.dt = self.docs['dumpTimes'][0]
        self.nDt  = self.docs['params']['nStat']
        self.N  = self.docs['params']['ngrdStat']
        self.LStat = self.docs['params']['domainLength']
        posf = np.linspace(-self.LStat*0.5,0.5*self.LStat,self.N+1) # space vector
        dy = abs(posf[1]-posf[0])
        self.yStat = 0.5*dy+posf[:-1]
        self.tStat = np.linspace(0,self.tEnd, self.nDt)
        self.nDmp = int(self.docs['params']['tEnd']/self.docs['dumpTimes'][0])
        self.tDmp = np.linspace(0,self.tEnd,self.nDmp)

    def createGasPhase(self):
        self.gas = self.gasPhase(self)

    def createParticlePhase(self):
        self.part = self.particlePhase(self)

    class gasPhase:
        def __init__(self, supInst):
            self.supInst = supInst
            self.visc    = supInst.docs['params']['kvisc0']
            self.rho     = supInst.docs['params']['rho0']
            self.Prod    = supInst.docs['params']['Prod']
            self.N  = supInst.docs['params']['ngrdStat']
            self.LStat = supInst.docs['params']['domainLength']
            posf = np.linspace(-self.LStat*0.5,0.5*self.LStat,self.N+1) # space vector
            dy = abs(posf[1]-posf[0])
            self.yStat = 0.5*dy+posf[:-1]
            self.nCPU = len(os.listdir(supInst.drc))

        def postProc(self, nRlz):
            self.d_diff = np.zeros([self.yStat.size,self.supInst.nDt])
            self.k1     = np.zeros([self.yStat.size,self.supInst.nDt])
            for j in range(nRlz):            
                iCPU = int(self.nCPU*j/nRlz)+1
                path = self.supInst.drc + 'CPU_'+ str(iCPU)+'/data/data_'+ str(j+1)+'/data/data-average-output_'
                ext = '.dat'
                for i in range(0,self.supInst.nDt):
                    if i < 10:
                        filename = path+ '0000' + str(i) + ext
                    elif i < 100:
                        filename = path+ '000'+ str(i) + ext
                    elif i < 1000:
                        filename = path+ '00'+ str(i) + ext
                    elif i < 10000:
                        filename = path+ '0' + str(i) + ext
                    else:
                        filename = path+ str(i) + ext
                    #print(filename)
                    data = np.genfromtxt(filename)                
                    self.d_diff[:,i] += data[:,9] + data[:,11]
                    self.k1[:,i]      += data[:,7]
                    
            self.d_diff   /= nRlz               
            self.k1        /= nRlz
            self.d_diff_t  = np.zeros(self.supInst.nDt)
            self.k1_t  = np.zeros(self.supInst.nDt)
            for l in range(0,self.supInst.nDt):
                self.d_diff_t[l] = np.average(self.d_diff[:,l])
                self.k1_t[l] = np.average(self.k1[:,l])
            self.d  = -1*np.average(self.d_diff_t)              
            self.k2 = np.average(self.k1_t)              
    
        def energSpect(self, nCells, nRlz):
            # Frequency
            #yFFT = np.linspace(-0.5*self.LStat, 0.5*self.LStat, 4001) 
            yFFT = np.linspace(-np.pi, np.pi, nCells)
            yFFT = yFFT[:-1]+0.5*(yFFT[1]-yFFT[0])
            halfN = int(0.5*yFFT.size)
            self.frq = np.arange(halfN)#/self.LStat
            Ex = np.zeros(halfN)
            Ey = np.zeros(halfN)
            Ez = np.zeros(halfN)
            self.fname = 'ddd'
            self.E = np.zeros([halfN,3])
            nDmp = int(self.supInst.docs['params']['tEnd']/self.supInst.docs['dumpTimes'][0])
            for i in range(nRlz):
                iCPU = int(self.nCPU*i/nRlz)+1
                path = self.supInst.drc + 'CPU_'+ str(iCPU)+'/data/data_'+ str(i+1)+'/data/data_00000/dmp_'
                ext_path = '.dat'
                for j in range(0,nDmp):
                    if j < 10:
                        filename = path+ '0000' + str(j) + ext_path
                    elif j < 100:
                        filename = path+ '000'+ str(j) + ext_path
                    elif j < 1000:
                        filename = path+ '00'+ str(j) + ext_path
                    elif j < 10000:
                        filename = path+ '0' + str(j) + ext_path
                    else:
                        filename = path+ str(j) + ext_path
                    self.fname = filename  
                    tmpEx, tmpEy, tmpEz = self.FFT_TKE(filename, yFFT)
                    Ex += tmpEx
                    Ey += tmpEy
                    Ez += tmpEz
            self.E[:,0] = Ex/(nDmp*self.nCPU)
            self.E[:,1] = Ey/(nDmp*self.nCPU)
            self.E[:,2] = Ez/(nDmp*self.nCPU)
            self.u2 = 2*np.trapz(self.E[:,0])
            self.v2 = 2*np.trapz(self.E[:,1])
            self.w2 = 2*np.trapz(self.E[:,2])
            self.k = np.trapz(self.E[:,0]+self.E[:,1]+self.E[:,2])
            
        def FFT_TKE(self, filename, y):
            data      = np.genfromtxt(filename)
            N         = y.size
            halfN     = int(0.5*N)
            yRlz      = data[:,0]
            uRlz      = data[:,4]
            vRlz      = data[:,5]
            wRlz      = data[:,6]
            u         = np.interp(y,yRlz,uRlz)
            v         = np.interp(y,yRlz,vRlz)
            w         = np.interp(y,yRlz,wRlz)
            # FFT for each velocity component
            U         = np.fft.fft(u)/N # fft computing and normalization
            U         = np.sqrt(U*np.conj(U))
            Ex        = np.power(U[:halfN],2)
            Ex[0]    /= 2.0
            V         = np.fft.fft(v)/N # fft computing and normalization
            V         = np.sqrt(V*np.conj(V))
            Ey        = np.power(V[:halfN],2)
            Ey[0]    /= 2.0
            W         = np.fft.fft(w)/N # fft computing and normalization
            W         = np.sqrt(W*np.conj(W))
            Ez        = np.power(W[:halfN],2)
            Ez[0]    /= 2.0
            return np.real(Ex), np.real(Ey), np.real(Ez)


        def outputSpectr(self, filename):
            eta = np.power(np.power(self.visc,3)/self.d, 0.25)
            #Open output file
            file = open(filename,"w")
            #Write header
            file.write("# Energy spectrum of homogeneous shear turbulence")
            header = ['Grid points: ' + str(self.supInst.N), 'Kin. viscosity: '+str(self.visc), 'Total kin. energy: '+str(format(self.k, '.5f')), 'Turbulent dissipation: '+str(format(self.d,'.5f')), 'Kolmogorov length scale: '+str(format(eta,'.5f'))]
            for n,value in enumerate(header,1):
                file.write("\n# "+ value)
            file.write("\n#")
            dataname = ['f', 'Ex(f)', 'Ey(f)', 'Ez(f)', 'E(f)']
            for n,value  in enumerate(dataname,1):
                if n == 1:
                    file.write(("["+str(n-1)+",:]: "+ value).ljust(17))
                else:
                    file.write(("["+str(n-1)+",:]: "+ value).ljust(18))
            file.write("\n")
            data = np.column_stack((self.frq[1:]*eta, self.E[1:,0], self.E[1:,1], self.E[1:,2], self.E[1:,0]+self.E[1:,1]+self.E[1:,2]))
            n = data.shape[0]
            m = data.shape[1]
            for i in range(0,n):
                for j in range(0,m):
                    file.write(str("{:011.8f}".format(data[i,j])).ljust(18))
                file.write("\n")
            file.close()

        def flowCharac(self):
            q = np.power(2/3*self.k,0.5)
            self.dictFlowCharac = {'k': self.k,  'epsilon': self.d, 'uPrime': q}
            self.dictFlowCharac['eta'] = np.power(np.power(self.visc,3)/self.d, 0.25)
            self.dictFlowCharac['tau_eta'] = np.power(self.visc/self.d, 0.5)
            self.dictFlowCharac['L_int'] = np.power(q,3)/self.d
            self.dictFlowCharac['L'] = np.power(self.k,1.5)/self.d
            self.dictFlowCharac['L_taylor'] = np.sqrt((15*self.visc*q**2)/self.d)
            self.dictFlowCharac['Re_taylor'] = q*self.dictFlowCharac['L_taylor']/self.visc
            print(self.dictFlowCharac)
            return


    class particlePhase:
        def __init__(self, supInst):
            self.supInst = supInst
            self.Npart       = supInst.docs['params']['Nparticle']*supInst.docs['initParticleParams'][0][0]
            self.rho     = supInst.docs['params']['DensiParticle']
            self.diam    = supInst.docs['params']['Dparticle']

            self.nCPU = len(os.listdir(supInst.drc))

        def postProc(self, nRlz):
            self.kp = np.zeros(self.supInst.nDt)
            for j in range(nRlz):            
                iCPU = int(self.nCPU*j/nRlz)+1
                path = self.supInst.drc + 'CPU_'+ str(iCPU)+'/data/data_'+ str(j+1)+'/particle_data/particle_dmp_'
                ext = '.dat'
                for i in range(0,self.supInst.nDt):
                    if i < 10:
                        filename = path+ '0000' + str(i) + ext
                    elif i < 100:
                        filename = path+ '000'+ str(i) + ext
                    elif i < 1000:
                        filename = path+ '00'+ str(i) + ext
                    elif i < 10000:
                        filename = path+ '0' + str(i) + ext
                    else:
                        filename = path+ str(i) + ext
                    #print(filename)
                    data = np.genfromtxt(filename)                
                    self.kp[i] += np.average(0.5*(np.power(data[:,6],2)+np.power(data[:,7],2)+np.power(data[:,8],2)))
                    
            self.kp   /= nRlz 
            self.kp   = np.average(self.kp)
 





















