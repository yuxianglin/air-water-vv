from math import pi,tanh,sqrt,exp,log
import numpy as np
from matplotlib  import pyplot


def piersonMoskovitz(f,f0,alpha=8.1e-3,beta=0.74,g=9.8):
    return (5.0/16.0)*Hs**2*(f0**4/f**5)*np.exp((-5.0/4.0)*(f0/f)**4)

def sigma(omega,omega0):
    sigmaReturn = omega.copy()
    sigmaReturn[:] = 0.07
    sigmaReturn[omega > omega0] = 0.09
    return sigmaReturn

def JONSWAP(f,f0,Hs,g,gamma):
    omega = 2.0*pi*f
    omega0 = 2.0*pi*f0
    alpha = 2.0*pi*0.0624*(1.094-0.01915*log(gamma))/(0.23+0.0336*gamma-0.0185/(1.9+gamma))
    r = np.exp(- (omega-omega0)**2/(2*sigma(omega,omega0)**2*omega0**2))
    return (alpha*Hs**2*omega0**4/omega**5)*np.exp(-(5.0/4.0)*(omega0/omega)**4)*gamma**r


def dispersion(w,d, niter = 1000, g = 9.81):
# This function calculates k from linear dispersion relation
# d = depth, w is cyclical frequency and niter is the number of solution iteration (default = 1000)
    Kd = w*sqrt(d/g)
    for jj in range(niter):
        #Kdn_1 = Kd
        Kd = w*w*d/g/np.tanh(Kd)
        #Kdn_1 /=0.01*Kd
        #Kdn_1 -= 100.
        #Kdn_1 = abs(Kdn_1)
        #try: Kdn_1 = mean(Kdn_1)
        #except: continue
    #print "Solution convergence for dispersion relation %s percent" % Kdn_1
    return(Kd/d)

class RandomWaves:
    """
    Generate approximate random wave solutions

    :param Tp: peak period [T]
    :param Hs: significant wave height [L]
    :param  d: depth [L]
    :param fp: frequency [1/T]
    :param bandFactor: width factor for band  around fp [-]
    :param N: number of frequency bins [-]
    :param mwl: mean water level [L]
    """
    def __init__(self,
                 Tp = 5.0,         #s peak period
                 Hs = 2.0,         #m significant wave height
                 d = 2.0,           #m depth
                 fp = 1.0/5.0,      #peak  frequency
                 bandFactor = 2.0, #controls width of band  around fp
                 N = 101,          #number of frequency bins
                 mwl = 0.0,        #mean water level
                 g = 9.81,         #accelerationof gravity
                 spec_fun = JONSWAP): #wave spectrum
        self.Tp = Tp
        self.Hs = Hs
        self.d = d
        self.fp = fp
        self.bandFactor = bandFactor
        self.N = N
        self.mwl = mwl
        self.g = g
        self.fmax = self.bandFactor*self.fp
        self.fmin = self.fp/self.bandFactor
        self.df = (self.fmax-self.fmin)/float(self.N-1)
        self.fi=np.zeros(self.N,'d')
        for i in range(self.N):
            self.fi[i] = self.fmin+self.df*i
        self.ki = dispersion(2.0*pi*self.fi,self.d,g=self.g)
        self.phi = 2.0*pi*np.random.random(self.fi.shape[0])
        #ai = np.sqrt((Si_J[1:]+Si_J[:-1])*(fi[1:]-fi[:-1]))
        fim_tmp = (0.5*(self.fi[1:]+self.fi[:-1])).tolist()
        self.fim = np.array([fim_tmp[0]-0.5*self.df]+fim_tmp+[fim_tmp[-1]+0.5*self.df])
        self.Si_Jm = spec_fun(self.fim,f0=self.fp,Hs=self.Hs,g=self.g,gamma=3.3)
        self.ai = np.sqrt((self.Si_Jm[1:]+self.Si_Jm[:-1])*(self.fim[1:]-self.fim[:-1]))
    def eta(self,x,t):
        return (self.ai*np.cos(2.0*pi*self.fi*t - self.ki*x + self.phi)).sum()
    
    def u(self,x,z,t):
        Z = z - self.mwl
        return (2.0*pi*self.fi*self.ai*np.cos(2.0*pi*self.fi*t-self.ki*x+self.phi)*
                np.cosh(self.ki*(self.d+Z))/np.sinh(self.ki*self.d)).sum()

    def w(self,x,z,t):
        Z = z - self.mwl
        return (2.0*pi*self.fi*self.ai*np.sin(2.0*pi*self.fi*t-self.ki*x+self.phi)*
                np.sinh(self.ki*(self.d+Z))/np.sinh(self.ki*self.d)).sum()
    
if __name__ == '__main__':
    Tp = 5.0 #s peak period
    Hs = 2.0 #m significant wave height
    d = Hs #m depth
    fp = 1.0/Tp #peak  frequency
    bandFactor = 2.0 #controls width of band  around fp
    N = 101 # number of frequency bins
    mwl = 0.0 #mean water level    Hs = 2.0
    waves = RandomWaves(Tp = Tp,
                        Hs = Hs,
                        d = d,
                        fp = fp,
                        bandFactor = bandFactor,
                        N = N,
                        mwl = mwl,
                        g = 9.81)
    waves = RandomWaves(Tp = 1.94,
                        Hs = 0.1,
                        d  = 1.0,
                        fp = 1.0/1.94,
                        bandFactor = 2.0,
                        N = 101,
                        mwl = 1.0,
                        g = 9.8)#shouldn't mwl = d always?
    fp = 1.0/1.94
    d = 1.0
    Hs = 0.1
    Tp = 1.94
    
    kp = dispersion(2.0*pi*fp,d)
    nn = 51
    zmax = d+1.8*Hs/2.0
    xmax = 10*zmax
    T=0.0
    L = 6.0*pi/kp
    T_end = 10.0*Tp
    nt = 10*21
    x = np.linspace(0,xmax,nn)
    z = np.linspace(0,zmax,nn)
    X,Z = np.meshgrid(x,z)
    U = np.zeros((nn,nn),'d')
    W = np.zeros((nn,nn),'d')
    for n,T in enumerate(np.linspace(0,T_end,nt)):
        for J,xi in enumerate(x):
            for I,zi in enumerate(z):
                z_surf = d+waves.eta(xi,T)
                if zi < z_surf:
                    UIJ = waves.u(xi,zi,T)
                    WIJ = waves.w(xi,zi,T)
                    U[I,J] =  UIJ
                    W[I,J] = WIJ
                else:
                    U[I,J] = 0.0
                    W[I,J] = 0.0
        speed = np.sqrt(U*U + W*W)
        hl = [d+waves.eta(xi,T) for xi in x]
        pyplot.clf()
        #pyplot.axis('equal')
        pyplot.plot(x,hl)
        pyplot.streamplot(x, z, U, W,density=(1,1), color=speed,linewidth=2.5*speed/speed.max())
        #pyplot.contourf(x,z,speed)
        fig = pyplot.gcf()
        fig.set_size_inches(16.0,16.0*zmax/xmax)
        pyplot.savefig('frame%4.4d.png' % n)
    
