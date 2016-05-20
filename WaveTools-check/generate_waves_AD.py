from proteus import WaveTools as WT
import numpy as np

Tp = 1.75
Hs = 0.112
mwl = 0.305
depth = 1.0

waves = WT.RandomWaves( Tp=Tp,
                        Hs = Hs,
                        mwl = mwl,
                        depth = depth,
                        bandFactor = 2.0,
                        N = 2001, 
                        waveDir = np.array([1,0,0]), 
                        g = np.array([0,-9.81,0]),
                        spectName = "JONSWAP",
                        spectral_params={"gamma": 3.3, "TMA": True,"depth": depth})


                                      
Tend = 1000                                      
tList = np.linspace(0.0,Tend , Tend /(Tp/50.))
etaList = tList.copy()

x = [0.,0.,0.]

ii = -1
for t in tList:
    ii +=1
    etaList[ii] = waves.eta(x,t)



a=np.array(zip(tList,etaList))

np.savetxt("waves3D.csv", a, delimiter=',')

freq_array = np.arange(0.01, 2.0, 0.01)
spectrum_true = WT.JONSWAP(f=freq_array,
                           f0=1.0/Tp,
                           Hs=Hs,
                           gamma=3.3,
                           TMA=True, 
                           depth = depth)
spec_array = np.array(zip(freq_array, spectrum_true))

np.savetxt("true_spectrum.csv", spec_array, delimiter=',')