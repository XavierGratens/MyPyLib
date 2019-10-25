import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack



"Simulation of the SdH vs B trace"
"Number of points"
NP=1750
"Frequency"
F=157
"Decay rate constant"
d = 25
i=np.arange(NP)
"Applied magnetic field"
B = (i+5)/100 
Rxx=np.cos(2*np.pi*F/B)*1/(np.sinh(d/B))+np.random.normal(size=NP)*0.020

"Interpolation of the data Rxx vs 1/B"
plt.subplot(3,1,2)
plt.xlabel('1/B ($T^{-1}$)')
plt.ylabel('$R_{xx}$ ($\Omega$)')
Data= np.empty((B.shape[0],2))
Data[:, 1]= Rxx[:]
Data[:, 0]= 1/B[:]
TT=Data.transpose()
indice=np.argsort(TT)[:1]
Data[:]=Data[indice]
N = 1100
xN = Data[:N, 0]
yN= Data[:N, 1]
maxCol = np.amax(xN)
minCol = np.amin(xN)
Cste1=2
x_interp=np.arange(1024*Cste1)
x_interp =minCol +(maxCol-minCol)*x_interp/(1024*Cste1)
y_interp=np.interp(x_interp,xN,yN)

"Make FFT of the Rxx vs 1/B data"
yf = scipy.fftpack.fft(y_interp)
xf = np.linspace(0.0,1024*Cste1/(maxCol-minCol), 1024*Cste1)



"Plot of the data"

"Plot of the Rxx vs B trace"
plt.subplot(3,1,1)
plt.xlabel('B (T)')
plt.ylabel('$R_{xx}$ ($\Omega$)')
plt.plot(B,Rxx)


"Plot of the Rxx vs 1/B trace"
plt.subplot(3,1,2)
plt.xlabel('1/B ($T^{-1}$)')
plt.ylabel('$R_{xx}$ ($\Omega$)')
plt.plot(x_interp,y_interp)


"Plot of the FFT trace"
plt.subplot(3,1,3)
plt.xlabel('Frequency (T)')
plt.ylabel('Magnitude')
Cste2=12
plt.plot(xf[:1024//Cste2],np.abs(yf[:1024//Cste2]))

plt.subplots_adjust(top=0.98,hspace=0.500)

plt.show()

