import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize


"Constants"

"Bohr magneton in cgs"
muB_cgs=9.2741E-21

"Boltzman´s constant in cgs"
kB_cgs=1.3807E-16
" Avogadro´s constant "
Na = 6.022E+23 

"Spin for Co2+"
Spin=3/2

"Landé Factors"
g = 2.263

"Molar mass for Zn(1-x)CoxO"
M_Zn =65.39
M_Co = 58.9332
M_O = 15.9994

"Import the data"
#Data = np.loadtxt("C:/....../Exp4KI.txt")
Data = np.loadtxt("Exp4KI.txt")

Data=Data
Data[:, 0] =Data[:,0]*10
Data[0,0] = 0.00001



"Modified brillouin function in cgs (emu/g). Two parameters fit, the concentration (x) and the tempearture (Teff)" 
def BF(h, x, Teff):
    return g*muB_cgs*Spin*x*Na/(M_Zn*(1-x)+M_Co*x+M_O)*((((2*Spin + 1)/(2*Spin))*1/(np.tanh((2*Spin + 1)*(g*muB_cgs*Spin*(h*1000)/(kB_cgs*Teff))/(2*Spin))) - (1/(2*Spin))*(1/(np.tanh((g*muB_cgs*Spin*(h*1000)/(kB_cgs*Teff))/(2*Spin))))))

def residualBF(p, h, y):
    return y - BF(h, *p)

             
p0 = [0.025,4]
popt, pcov = optimize.leastsq(residualBF, p0, args=(Data[:, 0], Data[:, 1]))
er=residualBF(popt, Data[:, 0], Data[:, 1])
mean_y = np.mean(Data[:, 1])
SStot = np.sum(np.square(Data[:, 1]-mean_y))
SSres = np.sum(np.square(er))
R_Square = 1-SSres/SStot
n=Data.shape[0]             
print('Modified Brillouin fit result', 'x = ', popt[0],', Teff = ', popt[1])
print('R_Square =', R_Square,'Radj_Square =', 1-(1-R_Square)*((n-1)/(n-2-1)))


i=np.arange(101)
xi=i/100*np.max(Data[:,0])
xi[0] = 1/100
yn = BF(xi, *popt)

fig1=plt.figure(1)
plt.xlabel('H (kOe)')
plt.ylabel('M (emu/g)')

plt.scatter(Data[:,0],Data[:,1], label = 'Exp.')


plt.plot(xi,yn, label='B. F. Fit')
plt.text(50, 0.6, '$\mathrm{Zn_{1-x}Co_xO}$', fontsize=10)
plt.text(50, 0.5, 'T = 4 K', fontsize=10)
plt.text(50, 0.4, '$x$ = %s' % np.round(popt[0],5) , fontsize=10)
plt.text(50, 0.3, '$T_{eff}$ = %s K' % np.round(popt[1],2) , fontsize=10)
plt.text(50, 0.2, '$R^{2}$ = %s ' % np.round(R_Square,5), fontsize=10)
plt.legend()

plt.show()

  
