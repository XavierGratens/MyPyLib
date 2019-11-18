import numpy as np
import matplotlib.pyplot as plt

"MEAN FIELD THEORY: Weiss model of a ferromagnet with the effect of external magnetic field" 


"kB is the Boltsman constant"
kB=1.3807E-23
"muB is the Bohr magneton"
muB=9.2741E-24
"g is the Landé factor"
g = 2
"S is the spin"
S = 7/2
"Msat is the magnetization saturation value: Msat = 1221482.9922947765 A/m for EuS"
Msat=1221482.99


"Numerically finding the roots of the equation: B_S(S,x) -(S+1)*g*muB/(3*kB*Tc)*[kB*T/(g*muB*S)*x-B])" 

"B_S(S,x) is the Brillouin function"


def B_S(S,y):
    return (((2*S + 1)/(2*S))*1/(np.tanh((2*S + 1)*y/(2*S))) - (1/(2*S))*(1/(np.tanh(y/(2*S)))))
    


"Dichotomy funtion"
def Mybisect(f,a,b,c,d,e,x0,x1,error=10e-12):
    x = (x0 + x1) / 2
    n = 0
    while (x1 - x0) / 2 > error:
        n += 1
        if f(a,b,c,d,e,x) == 0:
            return x
        elif f(a,b,c,d,e,x0) * f(a,b,c,d,e,x) < 0:
            x1 = x
        else:
            x0 = x
        x = (x0 + x1) / 2
    return B_S(b,x)
    
def MdeW(g,S,T,Tc,B,x):
    return (B_S(S,x))-((S+1)*g*muB/(3*kB*Tc)*(kB*T/(g*muB*S)*x-B))
    
"Calculation of the volumic magnetization as function of the temperature in SI units "
def mvsT(g,S,Tc,B,NP):
    Result = np.zeros(shape=(NP,2))
    for i in range(NP):
        Result[i,0]=i/NP*2*Tc
        Result[i,1]=Mybisect(MdeW,g,S,Result[i,0],Tc,B,0.001,100)*Msat
    return Result    


"Numerically calculation of the susceptibility in SI units"
def XvsT(g,S,Tc,B,DeltaB,NP):
    Result = np.zeros(shape=(NP,2))
    for i in range(NP):
        Result[i,0]=i/NP*2*Tc
        Result[i,1]=(Mybisect(MdeW,g,S,Result[i,0],Tc,B+DeltaB,0.001,100)*Msat-Mybisect(MdeW,g,S,Result[i,0],Tc,B-DeltaB,0.001,100)*Msat)*4*np.pi*1E-7/(2*DeltaB)
    return Result    



"Number of points"
numbP=1000
"External magnetic fields in Tesla"
H0=0.01
H1=0.035
H2=0.07

"Curie temperature in Kelvin"
Tc_exp = 15.8




"Calculation of the magnetization as function of the temperature for 3 values of the applied magnetic field"
M0 = mvsT(2,S,Tc_exp,H0,numbP)
M1 = mvsT(2,S,Tc_exp,H1,numbP)
M2 = mvsT(2,S,Tc_exp,H2,numbP)
"Calculation of the susceptibility as function of the temperature for 3 values of the applied magnetic field"
X0=XvsT(2,S,Tc_exp,H0,0.001,numbP)
X1=XvsT(2,S,Tc_exp,H1,0.001,numbP)
X2=XvsT(2,S,Tc_exp,H2,0.001,numbP)

"Plot M vs T"
plt.subplot(2,1,1)
plt.ylabel('Magnetization ($10^5$ A/m)')
plt.plot(M0[:,0],M0[:,1]/1E5, label='H =%s T'% H0)
plt.plot(M1[:,0],M1[:,1]/1E5,label='H =%s T'% H1 )
plt.plot(M2[:,0],M2[:,1]/1E5, label='H =%s T'% H2)
plt.text(2, 6, '$T_c$ = %s K' % Tc_exp , fontsize=10)
plt.legend()

"Plot X vs T"
plt.subplot(2,1,2)
plt.xlabel('T (K)')
plt.ylabel('Susceptibility')
plt.plot(X0[:,0],X0[:,1],  label='H =%s T'% H0)
plt.plot(X1[:,0],X1[:,1], label='H =%s T'% H1)
plt.plot(X2[:,0],X2[:,1], label='H =%s T'% H2)

plt.subplots_adjust(top=0.98,hspace=0.2500)
plt.legend()

plt.show()




   






