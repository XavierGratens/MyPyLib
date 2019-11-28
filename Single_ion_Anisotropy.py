import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as LA


"Constants"

"Bohr magneton in cgs"
muB_cgs=9.2741E-21
"Planck Constant in SI"
hp=6.6262E-34

"Boltzman´s constant in cgs"
kB_cgs=1.3807E-16
" Avogadro´s constant "
Na = 6.022E23 

"Spin for Co2+"
Spin=3/2

"Molar mass for Zn(1-x)CoxO"
M_Zn =65.39
M_Co = 58.9332
M_O = 15.9994

"Making matrix"
i=np.arange(Spin*2+1)
Sz = Spin-i
a0=np.arange(int(Spin*2+1))
mI=np.identity(int(2*Spin+1))
a0=np.zeros(shape=(int(Spin*2+1)))

"Matrix Sz"
mSz=np.diag(Sz)

Sp=np.sqrt(Spin*(Spin+1)-(Spin-i)*(Spin-i+1))
mSp=np.diag(Sp)
mSp=np.delete(mSp,0,0)
mSp=np.vstack((mSp,a0))
mSm=mSp.transpose()

"Matrix Sx" 
mSx=(mSp+mSm)/2
"Matrix Sy"
mSy=(mSp-mSm)/2*-1j

"Landé Factors"
g_para = 2.236
g_perp = 2.277

"Density of ZnO in g/cm^3"
d_ZnO = 5.61

"Axial term D in kelvin"
D =3.971

"Crystal Field Matrices"
Axial_M = mSz*mSz




"Zeeman energy"
def Zeeman(g_para,g_perp,theta,phi):
    ZZ=g_perp*(np.sin(np.pi*theta/180)*np.sin(np.pi*phi/180)*mSx+np.sin(np.pi*theta/180)*np.cos(np.pi*phi/180)*mSy)+g_para*np.cos(np.pi*theta/180)*mSz
    return ZZ



"Energy Levels in Kelvin, Field in kOe"
def EvsH_cgs(g_para,g_perp,theta,phi,D,Field,NP):
   v0=[0,0,0,0,0]
   for i in range(NP+2):
       v=LA.eigvals((muB_cgs/kB_cgs) * Zeeman(g_para,g_perp,theta,phi) * i /NP*Field*1000 + D*Axial_M)
       v=v.real
       v=np.sort(v)
       R=np.hstack((i /NP*Field,v))
       v0=np.vstack((v0,R))
   return v0[1:]    




"Magnetization calculated from energy levels in emu/g, x is the Co concentration and Temp is the temperature in Kelvin"        
def MvsH_cgs(E,x,Temp):
    Zero=[0,0]
    Em=E.transpose()
    Em1=Em[1:]
    Em2=Em1.transpose()
    H= Em[0].transpose()
    S = np.zeros(shape=(Em2.shape[0],1))
    mag = np.zeros(shape=(S.shape[0]-1,2))
    for i in range(Em2.shape[0]):
        S[i] = np.sum(np.exp(-Em2[i]/Temp))
    for i in range(S.shape[0]-1):
        mag[i,1] = (np.log(S[i+1])-np.log(S[i]))/((H[1]-H[0])*1000)*Temp*Na*kB_cgs*(x/(M_Zn*(1-x)+M_Co*x+M_O))
        mag[i,0] = H[i] +(H[1]-H[0])/2
    return np.vstack((Zero,mag))    


"Calculation of the magnetization in emu/g, Field is the maximum value of the magnetic field in kOe, NPF is the number of points, x is the Co concentration" 
"Temp is the temperature in Kelvin and Angle is the angle between the magnetic field and the c axis of the wurtzite crystal"

def MvsH_Bulk(Field,NPF,x, Temp, Angle):
    MF=np.zeros(shape=(NPF+2,2))
    MFF=np.zeros(shape=(NPF+1,2))
    Result_a = EvsH_cgs(g_para,g_perp,Angle,0,D,Field,NPF)
    M = MvsH_cgs(Result_a,x,Temp)
    "Interpolation"
    x_i=np.arange(NPF+1)
    x_i =x_i/NPF*Field
    y_i=np.interp(x_i,M[:,0],M[:,1])
    MFF[:,0] = x_i
    MFF[:,1] = y_i
    return MFF





fig1=plt.figure(1)
plt.xlabel('H (kOe)')
plt.ylabel('Energy (K)')

"Calculation of the energy levels as function of the magnetic field for H || c"
Result=EvsH_cgs(g_para,g_perp,0,0,D,120,1000)

"Plot of the energy levels as function of the field"
plt.plot(Result[:,0],Result[:,1])
plt.plot(Result[:,0],Result[:,2])
plt.plot(Result[:,0],Result[:,3])
plt.plot(Result[:,0],Result[:,4])

plt.text(4, 26, '$D$ = %s K' % D , fontsize=10)
plt.text(4, 22, r'$\theta$ = 0', fontsize=10)




fig2=plt.figure(2)

plt.xlabel('H (kOe)')
plt.ylabel('Magnetization (emu/g)')

"Imput Parameters"

"Temperature"
T=0.5
"Cobalt concentration"
c = 0.01
"Calculation of the magnetization for theta = 0"
Mbulk=MvsH_Bulk(120,1000,c, T, 0)
"Calculation of the magnetization for theta = 90"
Mbulk1=MvsH_Bulk(120,1000,c, T, 90)

plt.text(80,(np.max(Mbulk1[:,1]))/2, "$Zn_{1-x}Co_{x}O$",fontsize=10)
plt.text(80, (np.max(Mbulk1[:,1]))/3, " x = %s " % c , fontsize=10)
plt.text(80, (np.max(Mbulk1[:,1]))/4, "$T$ = %s K" % T , fontsize=10)

"Plot of the magnetization as function of the magnetic field"
plt.plot(Mbulk[:,0],Mbulk[:,1],label='H || c')
plt.plot(Mbulk1[:,0],Mbulk1[:,1], label='H $\perp$ c')
plt.legend()
plt.show()

