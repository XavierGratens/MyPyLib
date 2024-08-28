import numpy as np
import matplotlib.pyplot as plt


eps=(-1+np.power(-3+0j,1/2))/2

eps2=np.power((-1+np.power(-3+0j,1/2))/2,2)


def Solution(a,b,c,d):
    Delta0=np.power(b,2)-3*a*c
    Delta1=2*np.power(b,3)-9*a*b*c+27*np.power(a,2)*d
    Cp=np.power((Delta1+np.power(np.power(Delta1,2)-4*np.power(Delta0,3)+0j,1/2))*1/2,1/3)
    Cm=np.power((Delta1-np.power(np.power(Delta1,2)-4*np.power(Delta0,3),1/2))*1/2,1/3)+0j
    x1=-1/(3*a)*(b+Cp+Delta0/Cp)
    x2=-1/(3*a)*(b+eps*Cp+Delta0/(eps*Cp))
    x3=-1/(3*a)*(b+eps2*Cp+Delta0/(eps2*Cp))
    x1=-1/(3*a)*(b+Cm+Delta0/Cm)
    x2=-1/(3*a)*(b+eps*Cm+Delta0/(eps*Cm))
    x3=-1/(3*a)*(b+eps2*Cm+Delta0/(eps2*Cm))
    return x1,x2,x3






#print(np.power(-3+0j,1/2))

#R1=Solution(2,3,-11,-6)

#print(R1)
def Solution2(a,b,c,d):
    Delta0=np.power(b,2)-3*a*c
    Delta1=2*np.power(b,3)-9*a*b*c+27*np.power(a,2)*d
    Cp=np.power((Delta1+np.power(np.power(Delta1,2)-4*np.power(Delta0,3)+0j,1/2))*1/2,1/3)
    Cm=np.power((Delta1-np.power(np.power(Delta1,2)-4*np.power(Delta0,3)+0j,1/2))*1/2,1/3)+0j
    x1=-b/(3*a)-1/(3*a)*Cp-1/(3*a)*Cm
    x2=-b/(3*a)+ (1+np.power(3,1/2)*1j)/(6*a)*Cp +(1-np.power(3,1/2)*1j)/(6*a)*Cm
    x3=-b/(3*a)+ (1-np.power(3,1/2)*1j)/(6*a)*Cp +(1+np.power(3,1/2)*1j)/(6*a)*Cm
    return x1,x2,x3

R2=Solution2(1,-6,11,-6)
print(R2)
