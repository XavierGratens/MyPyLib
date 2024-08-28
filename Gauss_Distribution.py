import numpy as np
import matplotlib.pyplot as plt


def gaussian(mu, sig, NP):
    Result = np.zeros(shape=(NP,2))
    for i in range(NP):
        x = mu-(NP/2-i)*0.01
        Result[i,0] = x
        Result[i,1] = np.exp((-np.power((x - mu)/sig, 2.))*1/2)*1/(sig*(np.power(2*np.pi,1/2)))
    return Result    

G=gaussian(10,0.1,100)



A=np.sum(G[:,1])*1/100
print(A)
fig1=plt.figure(1)
plt.plot(G[:,0],G[:,1])
plt.show()
