from scipy import misc
import numpy as np
import matplotlib.pyplot as plt

f = misc.imread('Fig1m.png')
plt.imshow(f)
plt.text(300, -50, 'Close the window to digitize', fontsize=10)
plt.show()




def Extract_A(Img):
    redImg=Img[:,:,0]
    Result=np.array([[0,0]])
    for i in range(redImg.shape[0]):
        print('Running Extract')
        for j in range (redImg.shape[1]):
            if redImg[i, j] == 0 :
                Result = np.append(Result,[[j, -i]],axis = 0)
    return Result

def Filterin(M,x1,x2,y1,y2):
    Result=np.array([[0,0]])
    for i in range(M.shape[0]):
        if all ([M[i,0] > x1, M[i,0] < x2, M[i,1] > y1, M[i,1] < y2]):
                Result = np.append(Result,[M[i]],axis = 0)
                print('Running Filter in')
    return Result[1:]

def Filterout(M,x1,x2,y1,y2):
    Result=np.array([[0,0]])
    for i in range(M.shape[0]):
        if any ([M[i,0] < x1, M[i,0] > x2, M[i,1] < y1, M[i,1] > y2]):
                Result = np.append(Result,[M[i]],axis = 0)
                print('Running Filter out')
    return Result[1:]

"Digitize the overall figure"
R=Extract_A(f)



fig2=plt.figure(2)
R1=np.zeros(shape=(R.shape[0],2))
R1[:,0] = (R[:,0]-131)/867*6.6
R1[:,1] = (R[:,1]+690)/640 
plt.scatter(R1[:,0],R1[:,1], s = 0.75)
plt.text((np.max(R1[:,0]-np.min(R1[:,0])))/2+np.min(R1[:,0]), 1.2, 'Close the window to delete bad points', horizontalalignment='center', fontsize=10)
plt.show()



"Filter data"
fig3=plt.figure(3)
R2=Filterin(R1,0.01,6.55,0.05,0.95)
R2=Filterout(R2,4.1,5.6,0.1,0.2)
R2=Filterout(R2,6.4,6.7,0.16,0.82)
R2=Filterout(R2,-0.1,0.1,0.16,0.82)

"Sort data"
R2T=R2.transpose()
indice=np.argsort(R2T)[:1]
R2[:]=R2[indice]

"Interpol the data"
maxCol = np.amax(R2[:,0])
minCol = np.amin(R2[:,0])
NP=100
x_i=np.arange(NP)
Result=np.zeros(shape=(NP,2))
x_i =minCol +(maxCol-minCol)/NP*x_i
y_i= np.interp(x_i,R2[:,0],R2[:,1])
Result[:,0]=x_i
Result[:,1]=y_i

plt.scatter(Result[:,0],Result[:,1], s = 0.75)
plt.show()

    
