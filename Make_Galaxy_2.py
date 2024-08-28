
import socket
import sys
import time
from pyqtgraph.Qt import QtGui, QtCore
import numpy as np
from scipy import optimize
import pyqtgraph as pg
import pyqtgraph.exporters
import matplotlib as plt
import pylab as pl
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
from PyQt5 import QtGui, QtCore
from PyQt5.QtCore import (QBasicTimer, Qt)
from PyQt5 import QtGui, QtCore
from PyQt5.QtCore import (QBasicTimer, Qt)
from PyQt5.QtGui import *
from PyQt5.QtWidgets import (QMainWindow, QWidget, QProgressBar, QPushButton, QLabel, QApplication, QVBoxLayout,
                             QHBoxLayout, QFileDialog, QSlider)


class Window(QtGui.QMainWindow):
    def __init__(self, parent=None):
        super(Window, self).__init__(parent)
        self.setGeometry(0, 30, 1690, 1000)

        #self. showMaximized()

        self.setWindowIcon(QIcon('xavier.jpg'))
        self.setWindowTitle('MY PLOT')
        self.setStyleSheet("background-color:white;")


        Color_Axis = '#000000'
        self.fontCss = {'font-family': "Times New Roman", 'font-size': '22.0pt', "color": Color_Axis}
        self.font = QFont("Times New Roman",22)

        self.Plot_A()
        self.MakeData()


    def Plot_A(self):

        pg.setConfigOption('foreground', (0,0,0))
        self.graph1 = pg.PlotWidget(self)

        self.graph1.setGeometry(240, 100, 816, 816)
        self.graph1.enableAutoScale()
        self.graph1.setBackground((255,255,255))
        self.graph1.setLabel('bottom', 'H',  'T')
        self.graph1.showAxis('right')
        self.graph1.getAxis('right').setStyle(showValues=False)
        self.graph1.getAxis('top').setStyle(showValues=False)
        self.graph1.showMaximized()
        self.graph1.setLabel('top')
        #self.graph1.setLabel('left', '<font> &chi;</sup> <sup>-1 </sup>;</font>', 'g/emu')
        #self.graph1.setLabel('left', '<font> &chi; <sup>-1</font>', 'g/emu')
        self.graph1.setLabel('left', 'M', 'emu/g')
        self.graph1.getAxis('bottom').setLabel(**self.fontCss)
        self.graph1.getAxis('bottom').tickFont = self.font
        self.graph1.getAxis("bottom").setStyle(tickTextOffset=20)
        self.graph1.getAxis("left").setStyle(tickTextOffset=20)
        #self.graph1.getAxis('right').setLabel(**self.fontCss)
        self.graph1.getAxis('left').setLabel(**self.fontCss)
        self.graph1.getAxis('left').setPen((0,0,0), width = 2)
        self.graph1.getAxis('right').setPen((0, 0, 0), width=2)
        self.graph1.getAxis('bottom').setPen((0, 0, 0), width=2)
        self.graph1.getAxis('top').setPen((0, 0, 0), width=2)
        self.graph1.getAxis('bottom').setHeight(100)
        self.graph1.getAxis('left').setWidth(100)
        #self.graph1.getAxis('bottom').setWidth(200)
        #self.graph1.getAxis('right').tickFont = self.font
        self.graph1.getAxis('left').tickFont = self.font
        self.Curve1 = pg.PlotDataItem(pen=pg.mkPen(color=(255, 0, 0), width=1), symbol='o', symbolBrush=(255,255, 255), symbolPen=(255,0, 0), symbolSize=5, name='red')
        self.Curve2 = pg.PlotDataItem(pen=pg.mkPen(color=(255, 0, 0), width=2), symbol='o', symbolBrush=(255, 255, 255),symbolPen=(0, 0, 0), symbolSize=5, name='red')
        #self.Curve3 = pg.PlotDataItem(pen=pg.mkPen(color=(0, 255, 0), width=0), symbol='o', symbolBrush=(255, 255, 255),symbolPen=(0, 0, 0), symbolSize=0, name='red')
        self.Curve3 = pg.PlotDataItem(pen=None, symbol='o', symbolBrush=(255, 255, 255),symbolPen=(0, 0, 0), symbolSize=2, name='red')
        self.p1 = self.graph1.plotItem
        self.graph1.setRange(xRange = [-4,4], yRange = [-4,4])



    def MakeData(self):
             nf=200
             nf1=200
             
             b =0.997306*1
             bb=0.997306
             aa=1.2851*0.92
             aa=1.2851*0.985
             
             e = np.sqrt(1-np.square(0.997306/(1.2851*1)))
             ee = np.sqrt(1-np.square(bb/aa))
             phi= 90
             phi1=0
             T = np.zeros((nf1, 2))
             T1 = np.zeros((nf1, 2))
             T2 = np.zeros((nf, 2))
             T3 = np.zeros((nf, 2))
             T4 = np.zeros((nf, 2))
             T_Ellipse = np.zeros((nf, 4))
             T2b = np.zeros((nf, 2))
             for k in range(0,nf1):
                 size=((k+0)*0.27*bb+bb*0)/5*1.2+np.exp(-k/100)*1.5
                 phi_k=k*1
                 T[k, 0] = np.cos(phi_k*2*np.pi/180)*size
                 T[k, 1] = np.sin(phi_k*2*np.pi/180)*size
                 T1[k, 0] = np.cos(phi_k*2*np.pi/180)*-size
                 T1[k, 1] = np.sin(phi_k*2*np.pi/180)*-size
             TSpiral=np.vstack((T,T1))    

             for t in range(0,30):    
                 for kk in range(0,nf):
                     phi=(t+10)*15
                     #phi=-97.2466*np.exp((t+1)/(-2.89383)) + -357.496*np.exp((t+1)/(-19.1445))+363.99+90*1
                     
                     #a_value=(t)*0.27*bb+bb
                     a_value=(t)*0.27*bb+bb*1
                     T2[kk, 0] = (a_value)/np.sqrt(1-np.square(ee*np.cos(kk*2*np.pi/180)))*(np.cos(kk*2*np.pi/180)*np.cos(phi*np.pi/180)-np.sin(kk*2*np.pi/180)*np.sin(phi*np.pi/180))
                     T2[kk, 1] = (a_value)/np.sqrt(1-np.square(ee*np.cos(kk*2*np.pi/180)))*(np.cos(kk*2*np.pi/180)*np.sin(phi*np.pi/180)+np.sin(kk*2*np.pi/180)*np.cos(phi*np.pi/180))

                   
                 T3=np.vstack((T3,T2))

             #Make simple Ellipse   
             #pp=0
             #for kk in range(0,nf):
                     #T2[k, 0] = np.cos(k*2*np.pi/180)*(b*t*0.1)/np.sqrt(1-np.square(e*np.cos(k*2*np.pi/180+(phi+t*-9.5+np.cos(t*100*np.pi/180)*0)*np.pi/180)))
                     #T2[k, 1] = np.sin(k*2*np.pi/180)*(b*t*0.1)/np.sqrt(1-np.square(e*np.cos(k*2*np.pi/180+(phi+t*-9.5+np.cos(t*100*np.pi/180)*0)*np.pi/180)))
                     #T4[kk, 0] = np.cos(kk*2*np.pi/180)*(bb*pp*0.5)/np.sqrt(1-np.square(ee*np.cos(kk*2*np.pi/180+(phi1+pp*-9.5*0+-1*pp)*np.pi/180)))
                     #T4[kk, 1] = np.sin(kk*2*np.pi/180)*(bb*pp*0.5)/np.sqrt(1-np.square(ee*np.cos(kk*2*np.pi/180+(phi1+pp*-9.5*0+-1*pp)*np.pi/180)))
              #      T4[kk, 0] = np.cos(kk*2*np.pi/180)*(bb*pp*0.5)/np.sqrt(1-np.square(ee*np.cos(kk*2*np.pi/180+(phi1+pp*-9.5*0+-2*pp)*np.pi/180)))
              #      T4[kk, 1] = np.sin(kk*2*np.pi/180)*(bb*pp*0.5)/np.sqrt(1-np.square(ee*np.cos(kk*2*np.pi/180+(phi1+pp*-9.5*0+-2*pp)*np.pi/180)))
                     
              #   T2b=np.vstack((T2b,T4))

             #SimpleEllipse with rotation  
             for kk in range(0,nf):
                     phi=+90-6.76519*1
                     nn=1
                     T_Ellipse[kk, 0] = (bb*(1+0.27*nn))/np.sqrt(1-np.square(ee*np.cos(kk*2*np.pi/180)))*(np.cos(kk*2*np.pi/180)*np.cos(phi*np.pi/180)-np.sin(kk*2*np.pi/180)*np.sin(phi*np.pi/180))
                     T_Ellipse[kk, 1] = (bb*(1+0.27*nn))/np.sqrt(1-np.square(ee*np.cos(kk*2*np.pi/180)))*(np.cos(kk*2*np.pi/180)*np.sin(phi*np.pi/180)+np.sin(kk*2*np.pi/180)*np.cos(phi*np.pi/180))
                     
             

    



             self.Curve1.setData(x=TSpiral[:, 0], y=(TSpiral[:, 1]))
             self.Curve2.setData(x=T_Ellipse[:, 0], y=(T_Ellipse[:, 1]))
             self.Curve3.setData(x=T3[:, 0], y=(T3[:, 1]))

             #self.Curve2.setData(x=xf, y=yf)
             #self.graph1.plotItem.addItem(self.Curve1)
             #self.graph1.plotItem.addItem(self.Curve2)
             self.graph1.plotItem.addItem(self.Curve3)
             
             

             #np.savetxt("C:/Xavier/Xavier 2019/ZnCoO Bruno/Sample_2%/Sample 2/Exp5K.txt", Exp, delimiter=' ')
             #np.savetxt("C:/Xavier/Xavier 2019/ZnCoO Bruno/Sample_2%/Sample 2/Fit5K.txt", Fit, delimiter=' ')
             #np.savetxt("C:/Xavier/Xavier 2019/ZnCoO_Adolf/1%/Exp4K.txt", Exp, delimiter=' ')
             #np.savetxt("C:/Xavier/Xavier 2019/ZnCoO_Adolf/1%/Fit4K.txt", Fit, delimiter=' ')
             print(TSpiral)
























def main():
    app = QtGui.QApplication(sys.argv)
    # app.setStyleSheet(qdarkstyle.load_stylesheet_pyqt5())
    # app.setStyleSheet(open("stylePySide.qss", 'r').read())
    # app.setStyleSheet(open("aqua.qss", 'r').read())
    #app.setStyleSheet(open("qdarkstyle\darkstyle.qss", 'r').read())
    #app.setStyleSheet(open("orange.qss", 'r').read())
    # app.setApplicationName('Simple Plot')
    window = Window()
    #window.showMaximized()
    window.show()
    app.exec_()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
