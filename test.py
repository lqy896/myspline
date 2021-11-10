#coding=utf-8
import matplotlib.pyplot as plt
import numpy as np
import os
import math



def main():
    
    x=np.arange(-1,5,0.1)
    y=2+2*x+5*x*x-x*x*x
    plt.plot(x,y,"r")
    x0=2
    y0=2+2*x0+5*x0*x0-x0*x0*x0
    plt.plot(x0,y0,"bo")
    k=2+10*x0-3*x0*x0
    x1=np.arange(x0,5,0.1)
    y1=(x1-x0)*k+y0
    plt.plot(x1,y1,"g")
 
    k1=-1/k
    x0=2
    y0=2+2*x0+5*x0*x0-x0*x0*x0
    x2=np.arange(x0,-5,-0.1)
    y2=(x2-x0)*k1+y0
    #plt.xlim(0, 50)
    #plt.xlim(0, 50)
    plt.plot(x2,y2,"y")
    #plt.axis('scaled')
    plt.axis('equal')


    plt.show()

if __name__ == '__main__':
    main()
