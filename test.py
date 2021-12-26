#coding=utf-8
import matplotlib.pyplot as plt
import numpy as np
import os
import math
import spline3 as spline3

def getXY3(x,a,b,c,d):
    y=a*x**3+b*x**2+c*x+d
    return y

def getIndex(s,arr):
    l=len(arr)
    for i in range(l):
        if(arr[i]>s):
            return i
    return None

#输入路径点，求得图像长度的方程
def testSD(x,y):
    s=[0]
    l=len(x)
    for i in range(l-1):
        s.append(math.sqrt((y[i+1]-y[i])**2+(x[i+1]-x[i])**2))
    #print(s)
    s=np.cumsum(s)
    #分别用3次多项式合成曲线，x y 都是相对于路长的函数
    #给出S值后，要先找出 属于那个区间的
    resx=spline3.mySolve3(s,x,False)
    resy=spline3.mySolve3(s,y,False)
    spline3.mySolve3(x,y,True)
    print(s)
    #把D距离等分
    i=1
    w=4
    """
    while i<l:
        for t in np.arange(s[i-1],s[i],0.1):
            x0=getXY3(t,resx[(i-1)*w+0],resx[(i-1)*w+1],resx[(i-1)*w+2],resx[(i-1)*w+3])
            y0=getXY3(t,resy[(i-1)*w+0],resy[(i-1)*w+1],resy[(i-1)*w+2],resy[(i-1)*w+3])
            plt.plot(x0,y0,"yo")
        i=i+1
    """
    d=2
    i=getIndex(d,s)
    print(i,"-----------")
    x0=getXY3(d,resx[(i-1)*w+0],resx[(i-1)*w+1],resx[(i-1)*w+2],resx[(i-1)*w+3])
    y0=getXY3(d,resy[(i-1)*w+0],resy[(i-1)*w+1],resy[(i-1)*w+2],resy[(i-1)*w+3])
    plt.plot(x0,y0,"ro")
    """
    return x0,y0

    i=0
    while i<l:
        plt.plot(x[i],y[i],"ro")
        i=i+1
    """   
    plt.show()


def main():
    #x=[1,2,3,4,5,8]
    #y=[2,4,5,8,9,10]
    x=[3, 4.5, 7, 9,10]
    y=[2.5, 4.2, 2, 1,3.6]

    testSD(x,y)

if __name__ == '__main__':
    main()
