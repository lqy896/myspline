#coding=utf-8
import matplotlib.pyplot as plt
import numpy as np
import os
import math

plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False
#使用 ax2+bx+c=y 的方式
#边界二阶导数为0 a1=0
#只需要8个未知数 8个方程组即可
#b1 c1 a2 b2 c2 a3 b3 c3
# x0，x3两个端点都有一个二次函数经过，可确定2个方程
# x1，x2两个中间点都有两个二次函数经过，可确定4个方程
# 中间点处必须连续，需要保证左右二次函数一阶导相等

def getData(r):
    data=[0 for v in range(r)]
    return data

def showPic2(x,y):
    parame = []
    yparam=[]
    n=len(x)
    w=3
    r=(n-1)*w
    i=1
    #i 与 i-1 组成方程
    while i<n:
        data=getData(r)
        data[(i-1)*w+0]=x[i]*x[i]
        data[(i-1)*w+1]=x[i]
        data[(i-1)*w+2]=1
        temp=data[1:]
        parame.append(temp)
        yparam.append(y[i])
        i=i+1
    i=0
    #i 与 i+1 组成方程
    while i<n-1:
        data=getData(r)
        data[i*w]=x[i]*x[i]
        data[i*w+1]=x[i]
        data[i*w+2]=1
        temp=data[1:]
        parame.append(temp)
        yparam.append(y[i])
        i=i+1
    #中间两两系数1阶段导数相同
    i=1
    while i<n-1:
        data=getData(r)
        data[(i-1)*w+0]=2*x[i]
        data[(i-1)*w+1]=1
        data[(i-1)*w+3]=-2*x[i]
        data[(i-1)*w+4]=-1
        temp=data[1:]
        parame.append(temp)
        yparam.append(0)
        i=i+1

    print(len(parame[0]))
    print(len(yparam))
    a = np.array(parame)
    b = np.array(yparam)

    res=np.linalg.solve(a,b)
    nres=[0]
    for r in res:
        nres.append(r)
    #结果中a1 默认为0
    #在计算差值点并绘图 差值区间工有w个
    i=0
    while i<n-1:
        if(x[i+1]>x[i]):
            xt,yt=getXY(x[i],x[i+1],0.01,nres[i*w],nres[i*w+1],nres[i*w+2])
        else:
            xt,yt=getXY(x[i],x[i+1],-0.01,nres[i*w],nres[i*w+1],nres[i*w+2])
        i=i+1
        plt.plot(xt,yt,'r')
    #print(a)
    #print(b)
    #print(res)
    #print(nres)
    i=0
    while i<n:
        plt.plot(x[i],y[i],'bo')
        i=i+1
    plt.title(u"任意点二次样条曲线差值算法")


def getXY(x0,x1,d,a,b,c):
    x=np.arange(x0,x1,d)
    y=a*x*x+b*x+c
    return x,y


def main():
    x = [3, 4.5, 7, 9,10, 15, 17, 20,28,33,40,56]
    y = [2.5, 4.2, 2, 1,3.6, 8.2, 2, 5,6,1,9,22]
    showPic2(x,y)
    plt.show()

if __name__ == '__main__':
    main()
