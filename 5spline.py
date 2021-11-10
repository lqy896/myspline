#coding=utf-8
import matplotlib.pyplot as plt
import numpy as np
import os
import math

"""
3次样条实现
本程序由 李勤宇 在参照了网络大神 https://www.cnblogs.com/yanshw/p/11194058.html
的基础上，优化了算法过程，实现了高性能样条曲线绘制算法，仅供参考，转发请注明出处
本人联系方式 15250035289 微信同号
"""

plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False
#使用 ax3+bx2+cx+d=y 的方式
#边界二阶导数为0 6a=0
#n+1 个点 需要4*n个方程
#a1,b1 c1 d1 a2 b2 c2 d2 a3 b3 c3 d3 a4 b4 c4 d4
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
    w=4
    r=(n-1)*w
    i=1
    #i 与 i-1 组成方程
    while i<n:
        data=getData(r)
        data[(i-1)*w+0]=x[i]**3
        data[(i-1)*w+1]=x[i]**2
        data[(i-1)*w+2]=x[i]
        data[(i-1)*w+3]=1
        parame.append(data)
        yparam.append(y[i])
        i=i+1
    i=0
    #i 与 i+1 组成方程
    while i<n-1:
        data=getData(r)
        data[i*w]=x[i]**3
        data[i*w+1]=x[i]**2
        data[i*w+2]=x[i]
        data[i*w+3]=1
        parame.append(data)
        yparam.append(y[i])
        i=i+1
    #中间两两系数1阶段导数相同
    i=1
    while i<n-1:
        data=getData(r)
        data[(i-1)*w+0]=3*x[i]**2
        data[(i-1)*w+1]=2*x[i]
        data[(i-1)*w+2]=1
        data[(i-1)*w+4]=-3*x[i]**2
        data[(i-1)*w+5]=-2*x[i]
        data[(i-1)*w+6]=-1
        parame.append(data)
        yparam.append(0)
        i=i+1

    #中间两两系数2阶段导数相同
    i=1
    while i<n-1:
        data=getData(r)
        data[(i-1)*w+0]=6*x[i]
        data[(i-1)*w+1]=2
        data[(i-1)*w+4]=-6*x[i]
        data[(i-1)*w+5]=-2
        parame.append(data)
        yparam.append(0)
        i=i+1

    #两边2阶导数为0
    data=getData(r)
    data[0]=6
    parame.append(data)
    yparam.append(0)

    data=getData(r)
    data[r-4]=6
    parame.append(data)
    yparam.append(0)

    print(len(parame[0]))
    print(len(yparam))
    a = np.array(parame)
    b = np.array(yparam)

    res=np.linalg.solve(a,b)
    #print(res)
  
    #结果中a1 默认为0
    #在计算差值点并绘图 差值区间工有w个
    i=0
    while i<n-1:
        if(x[i+1]>x[i]):
            xt,yt=getXY(x[i],x[i+1],0.01,res[i*w],res[i*w+1],res[i*w+2],res[i*w+3])
        else:
            xt,yt=getXY(x[i],x[i+1],-0.01,res[i*w],res[i*w+1],res[i*w+2],res[i*w+3])
        i=i+1
        plt.plot(xt,yt,'r')
    print(a)
    print(b)
    print(res)
    print(res)
    i=0
    while i<n:
        plt.plot(x[i],y[i],'bo')
        i=i+1
    plt.title(u"任意点三次样条曲线差值算法")


def getXY(x0,x1,l,a,b,c,d):
    x=np.arange(x0,x1,l)
    y=a*x**3+b*x**2+c*x+d
    return x,y

def main():
    """
    x = [3, 4.5, 7, 9]
    y = [2.5, 4.2, 2, 1]
    showPic2(x,y)
    x = [8, 15, 17, 20]
    y = [3.6, 8.2, 2, 5]
    showPic2(x,y)
    """
    x = [3, 4.5, 7, 9,10, 15, 17, 20,28,33,40,56]
    y = [2.5, 4.2, 2, 1,3.6, 8.2, 2, 5,6,1,9,22]
    #i=0
    #while
    showPic2(x,y)
    plt.show()

if __name__ == '__main__':
    main()
