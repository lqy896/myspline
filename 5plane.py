#coding=utf-8
import matplotlib.pyplot as plt
import numpy as np
import os
import math
"""
状态拟合方程
本程序由 李勤宇 实现，仅供参考，转发请注明出处
本人联系方式 15250035289 微信同号
"""

plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False
#使用 ax2+bx+c=y 的方式
#边界二阶导数为0 a1=0
#n+1 个点 需要3*n个方程
#b1 c1 a2 b2 c2 a3 b3 c3
# x0，x3两个端点都有一个二次函数经过，可确定2个方程
# x1，x2两个中间点都有两个二次函数经过，可确定4个方程
# 中间点处必须连续，需要保证左右二次函数一阶导相等

# 损失函数权重
KJ = 0.1
KT = 0.1
KD = 1.0
KLAT = 1.0
PathList=[]

#横向损失
def TransverselossD(t,a0,a1,a2,a3,a4,a5):
    d=a0+a1*t+a2*t**2+a3*t**3+a4*t**4+a5*t**5
    v=a1+2*a2*t+3*a3*t**2+4*a4*t**3+5*a5*t**4
    a=2*a2+6*a3*t+12*a4*t**2+20*a5*t**3
    aa=6*a3+24*a4*t+60*a5*t**2
    return KJ*d
    #return aa*KJ

#纵向损失
def TransverselossS(t,a0,a1,a2,a3,a4,a5,vc):
    d=a0+a1*t+a2*t**2+a3*t**3+a4*t**4+a5*t**5
    v=a1+2*a2*t+3*a3*t**2+4*a4*t**3+5*a5*t**4
    a=2*a2+6*a3*t+12*a4*t**2+20*a5*t**3
    aa=6*a3+24*a4*t+60*a5*t**2
    #return d*KD+KT*t+(v-vc)*(v-vc)*
    return KJ*d

#计算横向损失
def getDLoss(dt,a0,a1,a2,a3,a4,a5,d1):
    bc=0.01
    t=np.arange(0,dt+bc,bc)
    #y=a0+a1*t+a2*t**2+a3*t**3+a4*t**4+a5*t**5
    loss=0
    for i in t:
        loss=loss+TransverselossD(i,a0,a1,a2,a3,a4,a5)
    loss=loss+KT*dt+KD*d1**2
    return loss

#计算纵向损失
def getSLoss(dt,a0,a1,a2,a3,a4,a5,v,vc):
    bc=0.01
    t=np.arange(0,dt+bc,bc)
    y=a0+a1*t+a2*t**2+a3*t**3+a4*t**4+a5*t**5
    loss=0
    for i in t:
        loss=loss+TransverselossS(i,a0,a1,a2,a3,a4,a5,v)
    loss=loss+KT*dt+KD*(v-vc)**2
    #print(v,vc)
    return loss

def  getXY(dt,a0,a1,a2,a3,a4,a5):
    bc=0.01
    t=np.arange(0,dt+bc,bc)
    y=a0+a1*t+a2*t**2+a3*t**3+a4*t**4+a5*t**5
    return t,y


#拟合5次多项式
def plane5(x0,x1,c,dt,v,l,vc):
    a0=x0[0]
    a1=x0[1]
    a2=x0[2]/2
    parame=[[dt**3,dt**4,dt**5],[3*dt**2,4*dt**3,5*dt**4],[6*dt,12*dt**2,20*dt**3]]
    yparam=[x1[0]-a0-a1*dt-a2*dt**2,x1[1]-a1-2*a2*dt,x1[2]-2*a2]
    a = np.array(parame)
    b = np.array(yparam)
    res=np.linalg.solve(a,b)
    x,y=getXY(dt,a0,a1,a2,res[0],res[1],res[2])
    #print("系数0-5",a0,a1,a2,res[0],res[1],res[2])
    plt.plot(x,y,c)
    #返回损失
    if l==0:
        return getDLoss(dt,a0,a1,a2,res[0],res[1],res[2],x1[0])
    else:
        return getSLoss(dt,a0,a1,a2,res[0],res[1],res[2],v,vc)

    

class  Fpath:
    def __init__(self,t,d,s,v,loss):
        self.t=t
        self.s=s
        self.d=d
        self.v=v
        self.loss=loss
    

def  loopGet():
    #初始状态 最终状态 求轨迹方程 用5次多项式拟合

    xd0=[2,0,0]
    xs0=[1,0,0]

    dt=0.5
    dd=0.5
    ds=0.5

    d=xd0[0]
    s=xs0[0]
    
    t=0
    maxt=2
    maxv=3
    dv=0.5
    a=0
    maxd=1
    mind=-1
    maxs=3
    mins=1.5

    plt.title("5次轨迹规划")

    #计算相同时间下 横向、纵向可选轨迹，加入列表，最后求得相同时间下，损失最小的
    while t<maxt:
        #先遍历横向 从-1 到 1之间遍历
        t=t+dt
        d=maxd+dd
        while d>mind:
            d=d-dd
            loss=plane5(xd0,[d,0,0],'r',t,0,0,0)#0代表横向
            # 再遍历纵向 s的值无所谓，主要关心跟车速度 加速度
            s=mins
            v=dv
            while v<maxv:
                v=v+dv
                loss=loss+plane5(xs0,[s,v,0],'b',t,v,1,xs0[1])
                #两次loss 加
                fpath=Fpath(t,d,s,v,loss)
                PathList.append(fpath)   
    plt.title("5次轨迹规划")

    tmp=999999
    for ps in PathList:
        if ps.loss < tmp:
            tmp=ps.loss
            min=ps
    print("最小",min.loss)
    print("最佳时间",min.t)
    print("最佳纵向",min.s)
    print("最佳横向",min.d)
    print("最佳纵向速度",min.v)

    plane5(xd0,[min.d,0,0],'g',min.t,0,0,0)
    plane5(xs0,[min.s,v,0],'y',min.t,v,1,maxv)
 
    plt.show()

def main():
    loopGet()
   # x1=[3,2,1]
    #x2=[3,1,1]
    #x3=[3,2,1]
    #xt1=[4,3,0]
    #xt2=[6,3,0]
    #xt3=[8,3,0]
    #(x1,xt1,'r',1)
    #plane5(x2,xt2,'b',1)
    #plane5(x3,xt3,'g',1)
    plt.show()
    

if __name__ == '__main__':
    main()
