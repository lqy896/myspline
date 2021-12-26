#coding=utf-8
import matplotlib.pyplot as plt
import numpy as np
import os
import math

from numpy.lib.arraysetops import setdiff1d
import spline3 as spline3

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

#全局路径点
xp=[3, 4.5, 7, 9,10]
yp=[2.5, 4.2, 2, 1,3.6]

xishu=[]
w=4

#计算Y值 3 次函数的路点
def getY(x,r,i):
        yp=r[i*w+0]*x**3+r[i*w+1]*x**2+r[i*w+2]*x+r[i*w+3]
        return yp

def getfrent(t,res):
        s=res[0]+res[1]*t+res[2]*t**2+res[3]*t**3+res[4]*t**4+res[5]*t**5
        v=res[1]+2*res[2]*t+3*res[3]*t**2+4*res[4]*t**3+5*res[5]*t**4
        a=2*res[2]+6*res[3]*t+12*res[4]*t**2+20*res[5]*t**3
        return s,v,a

def getXY3(x,a,b,c,d):
    y=a*x**3+b*x**2+c*x+d
    return y

#输入路径点，求得图像长度的方程
def getSDXY(x,y):
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
    #spline3.mySolve3(x,y)
    #x0=getXY3(t,resx[(i-1)*w+0],resx[(i-1)*w+1],resx[(i-1)*w+2],resx[(i-1)*w+3])
    #y0=getXY3(t,resy[(i-1)*w+0],resy[(i-1)*w+1],resy[(i-1)*w+2],resy[(i-1)*w+3])
    return s,resx,resy



#横向损失
def getDDloss(t,d,dt,a0,a1,a2,a3,a4,a5):
    #d=a0+a1*t+a2*t**2+a3*t**3+a4*t**4+a5*t**5
    #v=a1+2*a2*t+3*a3*t**2+4*a4*t**3+5*a5*t**4
    #a=2*a2+6*a3*t+12*a4*t**2+20*a5*t**3
    #aa=6*a3+24*a4*t+60*a5*t**2
    d_ddd = [6*a3+24*a4*t1+60*a5*t1**2 for t1 in np.arange(0,t,dt)]
    Jp = sum(np.power(d_ddd, 2))  # square of jerk
    return KJ * Jp + KT * t + KD * d**2

#纵向损失
def getSSloss(t,v,vc,dt,a0,a1,a2,a3,a4):
    d=a0+a1*t+a2*t**2+a3*t**3+a4*t**4
    v=a1+2*a2*t+3*a3*t**2+4*a4*t**3
    a=2*a2+6*a3*t+12*a4*t**2
    aa=6*a3+24*a4*t
    s_ddd= [6*a3+24*a4*t1 for t1 in np.arange(0,t,dt)]
    Js = sum(np.power(s_ddd, 2))  # square of jerk
    ds=(vc-v)**2
    return KJ * Js + KT * t + KD * ds

#loss=loss+KT*dt+KD*d1**2
#loss=loss+KT*dt+KD*(v-vc)**2

#5 次函数
def  getXY5(dt,a0,a1,a2,a3,a4,a5,color):
    bc=0.2
    t=np.arange(0,dt+bc,bc)
    d=a0+a1*t+a2*t**2+a3*t**3+a4*t**4+a5*t**5
    d1=a0+a1*dt+a2*dt**2+a3*dt**3+a4*dt**4+a5*dt**5
    dv1=a1+2*a2*dt+3*a3*dt**2+4*a4*dt**3+5*a5*dt**4
    da1=2*a2+6*a3*dt+12*a4*dt**2+20*a5*dt**3
     #绘制横向曲线
    #plt.plot(t,y,color)
    return d,d1,dv1,da1

#4 次函数
def  getXY4(dt,a0,a1,a2,a3,a4,color):
    bc=0.2
    t=np.arange(0,dt+bc,bc)
    s=a0+a1*t+a2*t**2+a3*t**3+a4*t**4
    s1=a0+a1*dt+a2*dt**2+a3*dt**3+a4*dt**4
    sd1=a1+2*a2*dt+3*a3*dt**2+4*a4*dt**3
    sa1=2*a2+6*a3*dt+12*a4*dt**2
     #绘制纵向曲线
    #plt.plot(t,y,color)
    return s,s1,sd1,sa1

#输入最好状态，计算系数 d s 系数
def plane5_best3(x0,x1,dt):
    a0=x0[0]
    a1=x0[1]
    a2=x0[2]/2
    parame=[[dt**3,dt**4,dt**5],[3*dt**2,4*dt**3,5*dt**4],[6*dt,12*dt**2,20*dt**3]]
    yparam=[x1[0]-a0-a1*dt-a2*dt**2,x1[1]-a1-2*a2*dt,x1[2]-2*a2]
    a = np.array(parame)
    b = np.array(yparam)
    res=np.linalg.solve(a,b)
    #print(res)
    return [a0,a1,a2,res[0],res[1],res[2]]


#x1 时刻 x[0] 不作为控制量时 用这个求解，求解的是4次方程
def plane5_best2(x0,x1,dt): 
    a0=x0[0]
    a1=x0[1]
    a2=x0[2]/2
    parame=[[3*dt**2,4*dt**3],[6*dt,12*dt**2]]
    yparam=[x1[1]-a1-2*a2*dt,x1[2]-2*a2]
    a = np.array(parame)
    b = np.array(yparam)
    res=np.linalg.solve(a,b)
    #print(res)
    return [a0,a1,a2,res[0],res[1]]


    

class  Fpath:
    def __init__(self,t,d,v,loss,ddr,ssr):
        self.t=t
        self.d=d
        self.v=v
        self.loss=loss
        self.ssr=ssr
        self.ddr=ddr
    

def loopGet(xd0,xs0):
    #初始状态 最终状态 求轨迹方程 用5次多项式拟合
    dt=0.2
    dd=0.2
    t=0
    maxt=2
    maxv=3
    dv=0.2
    maxd=1
    mind=-1
    vc=xs0[1]
    PathList=[]
    #计算相同时间下 横向、纵向可选轨迹，加入列表，最后求得相同时间下，损失最小的
    while t<maxt:
        #先遍历横向
        t=t+dt
        d=maxd
        while d>=mind:
            d=d-dd
            r5=plane5_best3(xd0,[d,0,0],t)
            #绘制图形
            getXY5(t,r5[0],r5[1],r5[2],r5[3],r5[4],r5[5],"r")
            #计算横向损失
            loss0=getDDloss(t,d,dt,r5[0],r5[1],r5[2],r5[3],r5[4],r5[5])
            # 再遍历纵向 s的值无所谓，主要关心跟车速度 加速度
            v=0
            while v<=maxv:
                v=v+dv
                r4=plane5_best2(xs0,[0,v,0],t)#与s无关
                #绘制图形
                getXY4(t,r4[0],r4[1],r4[2],r4[3],r4[4],"b")
                #计算纵向向损失
                loss1=getSSloss(t,v,vc,dt,r4[0],r4[1],r4[2],r4[3],r4[4])
                fp=Fpath(t,d,v,loss0+loss1,r5,r4)#t,d,v,loss,ssr,ddr
                PathList.append(fp)
    tmp=9999999
    for fp in PathList:
        if fp.loss < tmp:
            tmp=fp.loss
            min=fp
    print("最小",min.loss)
    print("最佳时间",min.t)
    print("最佳横向",min.d)
    print("最佳纵向速度",min.v)
    #绘制最佳
    d,d1,dv1,da1=getXY5(min.t,min.ddr[0],min.ddr[1],min.ddr[2],min.ddr[3],min.ddr[4],min.ddr[5],"y")
    s,s1,sv1,sa1=getXY4(min.t,min.ssr[0],min.ssr[1],min.ssr[2],min.ssr[3],min.ssr[4],"g")
    return d,d1,dv1,da1,s,s1,sv1,sa1

#获取全部路径下 再次计算的 x y 值
def getRange(xp,yp,res):
        l=len(xp)
        xarr=[]
        yarr=[]
        i=1
        d=0.1
        while i<l:
            x=np.arange(xp[i-1],xp[i],d)
            y=[]
            for xx in x:
                y.append(getY(xx,res,i-1))
            yarr.append(y)
            xarr.append(x)
            i=i+1
        return xarr,yarr

#计算各个点的s值 与它对应的坐标也是完全同形状数组
def  getS(xarr,yarr):
    a=[2,5,8,9]
    dxarr=[]
    dyarr=[]
    l=len(xarr)
    i=0
    while i<l:
        dx=np.diff(xarr[i])
        dy=np.diff(yarr[i])
        dx1=[0]
        dy1=[0]
        dx1.extend(dx)
        dy1.extend(dy)
        dxarr.append(dx1)
        dyarr.append(dy1)
        i=i+1
    #计算两两点dx dy距离
    l=len(dxarr)
    i=0
    sarr=[]
    while i<l:
        ll=len(dxarr[i])
        j=0
        ss=[]
        summ=0#求和
        while j<ll:
            s=dxarr[i][j]**2+dyarr[i][j]**2
            s=math.sqrt(s)
            summ=summ+s
            ss.append(summ)
            j=j+1
        sarr.append(ss)
        i=i+1
    #每一段的S距离 对应每一个点
    #print(sarr)
    #把下一个点的距离加上最后一个点的距离
    l=len(sarr)
    i=1
    nsarr=[]
    nsarr.append(sarr[0])
    while i<l:
        ll=len(sarr[i])
        d=sarr[i-1][-1:][0]#累加最后一个元素
        #print(d)
        j=0
        while j<ll:
            sarr[i][j]=sarr[i][j]+d#增加上一段最后值
            j=j+1
        nsarr.append(sarr[i])
        i=i+1
    #print(nsarr)
    return nsarr

def getIndex(s,arr):
    l=len(arr)
    for i in range(l):
        if(arr[i]>s):
            return i
    return None

#给一个S值，找到 最节点的点
def lookAt(s,ress,resx,resy):
    w=4#4个系数
    i=getIndex(s,ress)
    x0=getXY3(s,resx[(i-1)*w+0],resx[(i-1)*w+1],resx[(i-1)*w+2],resx[(i-1)*w+3])
    y0=getXY3(s,resy[(i-1)*w+0],resy[(i-1)*w+1],resy[(i-1)*w+2],resy[(i-1)*w+3])
    return x0,y0,i



#给一个S 一个D值，转换为 笛卡尔位置
def  getDCEXY(s,d,sarr,xarr,yarr,xishu):
    i,j=lookAt(s,sarr)
    x0=xarr[i][j]
    y0=yarr[i][j]
    plt.plot(x0,y0,"yo")
    getSD_DCE(s,d,x0,y0,xishu[i*w+0],xishu[i*w+1],xishu[i*w+2])

#给出系数 #使用 ax3+bx2+cx+d=y 的方式,求这点的夹角，并转换为 x y
def  getSD_DCE(s,d,x0,y0,a3,a2,a1):
    d1=3*a3*x0**2+2*a2*x0+a1
    jd=math.atan(d1)
    #法向单位向量
    x1=-math.sin(jd)
    y1=math.cos(jd)
    x2=x0+d*x1
    y2=y0+d*y1
    plt.plot(x2,y2,"yo")


def main():
    #先拟合图像,得出的是函数图像的系数
    xishu=spline3.mySolve3(xp,yp,True)
    #在计算散点
    #xarr,yarr=getRange(xp,yp,xishu)
    #再计算S
    #sarr=getS(xarr,yarr)
    #传入S,D 求出笛卡尔的值
    #getDCEXY(2,-1,sarr,xarr,yarr,xishu)

    #拟合s 与路点的关系
    ress,resx,resy=getSDXY(xp,yp)
    print(ress)
    #给定一个S字，求出 对应的XY值
    #x0,y0=lookAt(6,ress,resx,resy)
    #plt.plot(x0,y0,"ro")

    #给出t0 状态，规划出最好的轨迹
    #循环显示在笛卡尔中
    #min.t,min.d,min.s,min.v
    #返回 最佳时间 横向 纵向 纵向速度

    xd0=[1,0,0]
    xs0=[0,0,0] 
    """
    darr,d1,dv1,da1,sarr,s1,sv1,sa1=loopGet(xd0,xs0)
    ll=len(darr)
    #print(best_td,best_d,best_ts,best_s,best_v)
    for i in range(ll):
        x0,y0,t=lookAt(sarr[i],ress,resx,resy)
        plt.plot(x0,y0,"ro")
        #转换到迪卡尔坐标系
        getSD_DCE(sarr[i],darr[i],x0,y0,xishu[4*(t-1)+0],xishu[4*(t-1)+1],xishu[4*(t-1)+2])

    """
    for j in range(250):
        i=0
        darr,d1,dv1,da1,sarr,s1,sv1,sa1=loopGet(xd0,xs0)
        ll=len(darr)
        for i in range(ll):
            x0,y0,t=lookAt(sarr[i],ress,resx,resy)
            plt.plot(x0,y0,"ro")
            #转换到迪卡尔坐标系
            getSD_DCE(sarr[i],darr[i],x0,y0,xishu[4*(t-1)+0],xishu[4*(t-1)+1],xishu[4*(t-1)+2])
        xd0=[d1,dv1,da1]
        xs0=[s1,sv1,sa1]
 
     #循环传入最佳S D

    plt.show()
    

if __name__ == '__main__':
    main()
