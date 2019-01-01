import numpy as np
import Variables as me

def INSUartOut(sins):
    att = q2att(sins.qnb)
    out_data.Att[0]=att.i/glv.deg
    out_data.Att[1]=att.j/glv.deg
    out_data.Att[2]=CC180toC360(att.k)/glv.deg
    out_data.Vn[0]=sins.vn.i
    out_data.Vn[1]=sins.vn.j
    out_data.Vn[2]=sins.vn.k

    deg = sins.pos.j/glv.deg
    out_data.Pos[0]=deg
    out_data.Pos[1]=sins.pos.j/glv.deg-deg
    deg=sins.pos.i/glv.deg
    out_data.Pos[2]=deg
    out_data.Pos[3]=sins.pos.i/glv.deg-deg
    out_data.Pos[4]=sins.pos.k
    #TIM未知量？
    timcnt2=TIM2.CNT
    #TIM未知量？
    time_consume = (timcnt2-timcnt1)*0.1
    if time_cosume1<time_consume:
        time_consume1 = time_consume
    if math.fmod(sins.tk,5)<ts:
        time_consume1 = 0.0

class CSINS(object):
    def __init__(self,qnb0,vn0,pos0,tk0):

    def SetTauGA(self,tauG,tauA):

    def Update(selfpwm,pvm,nSamples,ts):

    def Extrap(selfwm,vm,ts):

    def lever(self,dl):

    def etm(self):

class CMahony(object):



class CKalman(object):
    def __init__(self,nq,nr):
        if nq>me.MMD or nr>me.MMD:
            return

        self.kftk=0.0
        self.Ft=self.Pk = np.mat([[0.0]*nq]*nq)
        self.Hk=np.mat([[0.0]*nq]*nr)
        self.Fading=np.mat([[1.0]*nr]*nq)
        self.zfdafa=0.1
        self.Qt = self.Pmin=self.Xk = np.zeros(nq,dtype='float').reshape(-1,1)
        self.Xmax=self.Pmax=np.full(nq,me.INF).reshape(-1,1)
        self.Zk = np.zeros(nr,dtype='float').reshape(-1,1)
        self.Rt = np.full(nr,me.INF).reshape(-1,1)
        self.rts = np.ones(nr,dtype='float').reshape(-1,1)
        self.Zfd = np.zeros(nr,dtype='float').reshape(-1,1)
        self.Zfd0 = np.full(nr,me.INF).reshape(-1,1)
        self.Rmax = np.full(nr,me.INF).reshape(-1,1)
        self.Rmin = np.zeros(nr,dtype='float').reshape(-1,1)
        self.Rbeta = np.ones(nr,dtype='float').reshape(-1,1)
        # for i in range(0,nr):
        #     self.Rmaxcount[me.MMD]=0
        #     self.Rmaxcount0[i]=5
        self.Rmaxcount=0
        self.Rmaxcount0=np.eye(nr,dtype='float')*5
        self.FBTau = self.FBMax = np.array([[me.INF]]*nq)
        self.FBXk = self.FXTotal = np.array([[0.0]]*nq)
        self.measFlat = 0

    def TimeUpdate(self,kfts,fback):

    def SetMeasFlag(self,flag):

    def MeasUpdate(self,fading):

    def RAdaptive(self,i,r,Pr):

    def RPkFading(self,i):

    def XPConstrain(self):

    def Feedback(self,fbts):

class CSINSKF(CKalman):
    def __init__(self,nq,nr):
        super(nq,nr)
    def init(self,sins0,grade):

    def SetFt(self):

    def SetHk(self):

    def Update(self,pwm,pvm,nSamples,ts):

    def Feedback(self,fbts):
        #先调用父类的FeedBack方法
    def QtMarkovGA(self,tauG,sRG,tauA=I31,sRA=031):



class CSINSTDKF(CSINSKF):
    def __init__(self,nq,nr):
        super(nq,nr)
    def TDReset(self):
        self.iter=-2
        self.ifn=0
        #八进制的031
        self.meanfn = 031
    def TDUpdate(self,pwm,pvm,nSamples,ts,nStep):

class CSINSGPSOD(CSINSTDKF):
    def __init__(self,nq=15,nr=10):

    #这里有ifndef语句
    def init(self,sins0,grade):

    def SetMeas(self):

    def SetMeasGPS(self,pgps,vgps):

    def SetMeasOD(self,dSod,ts):

    def SetMeasYaw(self,ymag):
        ...
    ...




class CKFApp(CSINSTDKF):
    ...


class CQEAHRS(CKalman):
    def __init__(self, ts):
        super(7, 3)


