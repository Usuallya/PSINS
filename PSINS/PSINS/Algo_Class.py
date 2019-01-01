import numpy as np
from tool import *
import Parameter as pa

conefactors = [[2.0 / 3.0, 0.0, 0.0, 0.0],
               [9.0 / 20.0, 27.0 / 20.0, 0.0, 0.0],
               [54.0 / 105.0, 92.0 / 105.0, 214.0 / 105.0, 0.0],
               [250.0 / 504.0, 525.0 / 504.0, 650.0 / 504.0, 1375.0 / 504.0],
               [0.0, 0.0, 0.0, 0.0]]

EPS = float(2.220446049e-16)
INF = float(3.402823466e+30)

class CEarth(object):
    a = b = f = e = e2 = wie = sl = sl2 = sl4 = cl = tl = RMh = RNh = clRNh = f_RMh = f_RNh = f_clRNh = 0.0
    pos = np.array([[0.0]] * 3)
    vn = np.array([[0.0]] * 3)
    wnie = np.array([[0.0]] * 3)
    wnen = np.array([[0.0]] * 3)
    wnin = np.array([[0.0]] * 3)
    gn = np.array([[0.0]] * 3)
    gcc = np.array([[0.0]] * 3)

    def __init__(self, a0=pa.Re, f0=pa.f, g0=pa.g0):
        self.a = a0
        self.f = f0
        self.wie = pa.wie0
        self.b = (1 - self.f) * self.a
        self.e = sqrt(self.a * self.a - self.b * self.b) / self.a
        self.e2 = self.e * self.e
        self.gn = np.array([[0], [0], [-pa.g0]])

    def Update(self, pos, vn=np.array([[0.0]] * 3)):
        self.pos = pos
        self.vn = vn
        self.sl, self.cl = sin(pos[0, 0]), cos(pos[0, 0])
        self.tl = self.sl / self.cl
        sq = 1 - self.e2 * self.sl * self.sl
        sq2 = sqrt(sq)
        self.RMh = self.a * (1 - self.e2) / sq / sq2 + pos[2, 0]
        self.f_RMh = 1.0 / self.RMh
        self.RNh = self.a / sq2 + pos[2, 0]
        self.clRNh = self.cl * self.RNh
        self.f_RNh = 1.0 / self.RNh
        self.f_clRNh = 1.0 / self.clRNh
        self.wnie[0, 0], self.wnie[1, 0], self.wnie[2, 0] = 0.0, self.wie * self.cl, self.wie * self.sl
        self.wnen[0, 0], self.wnen[1, 0] = -vn[1, 0] * self.f_RMh, vn[0, 0] * self.f_RNh,
        self.wnen[2, 0] = self.wnen[1, 0] * self.tl
        self.wnin = self.wnie + self.wnen
        self.sl2 = self.sl * self.sl
        self.sl4 = self.sl2 * self.sl2
        self.gn[2, 0] = -(pa.g0 * (1 + 5.27094e-3 * self.sl2 + 2.32718e-5 * self.sl4) - 3.086e-6 * pos[2, 0])
        self.gcc = self.gn - cross((self.wnie + self.wnin), vn)

    def vn2pos(self, vn, ts=1.0):
        return np.array([[vn[1, 0] * self.f_RMh],
                         [vn[0, 0] * self.f_clRNh],
                         [vn[2, 0]]]) * ts


class CIMU(object):
    nSamples = prefirst = 0
    phim = np.array([[0.0]] * 3)
    dvbm = np.array([[0.0]] * 3)
    wm_1 = np.array([[0.0]] * 3)
    vm_1 = np.array([[0.0]] * 3)

    def __init__(self):
        self.prefirst = 1

    def Update(self, wm, vm, nSamples):
        i = 0
        pcf = conefactors[nSamples - 2]
        cm = np.array([[0.0]] * 3)
        sm = np.array([[0.0]] * 3)
        wmm = np.array([[0.0]] * 3)
        vmm = np.array([[0.0]] * 3)

        self.nSamples = nSamples
        if nSamples == 1:
            if self.prefirst == 1:
                self.wm_1 = wm[0]
                self.vm_1 = vm[0]
                self.prefirst = 0
            cm = 1.0 / 12 * self.wm_1
            self.wm_1 = wm[0]
            sm = 1.0 / 12 * self.vm_1
            self.vm_1 = vm[0]
        if nSamples > 1:
            self.prefirst = 1
        for i in range(nSamples - 1):
            cm += pcf[i] * wm[i]
            sm += pcf[i] * vm[i]
            wmm += wm[i]
            vmm += vm[i]

        wmm += wm[i]
        vmm += vm[i]
        self.phim = wmm + cross(cm, wm[i])
        self.dvbm = vmm + cross(1.0 / 2 * wmm, vmm) + (cross(cm, vm[i]) + cross(sm, wm[i]))


def TauVal(x):
    if x > 3.402823466e+30 / 2:
        x = 0.0
    else:
        x = 1.0 / x


def ISZero(x):
    if -EPS < x[0, 0] < EPS and -EPS < x[1, 0] < EPS and -EPS < x[2, 0] < EPS:
        return True
    else:
        return False


class CSINS(object):
    nts = 0.0
    tk = 0.0
    eth = CEarth()
    imu = CIMU()
    qnb = np.array([[0.0]] * 4)
    Cnb = np.mat([[0.0] * 3] * 3)
    Cnb0 = np.mat([[0.0] * 3] * 3)
    Cbn = np.mat([[0.0] * 3] * 3)
    Kg = np.mat([[0.0] * 3] * 3)
    Ka = np.mat([[0.0] * 3] * 3)
    wib = np.array([[0.0]] * 3)
    fb = np.array([[0.0]] * 3)
    fn = np.array([[0.0]] * 3)
    an = np.array([[0.0]] * 3)
    web = np.array([[0.0]] * 3)
    wnb = np.array([[0.0]] * 3)
    att = np.array([[0.0]] * 3)
    vn = np.array([[0.0]] * 3)
    vb = np.array([[0.0]] * 3)
    pos = np.array([[0.0]] * 3)
    eb = np.array([[0.0]] * 3)
    db = np.array([[0.0]] * 3)
    tauGyro = np.array([[0.0]] * 3)
    tauAcc = np.array([[0.0]] * 3)
    Maa = np.mat([[0.0] * 3] * 3)
    Mav = np.mat([[0.0] * 3] * 3)
    Map = np.mat([[0.0] * 3] * 3)
    Mva = np.mat([[0.0] * 3] * 3)
    Mvv = np.mat([[0.0] * 3] * 3)
    Mvp = np.mat([[0.0] * 3] * 3)
    Mpv = np.mat([[0.0] * 3] * 3)
    Mpp = np.mat([[0.0] * 3] * 3)
    vnL = np.array([[0.0]] * 3)
    posL = np.array([[0.0]] * 3)
    CW = np.mat([[0.0] * 3] * 3)
    MpvCnb = np.mat([[0.0] * 3] * 3)

    def __init__(self, Q, V, P, tk0=0):
        self.tk = tk0
        self.nts = 0.0
        self.qnb = Q
        self.vn = V
        self.pos = P
        self.eth.Update(P, V)
        self.Cnb = Q2M(self.qnb)
        self.att = M2A(self.Cnb)
        self.Cnb0 = self.Cnb
        self.Cbn = self.Cnb.T
        self.vb = self.Cbn * self.vn
        self.Kg = np.mat([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        self.Ka = np.mat([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        self.etm()
        self.lever()

    def SetTauGA(self, tauG, tauA):
        self.tauGyro = np.array([[TauVal(tauG[0, 0])], [TauVal(tauG[1, 0])], [TauVal(tauG[2, 0])]])
        self.tauAcc = np.array([[TauVal(tauA[0, 0])], [TauVal(tauA[1, 0])], [TauVal(tauA[2, 0])]])

    def etm(self):
        tl, secl, wN, wU, vE, vN = self.eth.tl, 1 / self.eth.cl, self.eth.wnie[1, 0], self.eth.wnie[2, 0], \
                                   self.vn[0, 0], self.vn[1, 0]
        secl2 = secl * secl
        f_RMh, f_RNh, f_clRNh = self.eth.f_RMh, self.eth.f_RNh, self.eth.f_clRNh
        f_RMh2, f_RNh2 = f_RMh * f_RMh, f_RNh * f_RNh
        Avn = askew(self.vn)
        Mp1 = np.mat([[0.0, 0.0, 0.0],
                      [-wU, 0.0, 0.0],
                      [wN, 0.0, 0.0]])
        Mp2 = np.mat([[0.0, 0.0, vN * f_RMh2],
                      [0.0, 0.0, -vE * f_RNh2],
                      [vE * secl2 * f_RNh, 0.0, -vE * tl * f_RNh2]])
        self.Maa = askew(-self.eth.wnin)
        self.Mav = np.mat([[0.0, -f_RMh, 0.0],
                           [f_RNh, 0.0, 0.0],
                           [tl * f_RNh, 0.0, 0.0]])
        self.Map = Mp1 + Mp2
        self.Mva = askew(self.fn)
        self.Mvv = Avn * self.Mav - askew(self.eth.wnie + self.eth.wnin)
        self.Mvp = Avn * (Mp1 + self.Map)
        scl = self.eth.sl * self.eth.cl
        self.Mvp[2, 0] = self.Mvp[2, 0] - pa.g0 * (5.27094e-3 * 2 * scl + 2.32718e-5 * 4 * self.eth.sl2 * scl)
        self.Mvp[2, 2] = self.Mvp[2, 2] + 3.086e-6
        self.Mpv = np.mat([[0.0, f_RMh, 0.0],
                           [f_clRNh, 0.0, 0.0],
                           [0.0, 0.0, 1.0]])
        self.Mpp = np.mat([[0.0, 0.0, -vN * f_RMh2],
                           [vE * tl * f_clRNh, 0.0, -vE * secl * f_RNh2],
                           [0.0, 0.0, 0.0]])

    def lever(self, dl=np.array([[0.0]] * 3)):
        self.Mpv = np.mat([[0.0, self.eth.f_RMh, 0.0],
                           [self.eth.f_clRNh, 0.0, 0.0],
                           [0.0, 0.0, 1.0]])
        self.CW = self.Cnb * askew(self.web)
        self.MpvCnb = self.Mpv * self.Cnb
        self.vnL = self.vn + self.CW * dl
        self.posL = self.pos + self.MpvCnb * dl

    def Update(self, wm, vm, nSamples, ts):
        self.nts = nSamples * ts
        self.tk += self.nts
        nts2 = self.nts / 2
        self.imu.Update(wm, vm, nSamples)
        self.imu.phim = self.Kg * self.imu.phim - self.eb * self.nts
        self.imu.dvbm = self.Ka * self.imu.dvbm - self.db * self.nts
        vn01 = self.vn + self.an * nts2
        pos01 = self.pos + self.eth.vn2pos(vn01, nts2)
        self.eth.Update(pos01, vn01)
        self.wib = self.imu.phim / self.nts
        self.fb = self.imu.dvbm / self.nts
        self.web = self.wib - self.Cbn * self.eth.wnie
        self.wnb = self.wib - QV_cross(Q_cross(self.qnb, rv2Q(self.imu.phim / 2)), self.eth.wnin)
        self.fn = QV_cross(self.qnb, self.fb)
        self.an = QV_cross(rv2Q(-self.eth.wnin * nts2), self.fn) + self.eth.gcc
        vn1 = self.vn + self.an * self.nts
        self.pos = self.pos + self.eth.vn2pos(self.vn + vn1, nts2)
        self.vn = vn1
        self.Cnb0 = self.Cnb
        self.qnb = Q_cross(Q_cross(rv2Q(-self.eth.wnin * self.nts), self.qnb), rv2Q(self.imu.phim))
        self.Cnb = Q2M(self.qnb)
        self.att = M2A(self.Cnb)
        self.Cbn = self.Cnb.T
        self.vb = self.Cbn * self.vn


class CKalman(object):
    def __init__(self, nq, nr):
        # val in CSINSTDKF
        self.iter = -2
        self.ifn = 0
        self.tdts = 0.0
        self.measRes = 1
        self.Fk = np.mat([[0.0] * nq] * nq)
        self.Pk1 = np.mat([[0.0] * nq] * nq)
        self.Pxz = np.mat([[0.0]] * nq)
        self.Hi = np.mat([0.0] * nq)
        self.Qk = np.array([[0.0]] * nr)
        self.Kk = np.array([[0.0]] * nr)
        self.tmeas = np.array([[0.0]] * nr)
        self.meanfn = np.array([[0.0]] * 3)

        # val in CSINSKF
        self.sins = CSINS(Q=np.array([[1.0], [0.0], [0.0], [0.0]]),
                          V=np.array([[0.0]] * 3),
                          P=np.array([[0.0]] * 3))

        # val in CKalman
        self.kftk = 0.0
        self.nq = nq
        self.nr = nr
        self.Ft = np.mat([[0.0] * nq] * nq)
        self.Pk = np.mat([[0.0] * nq] * nq)
        self.Hk = np.mat([[0.0] * nq] * nr)
        self.Qt = np.array([[0.0]] * nq)
        self.Pmin = np.array([[0.0]] * nq)
        self.Xk = np.array([[0.0]] * nq)
        self.Pmax = np.array([[INF]] * nq)
        self.Rt = np.array([[0.0]] * nr)
        self.Zk = np.array([[0.0]] * nr)
        self.rts = np.array([[1.0]] * nr)
        self.Rmax = np.array([[INF]] * nr)
        self.Rmin = np.array([[0.0]] * nr)
        self.Rb = np.array([[0.0]] * nr)
        self.Rbeta = np.array([[1.0]] * nr)
        self.FBTau = np.array([[INF]] * nq)
        self.FBMax = np.array([[INF]] * nq)
        self.FBXk = np.array([[0.0]] * nq)
        self.FBTotal = np.array([[0.0]] * nq)
        self.measflag = 0

        # val in CKFApp
        self.tmeas = 0.0
        self.measGPSVn = np.array([[0.0]] * 3)
        self.measGPSPos = np.array([[0.0]] * 3)

    def Init(self, sins0):
        self.sins = sins0
        self.Pmax = np.array(
            [[(10.0 * pa.deg) ** 2], [(10.0 * pa.deg) ** 2], [(30.0 * pa.deg) ** 2],
             [50 ** 2], [50 ** 2], [50 ** 2],
             [(1.0e4 / pa.Re) ** 2], [(1.0e4 / pa.Re) ** 2], [1.0e4 ** 2],
             [(10 * pa.dph) ** 2], [(10 * pa.dph) ** 2], [(10 * pa.dph) ** 2],
             [(10 * pa.mg) ** 2], [(10 * pa.mg) ** 2], [(10 * pa.mg) ** 2]
             ]
        )
        self.Pmin = np.array(
            [[(0.01 * pa.min) ** 2], [(0.01 * pa.min) ** 2], [(0.1 * pa.min) ** 2],
             [0.01 ** 2], [0.01 ** 2], [0.1 ** 2],
             [(1.0 / pa.Re) ** 2], [(1.0 / pa.Re) ** 2], [0.1 ** 2],
             [(0.001 * pa.dph) ** 2], [(0.001 * pa.dph) ** 2], [(0.001 * pa.dph) ** 2],
             [(10 * pa.ug) ** 2], [(10 * pa.ug) ** 2], [(10 * pa.ug) ** 2]
             ]
        )
        self.Pk = np.mat(np.diagflat(
            [(1.0 * pa.deg) ** 2, (1.0 * pa.deg) ** 2, (10.0 * pa.deg) ** 2,
             1.0, 1.0, 1.0,
             (100.0 / pa.Re) ** 2, (100.0 / pa.Re) ** 2, 100.0 ** 2,
             (1.0 * pa.dph) ** 2, (1.0 * pa.dph) ** 2, (1.0 * pa.dph) ** 2,
             (1.0 * pa.mg) ** 2, (1.0 * pa.mg) ** 2, (1.0 * pa.mg) ** 2
             ])
        )
        self.Qt = np.array(
            [[(0.001 * pa.dpsh) ** 2], [(0.001 * pa.dpsh) ** 2], [(0.001 * pa.dpsh) ** 2],
             [(10.0 * pa.ugpsHz) ** 2], [(10.0 * pa.ugpsHz) ** 2], [(10.0 * pa.ugpsHz) ** 2],
             [0.0], [0.0], [0.0],
             [(0.0 * pa.dphpsh) ** 2], [(0.0 * pa.dphpsh) ** 2], [(0.0 * pa.dphpsh) ** 2],
             [(0.0 * pa.ugpsh) ** 2], [(0.0 * pa.ugpsh) ** 2], [(0.0 * pa.ugpsh) ** 2]
             ]
        )
        self.Rt = np.array(
            [[0.2 ** 2], [0.2 ** 2], [0.6 ** 2],
             [(10.0 / pa.Re) ** 2], [(10.0 / pa.Re) ** 2], [30.0 ** 2]
             ]
        )
        self.FBTau = np.array(
            [1.0, 1.0, 10.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0]
        ).reshape(-1, 1)
        self.SetHk()

    def Init16(self, sins0):
        self.sins = sins0
        self.kftk = self.sins.tk
        self.measGPSVn = np.array([[0.0]] * 3)
        self.measGPSPos = np.array([[0.0]] * 3)
        self.Pmax = np.array(
            [[(10.0 * pa.deg) ** 2], [(10.0 * pa.deg) ** 2], [(30.0 * pa.deg) ** 2],
             [50 ** 2], [50 ** 2], [50 ** 2],
             [(1.0e4 / pa.Re) ** 2], [(1.0e4 / pa.Re) ** 2], [1.0e4 ** 2],
             [(1000 * pa.dph) ** 2], [(1000 * pa.dph) ** 2], [(1000 * pa.dph) ** 2],
             [(100 * pa.mg) ** 2], [(100 * pa.mg) ** 2], [(100 * pa.mg) ** 2]
             ]
        )
        self.Pmin = np.array(
            [[(1 * pa.min) ** 2], [(1 * pa.min) ** 2], [(1 * pa.min) ** 2],
             [0.01 ** 2], [0.01 ** 2], [0.1 ** 2],
             [(1.0 / pa.Re) ** 2], [(1.0 / pa.Re) ** 2], [1.0 ** 2],
             [(1.0 * pa.dph) ** 2], [(1.0 * pa.dph) ** 2], [(1.0 * pa.dph) ** 2],
             [(0.1 * pa.mg) ** 2], [(0.1 * pa.mg) ** 2], [(0.1 * pa.mg) ** 2]
             ]
        )
        self.Pk = np.mat(np.diagflat(
            [(1.0 * pa.deg) ** 2, (1.0 * pa.deg) ** 2, (30.0 * pa.deg) ** 2,
             1.0, 1.0, 1.0,
             (100.0 / pa.Re) ** 2, (100.0 / pa.Re) ** 2, 100.0 ** 2,
             (100.0 * pa.dph) ** 2, (100.0 * pa.dph) ** 2, (100.0 * pa.dph) ** 2,
             (10.0 * pa.mg) ** 2, (10.0 * pa.mg) ** 2, (10.0 * pa.mg) ** 2
             ])
        )
        self.Qt = np.array(
            [[(1.0 * pa.dpsh) ** 2], [(1.0 * pa.dpsh) ** 2], [(1.0 * pa.dpsh) ** 2],
             [(100.0 * pa.ugpsHz) ** 2], [(100.0 * pa.ugpsHz) ** 2], [(100.0 * pa.ugpsHz) ** 2],
             [0.0], [0.0], [0.0],
             [0.0], [0.0], [0.0],
             [0.0], [0.0], [0.0]
             ]
        )
        self.Rt = np.array(
            [[0.5 ** 2], [0.5 ** 2], [0.5 ** 2],
             [(10.0 / pa.Re) ** 2], [(10.0 / pa.Re) ** 2], [10.0 ** 2]
             ]
        )
        self.Rmax = self.Rt * 100
        self.Rmin = self.Rt * 0.01
        self.Rb = np.array([[0.9]] * self.nr)
        self.FBTau = np.array(
            [1.0, 1.0, 10.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0]
        ).reshape(-1, 1)
        self.FBMax = np.array(
            [INF, INF, INF, INF, INF, INF, INF, INF, INF, 3000.0 * pa.dph, 3000.0 * pa.dph, 3000.0 * pa.dph,
             50.0 * pa.mg, 50.0 * pa.mg, 50.0 * pa.mg]).reshape(-1, 1)

    def SetFt(self):
        self.sins.etm()
        self.Ft = np.mat(
            [
                [self.sins.Maa[0, 0], self.sins.Maa[0, 1], self.sins.Maa[0, 2], self.sins.Mav[0, 0],
                 self.sins.Mav[0, 1],
                 self.sins.Mav[0, 2], self.sins.Map[0, 0], self.sins.Map[0, 1], self.sins.Map[0, 2],
                 -self.sins.Cnb[0, 0],
                 -self.sins.Cnb[0, 1], -self.sins.Cnb[0, 2], 0.0, 0.0, 0.0],
                [self.sins.Maa[1, 0], self.sins.Maa[1, 1], self.sins.Maa[1, 2], self.sins.Mav[1, 0],
                 self.sins.Mav[1, 1],
                 self.sins.Mav[1, 2], self.sins.Map[1, 0], self.sins.Map[1, 1], self.sins.Map[1, 2],
                 -self.sins.Cnb[1, 0],
                 -self.sins.Cnb[1, 1], -self.sins.Cnb[1, 2], 0.0, 0.0, 0.0],
                [self.sins.Maa[2, 0], self.sins.Maa[2, 1], self.sins.Maa[2, 2], self.sins.Mav[2, 0],
                 self.sins.Mav[2, 1],
                 self.sins.Mav[2, 2], self.sins.Map[2, 0], self.sins.Map[2, 1], self.sins.Map[2, 2],
                 -self.sins.Cnb[2, 0],
                 -self.sins.Cnb[2, 1], -self.sins.Cnb[2, 2], 0.0, 0.0, 0.0],
                [self.sins.Mva[0, 0], self.sins.Mva[0, 1], self.sins.Mva[0, 2], self.sins.Mvv[0, 0],
                 self.sins.Mvv[0, 1],
                 self.sins.Mvv[0, 2], self.sins.Mvp[0, 0], self.sins.Mvp[0, 1], self.sins.Mvp[0, 2], 0.0, 0.0, 0.0,
                 -self.sins.Cnb[0, 0], -self.sins.Cnb[0, 1], -self.sins.Cnb[0, 2]],
                [self.sins.Mva[1, 0], self.sins.Mva[1, 1], self.sins.Mva[1, 2], self.sins.Mvv[1, 0],
                 self.sins.Mvv[1, 1],
                 self.sins.Mvv[1, 2], self.sins.Mvp[1, 0], self.sins.Mvp[1, 1], self.sins.Mvp[1, 2], 0.0, 0.0, 0.0,
                 -self.sins.Cnb[1, 0], -self.sins.Cnb[1, 1], -self.sins.Cnb[1, 2]],
                [self.sins.Mva[2, 0], self.sins.Mva[2, 1], self.sins.Mva[2, 2], self.sins.Mvv[2, 0],
                 self.sins.Mvv[2, 1],
                 self.sins.Mvv[2, 2], self.sins.Mvp[2, 0], self.sins.Mvp[2, 1], self.sins.Mvp[2, 2], 0.0, 0.0, 0.0,
                 -self.sins.Cnb[2, 0], -self.sins.Cnb[2, 1], -self.sins.Cnb[2, 2]],
                [0.0, 0.0, 0.0, self.sins.Mpv[0, 0], self.sins.Mpv[0, 1], self.sins.Mpv[0, 2], self.sins.Mpp[0, 0],
                 self.sins.Mpp[0, 1], self.sins.Mpp[0, 2], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, self.sins.Mpv[1, 0], self.sins.Mpv[1, 1], self.sins.Mpv[1, 2], self.sins.Mpp[1, 0],
                 self.sins.Mpp[1, 1], self.sins.Mpp[1, 2], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, self.sins.Mpv[2, 0], self.sins.Mpv[2, 1], self.sins.Mpv[2, 2], self.sins.Mpp[2, 0],
                 self.sins.Mpp[2, 1], self.sins.Mpp[2, 2], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, self.sins.tauGyro[0, 0], 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, self.sins.tauGyro[1, 0], 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, self.sins.tauGyro[2, 0], 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, self.sins.tauAcc[0, 0], 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, self.sins.tauAcc[1, 0], 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, self.sins.tauAcc[2, 0]]
            ]
        )

    def SetHk(self):
        self.Hk[0, 3] = self.Hk[1, 4] = self.Hk[2, 5] = 1.0
        self.Hk[3, 6] = self.Hk[4, 7] = self.Hk[5, 8] = 1.0

    def TimeUpdate(self, kfts, fback=1):
        self.kftk += kfts
        self.SetFt()
        Fk = np.eye(self.nq) + self.Ft * kfts
        self.Xk = Fk * self.Xk
        self.Pk = Fk * self.Pk * Fk.T
        self.Pk += self.Qt * kfts
        if (fback):
            self.Feedback(kfts)

    def SetMeasFlag(self, flag):
        if flag == 0:
            self.measflag = 0
        else:
            self.measflag = self.measflag | flag

    def MeasUpdate(self, fading=1):
        self.SetMeas()
        for i in range(self.nr):
            if self.measflag & (0x01 << i):
                Hi = self.Hk[i]
                Pxz = self.Pk * Hi.T
                Pz0 = (Hi * Pxz)[0, 0]
                r = self.Zk[i, 0] - (Hi * self.Xk)[0, 0]
                self.RAdaptive(i, r, Pz0)
                Pzz = Pz0 + self.Rt[i, 0] / self.rts[i, 0]
                Kk = Pxz * (1.0 / Pzz)
                self.Xk += Kk * r
                self.Pk -= Kk * Pxz.T
        if fading > 1.0:
            self.Pk *= fading
        self.PKConstrain()
        self.Pk = symmetry(self.Pk)
        self.SetMeasFlag(0)

    def RAdaptive(self, i, r, Pr):
        if self.Rb[i, 0] > EPS:
            rr = r * r - Pr
            if rr < self.Rmin[i, 0]:
                rr = self.Rmin[i, 0]
            if rr > self.Rmax[i, 0]:
                self.Rt[i, 0] = self.Rmax[i, 0]
            else:
                self.Rt[i, 0] = (1.0 - self.Rbeta[i, 0]) * self.Rt[i, 0] + self.Rbeta[i, 0] + rr
            self.Rbeta[i, 0] = self.Rbeta[i, 0] / (self.Rbeta[i, 0] + self.Rb[i, 0])

    def PKConstrain(self):
        i, nq1 = 0, self.nq + 1
        while i < 15:
            if self.Pmin[i, 0] > self.Pk[i, i] > EPS:
                self.Pk[i, i] = self.Pmin[i, 0]
            elif self.Pk[i, i] > self.Pmax[i, 0]:
                sqf = sqrt(self.Pmax[i, 0] / self.Pk[i, i]) * 0.5
                for k in range(self.nq):
                    self.Pk[i, k] *= sqf
                    self.Pk[k, i] *= sqf
            i += 1

    def Feedback(self, fbts):
        for i in range(self.nq):
            if self.FBTau[i, 0] < INF / 2:
                if fbts < self.FBTau[i, 0]:
                    afa = fbts / self.FBTau[i, 0]
                else:
                    afa = 1.0
                self.FBXk[i, 0] = self.Xk[i, 0] * afa
                if self.FBTotal[i, 0] + self.FBXk[i, 0] > self.FBMax[i, 0]:
                    self.FBXk[i, 0] = self.FBMax[i, 0] - self.FBTotal[i, 0]
                elif self.FBTotal[i, 0] + self.FBXk[i, 0] < -self.FBMax[i, 0]:
                    self.FBXk[i, 0] = -self.FBMax[i, 0] - self.FBTotal[i, 0]
                self.Xk[i, 0] -= self.FBXk[i, 0]
                self.FBTotal[i, 0] += self.FBXk[i, 0]
            else:
                self.FBXk[i, 0] = 0.0
        self.sins.qnb = Q_cross(rv2Q(np.array([[self.FBXk[0, 0]], [self.FBXk[1, 0]], [self.FBXk[2, 0]]])), self.sins.qnb)
        self.sins.vn -= np.array([[self.FBXk[3, 0]], [self.FBXk[4, 0]], [self.FBXk[5, 0]]])
        self.sins.pos -= np.array([[self.FBXk[6, 0]], [self.FBXk[7, 0]], [self.FBXk[8, 0]]])
        self.sins.eb += np.array([[self.FBXk[9, 0]], [self.FBXk[10, 0]], [self.FBXk[11, 0]]])
        self.sins.db += np.array([[self.FBXk[12, 0]], [self.FBXk[13, 0]], [self.FBXk[14, 0]]])

    def Update(self, wm, vm, nSamples, ts):
        self.sins.Update(wm, vm, nSamples, ts)
        self.TimeUpdate(self.sins.nts)
        self.kftk = self.sins.tk
        self.MeasUpdate()

    def TDUpdate(self, wm, vm, nSamples, ts, nStep):
        self.sins.Update(wm, vm, nSamples, ts)
        self.Feedback(self.sins.nts)

        self.measRes = 0
        if nStep <= 0:
            nStep = 2 * (self.nq + self.nr) + 3
        self.tdts += self.sins.nts
        self.kftk = self.sins.tk
        self.meanfn += self.sins.fn
        self.ifn += 1
        for i in range(nStep):
            if self.iter == -2:
                if self.ifn == 0:
                    break
                vtemp = self.meanfn * 1.0 / self.ifn
                self.meanfn = np.array([[0.0]] * 3)
                self.ifn = 0
                self.sins.fn = vtemp
                self.SetFt()
                self.sins.fn = vtemp
                self.SetMeas()
            elif self.iter == -1:
                self.Fk = np.eye(self.nq) + self.Ft * self.tdts
                self.Qk = self.Qt * self.tdts
                self.Xk = self.Fk * self.Xk
                self.tdts = 0.0
            elif self.iter < self.nq:
                row = self.iter
                for i in range(self.nq):
                    f = 0.0
                    for j in range(self.nq):
                        f += self.Fk[row, j] * self.Pk[j, i]
                    self.Pk1[row, i] = f
            elif self.iter < self.nq * 2:
                row = self.iter - self.nq
                for i in range(self.nq):
                    f = 0.0
                    for j in range(self.nq):
                        f += self.Pk1[row, j] * self.Fk.T[j, i]
                    self.Pk[row, i] = f
                if row == self.nq - 1:
                    self.Pk += self.Qk
            elif self.iter < 2 * (self.nq + self.nr):
                row = int((self.iter - 2 * self.Ft.shape[0]) / 2)
                flag = self.measflag & (0x01 << row)
                if flag:
                    if (self.iter - 2 * self.Ft.shape[0]) % 2 == 0:
                        self.Hi = self.Hk[row]
                        self.Pxz = self.Pk * self.Hi.T
                        Pzz = (self.Hi * self.Pxz)[0, 0] + self.Rt[row, 0] / self.rts[row, 0]
                        self.Kk = self.Pxz * 1.0 / Pzz
                    else:
                        self.measRes |= flag
                        r = self.Zk[row, 0] - (self.Hi * self.Xk)[0, 0]
                        self.RAdaptive(row, r, (self.Hi * self.Pxz)[0, 0])
                        self.Xk += self.Kk * r
                else:
                    nStep += 1
            elif self.iter >= 2 * (self.nq + self.nr):
                self.PKConstrain()
                self.Pk = symmetry(self.Pk)
                self.SetMeasFlag(0)
                self.iter = -3
            self.iter += 1

    def SetMeas(self, vnm=None, posm=None, tm=None):
        if vnm is not None and posm is not None and tm is not None:
            self.measGPSVn = vnm
            self.measGPSPos = posm
            self.tmeas = tm
        else:
            dt = self.sins.tk - self.tmeas - 0.01 * 0
            if fabs(dt) > 0.5:
                return
            if not ISZero(self.measGPSVn):
                self.Zk[0, 0] = self.sins.vn[0, 0] - self.measGPSVn[0, 0] - (self.sins.an * dt)[0, 0]
                self.Zk[1, 0] = self.sins.vn[1, 0] - self.measGPSVn[1, 0] - (self.sins.an * dt)[1, 0]
                self.Zk[2, 0] = self.sins.vn[2, 0] - self.measGPSVn[2, 0] - (self.sins.an * dt)[2, 0]
                self.SetMeasFlag(7)
            if not ISZero(self.measGPSPos):
                self.Zk[3, 0] = self.sins.pos[0, 0] - self.measGPSPos[0, 0] - \
                                self.sins.eth.vn2pos(self.sins.vn, dt)[0, 0]
                self.Zk[4, 0] = self.sins.pos[1, 0] - self.measGPSPos[1, 0] - \
                                self.sins.eth.vn2pos(self.sins.vn, dt)[1, 0]
                self.Zk[5, 0] = self.sins.pos[2, 0] - self.measGPSPos[2, 0] - \
                                self.sins.eth.vn2pos(self.sins.vn, dt)[2, 0]
                self.SetMeasFlag(56)
            if self.measflag != 0:
                self.SetHk()
                self.measGPSVn = np.array([[0.0]] * 3)
                self.measGPSPos = np.array([[0.0]] * 3)
