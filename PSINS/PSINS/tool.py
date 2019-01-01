from math import *
import numpy as np


def M2Q(Cnb):
    q0, q1, q2, q3, qq4 = 0.0, 0.0, 0.0, 0.0, 0.0
    if Cnb[0, 0] >= Cnb[1, 1] + Cnb[2, 2]:
        q1 = 0.5 * sqrt(1 + Cnb[0, 0] - Cnb[1, 1] - Cnb[2, 2])
        qq4 = 4 * q1
        q0 = (Cnb[2, 1] - Cnb[1, 2]) / qq4
        q2 = (Cnb[0, 1] + Cnb[1, 0]) / qq4
        q3 = (Cnb[0, 2] + Cnb[2, 0]) / qq4
    elif Cnb[1, 1] >= Cnb[0, 0] + Cnb[2, 2]:
        q2 = 0.5 * sqrt(1 - Cnb[0, 0] + Cnb[1, 1] - Cnb[2, 2])
        qq4 = 4 * q2
        q0 = (Cnb[0, 2] - Cnb[2, 0]) / qq4
        q1 = (Cnb[0, 1] - Cnb[1, 0]) / qq4
        q3 = (Cnb[1, 2] + Cnb[2, 1]) / qq4
    elif Cnb[2, 2] >= Cnb[0, 0] + Cnb[1, 1]:
        q3 = 0.5 * sqrt(1 - Cnb[0, 0] - Cnb[1, 1] + Cnb[2, 2])
        qq4 = 4 * q3
        q0 = (Cnb[1, 0] - Cnb[0, 1]) / qq4
        q1 = (Cnb[0, 2] + Cnb[2, 0]) / qq4
        q2 = (Cnb[1, 2] + Cnb[2, 1]) / qq4
    else:
        q0 = 0.5 * sqrt(1 + Cnb[0, 0] + Cnb[1, 1] + Cnb[2, 2])
        qq4 = 4 * q2
        q1 = (Cnb[2, 1] - Cnb[1, 2]) / qq4
        q2 = (Cnb[0, 2] - Cnb[2, 0]) / qq4
        q3 = (Cnb[1, 0] - Cnb[0, 1]) / qq4
    nq = sqrt(q0 ** 2 + q1 ** 2 + q2 ** 2 + q3 ** 2)
    q0 /= nq
    q1 /= nq
    q2 /= nq
    q3 /= nq
    return np.array([q0, q1, q2, q3]).reshape(4, 1)


def A2Q(Eular):
    pitch, roll, yaw = Eular[0, 0] / 2, Eular[1, 0] / 2, Eular[2, 0] / 2
    sp, sr, sy = sin(pitch), sin(roll), sin(yaw)
    cp, cr, cy = cos(pitch), cos(roll), cos(yaw)
    q0 = cp * cr * cy - sp * sr * sy
    q1 = sp * cr * cy - cp * sr * sy
    q2 = cp * sr * cy + sp * cr * sy
    q3 = cp * cr * sy + sp * sr * cy
    return np.array([q0, q1, q2, q3]).reshape(4, 1)


def Q2M(Q):
    q11, q12, q13, q14 = Q[0, 0] * Q[0, 0], Q[0, 0] * Q[1, 0], Q[0, 0] * Q[2, 0], Q[0, 0] * Q[3, 0]
    q22, q23, q24 = Q[1, 0] * Q[1, 0], Q[1, 0] * Q[2, 0], Q[1, 0] * Q[3, 0]
    q33, q34 = Q[2, 0] * Q[2, 0], Q[2, 0] * Q[3, 0]
    q44 = Q[3, 0] * Q[3, 0]

    return np.mat([[q11 + q22 - q33 - q44, 2 * (q23 - q14), 2 * (q24 + q13)],
                   [2 * (q23 + q14), q11 - q22 + q33 - q44, 2 * (q34 - q12)],
                   [2 * (q24 - q13), 2 * (q34 + q12), q11 - q22 - q33 + q44]])


def M2A(Cnb):
    i = asinEx(Cnb[2, 1])
    j = atan2Ex(-Cnb[2, 0], Cnb[2, 2])
    k = atan2Ex(-Cnb[0, 1], Cnb[1, 1])
    return np.array([[i], [j], [k]])


def p_range(val, min, max):
    if val > max:
        val = max
    elif val < min:
        val = min
    return val


def atan2Ex(y, x):
    if sign(y) == 0 and sign(x) == 0:
        return 0.0
    else:
        return atan2(y, x)


def sign(x):
    if x < 2.220446049e-16 * -1:
        return -1
    elif x > 2.220446049e-16:
        return 1
    else:
        return 0


def asinEx(x):
    return asin(p_range(x, -1, 1))


def askew(v):
    return np.mat([[0.0, -v[2, 0], v[1, 0]],
                   [v[2, 0], 0.0, -v[0, 0]],
                   [-v[1, 0], v[0, 0], 0.0]])


def cross(v1, v2):
    return np.array([[v1[1, 0] * v2[2, 0] - v1[2, 0] * v2[1, 0]],
                     [v1[2, 0] * v2[0, 0] - v1[0, 0] * v2[2, 0]],
                     [v1[0, 0] * v2[1, 0] - v1[1, 0] * v2[0, 0]]
                     ])


def Q_cross(Q1, Q2):
    q0 = Q1[0, 0] * Q2[0, 0] - Q1[1, 0] * Q2[1, 0] - Q1[2, 0] * Q2[2, 0] - Q1[3, 0] * Q2[3, 0]
    q1 = Q1[0, 0] * Q2[1, 0] + Q1[1, 0] * Q2[0, 0] + Q1[2, 0] * Q2[3, 0] - Q1[3, 0] * Q2[2, 0]
    q2 = Q1[0, 0] * Q2[2, 0] + Q1[2, 0] * Q2[0, 0] + Q1[3, 0] * Q2[1, 0] - Q1[1, 0] * Q2[3, 0]
    q3 = Q1[0, 0] * Q2[3, 0] + Q1[3, 0] * Q2[0, 0] + Q1[1, 0] * Q2[2, 0] - Q1[2, 0] * Q2[1, 0]
    return np.array([q0, q1, q2, q3]).reshape(4, 1)


def QV_cross(Q, v):
    q0 = -Q[1, 0] * v[0, 0] - Q[2, 0] * v[1, 0] - Q[3, 0] * v[2, 0]
    q1 = Q[0, 0] * v[0, 0] + Q[2, 0] * v[2, 0] - Q[3, 0] * v[1, 0]
    q2 = Q[0, 0] * v[1, 0] + Q[3, 0] * v[0, 0] - Q[1, 0] * v[2, 0]
    q3 = Q[0, 0] * v[2, 0] + Q[1, 0] * v[1, 0] - Q[2, 0] * v[0, 0]

    i = -q0 * Q[1, 0] + q1 * Q[0, 0] - q2 * Q[3, 0] + q3 * Q[2, 0]
    j = -q0 * Q[2, 0] + q2 * Q[0, 0] - q3 * Q[1, 0] + q1 * Q[3, 0]
    k = -q0 * Q[3, 0] + q3 * Q[0, 0] - q1 * Q[2, 0] + q2 * Q[1, 0]

    return np.array([i, j, k]).reshape(3, 1)


def rv2Q(rv):
    F1 = 2.0 * 1.0
    F2 = F1 * 2 * 2
    F3 = F2 * 2 * 3
    F4 = F3 * 2 * 4
    F5 = F4 * 2 * 5
    n2 = rv[0, 0] * rv[0, 0] + rv[1, 0] * rv[1, 0] + rv[2, 0] * rv[2, 0]
    if n2 < pi / 180 * pi / 180:
        n4 = n2 * n2
        c = 1.0 - n2 * (1.0 / F2) + n4 * (1.0 / F4)
        f = 0.5 - n2 * (1.0 / F3) + n4 * (1.0 / F5)
    else:
        n_2 = sqrt(n2) / 2.0
        c = cos(n_2)
        f = sin(n_2) / n_2 * 0.5
    return np.array([c, f * rv[0, 0], f * rv[1, 0], f * rv[2, 0]]).reshape(4, 1)


def symmetry(m):
    for i in range(m.shape[0]):
        for j in range(i + 1, m.shape[0]):
            m[i, j] = m[j, i] = (m[i, j] + m[j, i]) * 0.5
    return m
