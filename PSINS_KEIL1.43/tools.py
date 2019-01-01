def q2att(qnb):
    q11 = qnb.q0 * qnb.q0, q12 = qnb.q0 * qnb.q1, q13 = qnb.q0 * qnb.q2, q14 = qnb.q0 * qnb.q3, \
    att=CVect3()                                                                                                                                                 q44 = qnb.q3 * qnb.q3;
    att.i = asinEx(2 * (q34 + q12));
    att.j = atan2Ex(-2 * (q24 - q13), q11 - q22 - q33 + q44);
    att.k = atan2Ex(-2 * (q23 - q14), q11 - q22 + q33 - q44);
    return att;

def asinEx(x):
    return