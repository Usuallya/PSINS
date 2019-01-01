from math import *

Re, f, wie0, g0 = 6378137.0, 1.0 / 298.257, 7.2921151467e-5, 9.7803267714
e = sqrt(2 * f - f * f)
e2 = e * e
mg = 1.0e-3 * g0
ug = 1.0e-6 * g0
deg = pi / 180.0
min = deg / 60.0
sec = min / 60.0
ppm = 1.0e-6
hur = 3600.0
dps = deg / 1.0
dph = deg / hur
dpsh = deg / sqrt(hur)
dphpsh = dph / sqrt(hur)
ugpsHz = ug / sqrt(1.0)
ugpsh = ug / sqrt(hur)
mpsh = 1 / sqrt(hur)
mpspsh = 1 / 1 / sqrt(hur)
ppmpsh = ppm / sqrt(hur)
secpsh = sec / sqrt(hur)