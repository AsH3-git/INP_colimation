#!/usr/bin/env python3
from constants import *
import numpy
import time
from math import sin, cos, sqrt, pi, exp, atan2, hypot
from scipy.integrate import quad, nquad
from scipy.special import erf, i0
import ROOT

tex0 = 0.0001779/10
x0 = 0.0001686/10
tey0 = 2.23998e-06
y0 = 1.339297e-05/10
#tey0 = tex0
#y0 = x0

E0 = 2000.e6
w0 = h * c / 527e-9 / e
print(w0)

def theta_c(E_0, w_0):
    return me/E_0*sqrt(1.+4.*E_0*w_0/me**2)

def w_max(E_0, w_0, alpha=pi):
    lam = 4 * E_0 * w_0 / me**2 * sin(0.5*alpha)
    return E_0*lam/(1+lam)

def omega_theta(theta, E_0, w_0, alpha=pi):
    return w_max(E_0, w_0, alpha) / ( 1 + (theta / theta_c(E_0, w_0))**2)
    
def theta_omega(w, E_0, w_0, alpha=pi):
    return theta_c(E_0, w_0) * sqrt(w_max(E_0, w_0, alpha) / w - 1) 

def delta_ring(x, y, r):
    if x**2+y**2 == r**2:
        return 1.
    else:
        return 0

L = 50.
hole = 0.0005
N = int(1e6)
wmax = w_max(E0, w0)
print(wmax*1e-6, "MeV")

#'''
with open("tx.dat", "w") as f:
    for n in range(N):
        w = numpy.random.uniform(0, wmax)
        tg = theta_omega(w, E0, w0)
        phi = numpy.random.uniform(0, 2*pi)

#        tgx = tg * numpy.random.choice([-1,1])
        tgx = tg * numpy.cos(phi)
        tgx += numpy.random.normal(0, tex0)
        tgy = tg * numpy.sin(phi)
        tgy += numpy.random.normal(0, tey0)

        x = numpy.random.normal(0, x0)
        y = numpy.random.normal(0, y0)
        
        f.write("{}\t{}\t{}\t{}\t{}\n".format(w, x, tgx, y, tgy))
#'''

# [0]*(erf(x/[1] - 1)+1)*exp(-0.5*(x - [1])^2 / [2]^2)

W = 100e6

def pdf(r, phi, phi0, r0=0, sx=1, sy=1):
    x = (r*cos(phi) - r0*cos(phi0))/sx
    y = (r*sin(phi) - r0*sin(phi0))/sy
#    print(r, phi, phi0)
    return r*exp(-0.5*x**2 - 0.5*y**2) / 2/pi/sx/sy
    


sx = numpy.hypot(tex0, x0/L)
sy = numpy.hypot(tey0, y0/L)
# 1d convolved
'''
func = lambda t : t*exp(-0.5*(t**2+tc**2)/sx**2)*i0(tc*t/sx**2)
r = quad(func, a=0, b=hole/L)[0]
if numpy.isnan(r):
    r = 0
'''

# 100 точек за 5 секунд

#'''
with open("collimation.dat", "w") as f:
#with open("theta_c.dat", "w") as f:
#    for w in numpy.geomspace(wmax, 0.001*wmax, 20):
#    for w in numpy.linspace(0.001*wmax, wmax, 20):
    for tc in [0, *numpy.geomspace(1e-6, hypot(5*sx,hole/L), 50)]:
        t = time.time()
#        tc = theta_omega(w, E0, w0)
        w = omega_theta(tc, E0, w0)
        if tc < hole/L + 5*max(sx, sy):
#            r = nquad(pdf, ranges=[(0, hole/L), (0.5*pi, pi), (0.5*pi, pi)], args=[tc, sx, sy])[0]*16
#            r = nquad(pdf, ranges=[(0, hole/L), (0, pi/2), (0, pi)],
#                      args=[tc, sx, sy],
#                      opts={"epsrel" : 1.e-5, "limit" : 100, "points" : [],
#                            "weight" : None, "wvar" : None, "wopts" : None})[0]*4
            r = nquad(pdf, ranges=[(0, hole/L), (0, pi/2), (0, pi)],
                      args=[tc, sx, sy],
                      opts={"epsrel" : 1.e-4, "limit" : 5, "points" : [],
                            "weight" : None, "wvar" : None, "wopts" : None}, full_output=True)
            print(r)
#            print(r[2])
            r = r[0]*4
        else:
            r = 0
        f.write("{}\t{}\t{}\n".format(w, r, time.time()-t))
        f.flush()
#        f.write("{}\t{}\n".format(w, tc))
#'''

nt_mc = ROOT.TNtuple("nt_mc", "nt_mc", "w:x:tgx:y:tgy")
nt_mc.ReadFile("tx.dat")

nt_an = ROOT.TNtuple("nt_an", "nt_an", "w:rate:time")
#nt_an = ROOT.TNtuple("nt_an", "nt_an", "t:r")
nt_an.ReadFile("collimation.dat")


nt_mc.Draw("w", "", "")
#nt_mc.Draw("sqrt(tgx**2+tgy**2)", "abs(w-{})<2e6".format(wmax-W))
#nt_mc.Draw("abs(tgx)", "abs(w-{})<2e6".format(wmax-W))
#nt_an.Draw("r*{}:t".format(N*1.6), "", "lsame")



ax = ROOT.gPad.GetPrimitive('htemp').GetXaxis()
k = N / ax.GetNbins() * ax.GetXmax() / wmax

nt_mc.Draw("w", "(tgx*{}+x)**2 + (tgy*{}+y)**2 <= {}".format(L, L, hole**2), "same")

nt_an.Draw("rate*{}:w".format(k*0.16), "", "lsame")

#time.sleep(10)
input()


