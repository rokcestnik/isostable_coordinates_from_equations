from eq_phase_rec import *
from matplotlib import pyplot # for a default plot
from math import sin, cos

# the system
def dx(state):
	return state[1] - sin(state[1])*state[0]/2
def dy(state):
	return -state[0] + cos(state[0])*state[1]/2
ders = [dx,dy]

period = oscillator_period(ders)
# phase response curve
prc = oscillator_PRC(ders, [1,0], period, dph=0.05)

floquet = oscillator_floquet(ders,period)
# amplitude response curve
arc = oscillator_ARC(ders, [1,0], period, floquet)

# plot
pyplot.plot(prc[0], prc[1], c='b')
pyplot.plot(arc[0], arc[1], c='y')
pyplot.xticks([0,pi/2,pi,3*pi/2,2*pi],[0,r"$\pi/2$",r"$\pi$",r"$3\pi/2$",r"$2\pi$"])
pyplot.grid(linestyle=':')
pyplot.xlabel(r"$\varphi$")
pyplot.ylabel(r"$Z(\varphi)$")
pyplot.show()
