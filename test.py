from eq_phase_rec import *
from matplotlib import pyplot
from math import sin, cos

# the system
def dx(state):
	return state[1] - sin(state[1])*state[0]/2
	#return state[1]
def dy(state):
	return -state[0] + cos(state[0])*state[1]/2
	#return 0.5*(1-state[0]**2)*state[1]-state[0]
ders = [dx,dy]

# phase response curve
period = oscillator_period(ders)
prc = oscillator_PRC(ders, [1,0], period, dph=0.05)

# amplitude response curve
floquet = oscillator_floquet(ders,period)
arc = oscillator_ARC(ders, [1,0], period, floquet, stimulation=0.005, dph=0.05)

# plot
pyplot.plot(prc[0], prc[1], c='b', label=r"$Z(\varphi)$")
pyplot.plot(arc[0], arc[1], c='#DD9900', label=r"$A(\varphi)$")
pyplot.plot([0,2*pi],[0,0], c='k', alpha=0.4)
pyplot.legend()
pyplot.xticks([0,pi/2,pi,3*pi/2,2*pi],[0,r"$\pi/2$",r"$\pi$",r"$3\pi/2$",r"$2\pi$"])
pyplot.grid(linestyle=':')
pyplot.xlabel(r"$\varphi$")
pyplot.ylabel(r"$Z(\varphi), A(\varphi)$")
pyplot.show()
