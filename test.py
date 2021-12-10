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
prc = oscillator_PRC(ders, [1,0], period)

# amplitude response curve
floquet = oscillator_floquet(ders,period)
arc = oscillator_ARC(ders, [1,0], period, floquet)

# plot
"""pyplot.plot(prc[0], prc[1], c='b', label=r"$Z(\varphi)$")
pyplot.plot(arc[0], arc[1], c='#DD9900', label=r"$A(\varphi)$")
pyplot.plot([0,2*pi],[0,0], c='k', alpha=0.4)
pyplot.legend()
pyplot.xticks([0,pi/2,pi,3*pi/2,2*pi],[0,r"$\pi/2$",r"$\pi$",r"$3\pi/2$",r"$2\pi$"])
pyplot.grid(linestyle=':')
pyplot.xlabel(r"$\varphi$")
pyplot.ylabel(r"$Z(\varphi), A(\varphi)$")
pyplot.show()"""

# limit cycle
"""
limitc = sample_limit_cycle(ders, 250, period)
iso_out = sample_local_isostable(ders, 250, period, floquet, 1)
iso_in = sample_local_isostable(ders, 250, period, floquet, -1)

limitc_phases = [oscillator_phase(limitc[i],ders,period) for i in range(len(limitc))]
iso_out_phases = [oscillator_phase(iso_out[i],ders,period) for i in range(len(iso_out))]
iso_in_phases = [oscillator_phase(iso_in[i],ders,period) for i in range(len(iso_in))]

limitc_ampls = [oscillator_amplitude(limitc[i], ders, period, floquet, limitc[0]) for i in range(len(limitc))]
iso_out_ampls = [oscillator_amplitude(iso_out[i],ders, period, floquet, limitc[0]) for i in range(len(iso_out))]
iso_in_ampls = [oscillator_amplitude(iso_in[i], ders, period, floquet, limitc[0]) for i in range(len(iso_in))]
"""

"""pyplot.plot(limitc_phases)
pyplot.plot(iso_out_phases)
pyplot.plot(iso_in_phases)
pyplot.show()

pyplot.plot(limitc_ampls)
pyplot.plot(iso_out_ampls)
pyplot.plot(iso_in_ampls)
pyplot.show()

pyplot.plot([limitc[i][0] for i in range(len(limitc))], [limitc[i][1] for i in range(len(limitc))])
pyplot.plot([iso_out[i][0] for i in range(len(iso_out))], [iso_out[i][1] for i in range(len(iso_out))])
pyplot.plot([iso_in[i][0] for i in range(len(iso_in))], [iso_in[i][1] for i in range(len(iso_in))])
pyplot.show()"""

local_iso_in = sample_local_isostable(ders, 250, period, floquet, -1)
local_iso_out = sample_local_isostable(ders, 250, period, floquet, 1)
isostables = oscillator_isostables(ders, local_iso_in, local_iso_out, period, floquet)
for iso in isostables:
	pyplot.plot([iso[i][0] for i in range(len(iso))],[iso[i][1] for i in range(len(iso))])
pyplot.show()
