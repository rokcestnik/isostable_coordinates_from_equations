from eq_phase_rec import *
from matplotlib import pyplot
from math import sin, cos

# the system
def dx(state):
	#return 1/6*state[0]+1*state[1]-(state[0]**2+state[1]**2)*(1/6*state[0]+0.7/6*state[1]) # SL (doesnt work very well)
	#return state[1] - sin(state[1])*state[0]/2 # sincos
	#return state[1] - state[0]/2 # cauchy
	#return state[1] # vdp
	return 1.5+state[0]**2*state[1]-3.5*state[0]-state[0] # brusselator
	#return 1*(1-state[0]**2)*state[0]-8*state[1] # rayleigh
def dy(state):
	#return 1/6*state[1]-1*state[0]-(state[0]**2+state[1]**2)*(1/6*state[1]-0.7/6*state[0]) # SL (doesnt work very well)
	#return -state[0] + cos(state[0])*state[1]/2 # sincos
	#return -state[0] + state[1]/(1+state[0]**2) # cauchy
	#return -state[0] + 0.5*(1-state[0]**2)*state[1] # vdp
	return 3.5*state[0]-state[0]**2*state[1] # brusselator
	#return state[0] # rayleigh
ders = [dx,dy]
thresh = 0.0 # most everything
thresh = 1.6 # brusellator

state = [1,1]
states = []
for t in range(100000):
	state = one_step_integrator(state, ders, 0.005)
	states.append(state)
pyplot.plot([states[i][0] for i in range(100000)])
pyplot.show()

# phase response curve
period = oscillator_period(ders, thr=thresh)
prc = oscillator_PRC(ders, [1,0], period, thr=thresh)
# amplitude response curve
floquet = oscillator_floquet(ders, period, thr=thresh)
arc = oscillator_ARC(ders, [1,0], period, floquet, thr=thresh)

# plot
pyplot.plot(prc[0], prc[1], c='b', label=r"$Z(\varphi)$")
pyplot.plot(arc[0], arc[1], c='#DD9900', label=r"$A(\varphi)$")
pyplot.plot([0,2*pi],[0,0], c='k', alpha=0.4)
pyplot.legend()
pyplot.xticks([0,pi/2,pi,3*pi/2,2*pi],[0,r"$\pi/2$",r"$\pi$",r"$3\pi/2$",r"$2\pi$"])
pyplot.grid(linestyle=':')
pyplot.xlabel(r"$\varphi$")
pyplot.ylabel(r"$Z(\varphi), A(\varphi)$")
pyplot.savefig("prc_arc.png", dpi=300)
pyplot.show()


# isospace
limitc = sample_limit_cycle(ders, 250, period, thr=thresh)
#isochrons, isostables = oscillator_isochrons_isostables(ders, initial_shift=0.025, isochron_resolution=10, number_of_isostables=10, amplitude_unit=0.015, amplitude_domain=10*0.015, dt=0.003) # SL (doesnt work very well)
#isochrons, isostables = oscillator_isochrons_isostables(ders, isochron_resolution=10, number_of_isostables=9, amplitude_unit=0.125, amplitude_domain=9*0.125, dt=0.003) # sincos
#isochrons, isostables = oscillator_isochrons_isostables(ders, isochron_resolution=10, number_of_isostables=7, amplitude_unit=0.07, amplitude_domain=7*0.07, dt=0.003) # cauchy
#isochrons, isostables = oscillator_isochrons_isostables(ders, isochron_resolution=10, number_of_isostables=5, amplitude_unit=0.06, amplitude_domain=5*0.06, dt=0.003) # vdp
isochrons, isostables = oscillator_isochrons_isostables(ders, isochron_resolution=10, number_of_isostables=6, amplitude_unit=0.01, amplitude_domain=6*0.01, thr=thresh, dt=0.003) # brusellator
#isochrons, isostables = oscillator_isochrons_isostables(ders, isochron_resolution=10, number_of_isostables=5, amplitude_unit=0.01, amplitude_domain=5*0.01, dt=0.003) # rayleigh


# plot
for iso in isostables:
	pyplot.plot([iso[i][0] for i in range(len(iso))],[iso[i][1] for i in range(len(iso))], c='#888888')
for iso in isochrons:
	pyplot.plot([iso[i][0] for i in range(len(iso))],[iso[i][1] for i in range(len(iso))], c='#CC8800')
pyplot.plot([limitc[i][0] for i in range(len(limitc))],[limitc[i][1] for i in range(len(limitc))], lw=2, c='#000000')
pyplot.xlabel(r"$x$")
pyplot.ylabel(r"$y$")
pyplot.axes().set_aspect('equal')
pyplot.savefig("isospace.png", dpi=300)
pyplot.show()
