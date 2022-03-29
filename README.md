Reconstruct the [isostable coordinates](http://www.scholarpedia.org/article/Isochron)
======================================================================================

[iso_rec_eq.py](iso_rec_eq.py) is a standalone set of Python functions that numerically estimate the isostable coordinates (phase and amplitude) from dynamical equations of an arbitrary 2D limit cycle system. All functions are documented with docstrings, run: "`pydoc iso_rec_eq`" in the same directory to see its compilation.

**Useful functions:** <br />
	- `oscillator_period()` - estimates the natural period, <br />
	- `oscillator_floquet()` - estimates the floquet exponent, <br /> 
	- `oscillator_phase()` - estimates the asymptotic phase of a chosen state, <br />
	- `oscillator_amplitude()` - estimates the isostable amplitude of a chosen state, <br />
	- `oscillator_PRC()` - estimates the phase response curve, <br />
	- `oscillator_ARC()` - estimates the isostable amplitude response curve, <br />
	- `oscillator_isochrons()` - estimates isochrons, <br />
	- `oscillator_isostables()` - estimates isostables, <br />
	- `oscillator_isochrons_isostables()` - estimates both isochrons and isostables (sharing computations). <br />
 
An example program [example.py](example.py) computes the phase/amplitude response and isochrons/isostables for several examples oscillators.
Minimal code:
```python
from iso_rec_eq import *
from matplotlib import pyplot
from math import sin, cos

# define the system dynamics
def dx(state):
	return state[1] - sin(state[1])*state[0]/2 
def dy(state):
	return -state[0] + cos(state[0])*state[1]/2 
ders = [dx,dy]

# estimate isospace
limitc = sample_limit_cycle(ders, sampling=250)
isochrons, isostables = oscillator_isochrons_isostables(ders)

# plot
pyplot.axes().set_aspect('equal')
for iso in isostables:
	pyplot.plot([iso[i][0] for i in range(len(iso))],[iso[i][1] for i in range(len(iso))], c='#888888')
for iso in isochrons:
	pyplot.plot([iso[i][0] for i in range(len(iso))],[iso[i][1] for i in range(len(iso))], c='#CC8800')
pyplot.plot([limitc[i][0] for i in range(len(limitc))],[limitc[i][1] for i in range(len(limitc))], lw=2, c='#000000')
pyplot.xlabel(r"$x$")
pyplot.ylabel(r"$y$")
pyplot.savefig("isospace.png", dpi=300)
pyplot.show()
```

![isospace](https://user-images.githubusercontent.com/29438168/160646962-6cdc0fdf-5008-4445-ae9b-bcd7c09363f0.png)
