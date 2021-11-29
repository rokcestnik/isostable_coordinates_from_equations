import random


def one_step_integrator(state, ders, ps, dt):
	"""RK4 integrates state with derivative and input for one step of dt
	
	:param state: state of the variables
	:param ders: derivative functions
	:param ps: perturbations [p1,p2,...] (default 0)
	:param dt: time step
	:return: state after one integration step"""
	D = len(state)
	# 1
	k1 = [ders[i](state,ps) for i in range(D)]
	# 2
	state2 = [state[i]+k1[i]*dt/2.0 for i in range(D)]
	k2 = [ders[i](state2,ps) for i in range(D)]
	# 3
	state3 = [state[i]+k2[i]*dt/2.0 for i in range(D)] 
	k3 = [ders[i](state3,ps) for i in range(D)]
	# 4
	state4 = [state[i]+k3[i]*dt for i in range(D)] 
	k4 = [ders[i](state4,ps) for i in range(D)]
	# put togeather
	statef = [state[i] + (k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0*dt for i in range(D)]
	return statef


def generate_signal(ders, n, sampling, initial_state=None, number_of_variables=1, number_of_perturbations=1, warmup_time=1000.0, tau=3.0, eps=0.5, dt=0.01):
	"""generates signal for the oscillator driven by correlated noise from dynamical equations
	
	:param ders: a list of state variable derivatives
	:param n: length of time series
	:param sampling: sampling rate
	:param initial_state: initial state
	:param number_of_variables: number of variables returned, not including the input (default 1)
	:param number_of_perturbations: number of perturbations (default 1)
	:param warmup_time: time for relaxing to the stationary regime (default 1000)
	:param tau: noise correlation time (default 3.0)
	:param eps: noise strength (default 0.5)
	:param dt: time step (default 0.01)
	:return: time series of the signal and driving noise"""
	# initial conditions
	if(initial_state==None):
		state = [random.gauss(0,1) for i in range(len(ders))]
	else:
		state = initial_state
	res_s = [[] for i in range(number_of_variables)] # resulting signal
	res_p = [[] for j in range(number_of_perturbations)] # resulting perturbations
	ps = [0 for p in range(number_of_perturbations)]
	# warmup
	for i in range(round(warmup_time/dt)):
		for p in range(number_of_perturbations):
			ps[p] = ps[p] - (ps[p]/tau - eps*sqrt(2/tau)*random.gauss(0,1)/sqrt(dt))*dt 
		state = one_step_integrator(state, ders, ps, dt)
	# real integration
	for i in range(n*sampling):
		for p in range(number_of_perturbations):
			ps[p] = ps[p] - (ps[p]/tau - eps*sqrt(2/tau)*random.gauss(0,1)/sqrt(dt))*dt 
		state = one_step_integrator(state, ders, ps, dt)
		for c in range(number_of_variables):
			res_s[c].append(state[c])
		for p in range(number_of_perturbations):
			res_p[p].append(ps[p])
	# sampling
	res_s = [res_s[i][::sampling] for i in range(number_of_variables)]
	res_p = [res_p[i][::sampling] for i in range(number_of_perturbations)]
	return res_s, res_p
