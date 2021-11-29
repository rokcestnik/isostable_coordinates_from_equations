from math import pi, floor, sqrt


###################################################
#####   PHASE RECONSTRUCTION FROM EQUATIONS   #####
###################################################


def one_step_integrator(state, ders, dt):
	"""RK4 integrates state with derivative for one step of dt
	
	:param state: state of the variables
	:param ders: derivative functions
	:param dt: time step
	:return: state after one integration step"""
	D = len(state)
	# 1
	k1 = [ders[i](state) for i in range(D)]
	# 2
	state2 = [state[i]+k1[i]*dt/2.0 for i in range(D)]
	k2 = [ders[i](state2) for i in range(D)]
	# 3
	state3 = [state[i]+k2[i]*dt/2.0 for i in range(D)] 
	k3 = [ders[i](state3) for i in range(D)]
	# 4
	state4 = [state[i]+k3[i]*dt for i in range(D)] 
	k4 = [ders[i](state4) for i in range(D)]
	# put togeather
	statef = [state[i] + (k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0*dt for i in range(D)]
	return statef


def integrate_period(state, ders, period, dt):
	"""integrates state with derivative for one period
	
	:param state: state of the variables
	:param ders: derivative functions
	:param period: period of integration
	:param dt: time step
	:return: state after one period"""
	for i in range(floor(period/dt)):
		state = one_step_integrator(state, ders, dt)
	state = one_step_integrator(state, ders, period-floor(period/dt)*dt)
	return state


def oscillator_period(ders, initial_state=None, warmup_time=1500.0, thr=0.0, dt=0.005):
	"""calculates the natural period of the oscillator from dynamical equations
	
	:param ders: a list of state variable derivatives
	:param initial_state: initial state (default None)
	:param warmup_time: time for relaxing to the stable orbit (default 1000)
	:param thr: threshold for determining period (default 0.0)
	:param dt: time step (default 0.005)
	:return: natural period"""
	# initial conditions
	if(initial_state==None):
		state = [0.5 for i in range(len(ders))]
	else:
		state = initial_state
	# warmup
	for i in range(round(warmup_time/dt)):
		state = one_step_integrator(state, ders, dt)
	# integration up to x = thr
	xh = state[0]
	while((state[0] > thr and xh < thr) == False):
		xh = state[0]
		state = one_step_integrator(state, ders, dt)
	# Henon trick
	dt_beggining = 1.0/ders[0](state)*(state[0]-thr)
	# spoil condition and go again to x = thr (still counting time)
	xh = state[0]
	time = 0
	while((state[0] > thr and xh < thr) == False):
		xh = state[0]
		state = one_step_integrator(state, ders, dt)
		time = time + dt
	# another Henon trick
	dt_end = 1.0/ders[0](state)*(state[0]-thr)
	return time + dt_beggining - dt_end


def oscillator_phase(state, ders, period, warmup_periods=5, thr=0.0, dt=0.005):
	"""calculates the asymptotic phase of the oscillator from dynamical equations
	
	:param state: state of the system
	:param ders: a list of state variable derivatives
	:param period: oscillator period
	:param warmup_periods: how many periods to wait for evaluating the asymptotic phase shift (default 5)
	:param thr: threshold determining zero phase (default 0.0)
	:param dt: time step (default 0.005)
	:return: asymptotic phase of state"""
	# integrate for some periods (relax to the limit cycle)
	for p in range(warmup_periods):
		state = integrate_period(state, ders, period, dt)
	# go to x = thr (counting time)
	time = 0
	xh = state[0]
	while((state[0] > thr and xh < thr) == False):
		xh = state[0]
		state = one_step_integrator(state, ders, dt)
		time = time + dt
	# Henon trick
	dt_end = 1.0/ders[0](state)*(state[0]-thr)
	return 2*pi*(1-(time-dt_end)/period)


def oscillator_PRC(ders, direction, initial_state=None, initial_warmup_periods=10, stimulation=0.05, warmup_periods=3, dph=0.1, thr=0.0, dt=0.005):
	"""calculates the phase response curve from dynamical equations
	
	:param ders: a list of state variable derivatives
	:param direction: direction in which the response is probed
	:param initial_state: initial state (default None)
	:param initial_warmup_periods: time for relaxing to the stable orbit (default 10)
	:param stimulation: strength of the stimulation (default 0.05)
	:param warmup_periods: how many periods to wait for evaluating the asymptotic phase shift (default 3)
	:param dph: phase resolution (default 0.1)
	:param thr: threshold for determining period (default 0.0)
	:param dt: time step (default 0.005)
	:return: the phase response curve"""
	period = oscillator_period(ders)
	PRC = [[dph*i for i in range(floor(2*pi/dph))],[0 for i in range(floor(2*pi/dph))]] # PRC list
	# initial conditions
	if(initial_state==None):
		state = [0.5 for i in range(len(ders))]
	else:
		state = initial_state
	# normalize the direction
	norm = sqrt(sum(direction[i]**2 for i in range(len(direction))))
	direction = [direction[i]/norm for i in range(len(direction))]
	# warmup
	for p in range(initial_warmup_periods):
		state = integrate_period(state, ders, period, dt)
	# integration up to x = thr (zero phase)
	xh = state[0]
	while((state[0] > thr and xh < thr) == False):
		xh = state[0]
		state = one_step_integrator(state, ders, dt)
	# Henon trick
	dt_beggining = 1.0/ders[0](state)*(state[0]-thr)
	# get limit cycle states
	limitc = []
	for ph in PRC[0]:
		limitc.append(state)
		state = integrate_period(state, ders, dph/(2*pi)*period, dt)
	# shift the states and evaluate the phase difference
	for s in range(len(PRC[0])):
		state_stim = [limitc[s][i]+direction[i]*stimulation for i in range(len(state))] # shift the state
		PRC[1][s] = (oscillator_phase(state_stim, ders, period)-PRC[0][s])/stimulation
	return PRC

