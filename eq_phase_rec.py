from math import pi, floor, sqrt, exp, log


###################################################
#####   PHASE RECONSTRUCTION FROM EQUATIONS   #####
###################################################


def distance(state1, state2):
	"""the euclidian distance between two states"""
	return sqrt(sum((state1[i]-state2[i])**2 for i in range(len(state1))))


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


def oscillator_floquet(ders, period, initial_state=None, warmup_time=1500.0, shift=0.05, dph=0.1, thr=0.0, dt=0.005):
	"""calculates the floquet exponent of the oscillator from dynamical equations
	
	:param ders: a list of state variable derivatives
	:param period: oscillator period
	:param initial_state: initial state (default None)
	:param warmup_time: time for relaxing to the stable orbit (default 1000)
	:param shift: proportional shift of the points from the limit cycle (default 0.01)
	:param dph: phase resolution (default 0.1)
	:param thr: threshold for determining period (default 0.0)
	:param dt: time step (default 0.005)
	:return: floquet exponent"""
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
	dt_over = 1.0/ders[0](state)*(state[0]-thr)
	state = one_step_integrator(state, ders, -dt_over)
	# get limit cycle states
	limitc = []
	for ph in [dph*i for i in range(floor(2*pi/dph))]:
		limitc.append(state)
		state = integrate_period(state, ders, dph/(2*pi)*period, dt)
	# average
	lc_avg = [sum([limitc[i][s] for i in range(len(limitc))])/len(limitc) for s in range(len(ders))]
	# measure the exponent
	exponent_sum = 0
	exponent_measures = 0
	for i in range(len(limitc)):
		# state in
		state_out = [lc_avg[s]+(limitc[i][s]-lc_avg[s])*(1+shift) for s in range(len(ders))]
		state_1 = integrate_period(state_out, ders, period, dt)
		state_2 = integrate_period(state_1, ders, period, dt)
		d1 = distance(state_out, state_1)
		d2 = distance(state_1, state_2)
		exponent_sum += log(d1/d2)/period
		exponent_measures += 1
		# state out
		state_in = [lc_avg[s]+(limitc[i][s]-lc_avg[s])*(1-shift) for s in range(len(ders))]
		state_1 = integrate_period(state_in, ders, period, dt)
		state_2 = integrate_period(state_1, ders, period, dt)
		d1 = distance(state_in, state_1)
		d2 = distance(state_1, state_2)
		exponent_sum += log(d1/d2)/period
		exponent_measures += 1
	return exponent_sum/exponent_measures


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


def inside_limit_cycle(state, ders, period, dt=0.005):
	"""determines whether the state is inside or outside the limit cycle based on the standard deviation of the trajectories forward and backward in time
	
	:param state: state of the system
	:param ders: a list of state variable derivatives
	:param period: oscillator period
	:param dt: time step (default 0.005)
	:return: true if inside limit cycle, false otherwise"""
	state0 = state.copy()
	# forward trajectory
	forward = []
	for i in range(floor(period/dt)):
		state = one_step_integrator(state, ders, dt)
		forward.append(state)
	lc_avg_f = [sum([forward[i][s] for i in range(len(forward))])/len(forward) for s in range(len(ders))]
	lc_std_f = [sum([(forward[i][s]-lc_avg_f[s])**2 for i in range(len(forward))])/len(forward) for s in range(len(ders))]
	std_f = sum(lc_std_f)
	# backward trajectory
	state = state0
	backward = []
	for i in range(floor(period/dt)):
		state = one_step_integrator(state, ders, -dt)
		backward.append(state)
	lc_avg_b = [sum([backward[i][s] for i in range(len(backward))])/len(backward) for s in range(len(ders))]
	lc_std_b = [sum([(backward[i][s]-lc_avg_b[s])**2 for i in range(len(backward))])/len(backward) for s in range(len(ders))]
	std_b = sum(lc_std_b)
	if(std_f > std_b):
		return True
	return False


def oscillator_amplitude(state, ders, period, floquet, zero_phase_lc, thr=0.0, dt=0.005):
	"""calculates the isostable amplitude of the oscillator from dynamical equations
	
	:param state: state of the system
	:param ders: a list of state variable derivatives
	:param period: oscillator period
	:param floquet: floquet exponent
	:param zero_phase_lc: zero phase limit cycle state
	:param thr: threshold determining zero phase (default 0.0)
	:param dt: time step (default 0.005)
	:return: isostable amplitude of state"""
	phase = oscillator_phase(state, ders, period)
	# evolve to 0 isochron
	state = integrate_period(state, ders, (1-phase/(2*pi))*period, dt)
	# amplitude sign
	if(inside_limit_cycle(state, ders, period)):
		sign = -1
	else:
		sign = 1
	return 0.2*sign*distance(state,zero_phase_lc)*exp(floquet*(1-phase/(2*pi))*period) # multiplied with an arbitrary constant 0.2


def oscillator_PRC(ders, direction, period, initial_state=None, initial_warmup_periods=10, stimulation=0.05, warmup_periods=3, dph=0.1, thr=0.0, dt=0.005):
	"""calculates the phase response curve from dynamical equations
	
	:param ders: a list of state variable derivatives
	:param direction: direction in which the response is probed
	:param period: oscillator period
	:param initial_state: initial state (default None)
	:param initial_warmup_periods: time for relaxing to the stable orbit (default 10)
	:param stimulation: strength of the stimulation (default 0.05)
	:param warmup_periods: how many periods to wait for evaluating the asymptotic phase shift (default 3)
	:param dph: phase resolution (default 0.1)
	:param thr: threshold for determining period (default 0.0)
	:param dt: time step (default 0.005)
	:return: the phase response curve"""
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
	dt_over = 1.0/ders[0](state)*(state[0]-thr)
	state = one_step_integrator(state, ders, -dt_over)
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


def oscillator_ARC(ders, direction, period, floquet, initial_state=None, initial_warmup_periods=15, stimulation=0.05, dph=0.1, thr=0.0, dt=0.005):
	"""calculates the amplitude response curve from dynamical equations
	
	:param ders: a list of state variable derivatives
	:param direction: direction in which the response is probed
	:param period: oscillator period
	:param floquet: floquet exponent
	:param initial_state: initial state (default None)
	:param initial_warmup_periods: time for relaxing to the stable orbit (default 15)
	:param stimulation: strength of the stimulation (default 0.05)
	:param dph: phase resolution (default 0.1)
	:param thr: threshold for determining period (default 0.0)
	:param dt: time step (default 0.005)
	:return: the phase response curve"""
	ARC = [[dph*i for i in range(floor(2*pi/dph))],[0 for i in range(floor(2*pi/dph))]] # ARC list
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
	dt_over = 1.0/ders[0](state)*(state[0]-thr)
	state = one_step_integrator(state, ders, -dt_over)
	zero_phase_lc = state.copy()
	# get limit cycle states
	limitc = []
	for ph in ARC[0]:
		limitc.append(state)
		state = integrate_period(state, ders, dph/(2*pi)*period, dt)
	# shift the states and evaluate the amplitude difference
	for s in range(len(ARC[0])):
		state_stim = [limitc[s][i]+direction[i]*stimulation for i in range(len(state))] # shift the state
		ARC[1][s] = oscillator_amplitude(state_stim, ders, period, floquet, zero_phase_lc)/stimulation
	return ARC

