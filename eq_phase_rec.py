from math import pi, floor, sqrt, exp, log


###################################################
#####   PHASE RECONSTRUCTION FROM EQUATIONS   #####
###################################################


def distance(state1, state2):
	"""the euclidian distance between two states"""
	return sqrt(sum((state1[i]-state2[i])**2 for i in range(len(state1))))


def set_initial_state(initial_state, ders):
	if(initial_state==None):
		return [0.5 for i in range(len(ders))]
	return initial_state


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


def integrate_period(state, ders, period, dt=0.005):
	"""integrates state with derivative for one period
	
	:param state: state of the variables
	:param ders: derivative functions
	:param period: period of integration
	:param dt: time step (default 0.005)
	:return: state after one period"""
	for i in range(floor(period/dt)):
		state = one_step_integrator(state, ders, dt)
	state = one_step_integrator(state, ders, period-floor(period/dt)*dt)
	return state


def integrate_up_to_thr(state, ders, thr, dt=0.005):
	"""integrates state with derivative for up to x = thr
	
	:param state: state of the variables
	:param ders: derivative functions
	:param thr: threshold
	:param dt: time step (default 0.005)
	:return: state at threshold crossing"""
	xh = state[0]
	while((state[0] > thr and xh < thr) == False):
		xh = state[0]
		state = one_step_integrator(state, ders, dt)
	# Henon trick
	dt_over = 1.0/ders[0](state)*(state[0]-thr)
	state = one_step_integrator(state, ders, -dt_over)
	return state


def time_up_to_thr(state, ders, thr, dt=0.005):
	"""counts time up to x = thr
	
	:param state: state of the variables
	:param ders: derivative functions
	:param thr: threshold
	:param dt: time step (default 0.005)
	:return: time until threshold crossing"""
	time = 0
	xh = state[0]
	while((state[0] > thr and xh < thr) == False):
		xh = state[0]
		state = one_step_integrator(state, ders, dt)
		time += dt
	# Henon trick
	dt_over = 1.0/ders[0](state)*(state[0]-thr)
	return time-dt_over


def oscillator_period(ders, initial_state=None, warmup_time=1500.0, thr=0.0, dt=0.005):
	"""calculates the natural period of the oscillator from dynamical equations
	
	:param ders: a list of state variable derivatives
	:param initial_state: initial state (default None)
	:param warmup_time: time for relaxing to the stable orbit (default 1000)
	:param thr: threshold for determining period (default 0.0)
	:param dt: time step (default 0.005)
	:return: natural period"""
	# initial conditions
	state = set_initial_state(initial_state, ders)
	# warmup
	state = integrate_period(state, ders, warmup_time, dt)
	# integration up to x = thr
	state = integrate_up_to_thr(state, ders, thr, dt)
	return time_up_to_thr(state, ders, thr, dt)


def sample_limit_cycle(ders, sampling, period=None, initial_state=None, warmup_time=1500.0, thr=0.0, dt=0.005):
	"""samples the limit cycle
	
	:param ders: a list of state variable derivatives
	:param sampling: the number of samples
	:param period: oscillator period (default None - gets calculated)
	:param initial_state: initial state (default None)
	:param warmup_time: time for relaxing to the stable orbit (default 1500)
	:param thr: threshold for determining period (default 0.0)
	:param dt: time step (default 0.005)
	:return: sampled limit cycle"""
	# initial conditions
	state = set_initial_state(initial_state, ders)
	# warmup
	state = integrate_period(state, ders, warmup_time, dt)
	# if period == None, then calculate it
	if(period == None): 
		state = integrate_up_to_thr(state, ders, thr, dt)
		period = time_up_to_thr(state, ders, thr, dt)
	# integration up to x = thr
	state = integrate_up_to_thr(state, ders, thr, dt)
	# now sample the limit cycle
	limit_cycle = []
	for i in range(sampling):
		limit_cycle.append(state)
		state = integrate_period(state, ders, period/sampling, dt)
	limit_cycle.append(limit_cycle[0]) # close the curve
	return limit_cycle


def oscillator_floquet(ders, period, initial_state=None, warmup_time=1500.0, shift=0.01, sampling=25, thr=0.0, dt=0.005):
	"""calculates the floquet exponent of the oscillator from dynamical equations
	
	:param ders: a list of state variable derivatives
	:param period: oscillator period
	:param initial_state: initial state (default None)
	:param warmup_time: time for relaxing to the stable orbit (default 1500)
	:param shift: proportional shift of the points from the limit cycle (default 0.01)
	:param sampling: limit cycle sampling (default 25)
	:param thr: threshold for determining period (default 0.0)
	:param dt: time step (default 0.005)
	:return: floquet exponent"""
	# sample limit cycle
	limitc = sample_limit_cycle(ders, sampling, period, initial_state, warmup_time, thr, dt)
	# average
	lc_avg = [sum([limitc[i][s] for i in range(len(limitc))])/len(limitc) for s in range(len(ders))]
	# measure the exponent
	exponent_sum = 0
	exponent_measures = 0
	for i in range(len(limitc)):
		# state out
		state_out = [lc_avg[s]+(limitc[i][s]-lc_avg[s])*(1+shift) for s in range(len(ders))]
		state_1 = integrate_period(state_out, ders, period, dt)
		state_2 = integrate_period(state_1, ders, period, dt)
		d1 = distance(state_out, state_1)
		d2 = distance(state_1, state_2)
		exponent_sum += log(d1/d2)/period
		exponent_measures += 1
		# state in
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
	state = integrate_period(state, ders, warmup_periods*period, dt)
	# go to x = thr (counting time)
	time = time_up_to_thr(state, ders, thr, dt)
	return 2*pi*(1-time/period)


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


def oscillator_amplitude(state, ders, period, floquet, zero_phase_lc, phase_warmup_periods=5, thr=0.0, dt=0.005):
	"""calculates the isostable amplitude of the oscillator from dynamical equations
	
	:param state: state of the system
	:param ders: a list of state variable derivatives
	:param period: oscillator period
	:param floquet: floquet exponent
	:param zero_phase_lc: zero phase limit cycle state
	:param phase_warmup_periods: how many periods to wait for evaluating the asymptotic phase shift (default 5)
	:param thr: threshold determining zero phase (default 0.0)
	:param dt: time step (default 0.005)
	:return: isostable amplitude of state"""
	# get phase
	phase = oscillator_phase(state, ders, period, phase_warmup_periods)
	# calculate time to evolve to zero isochron
	time = (1-phase/(2*pi))*period
	# evolve to 0 isochron
	state = integrate_period(state, ders, time, dt)
	# amplitude sign
	if(inside_limit_cycle(state, ders, period)):
		sign = -1
	else:
		sign = 1
	return 0.2*sign*distance(state,zero_phase_lc)*exp(floquet*time) # multiplied with an arbitrary constant 0.2


def oscillator_PRC(ders, direction, period, initial_state=None, initial_warmup_periods=10, stimulation=0.05, warmup_periods=3, sampling=100, thr=0.0, dt=0.005):
	"""calculates the phase response curve from dynamical equations
	
	:param ders: a list of state variable derivatives
	:param direction: direction in which the response is probed
	:param period: oscillator period
	:param initial_state: initial state (default None)
	:param initial_warmup_periods: time for relaxing to the stable orbit (default 10)
	:param stimulation: strength of the stimulation (default 0.05)
	:param warmup_periods: how many periods to wait for evaluating the asymptotic phase shift (default 3)
	:param sampling: phase sampling (default 100)
	:param thr: threshold for determining period (default 0.0)
	:param dt: time step (default 0.005)
	:return: the phase response curve"""
	PRC = [[2*pi/sampling*i for i in range(sampling)],[0 for i in range(sampling)]] # PRC list
	# normalize the direction
	norm = sqrt(sum(direction[i]**2 for i in range(len(direction))))
	direction = [direction[i]/norm for i in range(len(direction))]
	# sample limit cycle
	limitc = sample_limit_cycle(ders, sampling, period, initial_state, initial_warmup_periods*period, thr, dt)
	# shift the states and evaluate the phase difference
	for s in range(len(PRC[0])):
		state_stim = [limitc[s][i]+direction[i]*stimulation for i in range(len(limitc[s]))] # shift the state
		PRC[1][s] = (oscillator_phase(state_stim, ders, period)-PRC[0][s])/stimulation
	return PRC


def oscillator_ARC(ders, direction, period, floquet, initial_state=None, initial_warmup_periods=15, stimulation=0.05, sampling=100, thr=0.0, dt=0.005):
	"""calculates the amplitude response curve from dynamical equations
	
	:param ders: a list of state variable derivatives
	:param direction: direction in which the response is probed
	:param period: oscillator period
	:param floquet: floquet exponent
	:param initial_state: initial state (default None)
	:param initial_warmup_periods: time for relaxing to the stable orbit (default 15)
	:param stimulation: strength of the stimulation (default 0.05)
	:param sampling: phase sampling (default 100)
	:param thr: threshold for determining period (default 0.0)
	:param dt: time step (default 0.005)
	:return: the amplitude response curve"""
	ARC = [[2*pi/sampling*i for i in range(sampling)],[0 for i in range(sampling)]] # ARC list
	# normalize the direction
	norm = sqrt(sum(direction[i]**2 for i in range(len(direction))))
	direction = [direction[i]/norm for i in range(len(direction))]
	# sample limit cycle
	limitc = sample_limit_cycle(ders, sampling, period, initial_state, initial_warmup_periods*period, thr, dt)
	# shift the states and evaluate the amplitude difference
	for s in range(len(ARC[0])):
		state_stim = [limitc[s][i]+direction[i]*stimulation for i in range(len(limitc[s]))] # shift the state
		ARC[1][s] = oscillator_amplitude(state_stim, ders, period, floquet, limitc[0], phase_warmup_periods=1)/stimulation
	return ARC


def sample_local_isostable(ders, sampling, period, floquet, sign, shift=0.005, initial_state=None, warmup_periods=10, thr=0.0, dt=0.005):
	"""samples the isostable close to the limit cycle
	
	:param ders: a list of state variable derivatives
	:param sampling: the number of samples (typically should be 200 or more to work properly)
	:param period: oscillator period
	:param floquet: floquet exponent
	:param sign: whether the isostable is inside or outside the limit cycle (1 means out, -1 in)
	:param shift: proportional shift of the points from the limit cycle (default 0.005)
	:param initial_state: initial state (default None)
	:param warmup_periods: number of periods to relax to the limit cycle (default 10)
	:param thr: threshold for determining period (default 0.0)
	:param dt: time step (default 0.005)
	:return: sampled isostable close to limit cycle"""
	# sample limit cycle
	limitc = sample_limit_cycle(ders, sampling, period, initial_state, warmup_periods*period, thr, dt)
	# average
	lc_avg = [sum([limitc[i][s] for i in range(len(limitc))])/len(limitc) for s in range(len(ders))]
	# shift the first point
	state_sh = [lc_avg[s]+(limitc[0][s]-lc_avg[s])*(1+sign*shift) for s in range(len(ders))]
	# estimate the phase and integrate the phase appropriately so its aligned with the point on the limit cycle
	phase = oscillator_phase(state_sh, ders, period)
	if(phase > pi):
		phase = phase-2*pi
	state_sh = integrate_period(state_sh, ders, -phase/(2*pi)*period, -phase/abs(phase)*dt)
	# estimate the amplitude for reference throughout the sampling
	ampl0 = oscillator_amplitude(state_sh, ders, period, floquet, limitc[0])
	# now integrate the state, each time adjusting for amplitude decay
	iso = []
	for s in range(sampling):
		iso.append(state_sh)
		state_sh = integrate_period(state_sh, ders, period/sampling, dt)
		# adjust for amplitude decay
		factor = exp(floquet*period/sampling)
		state_sh = [limitc[s][i]+(state_sh[i]-limitc[s][i])*factor for i in range(len(ders))]
		# additionally adjust for any small phase shift due to higher order effects (if it were infinitesimal this was not needed)
		phase_diff = 2*pi*(s+1)/sampling-oscillator_phase(state_sh, ders, period)
		state_sh = integrate_period(state_sh, ders, phase_diff/(2*pi)*period, phase_diff/abs(phase_diff)*dt)
		# and again adjust amplitude by measuring it
		ampl = oscillator_amplitude(state_sh, ders, period, floquet, limitc[0])
		state_sh = [limitc[s][i]+(state_sh[i]-limitc[s][i])*(ampl0/ampl) for i in range(len(ders))]
	iso.append(iso[0]) # close the curve
	return iso


def oscillator_isostables(ders, local_iso_in, local_iso_out, period, floquet, number_of_periods=5, isostables_per_period=2, dt=0.005):
	"""estimates the isostables both inside and outside of the limit cycle
	
	:param ders: a list of state variable derivatives
	:param local_iso_in: sampling of a local isostable inside the limit cycle
	:param local_iso_out: sampling of a local isostable outside the limit cycle
	:param period: oscillator period
	:param floquet: floquet exponent
	:param number_of_periods: how many periods to evolve (default 5)
	:param isostables_per_period: how many isostables to plot per period (default 2)
	:param dt: time step (default 0.005)
	:return: isostables of the limit cycle"""
	# initialize running isostables
	iso_in = local_iso_in.copy()
	iso_out = local_iso_out.copy()
	# declare the array of isostables
	isos = [iso_in, iso_out]
	# evolve the isostables and save them
	evolve_time = period/isostables_per_period
	for e in range(int(number_of_periods*isostables_per_period)):
		for i in range(len(iso_in)):
			iso_in[i] = integrate_period(iso_in[i], ders, -evolve_time, -dt)
			iso_out[i] = integrate_period(iso_out[i], ders, -evolve_time, -dt)
		isos.append(iso_in.copy())
		isos.append(iso_out.copy())
	return isos
