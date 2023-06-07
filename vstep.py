def vcalc(odotx, mfdF, mfdFprev, deltatn, pmass):

	vnext = odotx + ((mfdF+mfdFprev)/pmass)*deltatn/2.0

	return vnext