def movex(x, odotx, acel, deltatn):

	xnext = x + deltatn*odotx + (acel*(deltatn**2.0))/2.0

	return xnext