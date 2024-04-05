def movex(x, odotx, acel, deltatn):

	xnext = x + deltatn*odotx + (acel*(deltatn**2.0))/2.0 	#x(t+dt)=x(t)+dt*v(t)+1/2*a(t)*dt^2

	return xnext