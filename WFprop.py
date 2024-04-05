import numpy as np

def stepWF (deltate, cr, ci, H):

	odotcr = np.dot(H, ci)
	cr = cr + (deltate/2.0)*odotcr	#cr(t+dt/2)=cr(t)+dt/2*H*ci(t)
	odotci = -1.0*np.dot(H, cr)
	ci = ci + deltate*odotci		#ci(t+dt)=ci(t)-dt*H*cr(t+dt/2)
	odotcr = np.dot(H, ci)
	cr = cr + (deltate/2.0)*odotcr	#cr(t+dt)=cr(t+dt/2)+dt/2*H*ci(t+dt)

	return cr, ci