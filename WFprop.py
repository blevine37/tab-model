def stepWF (deltate, cr, ci, H):

	import numpy as np 

	odotcr = np.dot(H, ci)
	cr = cr + (deltate/2.0)*odotcr
	odotci = -1.0*np.dot(H, cr)
	ci = ci + deltate*odotci
	odotcr = np.dot(H, ci)
	cr = cr + (deltate/2.0)*odotcr

	return cr, ci