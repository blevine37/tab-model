# Calculate dH/dR
def dHcalc(dimH,ndof,x, w1, w2, c, delta):
	import numpy as np
	dH=np.zeros((ndof, dimH, dimH))

	dH[0, 0, 0] = -1.0*w1
	
	i = 1
	while i < dimH:
		dH[0, i, i] = w2
		dH[1, i, 0] = c
		dH[1, 0, i] = c
		i = i + 1
	pass

	return dH
