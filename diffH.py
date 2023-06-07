# Calculate dH/dR
def dHcalc(dimH,x1, x2, w1, w2, c, delta):
	import numpy as np
	
	dH1 = np.zeros((dimH, dimH))
	dH2 = np.zeros((dimH, dimH))

	dH1[0][0] = -1.0*w1
	
	i = 1
	while i < dimH:
		dH1[i][i] = w2
		dH2[i][0] = c
		dH2[0][i] = c
		i = i + 1
	pass

	return dH1, dH2
