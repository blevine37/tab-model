# Calculate Ehrenfest energy
def effE (dimH, ct, x1, x2, w1, w2, c, delta):
	
	import numpy as np 

	ctbra = np.transpose(np.conjugate(ct))

	H = buildH(dimH, x1, x2, w1, w2, c, delta)

	cnorm = np.dot(ctbra, ct)

	energy = np.dot(ctbra, np.dot(H, ct))/cnorm

	energy = energy.real

	return energy
