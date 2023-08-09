# Calculate Ehrenfest energy
def effE (dimH, ct, x1, A, B, C):
	
	import numpy as np 

	ctbra = np.transpose(np.conjugate(ct))

	H = buildH(dimH, x1, A, B, C)

	cnorm = np.dot(ctbra, ct)

	energy = np.dot(ctbra, np.dot(H, ct))/cnorm

	energy = energy.real

	return energy
