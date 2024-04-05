#Build the Hamiltonian for at position (x1)
import numpy as np

def buildH(dimH, x1, A, B, C):

	#initialize Hamiltonian matrix
	H = np.zeros((dimH,dimH))

	H[0,0] = A   #diabatic potential
	H[1,1] = -A
	
	if x1 < 0:
		H[0,1] = B * np.exp(C * x1)
		H[1,0] = H[0,1]
	else:
		H[0,1] = B * (2 - np.exp(-C * x1))
		H[1,0] = H[0,1]

	return H
