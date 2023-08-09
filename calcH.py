def buildH(dimH, x1, A, B, C):
#Build the Hamiltonian for at position (x1,x2)
	
	import numpy as np

	#initialize Hamiltonian matrix
	H = np.zeros((dimH,dimH))

    ########################################################
	H[0][0] = -1.0*w1*x1   #diabatic potential
	
	i = 1

	while i < dimH:
		H[i][i] = w2*x1 - (i-1)*delta  #diabatic potential
		H[i][0] = c*x2		#diabatic coupling
		H[0][i] = H[i][0]
		i = i + 1
	pass
		
	i = 5

	while i < dimH:
		H[i][i] = H[i][i] - 0.08
		i = i + 1
	pass

	#######################################################

	return H
