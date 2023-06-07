def eigsort (dimH,w,VR):
#Sorts eigenvalues and eigenvectors of a matrix from lowest energy to highest energy	
	
	import numpy as np
	
	vtempl = np.zeros((dimH))
	vtemph = np.zeros((dimH))
	
	k = 0
	while k < dimH:
		j = k+1
		while j < dimH:
			if (w[k] > w[j]):
				templ = w[j]
				temph = w[k]
				
				w[k] = templ
				w[j] = temph
				
				i = 0
				while i < dimH:
					vtempl[i] = V[i,j]
					vtemph[i] = V[i,k]
					
					V[i,k] = vtempl[i]
					V[i,j] = vtemph[i]
					i = i + 1
				pass
			pass
			j = j+1
		pass
		k = k+1
	pass

	return w,VR
