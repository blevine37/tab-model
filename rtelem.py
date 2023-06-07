def minelem (eseg,rank,nzthresh): # For a given sub-block of the density matrix
	"""Returns the element with fastest coherence loss"""

	# Standard library imports ====================================
	import numpy as np
	import sys

	leave = 'no'

	dimH = rank

	i = 0

	p = 0 	# first coordinate of minimum element of eseg
	q = 0 	# second coordinate of min. element

	check = 1.1 	# placeholder value for scanning for minimum element of eseg
	while i < dimH - 1:
		j = i + 1
		while j < dimH:
			if (abs(eseg[i][j]) <= nzthresh):
				pass
			else:
				if (check > eseg[i][j]):
					check = eseg[i][j]
					p = i
					q = j
				pass
			pass
			j = j + 1
		pass
		i = i + 1
	pass

	if (p == q):
		leave = 'yes'
	pass

	return p, q, leave


