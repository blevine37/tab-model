def rstates(rank,eseg,bcore,nzthresh):   # comments
	"""I really should define this better"""

	import numpy as np
	import sys
	import random

	icstates = [] 	# list of unique states that have been found to have 0 elements
			# in eseg in line with the current bcore element

	p = bcore[0]
	q = bcore[1]
	sum = []

	dimH = rank

	# Though horribly inefficient, repeating the original algorithm to avoid errors, then appending this
	        # loop over the row (changes column index)
	rowstates = []
	j = 0
	while j < dimH:
		if (j != p):
			if (eseg[p][j] <= nzthresh):
				rowstates.append(j)
			pass
		pass
		j = j + 1
	pass

        # loop over the column (changes row index)
	colstates = []
	i = 0
	while i < dimH:
		if (i != q):
			if (eseg[i][q] <= nzthresh):
				colstates.append(i)
			pass
		pass
		i = i + 1
	pass

	sum = rowstates + colstates


	i = 0
	while i < dimH -1:
		j = i + 1
		while j < dimH:
			if (eseg[i][j] <= nzthresh):
				if (i == p or i == q):
					if (j != p and j != q):
						sum.append(j)
					pass
				else:
					if (j == p or j == q):
						sum.append(i)
					else:
						k = 0
						flag = 'none'
						check2 = len(sum)
						while k < check2:
							if (i == sum[k]):
								flag = 'i'
							pass
							if (j == sum[k]):
								flag = 'j'
							k = k + 1
						pass
						if (flag == 'none'):	
							check = random.random()
							if (check <= 0.5):
								sum.append(i)
							else:
								sum.append(j)
							pass
						elif (flag == 'i'):
							sum.append(i)
						elif (flag == 'j'):
							sum.append(j)
						pass
					pass
				pass
			pass
			j = j + 1
		pass
		i = i + 1
	pass

	temp3 = sum.sort()
	usum = []

	usum.append(sum[0])
	i = 1
	while i < len(sum):
		j = 0
		add = 0
		while j < len(usum):
			if (usum[j] == sum[i]):
				add = 1
			pass
			j = j + 1
		pass
		if (add == 0):
			usum.append(sum[i])
		pass
		i = i + 1
	pass

#	print 'sum in sremove'
#	print sum
#	print 'usum'
#	print usum
	icstates = usum[:]

#	print 'icstates in sremove'
#	print icstates

	return icstates

