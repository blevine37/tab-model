def gcollapse(dimH,deltatn,aforce1,poparray,dcp1,nzthresh,errortol,npthresh,pehrptol,odotrho,tolodotrho,nta,dtw,zpop,dgscale): # Determines which coherent sub-block 
	"""of the electronic density matrix the WF collapses into"""

	# Standard library imports =====================================
	import numpy as np
	import math
	import sys
	import random
	from scipy.optimize import lsq_linear

	# Custom function imports ======================================
	from rtelem import minelem
	from sremove import rstates


	# Rules ========================================================
	# i, j, k, l, m, and n are all reserved for integer incrementing


	# General setup ================================================
	npop = np.zeros((dimH)) 	# Stores output electronic populations

	invtau = np.zeros((dimH,dimH)) 	# Stores inverse of state-pairwise decoherence times

	i = 0
	while i < dimH:
		j = i + 1
		while j < dimH:
			invtau[i][j] = ((aforce1[i]-aforce1[j])**(2.0)/(8.0*dcp1))**(0.50)
			invtau[j][i] = invtau[i][j]
			j = j + 1
		i = i + 1
	pass

	# computing the weight function

#	w = np.zeros((dimH,nta+1))
#	zpoptime = np.zeros((dimH))

#	i = 0
#	while i < dimH:
#		if (odotrho[i] >= tolodotrho):
#			zpoptime[i] = poparray[i]/odotrho[i]
#			k = 0
#			sum = 0.0
#			while k < nta:
#				if ((sum + odotrho[i]*dtw) <= poparray[i]):
#					w[i][k] = odotrho[i]*dtw/poparray[i]
#				else: 
#					if (sum < poparray[i]):
#						w[i][k] = (poparray[i] - sum)/poparray[i]
#					else:
#						w[i][k] = 0.0
#					pass
#				pass
#				sum = sum + odotrho[i]*dtw
#				k = k + 1
#			pass
#			if (sum < poparray[i]):
#				w[i][-1] = (poparray[i] - sum)/poparray[i]
#			pass
#		i = i + 1
#	pass				

#	Not super efficient, could put sooner, but here is the rank reduction 
#	to working with only populated electronic states -> add a tolerance to pass
#	from namd-main.py

	vstates = []

	i = 0
	while i < dimH:
		if (poparray[i] >= zpop):
			vstates.append(i)
		i = i + 1
	pass

	rank = len(vstates)

#	print 'odotrho'
#	print odotrho
#	print 'w'
#	print w

#	sys.exit()
	# constructing vectorized target density matrix
	# This target is rho**(d) in the original python manuscript
	# However it is organized into a column vector of the diagonal and
	# lower triangular elements so it can be used as the target in a 
	# library linear least squares optimization

	eseg = np.zeros((rank,rank))
	nblock = 0

	vtarget = []
	i = 0
	while i < rank:
		j = i
		while j < rank:
# lots of logic in order to determine how integration will work based on current rules
			if (i == j):
				velem = dgscale*poparray[vstates[i]]
				eseg[i][i] = 1.0
#			elif (odotrho[vstates[i]] < tolodotrho and odotrho[vstates[j]] < tolodotrho):
#		                velem = ((poparray[vstates[i]]*poparray[vstates[j]])**(0.5))*(1.0+((math.exp(-1.0*((deltatn+dtw*nta)**(2.0))*invtau[vstates[i]][vstates[j]]*invtau[vstates[i]][vstates[j]])-math.exp(-1.0*((dtw*nta)**(2.0))*invtau[vstates[i]][vstates[j]]*invtau[vstates[i]][vstates[j]]))/math.exp(-1.0*((dtw*nta)**(2.0))*invtau[vstates[i]][vstates[j]]*invtau[vstates[i]][vstates[j]])))
#			elif (odotrho[vstates[i]] >= tolodotrho and odotrho[vstates[j]] < tolodotrho):
#				k = 0
#				nsum = 0.0
#				dsum = 0.0
#				while k <= nta:
#					nsum = nsum + w[vstates[i]][k]*(math.exp(-1.0*((deltatn+dtw*k)**(2.0))*invtau[vstates[i]][vstates[j]]*invtau[vstates[i]][vstates[j]])-math.exp(-1.0*((dtw*k)**(2.0))*invtau[vstates[i]][vstates[j]]*invtau[vstates[i]][vstates[j]]))
#					dsum = dsum + w[vstates[i]][k]*math.exp(-1.0*((dtw*k)**(2.0))*invtau[vstates[i]][vstates[j]]*invtau[vstates[i]][vstates[j]])
#					k = k + 1
#				velem = ((poparray[vstates[i]]*poparray[vstates[j]])**(0.5))*(1.0+nsum/dsum)
#			elif (odotrho[vstates[j]] >= tolodotrho and odotrho[vstates[i]] < tolodotrho):
 #                               k = 0
  #                              nsum = 0.0
   #                             dsum = 0.0
    #                            while k <= nta:
     #                                   nsum = nsum + w[vstates[j]][k]*(math.exp(-1.0*((deltatn+dtw*k)**(2.0))*invtau[vstates[i]][vstates[j]]*invtau[vstates[i]][vstates[j]])-math.exp(-1.0*((dtw*k)**(2.0))*invtau[vstates[i]][vstates[j]]*invtau[vstates[i]][vstates[j]]))
      #                                  dsum = dsum + w[vstates[j]][k]*math.exp(-1.0*((dtw*k)**(2.0))*invtau[vstates[i]][vstates[j]]*invtau[vstates[i]][vstates[j]])
       #                                 k = k + 1
        #                        velem = ((poparray[vstates[i]]*poparray[vstates[j]])**(0.5))*(1.0+nsum/dsum)
	#		elif (odotrho[vstates[i]] >= tolodotrho and odotrho[vstates[j]] >= tolodotrho): # clean up notation and just use derivatives
	#			if (zpoptime[vstates[i]] > zpoptime[vstates[j]]):
	 #                               k = 0
	  #                              nsum = 0.0
	   #                             dsum = 0.0
	    #                            while k <= nta:
	     #                                   nsum = nsum + w[vstates[j]][k]*(math.exp(-1.0*((deltatn+dtw*k)**(2.0))*invtau[vstates[i]][vstates[j]]*invtau[vstates[i]][vstates[j]])-math.exp(-1.0*((dtw*k)**(2.0))*invtau[vstates[i]][vstates[j]]*invtau[vstates[i]][vstates[j]]))
	      #                                  dsum = dsum + w[vstates[j]][k]*math.exp(-1.0*((dtw*k)**(2.0))*invtau[vstates[i]][vstates[j]]*invtau[vstates[i]][vstates[j]])
	       #                                 k = k + 1
	        #                        velem = ((poparray[vstates[i]]*poparray[vstates[j]])**(0.5))*(1.0+nsum/dsum)
#				else:
#                                        k = 0
 #                                       nsum = 0.0
  #                                      dsum = 0.0
   #                                     while k <= nta:
#	                                        nsum = nsum + w[vstates[i]][k]*(math.exp(-1.0*((deltatn+dtw*k)**(2.0))*invtau[vstates[i]][vstates[j]]*invtau[vstates[i]][vstates[j]])-math.exp(-1.0*((dtw*k)**(2.0))*invtau[vstates[i]][vstates[j]]*invtau[vstates[i]][vstates[j]]))
#	                                        dsum = dsum + w[vstates[i]][k]*math.exp(-1.0*((dtw*k)**(2.0))*invtau[vstates[i]][vstates[j]]*invtau[vstates[i]][vstates[j]])
#	                                        k = k + 1
#	                                velem = ((poparray[vstates[i]]*poparray[vstates[j]])**(0.5))*(1.0+nsum/dsum)
			else:
				velem = ((poparray[vstates[i]]*poparray[vstates[j]])**(0.5))*math.exp(-1.0*deltatn*invtau[vstates[i]][vstates[j]])
				eseg[i][j] = math.exp(-1.0*deltatn*invtau[i][j])
				eseg[j][i] = eseg[i][j]
			pass
			vtarget.append(velem)
			j = j + 1
		pass
		i = i + 1
	pass


	iter = 0        # Used to track what column is sent to minelem
	bcore = []      # block coordinates for fastest decaying element
	listbank = []

	p, q, leave = minelem(eseg,rank,nzthresh)

	bcore.append(p)
	bcore.append(q)

#       print 'bcore'
#       print bcore

#       print 'pre-algorithm check before individual element stuff is made'
#       sys.exit()

	while leave == 'no':

                # Finding the largest possible coherent sub-block
		i = 0
		bstates = []
		while i < rank:         # starting by adding all possible states
			bstates.append(i)
			i = i + 1
		pass

                # checking for zeroes in eseg such that blocks can be removed from the coherent block
		icstates = []   # list of states with zero elemnts in rows or columns shared by the
                                # coherence indicated by bcore

		if (nblock > 0):
			icstates = rstates(rank,eseg,bcore,nzthresh)    # function identifying the icoherent states
		pass

#               print 'icstates', icstates

		i = 0
		if (nblock == 0):
			pass
		else:
#                       print('is this getting killed by poor loop navigation?')
			while i < len(icstates):
				j = -(i + 1)
				temp = bstates.pop(icstates[j])
				i = i + 1
			pass
		pass

		listbank.append(bstates[:])

#               print 'bstates for current block', bstates

                # projecting out the current block from eseg
		iter = len(bstates)
		i = 0
		temp4 = eseg[bcore[0]][bcore[1]]
		while i < iter:
			j = 0
			while j < iter:
				eseg[bstates[i]][bstates[j]] = eseg[bstates[i]][bstates[j]] - temp4
				j = j + 1
			pass
			i = i + 1
		pass

#               print 'updated eseg'
#               print eseg

		nblock = nblock + 1
		if (nblock > 500):
			print('nblock blew up')
			print('Pnc matrix')
#                        print(Pnc)
			print('current eseg')
			print(eseg)
			sys.exit()
		pass

		bcore = []      # block coordinates for fastest decaying elements

		p, q, leave = minelem(eseg,rank,nzthresh)

		bcore.append(p)
		bcore.append(q)

	pass

	i = 0
	while i < rank:
		if (abs(poparray[vstates[i]]*eseg[i][i]) > nzthresh):
			listbank.append([])
			listbank[-1].append(i)
		pass
		i = i + 1
	pass


#	print 'vstates', vstates
#	print 'listbank'
#	print listbank




#	print 'vtarget'
#	print vtarget

#	sys.exit()

	# generating the vectorized set of coherent block density matrices
	# for the TAB wave function collapse step
	# These are determined using the current time step electronic populations,
	# and therefore must be recomputed at each collapse step. (Unlike the list of 
	# states populated in each block [listbank])

	A = [] 	#stores the block density matrices in form of A[block-ID][element]
	aindex = 0


	i = 0
	while i < len(listbank):
		k = 0
		A.append([])
		while k < rank:
			l = k
			while l < rank:
				if (k in listbank[i] and l in listbank[i]):
					m = 0
					popsum = 0.0
					while m < len(listbank[i]):
						popsum = popsum + poparray[vstates[listbank[i][m]]]
						m = m + 1
					pass
					if (popsum <= nzthresh):
						aelem = 0.0
					else:
						if (k == l):
							aelem = dgscale*((poparray[vstates[k]]*poparray[vstates[l]])**(0.5))/popsum
						else:
							aelem = ((poparray[vstates[k]]*poparray[vstates[l]])**(0.5))/popsum
						pass
					pass
				else:
					aelem = 0.0
				A[aindex].append(aelem)
				l = l + 1
			pass
			k = k + 1
		pass
		i = i + 1
		aindex = aindex + 1
	pass

#	Test to see if there is linear dependence, and if not will do the real
#	linear least squares optimization
	At = np.transpose(A)

#	ehrv = []
#	ehrv.append([])

#	i = 0
#	while i < len(A[-1]):
#		ehrv[0].append(A[-1][i])
#		i = i + 1
#	pass

#	ehrvt = np.transpose(ehrv)
#	ehrw = lsq_linear(ehrvt, vtarget, bounds=(0.0, 2.0), method='bvls')

#	if (abs(ehrw.x[0]-1.0) <= pehrptol):
#		i = 0
#		return poparray
#	pass


#	print 'A[dimH-1]'
#	print A[dimH-1]

	# scipy linear least squares optimization
#	At = np.transpose(A)
	optw = lsq_linear(At, vtarget, bounds=(0.0, 2.0), method='bvls')

	# Error analysis of the linear least squares target wave function



	# Collapsing the wave function

	i = 0
	ptotal = 0.0
	while i < len(optw.x):
		ptotal = ptotal + optw.x[i]
		i = i + 1
	pass

#	print 'optw'
#	print optw.x

	check2 = random.random()
	check = check2*ptotal

	track = 0
	psum = 0.0
	i = 0
	j = 0
	while j < 2:
		psum = psum + optw.x[i]
		if (check <= psum):
			j = 3
			track = i
		else:
			i = i + 1
		pass
	pass

	# track is the collapsed into density matrix according to A
	# constructing npop from the density matrix

	if (track == 0):
		return poparray
	pass

	k = 0
	index = 0
	while k < rank:
		l = k
		while l < rank:
			if (k == l):
				npop[vstates[k]] = A[track][index]/dgscale
			index = index + 1
			l = l + 1
		k = k + 1
	pass

#	print 'npop', npop
#	print 'vectorized target'
#	print A[track]


	return npop

