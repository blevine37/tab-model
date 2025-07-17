#Using linear least squares to obtain collapse weights
#Dynamics across a two-state 'coupled harmonic oscillator' Conical Intersection
#with TAB
#-----------------Import python packages----------------------------------------
import numpy as np 		#Python matrix operation package
import sys			#use sys.exit() to stop program
import random 			#Random number generator
import math

#-----------------Import custom functions----------------------------------------
from calcH import buildH		#Construct the Hamiltonian matrix for a give x1, x2
from diffH import dHcalc
from mainout import writemain 		#Writes out select quantities of interest to an output file
from WFprop import stepWF		#Propogate WF forward in t
from movepos import movex		#Step positions forward
from vstep import vcalc			#Step velocity forward
from efF import calEff			#Calculate Ehrenfest forces
from hwrsort import eigsort		#Sorts eigens of a matrix in ascending order
from cgauss import gcollapse		#Collapses after coherence is lost into a pure state
		
#----------------Creating log file for output-----------------------------------
outp = (open('run.log', 'w'))
outs = str('All packages and functions loaded successfully \n')
outp.write(outs)

#------------------User-adjustable parameters-----------------------------------
# Threshold for considering numbers numerically zero
nzthresh = 1.0e-10

# Tolerance for cumulative errors in collapse probabilities
errortol = 1.0e-6

# Tolerance for negative probabilities
npthresh = 1.0e-7

# 
pehrptol = 1.0e-5

# Tolerance for a positive rate of population change for an individual electronic state
tolodotrho = 1.0e-5

# Scales diagonal elements in linear-least squares in order to emphasize relative population
# conservation through the collapse of the wave function
dgscale = 1.0e+5

# Discrete time step for integrating over simulation history
dtw = 0.010

# How many numerical steps in the integration over simulation history
nta = 600

#
zpop = 1.0e-6

# number of total states in a single trajectory
dimH = 9

#number of nuclear degrees of freedom
ndof = 2

# (absolute) slope of diabatic1 potential along x1 direction
w1 = 0.25

# slope of diabatic2-dimH potential along x1 direction
w2 = 0.25

# linear coupling constant
c = 0.025

# spacing distance
delta = 0.01

#nuclear simulation time step
deltatn = 0.05

#half the number of electronic time steps during a nuclear time step
hnstepe = 50

#electronic time step
deltate = deltatn/(2.0*hnstepe)

#Number of trajectories to be run in a calculation
trajnum = 1000
#Maximum number of nuclear time steps within a simulation
tstepmax = 6000

#Decoherence correction parameter in each direction
dcp = [6.0, 6.0]
rescale='old'

reverse=int(0) #reverse the velocity for frustrated hops; 0 for no, 1 for yes
#Particle mass (Nuclear mass)
pmass = 1845

# Large number checked against to ensure infinite loops
qkill = 100000

# Number of time steps between output writes
twrite = 10

#----------------------Initial Conditions-------------------------------
#Real part of time dependent WF ct
cr = np.zeros((dimH))
#Imaginary part of time dependent WF ct
ci = np.zeros((dimH))
#The Hamiltonian
H = np.zeros((dimH,dimH))

# ===================================================================
# Trajectory Initilization

# Initial population only on state-0
intpop = np.zeros((dimH))
intpop[0] = 1.000
x= np.zeros((ndof))	#Initial position vector
odotx= np.zeros((ndof))	#Initial velocity vector
#loops over trajectories-k
k = 101
while k <= trajnum:	
	#-----------Initial Conditions for trajectory-k -----------------
	#Initial conditions, buildH, and diffH are still defined manually along each dof
	np.random.seed(k)
	random.seed(k)	#Set the random seed for reproducibility
	x[0] = np.random.normal(0.0,0.204)-1.0	#Initial paritcle position on x1-direction
	x[1] = np.random.normal(0.0,0.204)	#Initial particle position on x2-direction
	odotx[0] = np.random.normal(10.0,2.451)/pmass	#Initial particle velocity on x1-direction
	odotx[1] = np.random.normal(10.0,2.451)/pmass	#Initial particle velocity on x2-direction
	KE = 0.5*pmass*(np.inner(odotx, odotx))	#Initial kinetic energy



	t = 0.0		#Initial simulation time

	#-----------Initilize the Time dependent wave function----------
	i = 0
	while i < dimH:
		cr[i] = 0
		ci[i] = 0
		i = i + 1
	pass
	cr[0] = 1.0
	#Total wave function
	ct = cr + 1j*ci 

	#Calculating the Hamiltonian matrix at initial positions
	H = buildH(dimH, x, w1, w2, c, delta)

	#------------Creating trajectory-k specific output files--------

	#Opening trajectory specific position output file
	line1 = str('pos')
	line2 = str(k)
	line3 = str('.dat')
	line = line1 + line2 + line3
	posout = (open(line, 'w'))

	#Heading for position output file of each trajectory
	line1 = str('t').rjust(8)
	heading = line1
	for linei in range(ndof):
		linei = str('x' + str(linei+1)).rjust(20)
		heading+=linei
	posout.write(heading + '\n')

	#Opening trajectory specific velocity output file
	line1 = str('vel')
	line2 = str(k)
	line3 = str('.dat')
	line = line1 + line2 + line3
	velout = (open(line, 'w'))

	#Heading for velocity output file of each trajectory
	line1 = str('t').rjust(8)
	for linei in range(ndof):
		linei = str('v' + str(linei+1)).rjust(20)
		heading+=linei
	velout.write(heading + '\n')

	#Opening trajactory specific energy and norm output file
	line1 = str('ene')
	line2 = str(k)
	line3 = str('.dat')
	line = line1 + line2 + line3
	eneout = (open(line, 'w'))

	#Heading for energy and norm output file
	line1 = str('t').rjust(8)
	line2 = str('Mean-Field Energy').rjust(20)
	line2p = str('PE Difference').rjust(20)
	line3 = str('Total Energy').rjust(20)
	line4 = str('norm').rjust(20)
	heading = line1 + line2 + line2p + line3 + line4 + '\n'
	eneout.write(heading)

	#Opening trajectory specific state population output file
	line1 = str('pop')
	line2 = str(k)
	line3 = str('.dat')
	line = line1 + line2 + line3
	popout = (open(line, 'w'))

	#Heading for population output file
	i = 0
	heading1 = []
	line1 = str('t').rjust(8)
	while i < dimH:
		temp1 = str('state ')
		temp2 = str(i)
		temp3 = str(' pop')
		temp4 = temp1+temp2+temp3
		heading1.append(temp4.rjust(20))
		i = i + 1
	pass
	line2 = ''.join(heading1)
	line3 = str('poptot').rjust(20)
	heading = line1 + line2 + line3 + '\n'
	popout.write(heading)
		
	#Opening trajectory specific diabatic population output file
	line1 = str('dpop')
	line2 = str(k)
	line3 = str('.dat')
	line = line1 + line2 + line3
	dpopout = (open(line, 'w'))
	
	#Heading for diabaticpopulation output file
	i = 0
	heading1 = []
	line1 = str('t').rjust(8)
	while i < dimH:
		temp1 = str('state ')
		temp2 = str(i)
		temp3 = str(' pop')
		temp4 = temp1+temp2+temp3
		heading1.append(temp4.rjust(20))
		i = i + 1
	pass
	line2 = ''.join(heading1)
	line3 = str('poptot').rjust(20)
	heading = line1 + line2 + line3 + '\n'
	dpopout.write(heading)

	#Writing out t = 0 outputs
	null = writemain(t,dimH,ndof,x,ct,odotx,H,posout,velout,eneout,popout,dpopout,outp,pmass)

	#----------Begin the simulation---------------------------------
	#Compute the Ehrenfest forces
	mfdF = calEff(dimH, ndof, ct, x, w1, w2, c, delta)
	#----changed upto here
	amp = np.zeros((dimH),dtype=complex)
	poparray = np.zeros((dimH))
	oldpop = np.zeros((dimH))

	ccon = np.conjugate(ct)
	ccont = np.transpose(ccon)
	cnorm = np.dot(ccont,ct)
	norm2ct = (cnorm.real)**(0.50)

	w, VR = np.linalg.eigh(H)
	sw, sVR = eigsort(dimH,w,VR)
	tsVR = np.transpose(sVR)
	
	temp1 = np.zeros((1),dtype=complex)
	
	i = 0
	while i < dimH:
		amp[i] = np.dot(tsVR[i,:],ct)/norm2ct	
		temp1[0] = amp[i]
		temp2 = np.conjugate(temp1)
		temp3 = np.transpose(temp2)
		temp4 = np.dot(temp3,temp1)
		poparray[i] = temp4.real
		i = i + 1
	pass
	i = 0
	while i < dimH:
		oldpop[i] = poparray[i]
		i = i + 1
	pass

	oldforce = np.zeros((ndof, dimH))
	dH = dHcalc(dimH,ndof,x,w1,w2,c,delta)   #diratives of diabatic Hamiltonian
	
	for j in range(ndof):
		for i in range(dimH):
			oldforce[j,i]= -np.dot(tsVR[i,:],np.dot(dH[j,:,:],sVR[:,i]))


	#Time steps count
	n = 1
	
	while n <= tstepmax:
		t = t + deltatn	#step time forward

		#Propagate WF through the first half of a time step using the H(t)
		i = 0
		while i < hnstepe:
			cr, ci = stepWF(deltate, cr, ci, H)	#Forward WF one electronic time step
			i = i + 1
		pass

		ct = cr+1j*ci

		#Step positions forward in time
		acel = mfdF/pmass

		
		x = movex(x, odotx, acel, deltatn)

		#Calculte H at the new position
		H = buildH(dimH, x, w1, w2, c, delta)

		#Propagate WF through the other half time step using H(t+dt)
		i = 0
		while i < hnstepe:
			cr, ci = stepWF(deltate, cr, ci, H)	#Forward WF one electronic time step
			i = i + 1
		pass

		ct = cr+1j*ci
		ccon = np.conjugate(ct)
		ccont = np.transpose(ccon)
		cnorm = np.dot(ccont,ct)
		norm2ct = (cnorm.real)**(0.50)

		#Storing forces from the last step
		mfdFprev = mfdF

		#Calculate Ehrenfest forces
		mfdF = calEff(dimH, ndof, ct, x, w1, w2, c, delta)
		#Step velocities forward in time
		odotx = vcalc(odotx, mfdF, mfdFprev, deltatn, pmass)
		

		#-------------- TAB Starts from Here ---------------------------------------
		poparray = np.zeros((dimH))	#array holding state populations
		
		ampdir = np.zeros((dimH),dtype=complex)	#Stores amplitude directions for each state
		amp = np.zeros((dimH),dtype=complex)		#Stores amplitudes for each state
		KE = 0.5*pmass*(np.inner(odotx, odotx))

		w, VR = np.linalg.eigh(H)
		sw, sVR = eigsort(dimH,w,VR)
		tsVR = np.transpose(sVR)

		Estates = np.zeros((dimH))
		Estates = sw

		temp1 = np.zeros((1),dtype=complex)

		i = 0
		while i < dimH:
			amp[i] = np.dot(tsVR[i,:],ct)/norm2ct
			temp1[0] = amp[i]
			temp2 = np.conjugate(temp1)
			temp3 = np.transpose(temp2)
			temp4 = np.dot(temp3,temp1)
			poparray[i] = temp4.real
			
			if (poparray[i] < nzthresh):
				ampdir[i] = 1.0
			else:
				ampdir[i] = amp[i]/(poparray[i]**(0.5))
			pass
			i = i + 1
		pass
		
		EMF = np.dot(ccont,np.dot(H,ct))/cnorm
		rEMF = EMF.real
		roldEMF = rEMF
		
		dH = dHcalc(dimH,ndof,x,w1,w2,c,delta)   #diratives of diabatic Hamiltonian
		newforce=np.zeros((ndof,dimH))        #Adiabatic State Force 

		for j in range(ndof):
			for i in range(dimH):
				newforce[j,i]= -np.dot(tsVR[i,:],np.dot(dH[j,:,:],sVR[:,i]))

		
		aforce = np.zeros((ndof,dimH))
		aforce = (oldforce + newforce)/2.0	#Average force for the current time step
		odotrho = np.zeros((dimH))
		i = 0
		while i < dimH:
			odotrho[i] = (poparray[i] - oldpop[i])/deltatn
			i = i + 1
		pass
		
		#----------- New Collapse Routine Goes Here ---------------------
		npop = np.zeros((dimH))
		npop, track = gcollapse(dimH,ndof,deltatn,aforce,poparray,dcp,nzthresh,errortol,npthresh,pehrptol,odotrho,tolodotrho,nta,dtw,zpop,dgscale)
		
		oldforce = np.zeros((ndof,dimH))
		oldforce = newforce

		poparray = npop

		namp = np.zeros((dimH),dtype=complex)
		nct = np.zeros((dimH,1),dtype=complex)
		
		i = 0
		while i < dimH:
			namp[i] = ampdir[i]*(npop[i]**(0.5))*norm2ct
			i = i+1
		pass
		#print ('namp',namp)
		#print ('namp norm', np.dot(namp, np.transpose(np.conj(namp))))
		#namp_residue= np.zeros((dimH,1),dtype=complex)
		i = 0
		while i < dimH:
			nct = nct + namp[i]*np.transpose([tsVR[i]])
			i = i+1
		pass
		EMF = (np.dot(np.conjugate(np.transpose(nct)),np.dot(H,nct))/np.dot(np.conjugate(np.transpose(nct)),nct)).item()
		if (abs(EMF.imag) > nzthresh):
			print('Warning: Mean-field energy has non-zero imaginary part, which is unexpected.')
			print('EMF:', EMF)
			sys.exit()
		rEMF = EMF.real
		a_f = np.dot((nct.flatten()),ct)
		vrescale = np.ones(ndof)/np.sqrt(ndof)  # Effective "NAC" vector for rescaling velocities. We normalize it later.
		deltav = 0.0
		if(track==0):
			vrescale = np.ones(ndof)
		else:
			print('old pop', oldpop)
			print('npop', npop)
			a_r = np.sqrt(np.maximum(0.0,1-abs(a_f)**2.0)) 
			print('a_f', a_f)
			print('a_r', a_r)
			if (abs(a_r) > nzthresh):
				print('nct', nct)
				print('ct', ct)
				ct_residue = (ct - a_f*nct.flatten()) / a_r
				print('ct_residue', ct_residue)
				#print('namp',namp)
				#print('amp_residue',amp_residue)
				F_i = calEff(dimH, ndof, ct, x, w1, w2, c, delta)
				print('F_i', F_i)
				F_f = calEff(dimH, ndof, nct.flatten(), x, w1, w2, c, delta)
				print('F_f', F_f)
				F_r = calEff(dimH, ndof, ct_residue, x, w1, w2, c, delta)
				vrescale = F_i - (abs(a_f)**2)*F_f - (abs(a_r)**2)*F_r
				vrescale = vrescale / np.linalg.norm(vrescale)
				print('vrescale', vrescale)
				


		#--Rescaling Kinetic Energy-------
		deltaEMF = rEMF - roldEMF
		nKE = KE - deltaEMF
		if not (track==0):
			if (nKE < 0.0):
				print ('frustrated hops are needed')
				odotx = ((-1)**reverse)*odotx	#reverse the velocity or not
				#restore the old wave function
			elif(rescale=='old'):
				scale = (nKE/KE)**0.50
				odotx = odotx*scale
				
				i = 0
				while i < dimH:	
					cr[i] = nct[i][0].real
					ci[i] = nct[i][0].imag
					i = i + 1
				pass
				
				ct = cr+1j*ci
				ccon = np.conjugate(ct)
				ccont = np.transpose(ccon)
				oldnorm = cnorm
				cnorm = np.dot(ccont,ct)
				norm2ct = (cnorm.real)**(0.50)

			else:
				tempv=np.dot(odotx, np.transpose(vrescale))
				print('tempv',tempv)
				print('deltaKE',-1*deltaEMF)
				discriminant = tempv**2 - (2*deltaEMF/pmass)
				if discriminant >= 0: #we pick the root of smaller absolute value
					if (tempv<0):
						deltav = -1*tempv - np.sqrt(discriminant)
					else:
						deltav = -1*tempv + np.sqrt(discriminant)
					
					print('deltav',deltav)
					print('odotx',odotx)
					print('vrescale',vrescale)
					oldKE=0.5*pmass*(np.inner(odotx, odotx))
					odotx = odotx + deltav*vrescale		
					print('odotx after rescale',odotx)
					newKE = 0.5*pmass*(np.inner(odotx, odotx))
					print('KE sanity check:', newKE-oldKE+ deltaEMF)
					i = 0
					while i < dimH:	
						cr[i] = nct[i][0].real
						ci[i] = nct[i][0].imag
						i = i + 1
					pass
					
					ct = cr+1j*ci
					ccon = np.conjugate(ct)
					ccont = np.transpose(ccon)
					oldnorm = cnorm
					cnorm = np.dot(ccont,ct)
					norm2ct = (cnorm.real)**(0.50)
				
				else:
					print("Warning: Negative discriminant encountered. Assuming frustrated hop..")
					print("Discriminant:", discriminant)
					odotx = ((-1)**reverse)*odotx	#reverse the velocity, or not

		
		#--Norm Conservation Check---------
#		if (abs(cnorm-oldnorm).real >= 1e-12 or abs(cnorm-oldnorm).imag >= 1e-12 ):
#			print ('oldnorm',oldnorm)
#			print ('new norm',cnorm)
#			print ('norm difference',abs(cnorm-oldnorm).real,abs(cnorm-oldnorm).imag)
#			sys.exit()
#		pass
		
		mfdF = calEff(dimH, ndof, ct, x, w1, w2, c, delta)
		
		

		i = 0
		while i < dimH:
			oldpop[i] = poparray[i]
			i = i + 1
		pass

		if (n%twrite == 0):
			null = writemain(t,dimH,ndof,x,ct,odotx,H,posout,velout,eneout,popout,dpopout,outp,pmass)
		pass
		
		n = n+1 #forward one time step
	pass
	
	lout1 = str('Program is finished calculating trajectory ')
	lout2 = str(k)
	lout = lout1 + lout2 + '\n'
	outp.write(lout)

	k = k + 1
pass

outs = str('Normal termination of program \n')
outp.write(outs)
