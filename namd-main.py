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

# (absolute) slope of diabatic1 potential along x1 direction
w1 = 0.25

# slope of diabatic2-dimH potential along x1 direction
w2 = 0.025

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

#Decoherence correction parameter on x1 direction
dcp1 = 6.0
#Decoherence correction parameter on x2 direction
dcp2 = 6.0

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

#loops over trajectories-k
k = 450
while k <= trajnum:	
	#-----------Initial Conditions for trajectory-k -----------------
	np.random.seed(k)
	x1 = np.random.normal(0.0,0.204)-1.0	#Initial paritcle position on x1-direction
	x2 = np.random.normal(0.0,0.204)	#Initial particle position on x2-direction
	odotx1 = np.random.normal(10.0,2.451)/pmass	#Initial particle velocity on x1-direction
	odotx2 = np.random.normal(10.0,2.451)/pmass	#Initial particle velocity on x2-direction
	KE = 0.5*pmass*(odotx1**2.0+odotx2**2.0)	#Initial kinetic energy
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
	H = buildH(dimH, x1, x2, w1, w2, c, delta)

	#------------Creating trajectory-k specific output files--------

	#Opening trajectory specific position output file
	line1 = str('pos')
	line2 = str(k)
	line3 = str('.dat')
	line = line1 + line2 + line3
	posout = (open(line, 'w'))

	#Heading for position output file of each trajectory
	line1 = str('t').rjust(8)
	line2 = str('x').rjust(20)
	line3 = str('y').rjust(20)
	heading = line1 + line2 + line3 + '\n'
	posout.write(heading)

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
	null = writemain(t,dimH,x1,x2,ct,odotx1,odotx2,H,posout,eneout,popout,dpopout,outp,pmass)

	#----------Begin the simulation---------------------------------
	#Compute the Ehrenfest forces
	mfdF1, mfdF2 = calEff(dimH,ct, x1, x2, w1, w2, c, delta)
	
	amp = np.zeros((dimH),dtype=np.complex)
	poparray = np.zeros((dimH))
	oldpop = np.zeros((dimH))

	ccon = np.conjugate(ct)
	ccont = np.transpose(ccon)
	cnorm = np.dot(ccont,ct)
	norm2ct = (cnorm.real)**(0.50)

	w, VR = np.linalg.eigh(H)
	sw, sVR = eigsort(dimH,w,VR)
	tsVR = np.transpose(sVR)
	
	temp1 = np.zeros((1),dtype=np.complex)
	
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

	oldforce1 = np.zeros((dimH))
	oldforce2 = np.zeros((dimH))
	dH1,dH2 = dHcalc(dimH,x1,x2,w1,w2,c,delta)   #diratives of diabatic Hamiltonian
	i = 0
	while i < dimH:
		force1 = -np.dot(tsVR[i,:],np.dot(dH1,sVR[:,i]))
		force2 = -np.dot(tsVR[i,:],np.dot(dH2,sVR[:,i]))
		oldforce1[i] = force1
		oldforce2[i] = force2
		i = i + 1
	pass

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
		acel1 = mfdF1/pmass
		acel2 = mfdF2/pmass
		
		x1 = movex(x1, odotx1, acel1, deltatn)
		x2 = movex(x2, odotx2, acel2, deltatn)

		#Calculte H at the new position
		H = buildH(dimH, x1, x2, w1, w2, c, delta)

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
		mfdFprev1, mfdFprev2 = mfdF1, mfdF2

		#Calculate Ehrenfest forces
		mfdF1, mfdF2 = calEff(dimH, ct, x1, x2, w1, w2, c, delta)

		#Step velocities forward in time
		odotx1 = vcalc(odotx1, mfdF1, mfdFprev1, deltatn, pmass)
		odotx2 = vcalc(odotx2, mfdF2, mfdFprev2, deltatn, pmass)
		

		#-------------- TAB Starts from Here ---------------------------------------
		poparray = np.zeros((dimH))	#array holding state populations
		
		ampdir = np.zeros((dimH),dtype=np.complex)	#Stores amplitude directions for each state
		amp = np.zeros((dimH),dtype=np.complex)		#Stores amplitudes for each state
		KE = 0.5*pmass*(odotx1**2.0+odotx2**2.0)

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
			
			if (poparray[i] == 0):
				ampdir[i] = 1.0
			else:
				ampdir[i] = amp[i]/(poparray[i]**(0.5))
			pass
			i = i + 1
		pass
		
		EMF = np.dot(ccont,np.dot(H,ct))/cnorm
		rEMF = EMF.real
		roldEMF = rEMF
		
		dH1,dH2 = dHcalc(dimH,x1,x2,w1,w2,c,delta)   #diratives of diabatic Hamiltonian
		newforce1=np.zeros((dimH))        #Adiabatic State Force along x1 direction
		newforce2=np.zeros((dimH))        #Adiabatic State Force along x2 direction
		i = 0 
		while i < dimH:
			force1 = -np.dot(tsVR[i,:],np.dot(dH1,sVR[:,i]))
			force2 = -np.dot(tsVR[i,:],np.dot(dH2,sVR[:,i]))
			newforce1[i] = force1
			newforce2[i] = force2
			i = i + 1
		pass
		
		aforce1 = np.zeros((dimH))
		aforce2 = np.zeros((dimH))
		
		odotrho = np.zeros((dimH))
		i = 0
		while i < dimH:
			aforce1[i] = (oldforce1[i] + newforce1[i])/2.0
			aforce2[i] = (oldforce2[i] + newforce2[i])/2.0
			odotrho[i] = (poparray[i] - oldpop[i])/deltatn
			i = i + 1
		pass
		
		#----------- New Collapse Routine Goes Here ---------------------
		npop = np.zeros((dimH))
		npop = gcollapse(dimH,deltatn,aforce1,aforce2,poparray,dcp1,dcp2,nzthresh,errortol,npthresh,pehrptol,odotrho,tolodotrho,nta,dtw,zpop,dgscale)
		
		oldforce1 = np.zeros((dimH))
		oldforce2 = np.zeros((dimH))
		i = 0 
		while i < dimH:
			oldforce1[i] = newforce1[i]
			oldforce2[i] = newforce2[i]
			i = i + 1
		pass

		poparray = npop
		
		namp = np.zeros((dimH),dtype=complex)
		nct = np.zeros((dimH,1),dtype=complex)
		
		i = 0
		while i < dimH:
			namp[i] = ampdir[i]*(npop[i]**(0.5))*norm2ct
			i = i+1
		pass

		i = 0
		while i < dimH:
			nct = nct + namp[i]*np.transpose([tsVR[i]])
			i = i+1
		pass
	
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
		
		#--Norm Conservation Check---------
#		if (abs(cnorm-oldnorm).real >= 1e-12 or abs(cnorm-oldnorm).imag >= 1e-12 ):
#			print ('oldnorm',oldnorm)
#			print ('new norm',cnorm)
#			print ('norm difference',abs(cnorm-oldnorm).real,abs(cnorm-oldnorm).imag)
#			sys.exit()
#		pass
		
		mfdF1, mfdF2 = calEff(dimH, ct, x1, x2, w1, w2, c, delta)
		EMF = np.dot(ccont,np.dot(H,ct))/cnorm
		rEMF = EMF.real
		
		#--Rescaling Kinetic Energy-------
		deltaKE = rEMF - roldEMF
		nKE = KE - deltaKE
	
		if (nKE < 0.0):
			print ('frustrated hops are needed')
			print ('aborting program')
			sys.exit()
		pass
		
		odotx2 = math.copysign(abs((2.0*nKE/pmass)-odotx1**2.0)**0.50,odotx2)	
		
		i = 0
		while i < dimH:
			oldpop[i] = poparray[i]
			i = i + 1
		pass

		if (n%twrite == 0):
			null = writemain(t,dimH,x1,x2,ct,odotx1,odotx2,H,posout,eneout,popout,dpopout,outp,pmass)
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
