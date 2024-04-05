#Using linear least squares to obtain collapse weights
#Dynamics across a two-state Tully model 3
#with TAB
#-----------------Import python packages----------------------------------------
import numpy as np		#Python matrix operation package
import sys				#use sys.exit() to stop program
import random			#Random number generator
import math

#-----------------Import custom functions----------------------------------------
from calcH import buildH		#Construct the Hamiltonian matrix for a given x1
from diffH import dHcalc		#Caculate  the derivative of the Hamiltonian matrix
from mainout import writemain 	#Write out select quantities of interest to an output file
from WFprop import stepWF		#Propogate WF forward in t
from movepos import movex		#Step positions forward
from vstep import vcalc			#Step velocity forward
from efF import calEff			#Calculate Ehrenfest forces
from cgauss import gcollapse	#Collapse after coherence is lost into a pure state
		
#----------------Creating log file for output-----------------------------------
outp = (open('run.log', 'w'))
outp.write('All packages and functions loaded successfully \n')

#----------------Creating Collapsed output file---------------------------------
outcollap = (open('collapse.log', 'w'))
outcollap.write('Collapse record at all trajectories \n')

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
dimH = 2

# parameters for Hamiltonian
A = 6.0e-4
B = 0.10
C = 0.90

#nuclear simulation time step
deltatn = 0.03

#half the number of electronic time steps during a nuclear time step
hnstepe = 50

#electronic time step
deltate = deltatn/(2.0*hnstepe)

#Number of trajectories to be run in a calculation
trajnum = 1

#Maximum number of nuclear time steps within a simulation
tstepmax = 120000

#Decoherence correction parameter on x1 direction
dcp1 = 6.0

#Particle mass (Nuclear mass)
pmass = 2000.0

# Large number checked against to ensure infinite loops
qkill = 100000

# Number of time steps between output writes
twrite = 10

#----------------------Initial Conditions-------------------------------
#Real part of time dependent WF ct
cr = np.zeros((dimH))
cr[0] = 1.0
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
k = 1
while k <= trajnum:	
	#-----------Initial Conditions for trajectory-k -----------------
	np.random.seed(1)
	#x1 = np.random.normal(-15.0,4.0)	#Initial paritcle position on x1-direction
	x1 = -15.0
	#odotx1 = np.random.normal(5.0,0.5)/pmass	#Initial particle velocity on x1-direction
	odotx1 = 5.0/pmass
	KE = 0.5*pmass*(odotx1**2.0)	#Initial kinetic energy
	t = 0.0		#Initial simulation time

	#-----------Initilize the Time dependent wave function----------
	ct = np.zeros((dimH),dtype=complex)
	ct[0] = 1.0
	
	#Calculating the Hamiltonian matrix at initial positions
	H = buildH(dimH, x1, A, B, C)

	#------------Creating trajectory-k specific output files--------

	#Opening trajectory specific position output file
	posout = (open('pos' + str(k) + '.dat', 'w'))

	#Heading for position output file of each trajectory
	posout.write('{:>8}{:>20}\n'.format('t','x'))

	#Opening trajactory specific energy and norm output file
	eneout = (open('ene' + str(k) + '.dat', 'w'))

	#Heading for energy and norm output file
	eneout.write('{:>8}{:>20}{:>20}{:>20}{:>20}\n'.format('t','Mean-Field Energy',
												'PE Difference','Total Energy','norm'))

	#Opening trajectory specific state population output file
	popout = (open('pop' + str(k) + '.dat', 'w'))

	#Heading for population output file
	heading = []
	for i in range(dimH):
		heading.append('{:>20}'.format('state ' + str(i) + ' pop'))
	jheading = ''.join(heading)
	popout.write('{:>8}{}{:>20}\n'.format('t', jheading, 'poptot'))
		
	#Opening trajectory specific diabatic population output file
	dpopout = (open('dpop' + str(k) + '.dat', 'w'))
	
	#Heading for diabaticpopulation output file
	heading = []
	for i in range(dimH):
		heading.append('{:>20}'.format('state ' + str(i) + ' pop'))
	jheading = ''.join(heading)
	dpopout.write('{:>8}{}{:>20}\n'.format('t', jheading, 'poptot'))

	#Writing out t = 0 outputs
	null = writemain(t, dimH, x1, ct, odotx1, H, posout, eneout, popout, dpopout, outp, pmass)

	#----------Begin the simulation---------------------------------
	#Compute the Ehrenfest forces
	mfdF1 = calEff(dimH, ct, x1, A, B, C)

	amp = np.zeros((dimH),dtype = complex)
	poparray = np.zeros((dimH))
	oldpop = np.zeros((dimH))

	normct = np.linalg.norm(ct)	#norm of complex vector ct

	w, VR = np.linalg.eigh(H)	#calculate the eigenvalues(sorted in ascending order)
								#and eigenvectors of H, solve TISE to get E_j and psi_j
	
	for i in range(dimH):
		amp[i] = np.dot(np.transpose(VR[:,i]), ct)/normct	#calculate the amplitude c_i of each state and 
															#normalize WF so that total population is 1
		poparray[i] = np.linalg.norm(amp[i])**2.0

	oldpop = poparray

	oldforce1 = np.zeros((dimH))
	dH1 = dHcalc(dimH, x1, A, B, C)

	for i in range(dimH):
		force1 = -np.dot(np.transpose(VR[:,i]),np.dot(dH1, VR[:,i]))	#force = -dE/dx = -<psi_i|dH/dx|psi_i>
		oldforce1[i] = force1
	
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
		
#continue from here
		x1 = movex(x1, odotx1, acel1, deltatn)

		#Calculte H at the new position
		H = buildH(dimH, x1, A, B, C)

		#Propagate WF through the other half time step using H(t+dt)
		i = 0
		while i < hnstepe:
			cr, ci = stepWF(deltate, cr, ci, H)	#Forward WF one electronic time step
			i = i + 1
		pass

		ct = cr+1j*ci
		ccont = np.transpose(np.conjugate(ct))
		cnorm = np.dot(ccont,ct)
		normct = np.linalg.norm(ct)

		#Storing forces from the last step
		mfdFprev1 = mfdF1
		#Calculate Ehrenfest forces
		mfdF1 = calEff(dimH, ct, x1, A, B, C)

		#Step velocities forward in time
		odotx1 = vcalc(odotx1, mfdF1, mfdFprev1, deltatn, pmass)
		

		#-------------- TAB Starts from Here ---------------------------------------
		poparray = np.zeros((dimH))	#array holding state populations
		
		ampdir = np.zeros((dimH),dtype=complex)	#Stores amplitude directions for each state
		amp = np.zeros((dimH),dtype=complex)		#Stores amplitudes for each state
		KE = 0.5*pmass*(odotx1**2.0)

		w, VR = np.linalg.eigh(H)
		tVR = np.transpose(VR)

		Estates = np.zeros((dimH))
		Estates = w

		temp1 = np.zeros((1),dtype=complex)

		i = 0
		while i < dimH:
			amp[i] = np.dot(tVR[i,:],ct)/normct
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
		
		EMF = np.dot(ccont,np.dot(H,ct))/normct**2.0
		rEMF = EMF.real
		roldEMF = rEMF
   
		dH1 = dHcalc(dimH, x1, A, B, C)
		newforce1=np.zeros((dimH))        #Adiabatic State Force along x1 direction
		i = 0 
		while i < dimH:
			force1 = -np.dot(tVR[i,:],np.dot(dH1,VR[:,i]))
			newforce1[i] = force1
			i = i + 1
		pass
		
		aforce1 = np.zeros((dimH))
		
		odotrho = np.zeros((dimH))
		i = 0
		while i < dimH:
			aforce1[i] = (oldforce1[i] + newforce1[i])/2.0
			odotrho[i] = (poparray[i] - oldpop[i])/deltatn
			i = i + 1
		pass
		
		#----------- New Collapse Routine Goes Here ---------------------
		npop = np.zeros((dimH))
		npop, track = gcollapse(dimH,deltatn,aforce1,poparray,dcp1,nzthresh,errortol,npthresh,pehrptol,odotrho,tolodotrho,nta,dtw,zpop,dgscale)
		
		oldforce1 = np.zeros((dimH))
		i = 0 
		while i < dimH:
			oldforce1[i] = newforce1[i]
			i = i + 1
		pass
		
		if (track!=0):
			outcollap.write("collapsed happened at " + str(k) + "\n")
			poparray = npop
		
			namp = np.zeros((dimH),dtype=complex)
			nct = np.zeros((dimH,1),dtype=complex)
		
			i = 0
			while i < dimH:
				namp[i] = ampdir[i]*(npop[i]**(0.5))*normct
				i = i+1
			pass

			i = 0
			while i < dimH:
				nct = nct + namp[i]*np.transpose([tVR[i]])
				i = i+1
			pass
	
			i = 0
			while i < dimH:	
				cr[i] = nct[i][0].real
				ci[i] = nct[i][0].imag
				i = i + 1
			pass
			
			oldct = ct
			ct = cr+1j*ci
			ccon = np.conjugate(ct)
			ccont = np.transpose(ccon)
			oldnorm = cnorm
			cnorm = np.dot(ccont,ct)
		
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
			odotx1 = math.copysign((2.0*nKE/pmass)**0.50,odotx1)
		else:
			odotx1 = odotx1
	
		#if track != 0 and nKE >= 0.0:
		#	odotx1 = math.copysign((2.0*nKE/pmass)**0.50,odotx1)  #NOT Reverse the momentum if collapsed
		#elif track != 0 and nKE < 0.0:
		#	odotx1 = odotx1
		#	ct = oldct

		mfdF1=calEff(dimH, ct, x1, A, B, C) #Calculate mean field force
		
		i = 0   #store the poparray for next odotrho
		while i < dimH:
			oldpop[i] = poparray[i]
			i = i + 1
		pass

		if (n%twrite == 0):
			null = writemain(t,dimH,x1,ct,odotx1,H,posout,eneout,popout,dpopout,outp,pmass)
		pass
		
		n = n+1 #forward one time step
	pass
	
	outp.write('Program is finished calculating trajectory ' + str(k) + '\n')

	k = k + 1
pass

outp.write('Normal termination of program \n')
