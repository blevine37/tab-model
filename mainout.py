#writing to output files
import numpy as np
import sys	

def writemain (t, dimH, x1, ct, odotx1, H, posout, eneout, popout, dpopout, outp, pmass):

	#Calculating MF energy 
	ccont = np.transpose(np.conjugate(ct))	#Complex transpose of ct
	norm = np.dot(ccont,ct) 	#norm^2 of the wf (used in calculating expectation values
	EMF = np.dot(ccont,np.dot(H, ct))/norm 	#Mean field energy <ct*|H|ct>/<ct*|ct>

	dp = 1e-6 	#Placeholder for double precision

	if (EMF.imag > dp):	#Checking if imaginary portion of EMF is large
		outp.write(str('imaginary mean field energy \n'))
		outp.write(str(EMF.imag) + '\n')
		sys.exit()
	pass

	EMFr = EMF.real 	#Real part of the MF energy
	KE = 0.5*pmass*(odotx1**2.0)	#Kinetic energy
	Etot = EMFr + KE 		#Total Energy

	if (norm.imag > dp): 	#Checking to see if imag part of wf norm is small
		outp.write(str('imaginary wave function norm \n'))
		outp.write(str(norm.imag) + '\n')
	pass

	sqrnorm = (norm.real)**(0.5) 	#Norm as Ben defines it sq(ct*ct)
	
	#---------------diabatic population---------------------
	dpoparray = np.zeros((dimH))
	for i in range(dimH):
		dpoparray[i] = np.linalg.norm(ct[i])**2.0
		
	dpoptot = sum(dpoparray)

	#---------------adiabatic population--------------------
	w,VR = np.linalg.eigh(H)
	dPE = np.amin(abs(EMFr-w))	#Energy difference between mean field energy and closest state energy

	amp = np.zeros((dimH),dtype = np.complex)
	poparray = np.zeros((dimH))
	for i in range(dimH):
		amp[i] = np.dot(np.transpose(VR[:,i]), ct)/sqrnorm
		poparray[i] = np.linalg.norm(amp[i])**2

	poptot = sum(poparray)

	# ======================================================
	# Formatting and writing the outputs
	posout.write('{:>8.4f}{:>20.10f}\n'.format(t, x1))

	eneout.write('{:>8.4f}{:>20.10f}{:>20.10f}{:>20.10f}{:>20.10f}\n'.format(t, EMFr, dPE, Etot, sqrnorm))
	
	lineout = []
	lineout = ['{:>20.10f}'.format(i) for i in poparray]
	popout.write('{:>8.4f}{}{:>20.10f}\n'.format(t, ''.join(lineout), poptot))
	
	lineout = []
	lineout = ['{:>20.10f}'.format(i) for i in dpoparray]
	dpopout.write('{:>8.4f}{}{:>20.10f}\n'.format(t, ''.join(lineout), dpoptot))
	
	return