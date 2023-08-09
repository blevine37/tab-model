def writemain (t,dimH,x1,ct,odotx1,H,posout,eneout,popout,dpopout,outp,pmass): #writing to output files
	import numpy as np
	import sys	
	from hwrsort import eigsort

	#Calculating MF energy 

	ccon = np.conjugate(ct) 	
	ccont = np.transpose(ccon)	#Complex transpose of ct

	norm = np.dot(ccont,ct) 	#norm of the wf (used in calculating expectation 
					#values

	EMF = np.dot(ccont,np.dot(H,ct))/norm 	#Mean field energy 

	dp = 1e-6 	#Placeholder for double precision

	if (EMF.imag > dp):	#Checking if imaginary portion of EMF is large
		outs = str('imaginary mean field energy \n')
		outp.write(outs)
		out1 = str(EMF.imag)
		outs = out1 + '\n'
		outp.write(outs)
		sys.exit()
	pass

	EMFr = EMF.real 	#Real part of the MF energy

	#Kinetic energy
	#KE = 0.5*pmass*(odotx1**2.0+odotx2**2.0)
	KE = 0.5*pmass*(odotx1**2.0)
	
	Etot = EMFr + KE 		#Total Energy

	inorm = norm.imag
	if (inorm > dp): 	#Checking to see if imag part of wf norm is small
		outs = str('imaginary wave function norm \n')
		outp.write(outs)
		out1 = str(inorm)
		outs = out1 + '\n'
		outp.write(outs)
		pass

	rnorm = norm.real 	#Real part of WF norm 
	sqrnorm = (rnorm)**(0.5) 	#Norm as Ben defines it sq(ct*ct)
	
	#---------------diabatic population---------------------
	dpoparray = np.zeros((dimH))
	i = 0
	while i < dimH:
		dpoparray[i] = np.dot(ccont[i],ct[i]).real
		i = i + 1
		
	dpoptot = sum(dpoparray)

	#---------------adiabatic population--------------------
	w,VR = np.linalg.eigh(H)
	sw,sVR = eigsort(dimH,w,VR)
	tsVR = np.transpose(sVR)

	dPE = np.amin(abs(EMFr-w))

	amp = np.zeros((dimH),dtype=np.complex)
	poparray = np.zeros((dimH))

	temp1 = np.zeros((1),dtype=np.complex)
	
	i = 0
	while i < dimH:
		amp[i] = np.dot(tsVR[i,:],ct)/sqrnorm
		temp1[0] = amp[i]
		temp2 = np.conjugate(temp1)
		temp3 = np.transpose(temp2)
		temp4 = np.dot(temp3,temp1)
		poparray[i] = temp4.real
		i = i + 1
	pass
	poptot = sum(poparray)

	# ======================================================
	# Formatting and writing the outputs

	line1 = format(t,'.4f').rjust(8)
	line2 = format(x1,'.10f').rjust(20)
	line3 = format(x2,'.10f').rjust(20)
	lineout = line1 + line2 + line3 + '\n'
	posout.write(lineout)

	line2 = format(EMFr,'.10f').rjust(20)
	line2p = format(dPE,'.10f').rjust(20)
	line3 = format(Etot,'.10f').rjust(20)
	line4 = format(sqrnorm,'.10f').rjust(20)
	lineout = line1 + line2 + line2p + line3 + line4 + '\n'
	eneout.write(lineout)
	
	i = 0
	lineout1 = []
	while i < dimH:
		lineTemp = format(poparray[i],'.10f').rjust(20)
		lineout1.append(lineTemp)
		i = i + 1
	pass
	line2 = ''.join(lineout1)
	line3 = format(poptot,'.10f').rjust(20)
	lineout = line1 + line2 + line3 + '\n'
	popout.write(lineout)
	
	i = 0
	lineout1 = []
	while i < dimH:
		lineTemp = format(dpoparray[i],'.10f').rjust(20)
		lineout1.append(lineTemp)
		i = i + 1
	pass
	line2 = ''.join(lineout1)
	line3 = format(dpoptot,'.10f').rjust(20)
	lineout = line1 + line2 + line3 + '\n'
	dpopout.write(lineout)
	
	return




