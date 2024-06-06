# Calculate Ehrenfest force
def calEff(dimH, ct, x1, w1, w2, delta, c):
	
	import numpy as np

	from diffH import dHcalc
	
	ctbra = np.transpose(np.conjugate(ct))
	
	dH1= dHcalc(dimH, x1, w1, w2, c)

	cnorm = np.dot(ctbra, ct)

	Efft1 = -np.dot(ctbra, np.dot(dH1, ct))/cnorm
	Eff1 = Efft1.real

	#Efft2 = -np.dot(ctbra, np.dot(dH2, ct))/cnorm
	#Eff2 = Efft2.real

	return Eff1


	
