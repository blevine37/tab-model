# Calculate Ehrenfest force
def calEff(dimH, ct, x1, x2, w1, w2, c, delta):
	
	import numpy as np

	from diffH import dHcalc
	
	ctbra = np.transpose(np.conjugate(ct))
	
	dH1, dH2 = dHcalc(dimH, x1, x2, w1, w2, c, delta)

	cnorm = np.dot(ctbra, ct)

	Efft1 = -np.dot(ctbra, np.dot(dH1, ct))/cnorm
	Eff1 = Efft1.real

	Efft2 = -np.dot(ctbra, np.dot(dH2, ct))/cnorm
	Eff2 = Efft2.real

	return Eff1, Eff2


	
