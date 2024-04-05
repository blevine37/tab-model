# Calculate Ehrenfest force
import numpy as np
from diffH import dHcalc

def calEff(dimH, ct, x1, A, B, C):
	
	ctbra = np.transpose(np.conjugate(ct))
	
	dH1= dHcalc(dimH, x1, A, B, C)

	cnorm = np.dot(ctbra, ct)

	Efft1 = -np.dot(ctbra, np.dot(dH1, ct))/cnorm
	Eff1 = Efft1.real

	return Eff1


	
