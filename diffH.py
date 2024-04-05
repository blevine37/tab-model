# Calculate dH/dR
import numpy as np

def dHcalc(dimH, x1, A, B, C):
	
	dH1 = np.zeros((dimH, dimH))

	# Define dH1	
	dH1[0,0] = 0
	dH1[1,1] = 0
	
	if x1 < 0:
		dH1[0,1] = B * C * np.exp(C * x1)
		dH1[1,0] = dH1[0,1]
	else:
		dH1[0,1] = B * C * np.exp(-C * x1)
		dH1[1,0] = dH1[0,1]
	
	return dH1
