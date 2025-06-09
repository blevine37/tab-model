#Calculate Ehrenfest force
import numpy as np
from diffH import dHcalc

def calEff(dimH, ndof, ct, x, w1, w2, c, delta):
    # build bra for later
    ctbra = ct.conj().T                     # shape (dimH,)

    # dH has shape (ndof, dimH, dimH)
    dH = dHcalc(dimH, ndof, x, w1, w2, c, delta)

    # normalization scalar
    cnorm = ctbra @ ct                      # same as np.dot(ctbra, ct)

    # 1) first mat‐vec: each slice dH[k] dot ct → shape (ndof, dimH)
    tmp = dH @ ct

    # 2) then dot each of those rows with ctbra → shape (ndof,)
    #    and scale / take real part
    Eff = - (tmp @ ctbra) / cnorm
    return Eff.real


	
