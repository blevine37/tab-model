import numpy as np


def h_linear_tab(x, params, dimH):
    """Current TAB Hamiltonian form."""
    w1 = params["w1"]
    w2 = params["w2"]
    c = params["c"]
    delta = params["delta"]
    epsil = params["epsil"]

    H = np.zeros((dimH, dimH))
    H[0, 0] = -1.0 * w1 * x[0]

    i = 1
    while i < dimH:
        H[i, i] = w2 * x[0] - (i - 1) * delta
        H[i, 0] = c * x[1]
        H[0, i] = H[i, 0]
        i = i + 1

    i = 5
    while i < dimH:
        H[i, i] = H[i, i] - epsil
        i = i + 1

    return H


def dh_linear_tab(x, params, dimH, ndof):
    """Derivative of current TAB Hamiltonian."""
    w1 = params["w1"]
    w2 = params["w2"]
    c = params["c"]

    dH = np.zeros((ndof, dimH, dimH))
    dH[0, 0, 0] = -1.0 * w1

    i = 1
    while i < dimH:
        dH[0, i, i] = w2
        dH[1, i, 0] = c
        dH[1, 0, i] = c
        i = i + 1

    return dH


def h_linear_tab_general(x, params, dimH):
    """Generalized TAB Hamiltonian: diag_dof and coup_dof set which nuclear DOF
    drives the diagonal and off-diagonal elements respectively."""
    w1 = params["w1"]
    w2 = params["w2"]
    c = params["c"]
    delta = params["delta"]
    epsil = params["epsil"]
    diag_dof = params["diag_dof"]
    coup_dof = params["coup_dof"]

    H = np.zeros((dimH, dimH))
    H[0, 0] = -1.0 * w1 * x[diag_dof]

    i = 1
    while i < dimH:
        H[i, i] = w2 * x[diag_dof] - (i - 1) * delta
        H[i, 0] = c * x[coup_dof]
        H[0, i] = H[i, 0]
        i = i + 1

    i = 5
    while i < dimH:
        H[i, i] = H[i, i] - epsil
        i = i + 1

    return H


def dh_linear_tab_general(x, params, dimH, ndof):
    """Derivative of generalized TAB Hamiltonian."""
    w1 = params["w1"]
    w2 = params["w2"]
    c = params["c"]
    diag_dof = params["diag_dof"]
    coup_dof = params["coup_dof"]

    dH = np.zeros((ndof, dimH, dimH))
    dH[diag_dof, 0, 0] = -1.0 * w1

    i = 1
    while i < dimH:
        dH[diag_dof, i, i] = w2
        dH[coup_dof, i, 0] = c
        dH[coup_dof, 0, i] = c
        i = i + 1

    return dH


# Example placeholder for non-(w1,w2,delta,...) models:
def h_lvc(x, params, dimH):
    w10a = params["w10a"] /27.211
    w6a = params["w6a"] /27.211
    w1 = params["w1"] /27.211
    w9a = params["w9a"] /27.211
    delta = params["delta"] /27.211
    lambda_ = params["lambda"] /27.211
    k6a1 = params["k6a1"] /27.211
    k6a2 = params["k6a2"] /27.211
    k11 = params["k11"] /27.211
    k12 = params["k12"] /27.211
    k9a1 = params["k9a1"] /27.211
    k9a2 = params["k9a2"] /27.211

    H = np.zeros((dimH, dimH))
    H[0, 0] = 0.5 * (w10a * (x[0] ** 2) + w6a*(x[1] ** 2) +w1 * (x[2] ** 2) + w9a * (x[3] ** 2)) -delta + k6a1 * x[1] + k11 * x[2] + k9a1 * x[3]
    H[1, 1] = 0.5 * (w10a * (x[0] ** 2) + w6a * (x[1] ** 2) + w1 * (x[2] ** 2) + w9a * (x[3] ** 2)) +delta + k6a2 * x[1] + k12 * x[2] + k9a2 * x[3]
    H[0, 1] = lambda_ * x[0]
    H[1, 0] = H[0, 1]
    return H


def dh_lvc(x, params, dimH, ndof):

    w10a = params["w10a"] /27.211
    w6a = params["w6a"] /27.211
    w1 = params["w1"] /27.211
    w9a = params["w9a"] /27.211
    delta = params["delta"] /27.211
    lambda_ = params["lambda"] /27.211
    k6a1 = params["k6a1"] /27.211
    k6a2 = params["k6a2"] /27.211
    k11 = params["k11"] /27.211
    k12 = params["k12"] /27.211
    k9a1 = params["k9a1"] /27.211
    k9a2 = params["k9a2"] /27.211

    dH = np.zeros((ndof, dimH, dimH))
    dH[0, 0, 0] = w10a * x[0] 
    dH[0, 1, 1] = dH[0, 0, 0]
    dH[0, 0, 1] = lambda_
    dH[0, 1, 0] = dH[0, 0, 1]
    dH[1, 0, 0] = w6a * x[1] + k6a1
    dH[1, 1, 1] = w6a * x[1] + k6a2
    dH[2, 0, 0] = w1 * x[2] + k11
    dH[2, 1, 1] = w1 * x[2] + k12
    dH[3, 0, 0] = w9a * x[3] + k9a1
    dH[3, 1, 1] = w9a * x[3] + k9a2

    return dH


def h_lvc_general(x, params, dimH):
    """Generalized matrix-based Hamiltonian builder for N-state LVC models.

    Param shapes:
      w:   (ndof,)
      E0:  (dimH,)
      kl:  (dimH, ndof)
      kq:  (dimH, ndof)
      kb:  (dimH, ndof, ndof)
      lam: (dimH, dimH, ndof)         symmetric in state indices
      kbo: (dimH, dimH, ndof, ndof)   symmetric in state indices
    """
    w = params["w"]
    E0 = params["E0"]
    kl = params["kl"]
    kq = params["kq"]
    kb = params["kb"]
    lam = params["lam"]
    kbo = params["kbo"]

    H = np.zeros((dimH, dimH))
    ho = np.sum(0.5 * w * x**2)

    for state in range(dimH):
        H[state, state] = ho + E0[state] + np.sum(kl[state] * x) + np.sum(kq[state] * x**2) + np.einsum('ij,i,j->', kb[state], x, x)

    for i in range(dimH):
        for j in range(i + 1, dimH):
            off_diag = np.sum(lam[i, j] * x) + np.einsum('kl,k,l->', kbo[i, j], x, x)
            H[i, j] = off_diag
            H[j, i] = off_diag
    return H


def dh_lvc_general(x, params, dimH, ndof):
    """Generalized derivative of N-state LVC Hamiltonian."""
    w = params["w"]
    kl = params["kl"]
    kq = params["kq"]
    kb = params["kb"]
    lam = params["lam"]
    kbo = params["kbo"]

    dH = np.zeros((ndof, dimH, dimH))

    for i in range(ndof):
        ho_grad = w[i] * x[i]
        for state in range(dimH):
            grad_kb = 2.0 * np.dot(kb[state, i], x)
            dH[i, state, state] = ho_grad + kl[state, i] + 2.0 * kq[state, i] * x[i] + grad_kb

        for s1 in range(dimH):
            for s2 in range(s1 + 1, dimH):
                grad_off_diag = lam[s1, s2, i] + 2.0 * np.dot(kbo[s1, s2, i], x)
                dH[i, s1, s2] = grad_off_diag
                dH[i, s2, s1] = grad_off_diag
    return dH


def build_pyr24_params():
    """Build parameter tensors for Pyrazine 24-mode bi-linear LVC model"""
    ev_to_au = 1.0 / 27.211
    ndof = 24
    
    w = np.array([0.1139, 0.0739, 0.1258, 0.1525, 0.1961, 0.3788, 0.0937, 0.1219, 0.0873, 0.1669, 0.1891, 0.3769, 0.0423, 0.1190, 0.1266, 0.1408, 0.1840, 0.3734, 0.1318, 0.1425, 0.1756, 0.3798, 0.0521, 0.0973])
    delta = 0.42300
    
    kl = np.zeros((2, ndof))
    # S1
    kl[0, 1] = 0.09806
    kl[0, 2] = 0.05033
    kl[0, 3] = 0.14521
    kl[0, 4] = -0.04448
    kl[0, 5] = -0.02473
    # S2
    kl[1, 1] = -0.13545
    kl[1, 2] = 0.17100
    kl[1, 3] = 0.03746
    kl[1, 4] = 0.01677
    kl[1, 5] = -0.01619
    
    kq = np.zeros((2, ndof))
    qk1_vals = {0: -0.01159, 6: -0.02252, 7: -0.01825, 8: -0.00741, 9: 0.05183, 10: -0.05733, 11: -0.00333, 12: 0.01145, 13: -0.02040, 14: -0.04819, 15: -0.00792, 16: -0.02429, 17: -0.00492, 18: -0.00277, 19: 0.03924, 20: 0.00992, 21: -0.00110, 22: -0.02176, 23: 0.00315}
    for k,v in qk1_vals.items(): kq[0,k] = v
        
    qk2_vals = {0: -0.01159, 6: -0.03445, 7: -0.00265, 8: -0.00385, 9: 0.04842, 10: -0.06332, 11: -0.00040, 12: -0.01459, 13: -0.00618, 14: -0.00840, 15: 0.00429, 16: -0.00734, 17: 0.00062, 18: -0.01179, 19: 0.04000, 20: 0.01246, 21: 0.00069, 22: -0.02214, 23: -0.00496}
    for k,v in qk2_vals.items(): kq[1,k] = v
        
    kb = np.zeros((2, ndof, ndof))
    def add_kb(s, i, j, val):
        kb[s, i, j] = val
        kb[s, j, i] = val 
        
    add_kb(0, 1, 2, 0.00108)
    add_kb(0, 1, 3, 0.00204)
    add_kb(0, 1, 4, 0.00135)
    add_kb(0, 1, 5, -0.00285)
    add_kb(0, 2, 3, -0.00474)
    add_kb(0, 2, 4, -0.00154)
    add_kb(0, 2, 5, -0.00163)
    add_kb(0, 3, 4, 0.00872)
    add_kb(0, 3, 5, 0.00474)
    add_kb(0, 4, 5, 0.00143)
    add_kb(0, 6, 7, -0.00049)
    add_kb(0, 8, 9, 0.01321)
    add_kb(0, 8, 10, -0.00717)
    add_kb(0, 8, 11, 0.00515)
    add_kb(0, 9, 10, -0.03942)
    add_kb(0, 9, 11, 0.00170)
    add_kb(0, 10, 11, -0.00204)
    add_kb(0, 12, 13, 0.00100)
    add_kb(0, 14, 15, 0.00525)
    add_kb(0, 14, 16, -0.00485)
    add_kb(0, 14, 17, -0.00326)
    add_kb(0, 15, 16, 0.00852)
    add_kb(0, 15, 17, 0.00888)
    add_kb(0, 16, 17, -0.00443)
    add_kb(0, 18, 19, 0.00016)
    add_kb(0, 18, 20, -0.00250)
    add_kb(0, 18, 21, 0.00357)
    add_kb(0, 19, 20, -0.00197)
    add_kb(0, 19, 21, -0.00355)
    add_kb(0, 20, 21, 0.00623)
    add_kb(0, 22, 23, -0.00624)
    
    add_kb(1, 1, 2, -0.00298)
    add_kb(1, 1, 3, 0.00189)
    add_kb(1, 1, 4, 0.00203)
    add_kb(1, 1, 5, -0.00128)
    add_kb(1, 2, 3, -0.00155)
    add_kb(1, 2, 4, -0.00311)
    add_kb(1, 2, 5, -0.00600)
    add_kb(1, 3, 4, 0.01194)
    add_kb(1, 3, 5, 0.00334)
    add_kb(1, 4, 5, 0.00713)
    add_kb(1, 6, 7, 0.00911)
    add_kb(1, 8, 9, -0.00661)
    add_kb(1, 8, 10, 0.00429)
    add_kb(1, 8, 11, -0.00246)
    add_kb(1, 9, 10, -0.03034)
    add_kb(1, 9, 11, -0.00185)
    add_kb(1, 10, 11, -0.00388)
    add_kb(1, 12, 13, -0.00091)
    add_kb(1, 14, 15, 0.00536)
    add_kb(1, 14, 16, -0.00097)
    add_kb(1, 14, 17, 0.00034)
    add_kb(1, 15, 16, 0.00209)
    add_kb(1, 15, 17, -0.00049)
    add_kb(1, 16, 17, 0.00346)
    add_kb(1, 18, 19, -0.00884)
    add_kb(1, 18, 20, 0.07000)
    add_kb(1, 18, 21, -0.01249)
    add_kb(1, 19, 20, -0.05000)
    add_kb(1, 19, 21, 0.00265)
    add_kb(1, 20, 21, -0.00422)
    add_kb(1, 22, 23, -0.00261)
    
    lam = np.zeros((2, 2, ndof))
    lam[0, 1, 0] = lam[1, 0, 0] = 0.20804

    kbo = np.zeros((2, 2, ndof, ndof))
    def add_kbo(i, j, val):
        kbo[0, 1, i, j] = kbo[0, 1, j, i] = val
        kbo[1, 0, i, j] = kbo[1, 0, j, i] = val
        
    add_kbo(0, 1, 0.01000)
    add_kbo(0, 2, 0.00553)
    add_kbo(0, 3, 0.00126)
    add_kbo(0, 4, 0.00799)
    add_kbo(0, 5, 0.00514)
    add_kbo(6, 8, -0.01372)
    add_kbo(6, 9, -0.00466)
    add_kbo(6, 10, 0.00329)
    add_kbo(6, 11, -0.00031)
    add_kbo(7, 8, 0.00598)
    add_kbo(7, 9, -0.00914)
    add_kbo(7, 10, 0.00961)
    add_kbo(7, 11, 0.00500)
    add_kbo(12, 14, -0.01056)
    add_kbo(12, 15, 0.00559)
    add_kbo(12, 16, 0.00401)
    add_kbo(12, 17, -0.00226)
    add_kbo(13, 14, -0.01200)
    add_kbo(13, 15, -0.00213)
    add_kbo(13, 16, 0.00328)
    add_kbo(13, 17, -0.00396)
    add_kbo(18, 22, 0.00118)
    add_kbo(19, 22, -0.00009)
    add_kbo(20, 22, -0.00285)
    add_kbo(21, 22, -0.00095)
    add_kbo(18, 23, 0.01281)
    add_kbo(19, 23, -0.01780)
    add_kbo(20, 23, 0.00134)
    add_kbo(21, 23, -0.00481)

    return {
        "w": w * ev_to_au,
        "E0": np.array([-delta, delta]) * ev_to_au,
        "kl": kl * ev_to_au,
        "kq": kq * ev_to_au,
        "kb": kb * ev_to_au,
        "lam": lam * ev_to_au,
        "kbo": kbo * ev_to_au
    }


def build_pyr4bath_params():
    """Build parameter tensors for Pyrazine 24-mode (4 system + 20 bath) linear LVC model."""
    ev_to_au = 1.0 / 27.211
    ndof = 24
    
    # Frequencies (v10a, v6a, v1, v9a, 1b...20b)
    w = np.array([
        0.09357, 0.0740, 0.1273, 0.1568,
        0.0400, 0.0589, 0.0778, 0.0968, 0.1157, 0.1347, 0.1536, 0.1726, 0.1915, 0.2105,
        0.2294, 0.2484, 0.2673, 0.2863, 0.3052, 0.3242, 0.3431, 0.3621, 0.3810, 0.4000
    ])
    
    delta = 0.46165
    
    kl = np.zeros((2, ndof))
    # State S1 linear coupling constants (k_i)
    kl[0] = np.array([
        0.0, -0.0964, 0.0470, 0.1594,
        0.0069, 0.0112, 0.0102, 0.0188, 0.0261, 0.0308, 0.0210, 0.0265, 0.0196, 0.0281,
        0.0284, 0.0361, 0.0560, 0.0433, 0.0625, 0.0717, 0.0782, 0.0780, 0.0269, 0.0306
    ])
    
    # State S2 linear coupling constants (k_i)
    kl[1] = np.array([
        0.0, 0.1194, 0.2012, 0.0484,
        -0.0069, -0.0112, -0.0102, -0.0188, -0.0261, -0.0308, -0.0210, -0.0265, -0.0196, -0.0281,
        -0.0284, -0.0361, -0.0560, -0.0433, -0.0625, -0.0717, -0.0782, -0.0780, -0.0269, -0.0306
    ])
    
    # Off-diagonal linear coupling
    lam = np.zeros((2, 2, ndof))
    lam[0, 1, 0] = lam[1, 0, 0] = 0.1825 # lambda on v10a

    # Zeros for higher order terms not present in the model
    kq = np.zeros((2, ndof))
    kb = np.zeros((2, ndof, ndof))
    kbo = np.zeros((2, 2, ndof, ndof))

    return {
        "w": w * ev_to_au,
        "E0": np.array([-delta, delta]) * ev_to_au,
        "kl": kl * ev_to_au,
        "kq": kq * ev_to_au,
        "kb": kb * ev_to_au,
        "lam": lam * ev_to_au,
        "kbo": kbo * ev_to_au
    }


def build_fulvene_params():
    """Build parameter tensors for Fulvene 30-mode 2-state LVC model.

    Source: MCTDH operator file (Quantics package, vcham).
    All energies in eV; converted to a.u. on return.

    Conventions:
      - Modes v1..v30 -> indices 0..29.
      - States S1, S2 -> indices 0, 1 (E1=0, E2=4.16142 eV).
      - kappa{s}_{m} = on-diagonal linear coupling on state s, mode m.
      - lambda1_2_{m} = off-diagonal linear coupling on mode m.
      - All quadratic, bilinear, and cubic sections of the op file are empty.
    """
    ev_to_au = 1.0 / 27.211
    ndof = 30

    # Frequencies (eV) for v1 .. v30
    w = np.array([
        0.02610, 0.04530, 0.06174, 0.07760, 0.08560, 0.08750, 0.09489, 0.09766, 0.10603, 0.11074,
        0.11154, 0.11242, 0.11739, 0.12811, 0.12909, 0.14578, 0.14629, 0.16851, 0.18011, 0.18285,
        0.19361, 0.20310, 0.20964, 0.22047, 0.41409, 0.42437, 0.42759, 0.42790, 0.43081, 0.43248
    ])

    # State energies (eV)
    E0 = np.array([0.00000, 4.16142])

    # On-diagonal linear couplings (kappa).  kappa{state}_{mode} -> kl[state-1, mode-1].
    kl = np.zeros((2, ndof))
    # State 1 (index 0)
    kl[0,  5] = -0.00799   # kappa1_6
    kl[0, 12] =  0.00608   # kappa1_13
    kl[0, 14] =  0.00487   # kappa1_15
    kl[0, 16] = -0.01463   # kappa1_17
    kl[0, 19] = -0.01265   # kappa1_20
    kl[0, 20] = -0.00096   # kappa1_21
    kl[0, 21] = -0.00421   # kappa1_22
    kl[0, 23] = -0.00090   # kappa1_24
    kl[0, 24] =  0.00464   # kappa1_25
    kl[0, 27] =  0.01287   # kappa1_28
    kl[0, 29] =  0.05083   # kappa1_30
    # State 2 (index 1)
    kl[1,  5] = -0.16463   # kappa2_6
    kl[1, 12] =  0.16902   # kappa2_13
    kl[1, 14] =  0.20870   # kappa2_15
    kl[1, 16] = -0.19129   # kappa2_17
    kl[1, 19] =  0.12386   # kappa2_20
    kl[1, 20] = -0.06968   # kappa2_21
    kl[1, 21] = -0.38302   # kappa2_22
    kl[1, 23] =  0.49293   # kappa2_24
    kl[1, 24] = -0.05766   # kappa2_25
    kl[1, 27] =  0.01640   # kappa2_28
    kl[1, 29] =  0.03264   # kappa2_30

    # Off-diagonal linear couplings (lambda).  lambda1_2_{mode} -> lam[0, 1, mode-1].
    lam = np.zeros((2, 2, ndof))
    lam_12 = np.zeros(ndof)
    lam_12[ 1] = -0.02235   # lambda1_2_2
    lam_12[ 8] = -0.14659   # lambda1_2_9
    lam_12[13] =  0.04999   # lambda1_2_14
    lam_12[15] =  0.02678   # lambda1_2_16
    lam_12[17] = -0.07546   # lambda1_2_18
    lam_12[18] = -0.14017   # lambda1_2_19
    lam_12[22] =  0.04228   # lambda1_2_23
    lam_12[25] =  0.00122   # lambda1_2_26
    lam_12[26] =  0.00416   # lambda1_2_27
    lam_12[28] = -0.00788   # lambda1_2_29
    lam[0, 1] = lam_12
    lam[1, 0] = lam_12

    # Higher-order terms not present in the op file
    kq = np.zeros((2, ndof))
    kb = np.zeros((2, ndof, ndof))
    kbo = np.zeros((2, 2, ndof, ndof))

    return {
        "w":   w   * ev_to_au,
        "E0":  E0  * ev_to_au,
        "kl":  kl  * ev_to_au,
        "kq":  kq  * ev_to_au,
        "kb":  kb  * ev_to_au,
        "lam": lam * ev_to_au,
        "kbo": kbo * ev_to_au,
    }


def build_dmabn_params():
    """Build parameter tensors for DMABN 57-mode 3-state LVC model.

    Source: MCTDH operator file (Quantics package, vcham).
    States: S1 (E1=0.0), S2 (E2=4.93397 eV), S3 (E3=5.31142 eV).
    Only kappa (on-diagonal linear) and lambda (off-diagonal linear) terms
    are present; quadratic, bilinear, and cubic sections are empty.

    Conventions:
      - Modes v1..v57 -> indices 0..56.
      - States S1, S2, S3 -> indices 0, 1, 2.
      - kappa{s}_{m} -> kl[s-1, m-1].
      - lambda{a}_{b}_{m} -> lam[a-1, b-1, m-1] (symmetric).
    """
    ev_to_au = 1.0 / 27.211
    ndof = 57
    dimH = 3

    # Frequencies (eV) for v1..v57
    w = _dmabn_omega_ev.copy()

    # State energies (eV)
    E0 = np.array([0.00000, 4.93397, 5.31142])

    # On-diagonal linear couplings (kappa).
    kl = np.zeros((dimH, ndof))
    kappa1 = {2: -0.00036, 3: -0.00082, 4: -0.00045, 6: 0.00032, 7: 0.00037,
              8: -0.00043, 10: -0.00049, 11: -0.00071, 12: -0.00038,
              14: 0.00018, 16: -0.00019, 18: -0.00044, 21: 0.00033,
              22: -0.00038, 24: -0.00026, 25: -0.00015, 28: 0.00014,
              30: -0.00033, 31: 0.00046, 35: -0.00066, 36: -0.00025,
              40: 0.00036, 42: 0.00014, 45: 0.00012, 47: -0.00014,
              48: 0.00011, 52: -0.00038, 53: -0.00018, 55: -0.00028,
              56: -0.00027}
    for m, v in kappa1.items(): kl[0, m] = v

    kappa2 = {0: 0.01981, 1: 0.00114, 2: 0.00679, 3: -0.00129, 4: -0.01080,
              5: -0.00051, 6: 0.00074, 7: 0.01364, 8: 0.01326, 9: 0.00025,
              10: -0.00067, 11: -0.02516, 12: -0.01572, 13: 0.00029,
              14: 0.00120, 15: -0.00017, 16: 0.00630, 17: 0.00095,
              18: -0.09430, 20: -0.00292, 21: -0.03331, 22: 0.04050,
              24: -0.00984, 25: -0.00036, 26: -0.01563, 27: 0.00299,
              28: -0.00022, 29: -0.01587, 30: 0.00569, 31: 0.07393,
              32: 0.00135, 33: -0.00099, 34: 0.00021, 35: 0.07246,
              36: 0.00079, 37: -0.00031, 38: 0.00674, 39: 0.01043,
              40: -0.00021, 41: -0.00014, 42: 0.02299, 43: 0.00436,
              44: -0.00038, 45: 0.10678, 46: -0.02123, 47: -0.00070,
              48: -0.01697, 49: 0.00086, 50: 0.00237, 51: -0.00024,
              52: -0.00153, 53: -0.02192, 54: 0.00756, 55: -0.00137,
              56: 0.03034}
    for m, v in kappa2.items(): kl[1, m] = v

    kappa3 = {0: 0.00797, 1: -0.00096, 2: -0.00024, 3: -0.00094, 4: -0.00628,
              5: -0.00062, 6: 0.00110, 7: 0.00920, 8: -0.00695, 9: -0.00019,
              10: -0.00010, 11: -0.01937, 12: -0.01093, 13: -0.00014,
              14: -0.00166, 16: -0.00475, 17: -0.00549, 18: -0.07216,
              19: 0.00024, 20: 0.00189, 21: -0.03457, 22: 0.04590,
              24: -0.01089, 25: -0.00013, 26: -0.01695, 27: 0.00394,
              28: 0.00141, 29: -0.02115, 30: -0.07638, 31: 0.04662,
              32: 0.00065, 33: 0.00011, 35: -0.00036, 36: -0.00018,
              37: -0.00029, 38: 0.00487, 39: -0.00584, 40: 0.00064,
              42: -0.01030, 43: 0.06696, 44: 0.00062, 45: -0.12549,
              46: -0.10953, 47: 0.00058, 48: 0.02129, 49: 0.00104,
              50: 0.00253, 51: -0.00031, 52: 0.00905, 53: -0.01492,
              54: 0.00526, 55: -0.00034, 56: 0.00257}
    for m, v in kappa3.items(): kl[2, m] = v

    # Off-diagonal linear couplings (lambda), symmetric in state indices.
    lam = np.zeros((dimH, dimH, ndof))

    lam_12_vals = {0: 0.00023, 1: 0.00716, 2: -0.00237, 3: 0.00655,
                   4: -0.00235, 5: 0.02356, 6: -0.01311, 7: -0.00098,
                   8: 0.00015, 9: -0.02099, 10: 0.04703, 11: -0.00047,
                   12: -0.00013, 13: -0.05253, 14: -0.00131, 15: 0.01382,
                   16: 0.00073, 17: -0.00052, 19: -0.00319, 22: -0.00066,
                   23: 0.00046, 24: -0.00050, 25: -0.03919, 26: 0.00331,
                   27: 0.01165, 28: -0.07418, 29: -0.00155, 30: 0.00043,
                   31: 0.00221, 32: -0.23453, 33: -0.26477, 34: 0.33483,
                   35: -0.00165, 36: 0.05187, 37: 0.02173, 38: -0.00409,
                   39: 0.01359, 40: 0.13492, 41: 0.06252, 42: -0.00181,
                   43: 0.00021, 44: 0.03665, 45: -0.00047, 46: 0.00018,
                   47: -0.01867, 48: 0.00063, 49: 0.00470, 50: -0.00202,
                   51: -0.00113, 53: -0.00052, 54: -0.00103, 55: 0.00272,
                   56: 0.00041}
    for m, v in lam_12_vals.items():
        lam[0, 1, m] = lam[1, 0, m] = v

    lam_13_vals = {0: -0.06569, 2: -0.02458, 4: 0.02188, 5: 0.00162,
                   6: -0.00029, 7: -0.02782, 8: -0.02332, 9: 0.00032,
                   10: -0.00013, 11: -0.00192, 12: 0.01311, 13: 0.00024,
                   14: -0.00651, 15: 0.00185, 16: -0.05622, 17: -0.00342,
                   18: 0.01769, 19: -0.00020, 20: 0.00219, 21: -0.02000,
                   22: 0.02829, 24: -0.08859, 25: -0.00015, 26: 0.03954,
                   27: -0.00879, 28: -0.00048, 29: 0.02818, 30: 0.08258,
                   31: -0.04090, 32: 0.00041, 33: 0.00041, 34: -0.00069,
                   35: -0.13269, 36: -0.00241, 37: 0.00037, 38: -0.01557,
                   39: 0.04191, 40: -0.00465, 41: 0.00085, 42: -0.04983,
                   43: 0.09946, 45: 0.16533, 46: 0.05411, 47: 0.00277,
                   48: 0.07859, 49: -0.00343, 50: -0.00810, 51: -0.00015,
                   52: 0.01106, 53: -0.00814, 54: 0.00304, 56: -0.00324}
    for m, v in lam_13_vals.items():
        lam[0, 2, m] = lam[2, 0, m] = v

    lam_23_vals = {1: 0.00195, 3: 0.00617, 4: -0.00123, 5: 0.01004,
                   6: -0.00077, 8: 0.00014, 9: -0.00338, 10: 0.01277,
                   11: 0.00013, 13: 0.02641, 14: 0.00099, 15: 0.00530,
                   16: 0.00031, 19: -0.00208, 23: -0.00278, 24: -0.00030,
                   25: -0.01234, 26: 0.00099, 27: 0.00252, 28: -0.01659,
                   29: -0.00026, 30: -0.00020, 31: 0.00057, 32: -0.05087,
                   33: -0.00962, 34: 0.03118, 35: -0.00021, 36: 0.01499,
                   37: 0.00439, 38: -0.00069, 39: 0.00299, 40: 0.02879,
                   41: 0.02474, 42: -0.00032, 43: 0.00058, 44: -0.09530,
                   45: -0.00035, 46: -0.00014, 47: -0.00565, 48: 0.00019,
                   49: 0.00185, 50: -0.00073, 51: -0.00199, 53: 0.00635,
                   54: 0.01709, 55: -0.00350}
    for m, v in lam_23_vals.items():
        lam[1, 2, m] = lam[2, 1, m] = v

    # Higher-order terms not present in the op file
    kq = np.zeros((dimH, ndof))
    kb = np.zeros((dimH, ndof, ndof))
    kbo = np.zeros((dimH, dimH, ndof, ndof))

    return {
        "w":   w   * ev_to_au,
        "E0":  E0  * ev_to_au,
        "kl":  kl  * ev_to_au,
        "kq":  kq  * ev_to_au,
        "kb":  kb  * ev_to_au,
        "lam": lam * ev_to_au,
        "kbo": kbo * ev_to_au,
    }


def _lvc_ic(pmass):
    """Wigner-like ICs for LVC models: zero mean, sqrt(0.5) width in mass-weighted coords."""
    ndof = len(pmass)
    return {
        "x_mean": [0.0] * ndof,
        "x_std":  [np.sqrt(0.5)] * ndof,
        "v_mean": [0.0] * ndof,
        "v_std":  [np.sqrt(0.5) / m for m in pmass],
    }


def _tab_ic(ndof, pmass):
    """Standard linear-TAB ICs: particle starts at x[0]=-1 moving toward crossing."""
    return {
        "x_mean": [-1.0] + [0.0] * (ndof - 1),
        "x_std":  [0.204] * ndof,
        "v_mean": [10.0 / pmass[i] for i in range(ndof)],
        "v_std":  [2.451 / pmass[i] for i in range(ndof)],
    }


_DCP_LVC = 0.5
_DCP_TAB = 6.0


_pyr4bath_pmass = list(27.211 / np.array([
    0.09357, 0.0740, 0.1273, 0.1568,
    0.0400, 0.0589, 0.0778, 0.0968, 0.1157, 0.1347, 0.1536, 0.1726, 0.1915, 0.2105,
    0.2294, 0.2484, 0.2673, 0.2863, 0.3052, 0.3242, 0.3431, 0.3621, 0.3810, 0.4000
]))

_pyr24_pmass = [238.9025, 368.2138, 216.3037, 178.4328, 138.7608, 71.8347, 290.4055,
                223.2240, 311.6953, 163.0377, 143.8974, 72.1969, 643.2861, 228.6639,
                214.9368, 193.2599, 147.8859, 72.8736, 206.4568, 190.9544, 154.9601,
                71.6456, 522.2841, 279.6608]

_fulvene_omega_ev = np.array([
    0.02610, 0.04530, 0.06174, 0.07760, 0.08560, 0.08750, 0.09489, 0.09766, 0.10603, 0.11074,
    0.11154, 0.11242, 0.11739, 0.12811, 0.12909, 0.14578, 0.14629, 0.16851, 0.18011, 0.18285,
    0.19361, 0.20310, 0.20964, 0.22047, 0.41409, 0.42437, 0.42759, 0.42790, 0.43081, 0.43248
])
_fulvene_pmass = list(27.211 / _fulvene_omega_ev)  # 1/omega in a.u. (mass-weighted normal coords)

_dmabn_omega_ev = np.array([
    0.00622, 0.00909, 0.01033, 0.01673, 0.02188, 0.02317, 0.03186, 0.03575, 0.04213, 0.05324,
    0.05919, 0.06109, 0.06234, 0.07081, 0.07103, 0.08237, 0.08386, 0.09314, 0.10163, 0.10377,
    0.10575, 0.12272, 0.12294, 0.12408, 0.12707, 0.13479, 0.14102, 0.14110, 0.14298, 0.14786,
    0.14935, 0.15535, 0.16070, 0.16529, 0.17258, 0.17568, 0.17782, 0.18109, 0.18219, 0.18324,
    0.18354, 0.18661, 0.18749, 0.19708, 0.20350, 0.21190, 0.29621, 0.37412, 0.37509, 0.38359,
    0.38373, 0.39291, 0.39407, 0.39931, 0.39943, 0.40299, 0.40309
])
_dmabn_pmass = list(27.211 / _dmabn_omega_ev)  # 1/omega in a.u. (mass-weighted normal coords)

_pyr4_pmass  = [290.809, 367.716, 213.755, 173.539]
_tab2_pmass  = [1845.0, 1845.0]
_tab3_pmass  = [1845.0, 1845.0, 1845.0]

MODELS = {
    "pyr4bath": {
        "dimH": 2,
        "ndof": 24,
        "pmass": _pyr4bath_pmass,
        "params": build_pyr4bath_params(),
        "h_builder": h_lvc_general,
        "dh_builder": dh_lvc_general,
        "ic": _lvc_ic(_pyr4bath_pmass),
        "dcp": [_DCP_LVC] * 24,
        "init_state": 1,
    },
    "fulvene": {
        "dimH": 2,
        "ndof": 30,
        "pmass": _fulvene_pmass,
        "params": build_fulvene_params(),
        "h_builder": h_lvc_general,
        "dh_builder": dh_lvc_general,
        "ic": _lvc_ic(_fulvene_pmass),
        "dcp": [_DCP_LVC] * 30,
        "init_state": 1,
    },
    "dmabn": {
        "dimH": 3,
        "ndof": 57,
        "pmass": _dmabn_pmass,
        "params": build_dmabn_params(),
        "h_builder": h_lvc_general,
        "dh_builder": dh_lvc_general,
        "ic": _lvc_ic(_dmabn_pmass),
        "dcp": [_DCP_LVC] * 57,
        "init_state": 2,  # Start in S3 (idx 2)
    },
    "pyr24": {
        "dimH": 2,
        "ndof": 24,
        "pmass": _pyr24_pmass,
        "params": build_pyr24_params(),
        "h_builder": h_lvc_general,
        "dh_builder": dh_lvc_general,
        "ic": _lvc_ic(_pyr24_pmass),
        "dcp": [_DCP_LVC] * 24,
        "init_state": 1,
    },
    "pyr4": {
        "dimH": 2,
        "ndof": 4,
        "pmass": _pyr4_pmass,
        "params": {"w10a": 0.09357, "w6a": 0.0740, "w1": 0.1273, "w9a": 0.1568,
                   "delta": 0.46165, "lambda": 0.1825, "k6a1": -0.0964, "k6a2": 0.1194,
                   "k11": 0.0470, "k12": 0.2012, "k9a1": 0.1594, "k9a2": 0.0484},  # all in eV
        "h_builder": h_lvc,
        "dh_builder": dh_lvc,
        "ic": _lvc_ic(_pyr4_pmass),
        "dcp": [_DCP_LVC] * 4,
        "init_state": 1,
    },
    "m3_10000_3d": {
        "dimH": 3,
        "ndof": 3,
        "pmass": _tab3_pmass,
        "params": {"delta": 0.01, "w1": 0.25, "w2": 0.025, "c": 0.025, "epsil": 0.0},
        "h_builder": h_linear_tab,
        "dh_builder": dh_linear_tab,
        "ic": _tab_ic(3, _tab3_pmass),
        "dcp": [_DCP_TAB] * 3,
        "init_state": 0,
    },
    "m3_x2coup": {
        # Coupling along x[2], diagonal driven by x[0], x[1] is spectator.
        # x[1] is given large initial KE (spectator reservoir) while x[2] (the
        # coupled direction) has small initial KE.  This makes p-rescaling
        # unphysical: it satisfies energy conservation by drawing from x[1] even
        # though that DOF does not participate in the hop.  h/gh correctly
        # restrict rescaling to x[2] and will encounter frustrated hops more
        # often, producing different population dynamics.
        "dimH": 3,
        "ndof": 3,
        "pmass": _tab3_pmass,
        "params": {"delta": 0.01, "w1": 0.25, "w2": 0.025, "c": 0.025, "epsil": 0.0,
                   "diag_dof": 0, "coup_dof": 1},
        "h_builder": h_linear_tab_general,
        "dh_builder": dh_linear_tab_general,
        "ic": {
            "x_mean": [-1.0, 0.0, 0.0],
            "x_std":  [0.204, 0.204, 0.204],
            # x[0]: reaction coordinate — same as standard linear TAB
            # x[1]: spectator with large KE reservoir
            # x[2]: coupled direction — small KE so hops along it are often frustrated
            "v_mean": [0.0 / _tab3_pmass[0], 0.0, 20.0],
            "v_std":  [2.451 / _tab3_pmass[0],
                       2.451  / _tab3_pmass[1],
                       2.451   / _tab3_pmass[2]],
        },
        "dcp": [_DCP_TAB] * 3,
        "init_state": 0,
    },
    "m9_5000": {
        "dimH": 9,
        "ndof": 2,
        "pmass": _tab2_pmass,
        "params": {"delta": 0.005, "w1": 0.25, "w2": 0.025, "c": 0.025, "epsil": 0.0},
        "h_builder": h_linear_tab,
        "dh_builder": dh_linear_tab,
        "ic": _tab_ic(2, _tab2_pmass),
        "dcp": [_DCP_TAB] * 2,
        "init_state": 0,
    },
    "m9_500_split": {
        "dimH": 9,
        "ndof": 2,
        "pmass": _tab2_pmass,
        "params": {"delta": 0.0005, "w1": 0.25, "w2": 0.025, "c": 0.025, "epsil": 0.08},
        "h_builder": h_linear_tab,
        "dh_builder": dh_linear_tab,
        "ic": _tab_ic(2, _tab2_pmass),
        "dcp": [_DCP_TAB] * 2,
        "init_state": 0,
    },
    "m9_10000_split_steep": {
        "dimH": 9,
        "ndof": 2,
        "pmass": _tab2_pmass,
        "params": {"delta": 0.01, "w1": 0.25, "w2": 0.25, "c": 0.025, "epsil": 0.08},
        "h_builder": h_linear_tab,
        "dh_builder": dh_linear_tab,
        "ic": _tab_ic(2, _tab2_pmass),
        "dcp": [_DCP_TAB] * 2,
        "init_state": 0,
    },
    "m3_10000": {
        "dimH": 3,
        "ndof": 2,
        "pmass": _tab2_pmass,
        "params": {"delta": 0.01, "w1": 0.25, "w2": 0.025, "c": 0.025, "epsil": 0.0},
        "h_builder": h_linear_tab,
        "dh_builder": dh_linear_tab,
        "ic": _tab_ic(2, _tab2_pmass),
        "dcp": [_DCP_TAB] * 2,
        "init_state": 0,
    },
    "m9_10000_split": {
        "dimH": 9,
        "ndof": 2,
        "pmass": _tab2_pmass,
        "params": {"delta": 0.01, "w1": 0.25, "w2": 0.025, "c": 0.025, "epsil": 0.08},
        "h_builder": h_linear_tab,
        "dh_builder": dh_linear_tab,
        "ic": _tab_ic(2, _tab2_pmass),
        "dcp": [_DCP_TAB] * 2,
        "init_state": 0,
    },
    "z_example": {
        "dimH": 3,
        "ndof": 3,
        "pmass": _tab3_pmass,
        "params": {"delta": 0.1, "w1": 0.25, "w2": 0.025, "c": 0.025, "epsil": 0.0},
        "h_builder": h_linear_tab,
        "dh_builder": dh_linear_tab,
        "ic": {
            "x_mean": [-1.0, 0.0, 0.0],
            "x_std":  [0.204, 0.204, 0.204],
            # x[0]: reaction coordinate — same as standard linear TAB
            # x[1]: spectator with large KE reservoir
            # x[2]: coupled direction — small KE so hops along it are often frustrated
            "v_mean": [10.0 / _tab3_pmass[0], 10.0/_tab3_pmass[1], 50.0/_tab3_pmass[2]],
            "v_std":  [2.451 / _tab3_pmass[0],
                       2.451  / _tab3_pmass[1],
                       2.451   / _tab3_pmass[2]],
        },
        "dcp": [_DCP_TAB] * 3,
        "init_state": 0,
    },
}


def get_model(model_name):
    if model_name not in MODELS:
        valid = ", ".join(sorted(MODELS.keys()))
        raise ValueError(f"Unknown model '{model_name}'. Valid models: {valid}")
    return MODELS[model_name]