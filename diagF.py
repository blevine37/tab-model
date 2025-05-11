def Force_diag(dH1, dH2, dimH):
    import numpy as np
    import jadoc #uses Joint Approximate Diagonalization under Orthogonality Constraints. Method and code from  https://arxiv.org/abs/2409.02005
    F_tensor = np.empty((2,dimH,dimH),dtype=complex)
    F_tensor[0] = dH1
    F_tensor[1] = dH2

    FvR = jadoc.PerformJADOC(F_tensor, iS=dimH)
    #make FvR real
    FvR = np.real(FvR)
    print (FvR)
    return FvR