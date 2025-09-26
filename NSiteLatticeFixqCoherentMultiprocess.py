import numpy as np
import sys
from multiprocessing import Pool, cpu_count
from functools import partial
sys.path.append(r"C:/Users/xylin/Documents/Researches/multiphonon-direct-detection/StructureFactor1D.py")
from StructureFactor1D import StructureFactor1D


def multi_processing_helper(q_value, nPhonon, 
                            omegalength=100, nLattice=200, max_omega = 150*1e-6): 
    sFacCalc = StructureFactor1D(nLattice = nLattice,q_value=q_value, couplingConst=28,)
    xlist, delta_omega = np.linspace(0, max_omega,omegalength, retstep=True) 
    s_tot = np.zeros(xlist.shape)
    sFacCalc.setup_DoS(max_omega, nBins=500)
    sFacCalc.init_recursion(print_message=False)
    s_diag = sFacCalc.get_binned_s_diag(xlist, delta_omega)
    s_offdiag = sFacCalc.get_binned_s_offdiag(xlist, delta_omega)
    # update the sum of s factor
    s_tot += s_diag 
    s_tot += s_offdiag
    for i in range(2, nPhonon+1):
        sFacCalc.update_s_factor(print_message=False)
        # get the binned s factor of order i
        s_diag = sFacCalc.get_binned_s_diag(xlist, delta_omega)
        s_offdiag = sFacCalc.get_binned_s_offdiag(xlist, delta_omega)
        # update the sum of s factor
        s_tot += s_diag 
        s_tot += s_offdiag
    return s_tot

phonon_func = partial(multi_processing_helper, omegalength=100, nLattice=200, max_omega = 150*1e-6)
nPhonon5func = partial(phonon_func, nPhonon = 3)
nPhonon10func = partial(phonon_func, nPhonon = 5)


if __name__ == "__main__":
    num = np.sqrt(2*26161e3*13.13e-6) # uses averaged omega (not inverse of averaged inverse)
    qlist1 = [num,num*1.5]
    qlist2 = [num*2,num*3]
    with Pool(cpu_count()) as pool:
        result5phonon = pool.map(nPhonon5func, qlist1)
        result10phonon = pool.map(nPhonon10func, qlist2)
    final_results = np.concatenate((result5phonon, result10phonon))
    print(final_results.shape)


