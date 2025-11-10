import numpy as np
import sys
from multiprocessing import Pool, cpu_count
from functools import partial
from datetime import datetime
sys.path.append(r"C:/Users/xylin/Documents/Researches/multiphonon-direct-detection/StructureFactor1D.py")
from StructureFactor1D import StructureFactor1D


def multi_processing_helper(q_value, nPhonon, 
                            omegalength, nLattice, max_omega): 
    sFacCalc = StructureFactor1D(nLattice = nLattice,q_value=q_value, couplingConst=28,)
    xlist, delta_omega = np.linspace(0, max_omega,omegalength, retstep=True) 
    s_tot = np.zeros(xlist.shape)
    s_diag_tot = np.zeros(xlist.shape)
    sFacCalc.setup_DoS(max_omega, nBins=1500)
    sFacCalc.init_recursion(print_message=False)
    s_diag = sFacCalc.get_binned_s_diag(xlist, delta_omega)
    s_offdiag = sFacCalc.get_binned_s_offdiag(xlist, delta_omega)
    # update the sum of s factor
    s_diag_tot += s_diag 
    s_tot += s_diag 
    s_tot += s_offdiag
    for i in range(2, nPhonon+1):
        sFacCalc.update_s_factor(print_message=False)
        # get the binned s factor of order i
        s_diag = sFacCalc.get_binned_s_diag(xlist, delta_omega)
        s_offdiag = sFacCalc.get_binned_s_offdiag(xlist, delta_omega)
        # update the sum of s factor
        s_diag_tot += s_diag 
        s_tot += s_diag 
        s_tot += s_offdiag
    return s_diag_tot, s_tot

omegalength=1200
max_omega = 250*1e-6
nLattice=120
phonon_func = partial(multi_processing_helper, omegalength=omegalength, nLattice=nLattice, max_omega = max_omega)
nPhonon10func = partial(phonon_func, nPhonon = 10)
nPhonon20func = partial(phonon_func, nPhonon = 20)

if __name__ == "__main__":
    omegalist = np.linspace(0, max_omega,omegalength, ) 
    qscale = np.sqrt(2*26161e3*13.13e-6) # uses averaged omega (not inverse of averaged inverse)
    qlenth = 1200
    qlist = np.linspace(0.0004, 2.5*qscale, qlenth) # 2.5 sqrt(2m omega_bar)
    print("Start time:", datetime.now())
    with Pool(cpu_count()//2) as pool:
        allresult10phonon = pool.map(nPhonon10func, qlist[:qlenth//2])
        allresult20phonon = pool.map(nPhonon20func, qlist[qlenth//2:])
    res10diag, res10 = zip(*allresult10phonon)
    res20diag, res20 = zip(*allresult20phonon)
    finalresdiag = np.concatenate((res10diag, res20diag))
    finalres = np.concatenate((res10, res20))
     # shape: (q length, omega length)
    print("End time:", datetime.now())
    print(finalres.shape)
    np.savez('testMultiprocess',s_tot=finalres, s_diag=finalresdiag, qlist=qlist, omegalist=omegalist)


