import numpy as np
from scipy.special import factorial
from scipy.stats import binned_statistic

class StructureFactor1D:
    def __init__(self, nLattice,q_value,mass=26161e3, omegaNaught=10e-6, latticeConst=2*np.pi/2.28, couplingConst=28, energyThreshold = 1e-6):
        self.nLattice = nLattice # number of lattice
        # self.nPhononMax = nPhononMax
        self.mass = mass # Silicon: 26161e3 keV
        self.omegaNaught = omegaNaught # Silicon: 10e-6 keV, maximum energy 20 meV
        self.latticeConst = latticeConst # Silicon: 2*np.pi/2.28 keV^-1
        self.couplingConst = couplingConst # A, default: 28 for silicon
        self.volume = nLattice*8/(latticeConst**3) # propotional to nLattice
        self.q_value = q_value
        #### omega list with energy cut ####
        self.omegaList = 2*self.omegaNaught*np.sin(np.pi*np.arange(1, self.nLattice)/self.nLattice) # num of energy levels: nLattice-1
        filter = np.argwhere(self.omegaList >= energyThreshold)
        self.omegaList = self.omegaList[filter].reshape(-1)
        self.nuList = np.arange(1, self.nLattice)[filter]
        ################################
        self.l_diff = np.arange(1, len(self.nuList))
        self.currIntS_diag = None # integrated diagonal structure factor (the integration bins are defined in setup_DoS)
        self.sOrder = None # number of phonons corresponding to currIntS (or the order of currIntS)

    ######## Public functions ########

    def setup_DoS(self,max_omega, nBins=500):
        '''
        Defines the bin used to calculate IntS
        '''
        self._DoS_max_omega = max_omega
        self._DoS_min_omega = 0 
        self._DoS_nBins = nBins
        self._DoS_x_list, self._DoS_delta_omega = np.linspace(0, max_omega,nBins, retstep=True) 

    def init_recursion(self):
        print('Processing n = 1...')
        self.sOrder = 1
        self._s_diag_init()
        # self.get_c_offdiag_init()

    def update_s_factor(self):
        '''
        Update the integrated diagonal s factor (both diagonal and off diagonal) using the recursive relation.
        '''
        self.sOrder += 1
        print(f'Processing n = {self.sOrder}...')
        newOmega, newS = self._calc_s_diag_rec()
        self.currIntS_diag = self._integrate_s_diag(newOmega, newS)

    def get_binned_s_diag(self, xlist, delta_omega):
        ''' 
        Gives the binned diagonal (incoherent) part of s factor: int_{each bin} d(omega') S(q, omega') / (bin size)
        xlist: x values used for plotting
        delta_omega: bin size, or step size of the xlist
        --------returns--------
        ylist: array with same shape as xlist 
        '''
        result = binned_statistic(self._DoS_x_list, self.currIntS_diag, statistic="sum",bins=len(xlist), range=(min(xlist), max(xlist)))
        ylist = result.statistic/delta_omega
        return ylist
    
    ######## Helper functions ########
        
    def _DebyeWallerConst(self):
        '''
        Get W(q) appearing in the debye waller factor
        '''
        return self.q_value**2 / (4*self.mass*self.nLattice) * np.sum(1/self.omegaList)

    def _s_diag_init(self):
        '''
        Get integrated diagonal s factor of order 1
        '''
        DWfactor =  np.exp(- 2 * self._DebyeWallerConst()) 
        s_values = 2 * np.pi / self.volume * self.nLattice * DWfactor * self.couplingConst**2 *  (self.q_value**2 /(2*self.mass*self.nLattice)) / self.omegaList
        self.currIntS_diag = self._integrate_s_diag(self.omegaList, s_values)
        # return self.currIntS_diag

    def _integrate_s_diag(self,allowed_omegas, s_diag):
        '''
        Helper function. Integrate the structure factor over omega in bins defined by setup_DoS (int_{omega}^{omega+Delta omega} S d(omega')) 
        --------returns--------
        intS: (DoS_nBins,) array, sum of s factor
        '''
        result = binned_statistic(allowed_omegas, s_diag, statistic="sum",bins=self._DoS_nBins, range=(self._DoS_min_omega, self._DoS_max_omega))
        # np.nan_to_num(result.statistic, nan=0.0)
        return result.statistic
    
    def _calc_s_diag_rec(self,):
        '''
        Helper function. 
        Use the recursive relation to calculate the structure factor of nPhonon-th order
        --------returns--------
        newOmega, newS: both are (a*len(omegaList),) array
        '''
        newIntS = (self.q_value**2 /(2*self.mass*self.nLattice))/self.sOrder * self.currIntS_diag[:, np.newaxis] / self.omegaList 
        newOmega = self._DoS_x_list[:, np.newaxis] + self.omegaList 
        return newOmega.T.reshape(-1,), newIntS.T.reshape(-1,)
    




########## The functions below are used for direct calculations (no recursive relation) ##############
    

    # def partition(self, max_range, max_sum):
    #     max_range = np.asarray(max_range, dtype = int).ravel()        
    #     if(max_range.size == 1):
    #         return np.arange(min(max_range[0],max_sum) + 1, dtype = int).reshape(-1,1)
    #     P = self.partition(max_range[1:], max_sum)
    #     # S[i] is the largest summand we can place in front of P[i]            
    #     S = np.minimum(max_sum - P.sum(axis = 1), max_range[0])
    #     offset, sz = 0, S.size
    #     out = np.empty(shape = (sz + S.sum(), P.shape[1]+1), dtype = int)
    #     out[:sz,0] = 0
    #     out[:sz,1:] = P
    #     for i in range(1, max_range[0]+1):
    #         ind, = np.nonzero(S)
    #         offset, sz = offset + sz, ind.size
    #         out[offset:offset+sz, 0] = i
    #         out[offset:offset+sz, 1:] = P[ind]
    #         S[ind] -= 1
    #     return out
    
    # def partitionwfilter(self, max_range, max_sum):
    #     arr = self.partition(max_range, max_sum)
    #     return arr[np.argwhere(np.sum(arr, axis=1) == max_sum)].reshape(-1, len(max_range))
    
    # def getJVector(self, nPhonon):
    #     '''
    #     For direct calculations (no recursive relation)
    #     Get all posible sets of j_nu given the total number of phonons. 
    #     nPhonon: int
    #     '''
    #     jVector = self.partitionwfilter(nPhonon*np.ones(len(self.omegaList)), nPhonon)
    #     return jVector
    
    # def get_diag(self, nPhonon):
    #     '''
    #     For direct calculations (no recursive relation)
    #     '''
    #     factor = 2 * np.pi / self.volume * np.exp(- 2 * self.getDebyeWallerConst()) 
    #     jVector = self.getJVector(nPhonon)
    #     prod_part = np.prod(1/(factorial(jVector)*(self.omegaList ** jVector)), axis=1)
    #     return factor * self.nLattice * self.couplingConst**2 * (self.q_value**2 /(2*self.mass*self.nLattice))**nPhonon * prod_part

    # def get_allowed_omega(self, nPhonon):
    #     '''
    #     For direct calculations (no recursive relation)
    #     '''
    #     jVector = self.getJVector(nPhonon)
    #     allowed_omega = np.sum(np.multiply(self.omegaList, jVector),axis=1)
    #     return allowed_omega

    # def get_binned_s_diag(self,nPhonon, xlist, delta_omega):
    #     allowed_omegas = self.get_allowed_omega(nPhonon)
    #     s_diag = self.get_diag(nPhonon)
    #     result = binned_statistic(allowed_omegas, s_diag, statistic="sum",bins=len(xlist), range=(min(xlist), max(xlist)))
    #     ylist = result.statistic/delta_omega
    #     return xlist, ylist

##################### End #######################



########### N-independent DW factor ###########

class StructureFactor1D_NIndepDW(StructureFactor1D):
    def getDebyeWallerConst(self):
        return self.q_value**2 / (4*self.mass*8.25e-6) # take W(q) = q^2/(4 m omega_DW) 
    



################# Sample code #######################
