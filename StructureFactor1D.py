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
        self.volume = nLattice*(latticeConst**3)/4 # propotional to nLattice
        self.q_value = q_value
        #### omega list with energy cut ####
        self.omegaList = 2*self.omegaNaught*np.sin(np.pi*np.arange(1, self.nLattice)/self.nLattice) # num of energy levels: nLattice-1
        filter = np.argwhere(self.omegaList >= energyThreshold)
        self.omegaList = self.omegaList[filter].reshape(-1)
        self.nuList = np.arange(1, self.nLattice)[filter]
        ################################
        self.l_diff = np.arange(1, nLattice)
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
        self.sOrder = 1
        print(f'Processing n = {self.sOrder}...')
        self._s_diag_init()
        self._c_offdiag_init()

    def update_s_factor(self):
        '''
        Update the integrated diagonal s factor (both diagonal and off diagonal) using the recursive relation.
        '''
        self.sOrder += 1
        print(f'Processing n = {self.sOrder}...')
        # update diagonal term
        newOmega, newS = self._calc_s_diag_rec()
        self.currIntS_diag = self._integrate_s_diag(newOmega, newS)
        # update off diagonal term
        newOmega, newC = self._calc_c_offdiag_rec()
        self.currIntC_offdiag = self._integrate_c_offdiag(newOmega, newC)

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
        return np.nan_to_num(ylist, nan=0.0) # cleaning up the nan
    
    def get_binned_s_offdiag(self, xlist, delta_omega):
        # sum over c_ll' to get s factor
        mult = (self.nLattice - self.l_diff)*np.real(self.currIntC_offdiag)
        sfactor_sum = 2* (self.couplingConst**2) * np.sum(mult, axis=1)
        result = binned_statistic(self._DoS_x_list, sfactor_sum, statistic="sum",bins=len(xlist), range=(min(xlist), max(xlist)))
        ylist = result.statistic/delta_omega
        return np.nan_to_num(ylist, nan=0.0) # cleaning up the nan
    
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

    def _c_offdiag_init(self,): # l' minus l
        '''
        Calculate the first c_ll' of the recursive relation. 
        currIntC_offdiag: 2d array with shape (# of nu, # of xi)
        '''
        DWfactor =  np.exp(- 2 * self._DebyeWallerConst()) 
        offdiag_factor = np.exp(1j*(self.q_value*self.latticeConst-2*np.pi* self.nuList/self.nLattice)* self.l_diff) 
        c_values =  2 * np.pi / self.volume * DWfactor *  (self.q_value**2 /(2*self.mass*self.nLattice)) * offdiag_factor/ self.omegaList[:,np.newaxis]
        self.currIntC_offdiag = self._integrate_c_offdiag(self.omegaList, c_values)

    def _integrate_c_offdiag(self,allowed_omegas, c_offdiag):
        '''
        Get the integrated version of c_ll'
        '''
        result = binned_statistic(allowed_omegas, c_offdiag.T, statistic="sum",bins=self._DoS_nBins, range=(self._DoS_min_omega, self._DoS_max_omega))
        return result.statistic.T

    def _calc_c_offdiag_rec(self,):
        exp_factor = np.exp(-2j*np.pi*self.nuList*self.l_diff/self.nLattice)
        rec_factor = exp_factor/self.omegaList[:, np.newaxis]
        newIntC = (self.q_value**2 /(2*self.mass*self.nLattice))/self.sOrder * rec_factor[:, np.newaxis, :]*self.currIntC_offdiag
        newOmega = self._DoS_x_list[:, np.newaxis] + self.omegaList 
        return newOmega.T.reshape(-1,), newIntC.reshape(-1,self.l_diff.shape[0])
    




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







################# Sample code #######################

# multplier = 1.8
# q_value  = multplier*np.sqrt(2*26161e3*13.13e-6) # uses averaged omega (not inverse of averaged inverse)

# nLattice = 200
# test = StructureFactor1D(nLattice = nLattice,q_value=q_value, couplingConst=28,)
# # color = ['r', 'g', 'b', 'y', 'c', 'tab:blue', 'tab:red', 'tab:yellow', "tab:brown", 'tab:purple']
# max_omega = 150*1e-6 # keV
# xlist, delta_omega = np.linspace(0, max_omega,100, retstep=True) 
# nPhonon = 10

# # impulse approx
# deltasq = q_value**2 * 13.13e-6/(2*test.mass)
# ia =  28**2 / (8/(test.latticeConst**3)) * np.sqrt(2*np.pi/deltasq) * np.exp(- (xlist-(q_value**2/(2*test.mass)))**2 / (2*deltasq)) # uses averaged omega (not inverse of averaged inverse)

# sum_diag = np.zeros(xlist.shape)
# sum_offdiag = np.zeros(xlist.shape)
# fig, ax = plt.subplots(figsize=(11,6), ncols=2, sharey=True)
# test.setup_DoS(max_omega, nBins=1000)
# test.init_recursion()
# # get the binned s factor of order 1
# s_diag = test.get_binned_s_diag(xlist, delta_omega)
# s_offdiag = test.get_binned_s_offdiag(xlist, delta_omega)
# # update the sum of s factor
# sum_diag += s_diag 
# sum_offdiag += s_offdiag
# # plot the current s factor
# print(f'Plotting n = 1...')
# ax[0].plot(xlist*1e6, s_diag,label=f'n = 1',  linestyle='--') 
# ax[1].plot(xlist*1e6, s_diag+s_offdiag,label=f'n = 1',) 
# for i in range(2,nPhonon+1):
#     test.update_s_factor()
#     # get the binned s factor of order i
#     s_diag = test.get_binned_s_diag(xlist, delta_omega)
#     s_offdiag = test.get_binned_s_offdiag(xlist, delta_omega)
#     # update the sum of s factor
#     sum_diag += s_diag 
#     sum_offdiag += s_offdiag
#     # plot the current s factor
#     print(f'Plotting n = {i}...')
#     ax[0].plot(xlist*1e6, s_diag,label=f'n = {i}',  linestyle='--') 
#     ax[1].plot(xlist*1e6, s_diag+s_offdiag,label=f'n = {i}',) 
#     # calculate the s factor of next order
# formatter = ScalarFormatter(useMathText=True)
# formatter.set_powerlimits((8,8))
# ax[0].yaxis.set_major_formatter(formatter)
# ax[1].yaxis.set_major_formatter(formatter)
# # plot the impulse approximations
# ax[0].plot(xlist*1e6, ia, linestyle=":", color="gray", label='Impulse')
# ax[1].plot(xlist*1e6, ia, linestyle=":", color="gray", label='Impulse')
# # plot the sum of s for all orders
# ax[0].plot(xlist*1e6,sum_diag, label="Total (incoherent)",linestyle="--",color="black",)
# ax[1].plot(xlist*1e6,sum_diag+sum_offdiag, label="Total",color="black",)
# # labels & titles
# txt1 = r"$\omega_{DW}$" #r"exp$(\frac{q^2}{4m\omega_{DW}})$"
# txt2 = r"$\sqrt{2m\bar\omega}$" 
# plt.ylim(-0.1e8, 7.2e8)
# plt.suptitle(f'N-independent DW, Recursive\nN = {nLattice}, {txt1} = {omega_dw*1e6:.2f} meV, q = {multplier}{txt2}')# , {txt} = {test.omegaNaught*1e6} meV
# ax[0].set_ylabel('Structure factor')
# ax[0].set_xlabel(r'$\omega$[meV]')
# ax[1].set_xlabel(r'$\omega$[meV]')
# ax[0].legend()
# ax[1].legend()
# ax[0].set_title('Incoherent Approximation')
# ax[1].set_title('Includes Coherent Part')
# plt.tight_layout()
# # plt.savefig(f"fullsfactorwrtomega_q{int(q_value)}_{nLattice}lattice_NindepDW_energyCut_recur_DoS", bbox_inches="tight")