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
        self.omegaList = 2*self.omegaNaught*np.sin(np.pi*np.arange(1, self.nLattice)/self.nLattice) # num of energy levels: nLattice-1
        filter = np.argwhere(self.omegaList >= energyThreshold)
        self.omegaList = self.omegaList[filter].reshape(-1)
        self.nuList = np.arange(1, self.nLattice)[filter]

    def getDebyeWallerConst(self):
        '''
        Get W(q) appearing in the debye waller factor
        '''
        return self.q_value**2 / (4*self.mass*self.nLattice) * np.sum(1/self.omegaList)
        

    def get_diag_ini(self):
        '''
        Get IntS of order 1
        '''
        print('Processing n = 1...')
        self.sOrder = 1
        DWfactor =  np.exp(- 2 * self.getDebyeWallerConst()) 
        s_values = 2 * np.pi / self.volume * self.nLattice * DWfactor * self.couplingConst**2 *  (self.q_value**2 /(2*self.mass*self.nLattice)) / self.omegaList
        self.currIntS = self.get_int_s_diag(np.row_stack((self.omegaList, s_values)))
        # return self.currIntS

    def get_diag_rec(self,):
        '''
        Use the recursive relation to calculate the structure factor of nPhonon-th order
        prevS: (2,a) array, first row: omega list 
                            second row: structure factor integrated over omega in bins defined by setup_DoS (int_{omega}^{omega+Delta omega} S d(omega'))
        --------returns--------
        newS: (2,a*len(omegaList)) array
        '''
        newIntS = (self.q_value**2 /(2*self.mass*self.nLattice))/self.sOrder * self.currIntS[1,:, np.newaxis] / self.omegaList 
        newOmega = self.currIntS[0,:, np.newaxis] + self.omegaList 
        return np.row_stack((newOmega.T.reshape(-1,), newIntS.T.reshape(-1,)))
    
    def setup_DoS(self,max_omega, nBins=500):
        '''
        Helper function. Defines the bin used to calculate int_S
        '''
        self._DoS_max_omega = max_omega
        self._DoS_min_omega = 0 
        self._DoS_nBins = nBins
        self._DoS_x_list, self._DoS_delta_omega = np.linspace(0, max_omega,nBins, retstep=True) 


    def get_int_s_diag(self,sFactor):
        '''
        Helper function. Integrate the structure factor over omega in bins defined by setup_DoS (int_{omega}^{omega+Delta omega} S d(omega')) 
        --------returns--------
        intS: (2,DoS_nBins) array
        '''
        allowed_omegas = sFactor[0,:]
        s_diag = sFactor[1,:]
        result = binned_statistic(allowed_omegas, s_diag, statistic="sum",bins=self._DoS_nBins, range=(self._DoS_min_omega, self._DoS_max_omega))
        # np.nan_to_num(result.statistic, nan=0.0)
        return np.row_stack((self._DoS_x_list,result.statistic))
    

    def update_s_diag(self):
        '''
        Update the integrated diagonal s factor using the recursive relation.
        '''
        self.sOrder += 1
        print(f'Processing n = {self.sOrder}...')
        sFactor = self.get_diag_rec()
        self.currIntS = self.get_int_s_diag(sFactor)
    

    def get_binned_s_diag(self, xlist, delta_omega):
        ''' 
        Call this function to get the binned s factor: int_{each bin} d(omega') S(q, omega') / (bin size)
        xlist: x values used for plotting
        delta_omega: bin size, or step size of the xlist
        --------returns--------
        ylist: array with same shape as xlist 
        '''
        allowed_omegas = self.currIntS[0,:]
        s_diag = self.currIntS[1,:]
        result = binned_statistic(allowed_omegas, s_diag, statistic="sum",bins=len(xlist), range=(min(xlist), max(xlist)))
        ylist = result.statistic/delta_omega
        return ylist


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

# omega_dw = 8.194e-6
# multplier = 2.5
# q_value  = multplier*np.sqrt(2*26161e3*13.124e-6) # uses averaged omega (not inverse of averaged inverse)


# nLattice = 1000
# test = StructureFactor1D_NIndepDW(nLattice = nLattice,q_value=q_value, couplingConst=28,)
# fmt = ['r-', 'g-', 'b-', 'y-', 'c-']
# max_omega = 150*1e-6 # keV
# xlist, delta_omega = np.linspace(0, max_omega,50, retstep=True) 
# nPhonon = 40

# # impulse approx
# deltasq = q_value**2 * 13.124e-6/(2*test.mass)
# ia =  28**2 / (8/(test.latticeConst**3)) * np.sqrt(2*np.pi/deltasq) * np.exp(- (xlist-(q_value**2/(2*test.mass)))**2 / (2*deltasq)) # uses averaged omega (not inverse of averaged inverse)

# sum = np.zeros(xlist.shape)
# fig, ax = plt.subplots(figsize=(8,6))
# test.setup_DoS(max_omega, nBins=2000)
# test.get_diag_ini()
# for i in range(1,nPhonon):
#     s_diag = test.get_binned_s_diag(xlist, delta_omega) # get the binned s factor of order i
#     s_diag = np.nan_to_num(s_diag, nan=0.0) # cleaning up the nan
#     sum += s_diag # update the sum of s factor
#     print(f'Plotting n = {i}...')
#     # plot the current s factor
#     ax.plot(xlist*1e6, s_diag,label=f'n = {i}') 
#     formatter = ScalarFormatter(useMathText=True)
#     formatter.set_powerlimits((8,8))
#     ax.yaxis.set_major_formatter(formatter)
#     # calculate the s factor of next order
#     test.update_s_diag()
# # plotting the last s factor
# print(f'Plotting n = {nPhonon}...')
# ax.plot(xlist*1e6, s_diag,label=f'n = {nPhonon}')
# formatter = ScalarFormatter(useMathText=True)
# formatter.set_powerlimits((8,8))
# ax.yaxis.set_major_formatter(formatter)
# # plot the impulse approximation
# ax.plot(xlist*1e6, ia, linestyle="--", color="gray", label='Impulse')
# # plot the sum of s for all orders
# ax.plot(xlist*1e6,sum, label="Total",color="black")
# # labels & titles
# txt1 = r"$\omega_{DW}$" #r"exp$(\frac{q^2}{4m\omega_{DW}})$"
# txt2 = r"$\sqrt{2m\bar\omega}$" 
# plt.ylim(-0.1e8, 7.2e8)
# plt.title(f'DoS + Incoherent approx, N-independent DW, Recursive\nN = {nLattice}, {txt1} = {omega_dw*1e6} meV, q = {multplier}{txt2}')# , {txt} = {test.omegaNaught*1e6} meV
# plt.ylabel('Structure factor')
# plt.xlabel(r'$\omega$[meV]')
# plt.legend()
# plt.savefig(f"sfactorwrtomega_q{int(q_value)}_{nLattice}lattice_NindepDW_energyCut_recur_DoS", bbox_inches="tight")