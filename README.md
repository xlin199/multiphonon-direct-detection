## Files

StructureFactor1D.py --
  1. Script calculating $S(\omega, q)$ as a function of $\omega$ for different $q$.
  2. Finished the incoherent approximation (see NSiteLatticeFixqRecursive.ipynb)

NSiteLatticeCoherentMultiprocess.py --
  1. Batch calculation of $S(\omega, q)$ for a list of $q$

TwoSite+NSiteLattice_Diagram.ipynb --
  1. Schematic diagrams and density of states plots for 2-site and N-site lattice models

TwoSiteLattice_fixq.ipynb --
  1. 2-site structure factor as a function of $\omega$ at fixed $q$ values 

NSiteLatticeStructureFactor.ipynb --
  1. N-site structure factor as a function of $\omega$ at fixed $q$ values 

NSiteLatticeFixq_2phonons.ipynb --
  1. N-site structure factor on $\omega_1$ vs $\omega_2$ plane, showing the dispersion relation of 2 phonons

NSiteLattice_2d.ipynb --
  1. N-site structure factor on ($q$, $\omega$) plane

ScatteringRate_dRdomega(dRdq, R_vs_mass)_massive(massless)Mediator.ipynb --
  1. (Partial) scattering rate and cross section for massive (massless) mediator, N-site lattice
	
ScatteringRate_dRdq_massiveMediator_3orders.ipynb --
  1. Compare dR/dq for a massive mediator computed with the first-, second-, and third-order structure factor separately (one panel each).

TwoSiteLattice_integratedS_3orders.ipynb --
  1. Compare structure factor integrated over $\omega$ fcomputed with the first-, second-, and third-order separately 



  
  ---
<!-- 
/archive/NSiteLatticeFixq.ipynb --
  1. Here focusing on plotting $S(\omega, q)$ as a function of $\omega$ for different $q$
  2. ~~Why doesn't it match with the impulse approximation?~~ used wrong DW factor. Fixed in newer notebooks. 
  3. Moved to the /archive folder

NSiteLatticeFixq_2phonons.ipynb
  1. Visualizes the phonon energy distribution for the 2-phonon structure factor, incoherent corrections included
  2. 2d plot ( $\omega_1$ vs $\omega_2$ ) with colorbar $S^{(2)}(\omega=\omega_1+\omega_2, q)$
  3. Matches the dispersion relation

NSiteLattice1D.ipynb -- 
  1. Here focusing on plotting $\int S(\omega, q) d\omega$
  2. Directly calculating the structure factor without any smearing of the phonon 
  distribution, and without using the recursion relation.
  3. Plots comparing the incoherent approx with the full coherent result - still useful! 
  4. Can do one version with energy cut and one without energy cut (N-dependent DW)
  
NSiteLattice2D.ipynb -- 
  1. Here focusing on plotting $S(\omega, q)$
  2. Still useful to have 2D plots of $S(\omega, q)$! [But first focus on the 1D plots]
  3. Currently old independent code from Momei -- would either need to perform checks against
  4. the code in NSiteLattice1D.ipynb, or can rewrite it given the latest results in calculating structure -->

---

<!-- TODO:
  1. Update notebooks for the N-site lattice 
       * $\int S(\omega, q) d\omega$ vs. q comparing incoherent and coherent
       * Plotting $S(\omega, q)$ as a function of $\omega$ for different $q$ 
          -- compare incoherent and coherent, compare impulse approx
       * 2D plot $S(\omega, q)$
       * Rate calculation (compare with single HO)
  . New notebook for reproducing single HO structure factor and rate
  1. Notebook for structure factor and rate for two-site lattice
       * $\int S(\omega, q) d\omega$ vs. q comparing incoherent and coherent
       * Plotting $S(\omega, q)$ as a function of $\omega$ for different $q$ 
          -- compare incoherent and coherent, compare impulse approx
       * 2D plot $S(\omega, q)$
       * Rate calculation (compare with single HO) -->

#### TODO
1. Rate calculation 
   - Sanity check: reproduce 2205.02250 Fig 15(a), silicon with 1 meV energy threshold (intermediate step: Fig 13 (c) $dR/d \omega$)


